import os
import random
import sys
from time import perf_counter
from pathlib import Path
import logging
import warnings

import pandas as pd
from pgmpy.estimators import MaximumLikelihoodEstimator
from pgmpy.inference import VariableElimination
from pgmpy.models import BayesianNetwork
from pgmpy.sampling import BayesianModelSampling
import numpy as np
from pgmpy.utils import get_example_model

pgmpy_logger = logging.getLogger("pgmpy")
pgmpy_logger.setLevel(logging.ERROR)
pgmpy_logger.propagate = False
warnings.filterwarnings("ignore", category=RuntimeWarning, module="networkx.utils.backends")
os.environ["NUMEXPR_MAX_THREADS"] = "32"

base_dir = Path(__file__).resolve().parent
build_dir = base_dir / "build"
sys.path.insert(0, str(build_dir))
print(f"Using Python executable: {sys.executable}")
EXPECTED_PYTHON = (3, 9)
if sys.version_info[:2] != EXPECTED_PYTHON:
    raise RuntimeError(
        "当前解释器是 Python "
        f"{sys.version_info.major}.{sys.version_info.minor}，"
        "但 decom_h.pyd 是按 Python 3.9 编译。"
        "请使用 D:\\Anaconda3\\python.exe 运行，"
        "或重新用当前 Python 版本重编译扩展。"
    )

os.environ["PATH"] = str(build_dir) + os.pathsep + r"D:\\vcpkg-master\\installed\\x64-windows\\bin" + os.pathsep + os.environ.get("PATH", "")

if hasattr(os, "add_dll_directory"):
    _dll_dirs = [
        os.add_dll_directory(str(build_dir)),
        os.add_dll_directory(r"D:\\vcpkg-master\\installed\\x64-windows\bin"),
    ]

import igraph as ig
try:
    import decom_h
except ImportError as exc:
    raise ImportError(
        "导入 decom_h 失败。常见原因是 Python 版本与 pyd ABI 不匹配或缺少依赖 DLL。"
        "当前项目构建配置为 Python 3.9。"
        "请优先使用 D:\\Anaconda3\\python.exe 运行该脚本。"
    ) from exc

def build_node_mappings(true_model):
    node_names = list(true_model.nodes())
    node_to_idx = {name: idx for idx, name in enumerate(node_names)}
    idx_to_node = {idx: name for name, idx in node_to_idx.items()}
    return node_names, node_to_idx, idx_to_node


def model_to_igraph(true_model, node_to_idx):
    ig_graph = ig.Graph(directed=True)
    ig_graph.add_vertices(len(node_to_idx))
    ig_graph.add_edges([(node_to_idx[u], node_to_idx[v]) for u, v in true_model.edges()])
    return ig_graph


def map_nodes_to_names(node_ids, idx_to_node):
    return [idx_to_node[idx] for idx in node_ids]


def induced_edges_as_names(true_model, node_set):
    return [(u, v) for u, v in true_model.edges() if u in node_set and v in node_set]


def run_experiment(r_repeat, sample_repeat, r_size, file_names, sample_sizes, seed=42):
    if r_repeat <= 0:
        raise ValueError("r_repeat must be positive")
    if sample_repeat <= 0:
        raise ValueError("sample_repeat must be positive")
    if r_size <= 0:
        raise ValueError("r_size must be positive")
    random.seed(seed)
    np.random.seed(seed)

    result_rows = []

    for file in file_names:
        true_model = get_example_model(file)
        node_names, node_to_idx, idx_to_node = build_node_mappings(true_model)
        ig_graph = model_to_igraph(true_model, node_to_idx)

        if r_size > len(node_names):
            raise ValueError(f"r_size={r_size} 超出节点数量范围(1, {len(node_names)})，file={file}")

        for sample_size in sample_sizes:
            r_case_rows = []

            for _ in range(r_repeat):
                r_names = random.sample(node_names, r_size)
                r_idx = [node_to_idx[node] for node in r_names]

                t_ang_pre = perf_counter()
                ang_nodes_idx = decom_h.dag_get_ancestors(ig_graph, r_idx)
                ang_nodes = map_nodes_to_names(ang_nodes_idx, idx_to_node)
                ang_model = BayesianNetwork()
                ang_model.add_nodes_from(ang_nodes)
                ang_model.add_edges_from(induced_edges_as_names(true_model, set(ang_nodes)))
                ang_pre_time = perf_counter() - t_ang_pre

                t_collapsed_pre = perf_counter()
                collapsed_nodes_idx = decom_h.dag_get_minimal_collapsible(ig_graph, r_idx)
                collapsed_nodes = map_nodes_to_names(collapsed_nodes_idx, idx_to_node)
                collapsed_model = BayesianNetwork()
                collapsed_model.add_nodes_from(collapsed_nodes)
                collapsed_model.add_edges_from(induced_edges_as_names(true_model, set(collapsed_nodes)))
                collapsed_pre_time = perf_counter() - t_collapsed_pre

                ang_sample_time = 0.0
                collapsed_sample_time = 0.0

                for _ in range(sample_repeat):
                    df = BayesianModelSampling(true_model).forward_sample(size=sample_size, show_progress=False)

                    t_ang_start = perf_counter()
                    ang_model.fit(df[ang_nodes], estimator=MaximumLikelihoodEstimator)
                    inference_ang = VariableElimination(ang_model)
                    _ = inference_ang.query(variables=r_names, evidence={}, show_progress=False)
                    ang_sample_time += perf_counter() - t_ang_start

                    t_collapsed_start = perf_counter()
                    collapsed_model.fit(df[collapsed_nodes], estimator=MaximumLikelihoodEstimator)
                    inference_collapsed = VariableElimination(collapsed_model)
                    _ = inference_collapsed.query(variables=r_names, evidence={}, show_progress=False)
                    collapsed_sample_time += perf_counter() - t_collapsed_start

                ang_total_time = ang_pre_time + ang_sample_time
                collapsed_total_time = collapsed_pre_time + collapsed_sample_time

                r_case_rows.append({
                    "ang_pre_time": ang_pre_time,
                    "collapsed_pre_time": collapsed_pre_time,
                    "ang_sample_time": ang_sample_time,
                    "collapsed_sample_time": collapsed_sample_time,
                    "ang_total_time": ang_total_time,
                    "collapsed_total_time": collapsed_total_time,
                    "ang_nodes": len(ang_nodes),
                    "collapsed_nodes": len(collapsed_nodes),
                })

            avg_ang_pre_time = sum(row["ang_pre_time"] for row in r_case_rows) / r_repeat
            avg_collapsed_pre_time = sum(row["collapsed_pre_time"] for row in r_case_rows) / r_repeat
            avg_ang_sample_time = sum(row["ang_sample_time"] for row in r_case_rows) / r_repeat
            avg_collapsed_sample_time = sum(row["collapsed_sample_time"] for row in r_case_rows) / r_repeat
            avg_ang_total_time = sum(row["ang_total_time"] for row in r_case_rows) / r_repeat
            avg_collapsed_total_time = sum(row["collapsed_total_time"] for row in r_case_rows) / r_repeat
            avg_ang_nodes = sum(row["ang_nodes"] for row in r_case_rows) / r_repeat
            avg_collapsed_nodes = sum(row["collapsed_nodes"] for row in r_case_rows) / r_repeat

            """
            print(
                f"file={file}, sample_size={sample_size}, r_repeat={r_repeat}, sample_repeat={sample_repeat}, "
                f"r_size={r_size}, ang_dim=({avg_ang_nodes:.1f}), ang_pre={avg_ang_pre_time:.4f}s, "
                f"ang_sample={avg_ang_sample_time:.4f}s, ang_total={avg_ang_total_time:.4f}s, "
                f"collapsed_dim=({avg_collapsed_nodes:.1f}), collapsed_pre={avg_collapsed_pre_time:.4f}s, "
                f"collapsed_sample={avg_collapsed_sample_time:.4f}s, collapsed_total={avg_collapsed_total_time:.4f}s"
            )
            """

            result_rows.append({
                "file": file,
                "sample_size": sample_size,
                "r_repeat": r_repeat,
                "sample_repeat": sample_repeat,
                "r_size": r_size,
                "ang_pre_time": avg_ang_pre_time,
                "collapsed_pre_time": avg_collapsed_pre_time,
                "ang_sample_time": avg_ang_sample_time,
                "collapsed_sample_time": avg_collapsed_sample_time,
                "ang_total_time": avg_ang_total_time,
                "collapsed_total_time": avg_collapsed_total_time,
                "ang_nodes": avg_ang_nodes,
                "collapsed_nodes": avg_collapsed_nodes,
            })

        print(f"completed {file} experiment")
    

    return pd.DataFrame(result_rows)


def main():
    # 用于快速切换的测试文件集合；需要做正式实验时可以换成更完整的网络列表。
    file_names = ["hailfinder", "hepar2", "andes", "munin"]
    # 每个样本量下重复采样的规模，代表数据集大小。
    sample_sizes = [500, 5000]
    # 对不同 R 组进行重复实验的次数，最后会对这部分取平均。
    r_repeat = 10
    # 同一个 R 下，重复抽样和学习/推理的次数。
    sample_repeat = 50
    # 每次随机抽取的 R 的节点数。
    r_size = 5
    # 全局随机种子，保证别人运行时可以复现相同的随机序列。
    seed = 20260410
    print(f"starting experiment with r_repeat={r_repeat}, sample_repeat={sample_repeat}, r_size={r_size}, ")

    result = run_experiment(
        r_repeat=r_repeat,
        sample_repeat=sample_repeat,
        r_size=r_size,
        file_names=file_names,
        sample_sizes=sample_sizes,
        seed=seed,
    )

    result_dir = base_dir / "result"
    result_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(result_dir / "collapsible_experiment.csv", index=False)
    print(result.head())


if __name__ == "__main__":
    main()