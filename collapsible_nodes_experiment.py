import os
import random
import sys
from time import perf_counter
from pathlib import Path
import logging
import warnings

import pandas as pd
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
        os.add_dll_directory(r"D:\\vcpkg-master\\installed\\x64-windows\\bin"),
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
    return node_names, node_to_idx


def model_to_igraph(true_model, node_to_idx):
    ig_graph = ig.Graph(directed=True)
    ig_graph.add_vertices(len(node_to_idx))
    ig_graph.add_edges([(node_to_idx[u], node_to_idx[v]) for u, v in true_model.edges()])
    return ig_graph


def run_experiment(r_repeat, r_size, file_names, seed=42):
    if r_repeat <= 0:
        raise ValueError("r_repeat must be positive")
    if r_size <= 0:
        raise ValueError("r_size must be positive")

    random.seed(seed)
    np.random.seed(seed)

    result_rows = []

    for file in file_names:
        true_model = get_example_model(file)
        node_names, node_to_idx = build_node_mappings(true_model)
        ig_graph = model_to_igraph(true_model, node_to_idx)

        if r_size > len(node_names):
            raise ValueError(f"r_size={r_size} 超出节点数量范围(1, {len(node_names)})，file={file}")

        ang_times = []
        collapsed_times = []
        ang_sizes = []
        collapsed_sizes = []

        for _ in range(r_repeat):
            r_names = random.sample(node_names, r_size)
            r_idx = [node_to_idx[node] for node in r_names]

            t_ang_start = perf_counter()
            ang_nodes_idx = decom_h.dag_get_ancestors(ig_graph, r_idx)
            ang_times.append(perf_counter() - t_ang_start)
            ang_sizes.append(len(ang_nodes_idx))

            t_collapsed_start = perf_counter()
            collapsed_nodes_idx = decom_h.dag_get_minimal_collapsible(ig_graph, r_idx)
            collapsed_times.append(perf_counter() - t_collapsed_start)
            collapsed_sizes.append(len(collapsed_nodes_idx))

        avg_ang_time = sum(ang_times) / r_repeat
        avg_collapsed_time = sum(collapsed_times) / r_repeat
        avg_ang_nodes = sum(ang_sizes) / r_repeat
        avg_collapsed_nodes = sum(collapsed_sizes) / r_repeat

        print(
            f"file={file}, r_repeat={r_repeat}, r_size={r_size}, "
            f"avg_ang_nodes={avg_ang_nodes:.1f}, avg_ang_time={avg_ang_time:.6f}s, "
            f"avg_collapsed_nodes={avg_collapsed_nodes:.1f}, avg_collapsed_time={avg_collapsed_time:.6f}s"
        )

        result_rows.append({
            "file": file,
            "network_nodes": len(node_names),
            "network_edges": true_model.number_of_edges(),
            "r_repeat": r_repeat,
            "r_size": r_size,
            "avg_ang_nodes": avg_ang_nodes,
            "avg_ang_time": avg_ang_time,
            "avg_collapsed_nodes": avg_collapsed_nodes,
            "avg_collapsed_time": avg_collapsed_time,
            "relative_submodel_size": avg_collapsed_nodes / avg_ang_nodes,
        })

    return pd.DataFrame(result_rows)


def main():
    # 用于快速切换的测试文件集合；需要做正式实验时可以换成更完整的网络列表。
    file_names = ["hailfinder", "hepar2", "andes", "munin"]
    # 对不同 R 组进行重复实验的次数，最后会对这部分取平均。
    r_repeat = 100
    # 每次随机抽取的 R 的节点数。
    r_size = 5
    # 全局随机种子，保证别人运行时可以复现相同的随机序列。
    seed = 20260410
    print(f"starting experiment with r_repeat={r_repeat}, r_size={r_size}, seed={seed}")

    result = run_experiment(
        r_repeat=r_repeat,
        r_size=r_size,
        file_names=file_names,
        seed=seed,
    )

    result_dir = base_dir / "result"
    result_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(result_dir / "collapsible_nodes_experiment.csv", index=False)
    print(result.head())


if __name__ == "__main__":
    main()