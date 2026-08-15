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


def extract_state_names(true_model, variables):
    state_names = {}
    for variable in variables:
        cpd = true_model.get_cpds(variable)
        if cpd is not None and hasattr(cpd, "state_names"):
            state_names[variable] = list(cpd.state_names[variable])
    return state_names


def inject_missing_values(df, miss_proportion, covariate_cols, entry_missing_rate=0.1):
    """Bernoulli row contamination + Bernoulli per-covariate missingness in contaminated rows."""
    if miss_proportion <= 0 or not covariate_cols:
        return df

    df_incomplete = df.copy()
    n_rows = len(df_incomplete)
    covariate_cols = list(covariate_cols)

    for row_pos in range(n_rows):
        # Row-level Bernoulli: whether this sample gets contaminated.
        if random.random() >= miss_proportion:
            continue

        # Cell-level Bernoulli within contaminated rows.
        for col in covariate_cols:
            if random.random() < entry_missing_rate:
                df_incomplete.iat[row_pos, df_incomplete.columns.get_loc(col)] = np.nan

    return df_incomplete


def compute_l1_distance(query_est, query_true):
    # Align by axis order first, then compare array entries.
    true_vars = list(query_true.variables)
    est_vars = list(query_est.variables)

    if set(true_vars) != set(est_vars):
        raise ValueError("query_est and query_true variables are not identical")

    perm = [est_vars.index(var) for var in true_vars]
    est_values = np.transpose(np.asarray(query_est.values), axes=perm)
    true_values = np.asarray(query_true.values)

    if est_values.shape != true_values.shape:
        raise ValueError(
            f"Aligned factor shapes differ: est={est_values.shape}, true={true_values.shape}"
        )

    return float(np.abs(est_values - true_values).sum())


def compute_js_divergence(query_est, query_true):
    # Keep consistent alignment with TVD to avoid factor-order mismatch.
    true_vars = list(query_true.variables)
    est_vars = list(query_est.variables)

    if set(true_vars) != set(est_vars):
        raise ValueError("query_est and query_true variables are not identical")

    perm = [est_vars.index(var) for var in true_vars]
    est_values = np.transpose(np.asarray(query_est.values), axes=perm)
    true_values = np.asarray(query_true.values)

    if est_values.shape != true_values.shape:
        raise ValueError(
            f"Aligned factor shapes differ: est={est_values.shape}, true={true_values.shape}"
        )

    p = est_values.astype(np.float64, copy=False)
    q = true_values.astype(np.float64, copy=False)
    m = 0.5 * (p + q)

    with np.errstate(divide="ignore", invalid="ignore"):
        kl_pm = np.where(p > 0, p * np.log(p / m), 0.0)
        kl_qm = np.where(q > 0, q * np.log(q / m), 0.0)

    return float(0.5 * np.sum(kl_pm) + 0.5 * np.sum(kl_qm))


def run_experiment(r_repeat, sample_repeat, Missing_proportions, r_size, file_names, sample_size=100, seed=42):
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

        for miss_proportion in Missing_proportions:
            r_case_rows = []

            for _ in range(r_repeat):
                r_names = random.sample(node_names, r_size)
                inference_true = VariableElimination(true_model)
                true_query = inference_true.query(variables=r_names, evidence={}, show_progress=False)

                r_idx = [node_to_idx[node] for node in r_names]


                ang_nodes_idx = decom_h.dag_get_ancestors(ig_graph, r_idx)
                ang_nodes = map_nodes_to_names(ang_nodes_idx, idx_to_node)
                ang_model = BayesianNetwork()
                ang_model.add_nodes_from(ang_nodes)
                ang_model.add_edges_from(induced_edges_as_names(true_model, set(ang_nodes)))
                ang_state_names = extract_state_names(true_model, ang_nodes)



                collapsed_nodes_idx = decom_h.cmdsa(ig_graph, r_idx)
                collapsed_nodes = map_nodes_to_names(collapsed_nodes_idx, idx_to_node)
                collapsed_model = BayesianNetwork()
                collapsed_model.add_nodes_from(collapsed_nodes)
                collapsed_model.add_edges_from(induced_edges_as_names(true_model, set(collapsed_nodes)))
                collapsed_state_names = extract_state_names(true_model, collapsed_nodes)

                covariate_cols = list(node_names)
                ang_js_values = []
                collapsed_js_values = []
                ang_complete_rows = []
                collapsed_complete_rows = []
                ang_errors = []
                collapsed_errors = []

                for _ in range(sample_repeat):
                    df = BayesianModelSampling(true_model).forward_sample(size=sample_size, show_progress=False)
                    df_incomplete = inject_missing_values(
                        df=df,
                        miss_proportion=miss_proportion,
                        covariate_cols=covariate_cols,
                        entry_missing_rate=0.01,  # 1% cell-level missingness within contaminated rows
                    )

                    # MLE cannot use missing values; keep complete-case rows for each method's variable subset.
                    ang_df_full = df_incomplete[ang_nodes].dropna(axis=0, how="any")
                    collapsed_df_full = df_incomplete[collapsed_nodes].dropna(axis=0, how="any")

                    ang_complete_rows.append(len(ang_df_full))
                    collapsed_complete_rows.append(len(collapsed_df_full))

                    if len(ang_df_full) > 0:
                        try:
                            ang_model.fit(
                                ang_df_full,
                                estimator=MaximumLikelihoodEstimator,
                                state_names=ang_state_names,
                            )
                            inference_ang = VariableElimination(ang_model)
                            ang_query = inference_ang.query(variables=r_names, evidence={}, show_progress=False)
                            ang_js_values.append(compute_js_divergence(ang_query, true_query))
                        except Exception as exc:
                            ang_errors.append(repr(exc))

                    if len(collapsed_df_full) > 0:
                        try:
                            collapsed_model.fit(
                                collapsed_df_full,
                                estimator=MaximumLikelihoodEstimator,
                                state_names=collapsed_state_names,
                            )
                            inference_collapsed = VariableElimination(collapsed_model)
                            collapsed_query = inference_collapsed.query(variables=r_names, evidence={}, show_progress=False)
                            collapsed_js_values.append(compute_js_divergence(collapsed_query, true_query))
                        except Exception as exc:
                            collapsed_errors.append(repr(exc))

                r_case_rows.append({
                    "ang_js": float(np.mean(ang_js_values)) if ang_js_values else np.nan,
                    "collapsed_js": float(np.mean(collapsed_js_values)) if collapsed_js_values else np.nan,
                    "ang_complete_rows": float(np.mean(ang_complete_rows)) if ang_complete_rows else 0.0,
                    "collapsed_complete_rows": float(np.mean(collapsed_complete_rows)) if collapsed_complete_rows else 0.0,
                    "ang_nodes": len(ang_nodes),
                    "collapsed_nodes": len(collapsed_nodes),
                    "ang_error": ang_errors[-1] if ang_errors else "",
                    "collapsed_error": collapsed_errors[-1] if collapsed_errors else "",
                })

            valid_ang = [row["ang_js"] for row in r_case_rows if not np.isnan(row["ang_js"])]
            valid_collapsed = [row["collapsed_js"] for row in r_case_rows if not np.isnan(row["collapsed_js"])]
            avg_ang_js = float(np.mean(valid_ang)) if valid_ang else np.nan
            avg_collapsed_js = float(np.mean(valid_collapsed)) if valid_collapsed else np.nan
            avg_ang_complete_rows = float(np.mean([row["ang_complete_rows"] for row in r_case_rows]))
            avg_collapsed_complete_rows = float(np.mean([row["collapsed_complete_rows"] for row in r_case_rows]))
            avg_ang_nodes = float(np.mean([row["ang_nodes"] for row in r_case_rows]))
            avg_collapsed_nodes = float(np.mean([row["collapsed_nodes"] for row in r_case_rows]))
            ang_error = next((row["ang_error"] for row in r_case_rows if row["ang_error"]), "")
            collapsed_error = next((row["collapsed_error"] for row in r_case_rows if row["collapsed_error"]), "")
            
            """
            print(
                f"file={file}, miss={miss_proportion:.2f}, "
                f"ang_js={avg_ang_js:.4f}, coll_js={avg_collapsed_js:.4f}, "
                f"ang_sample={avg_ang_complete_rows:.1f}, coll_sample={avg_collapsed_complete_rows:.1f}, "
                f"ang_dim=({avg_ang_nodes:.1f}), coll_dim=({avg_collapsed_nodes:.1f})"
            )
            """

            if not np.isnan(avg_ang_js) and not np.isnan(avg_collapsed_js):
                pass
            else:
                print(f"  ang_error={ang_error}")
                print(f"  collapsed_error={collapsed_error}")

            result_rows.append({
                "file": file,
                "sample_size": sample_size,
                "miss_proportion": miss_proportion,
                "r_repeat": r_repeat,
                "sample_repeat": sample_repeat,
                "r_size": r_size,
                "ang_js": avg_ang_js,
                "collapsed_js": avg_collapsed_js,
                "ang_complete_rows": avg_ang_complete_rows,
                "collapsed_complete_rows": avg_collapsed_complete_rows,
                "ang_nodes": avg_ang_nodes,
                "collapsed_nodes": avg_collapsed_nodes,
                "ang_error": ang_error,
                "collapsed_error": collapsed_error,
            })
            
        print(f"completed {file} experiment")
    return pd.DataFrame(result_rows)


def main():
    # 用于快速切换的测试文件集合；需要做正式实验时可以换成更完整的网络列表。
    file_names = ["hailfinder", "hepar2",  "andes", "munin"]
    # 每个样本量下重复采样的规模，代表数据集大小。
    Missing_proportions = [0.2,0.4,0.6,0.8]
    # 对不同 R 组进行重复实验的次数，最后会对这部分取平均。
    r_repeat = 10
    # 同一个 R 下，重复抽样和学习/推理的次数。
    sample_repeat = 50
    # 每次随机抽取的 R 的节点数。
    r_size = 5
    # 全局随机种子，保证别人运行时可以复现相同的随机序列。
    seed = 20260410

    result = run_experiment(
        r_repeat=r_repeat,
        sample_repeat=sample_repeat,
        r_size=r_size,
        file_names=file_names,
        Missing_proportions=Missing_proportions,
        seed=seed,
    )

    result_dir = base_dir / "result"
    result_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(result_dir / "miss_data_experiment.csv", index=False)
    print(result.head())


if __name__ == "__main__":
    main()