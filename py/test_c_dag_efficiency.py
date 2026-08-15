import os
import random
import sys
from itertools import combinations
from pathlib import Path
from time import perf_counter

import igraph as ig

from Graph_gererators import generate_connected_dag


base_dir = Path(__file__).resolve().parent
project_root = base_dir.parent
build_dir = project_root / "build"
sys.path.insert(0, str(build_dir))
os.environ["PATH"] = (
    str(build_dir)
    + os.pathsep
    + r"D:\\vcpkg-master\\installed\\x64-windows\bin"
    + os.pathsep
    + os.environ.get("PATH", "")
)

if hasattr(os, "add_dll_directory"):
    _dll_dirs = [
        os.add_dll_directory(str(build_dir)),
        os.add_dll_directory(r"D:\\vcpkg-master\\installed\\x64-windows\\bin"),
    ]

import decom_h


def nx_dag_to_igraph(g_nx):
    node_names = list(g_nx.nodes())
    name_to_index = {name: idx for idx, name in enumerate(node_names)}
    edges = [(name_to_index[u], name_to_index[v]) for (u, v) in g_nx.edges()]
    g_ig = ig.Graph(n=len(node_names), edges=edges, directed=True)
    return g_ig


def ancestor_moral_subgraph_with_mapping(g_dag, r_nodes):
    """Build moralized undirected ancestor subgraph of R and keep index mapping."""
    ancestors = set()
    for r in r_nodes:
        ancestors.update(g_dag.subcomponent(r, mode="IN"))

    local_to_global = sorted(ancestors)
    global_to_local = {v: i for i, v in enumerate(local_to_global)}

    # Directed ancestor-induced subgraph.
    g_an_r = g_dag.induced_subgraph(local_to_global)

    # Moralization: keep skeleton edges and connect co-parents.
    undirected_edges = set()
    for u, v in g_an_r.get_edgelist():
        a, b = (u, v) if u < v else (v, u)
        undirected_edges.add((a, b))

    for child in range(g_an_r.vcount()):
        parents = g_an_r.predecessors(child)
        if len(parents) < 2:
            continue
        for u, v in combinations(parents, 2):
            a, b = (u, v) if u < v else (v, u)
            undirected_edges.add((a, b))

    g_an_r_ug = ig.Graph(n=g_an_r.vcount(), edges=list(undirected_edges), directed=False)
    r_local = [global_to_local[r] for r in r_nodes]
    return g_an_r_ug, r_local, local_to_global


def run_one_case(n, p, r_size, seed):
    random.seed(seed)

    g_ig = generate_connected_dag(n, p)

    all_nodes = list(range(g_ig.vcount()))
    R = random.sample(all_nodes, r_size)

    # Python 版预处理
    t1 = perf_counter()
    g_an_r_ug, r_local, local_to_global = ancestor_moral_subgraph_with_mapping(g_ig, R)
    #print(f"python 预处理: |r_local|={r_local}, |local_to_global|={local_to_global}")
    h1_local = decom_h.get_minimal_collapsible(g_an_r_ug, r_local)
    ##print(f"python 预处理: h1_local={h1_local}")
    H_py = sorted([local_to_global[idx] for idx in h1_local])
    H_py = sorted(decom_h.dag_get_minimal_collapsible(g_ig, H_py))
    t2 = perf_counter()

    # C 版DAG版本（内部预处理，现已优化为批量边添加）
    t3 = perf_counter()
    H_c = sorted(decom_h.dag_get_minimal_collapsible(g_ig, R))
    t4 = perf_counter()

    same_h = H_py == H_c

    return {
            "m": g_ig.ecount(),
            "py_time": t2 - t1,
            "c_time": t4 - t3,
            "h_size": len(H_c),
            "same_h": same_h,
            "H_py": H_py,
            "H_c": H_c,
            "R": R,
            "g_ig": g_ig,
        }


if __name__ == "__main__":
    n_values = [1000, 2000]
    p_values = [0.05, 0.001]
    repeats = 20
    r_size = 10


    print("=== C实现 DAG Convex Hull 效率测试（Python vs C预处理对比）===")
    print(
        f"固定参数: n_values={n_values}, "
        f"p_values={p_values}, 每个配置重复{repeats}次, |R|={r_size}"
    )

    found_mismatch = False
    for n in n_values:
        if found_mismatch:
            break
        for p_idx, p in enumerate(p_values):
            if found_mismatch:
                break
            case_results = []
            for rep in range(repeats):
                seed = 20260318 + n * 10000 + p_idx * 100 + rep
                result = run_one_case(n=n, p=p, r_size=r_size, seed=seed)
                case_results.append(result)
                
                # 检查不一致
                if not result["same_h"]:
                    found_mismatch = True
                    print("\n" + "="*80)
                    print("❌ 发现不一致！")
                    print("="*80)
                    print(f"参数: n={n}, p={p}, seed={seed}, |R|={r_size}")
                    print(f"\nR (选中的{r_size}个节点): {result['R']}")
                    print(f"\nH_py (Python预处理结果): {result['H_py']}")
                    print(f"H_c  (C预处理结果): {result['H_c']}")
                    print(f"\n不一致节点:")
                    h_py_set = set(result['H_py'])
                    h_c_set = set(result['H_c'])
                    only_in_py = h_py_set - h_c_set
                    only_in_c = h_c_set - h_py_set
                    if only_in_py:
                        print(f"  仅在Python中: {sorted(only_in_py)}")
                    if only_in_c:
                        print(f"  仅在C中: {sorted(only_in_c)}")
                    
                    # 输出有向图的所有边
                    print(f"\n有向图 (共{n}个节点，{result['m']}条边):")
                    edges = result['g_ig'].get_edgelist()
                    print(f"  边列表: {edges}")
                    
                    print("\n终止后续测试循环。")
                    print("="*80 + "\n")
                    break

            if not found_mismatch:
                avg_m = sum(x["m"] for x in case_results) / repeats
                avg_py = sum(x["py_time"] for x in case_results) / repeats
                avg_c = sum(x["c_time"] for x in case_results) / repeats
                avg_h = sum(x["h_size"] for x in case_results) / repeats
                same_count = sum(1 for x in case_results if x["same_h"])

                print(
                    f"n={n}, p={p}, avg_m={avg_m:.1f}, "
                    f"py_time={avg_py:.6f}s, "
                    f"c_time={avg_c:.6f}s, "
                    f"avg_|H|={avg_h:.1f}, "
                    f"H_py==H_c: {same_count}/{repeats}"
                )
