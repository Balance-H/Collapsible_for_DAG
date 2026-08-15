//#define _CRTDBG_MAP_ALLOC //用于内存泄露检测
//#include <crtdbg.h>


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <igraph.h>
#include "decom_h.h"  // 你的函数声明头文件

namespace py = pybind11;

// 包装 components_forbidden
// C 函数声明（假设在 decom_h.h 中）
extern "C" igraph_error_t components_forbidden(
    const igraph_t *graph,
    igraph_vector_ptr_t *components,
    igraph_vector_ptr_t *boundaries,
    const igraph_vector_int_t *forbidden_vertices);

// 用于释放 igraph_vector_ptr_t 中每个元素的辅助函数
void free_vector_ptr(igraph_vector_ptr_t *vec_ptr) {
    for (long i = 0; i < igraph_vector_ptr_size(vec_ptr); ++i) {
        igraph_vector_int_t *v = static_cast<igraph_vector_int_t*>(VECTOR(*vec_ptr)[i]);
        igraph_vector_int_destroy(v);
        free(v);
    }
    igraph_vector_ptr_destroy(vec_ptr);
}

// C++ 封装接口，返回 Python 的 list[(component, boundary), ...]
py::list components_forbidden_wrapper(py::object graph_obj,
                             const std::vector<int> &forbidden_vertices = {}) {
    // 从 Python igraph 对象获取 igraph_t 指针
    py::object graph_capsule = graph_obj.attr("__graph_as_capsule")();
    igraph_t *graph = static_cast<igraph_t*>(PyCapsule_GetPointer(graph_capsule.ptr(), nullptr));
    if (!graph) throw std::runtime_error("Invalid igraph capsule");

    // 初始化禁忌节点向量
    igraph_vector_int_t forbidden_vec;
    igraph_vector_int_init(&forbidden_vec, 0);
    for (int v : forbidden_vertices) {
        igraph_vector_int_push_back(&forbidden_vec, v);
    }

    // 初始化存储组件和边界的向量
    igraph_vector_ptr_t components;
    igraph_vector_ptr_init(&components, 0);
    igraph_vector_ptr_t boundaries;
    igraph_vector_ptr_init(&boundaries, 0);

    // 调用 C 函数
    igraph_error_t err = components_forbidden(graph, &components, &boundaries, &forbidden_vec);
    igraph_vector_int_destroy(&forbidden_vec);

    if (err != IGRAPH_SUCCESS) {
        free_vector_ptr(&components);
        free_vector_ptr(&boundaries);
        throw std::runtime_error("components_forbidden failed");
    }

    // 转换为 Python list[(list[int], list[int]), ...]
    py::list result;
    igraph_integer_t n = igraph_vector_ptr_size(&components);
    for (igraph_integer_t i = 0; i < n; i++) {
        igraph_vector_int_t *comp = static_cast<igraph_vector_int_t*>(VECTOR(components)[i]);
        igraph_vector_int_t *bound = static_cast<igraph_vector_int_t*>(VECTOR(boundaries)[i]);

        py::list py_comp;
        for (long j = 0; j < igraph_vector_int_size(comp); ++j) {
            py_comp.append(VECTOR(*comp)[j]);
        }

        py::list py_bound;
        for (long j = 0; j < igraph_vector_int_size(bound); ++j) {
            py_bound.append(VECTOR(*bound)[j]);
        }

        result.append(py::make_tuple(py_comp, py_bound));
    }

    // 释放内存
    free_vector_ptr(&components);
    free_vector_ptr(&boundaries);

    return result;
}


py::list close_separator_wrapper(py::object graph_obj,
                                        int vertex,
                                        const std::vector<int> &forbidden_vertices = {}) {
    // 获取封装的图对象
    py::object graph_capsule = graph_obj.attr("__graph_as_capsule")();
    igraph_t *graph = static_cast<igraph_t*>(PyCapsule_GetPointer(graph_capsule.ptr(), nullptr));
    if (!graph) throw std::runtime_error("Invalid igraph capsule");

    // 初始化禁忌节点列表
    igraph_vector_int_t forbidden_vec;
    igraph_vector_int_init(&forbidden_vec, 0);
    for (auto v : forbidden_vertices) {
        igraph_vector_int_push_back(&forbidden_vec, v);
    }

    // 初始化 bound_b
    igraph_vector_int_t bound_b;
    igraph_vector_int_init(&bound_b, 0);

    igraph_error_t err = close_separator(graph, vertex,
                                                               &forbidden_vec, &bound_b);
    igraph_vector_int_destroy(&forbidden_vec);

    if (err != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&bound_b);
        throw std::runtime_error("close_separator_b failed");
    }

    // 转换 bound_b 为 Python 列表
    py::list py_bound_b;
    for (long i = 0; i < igraph_vector_int_size(&bound_b); ++i) {
        py_bound_b.append(VECTOR(bound_b)[i]);
    }

    igraph_vector_int_destroy(&bound_b);
    return py_bound_b;
}



extern "C" igraph_error_t get_minimal_collapsible(
    const igraph_t *graph,
    const igraph_vector_int_t *nodes,
    igraph_vector_int_t *H_out
);

extern "C" igraph_error_t dag_get_minimal_collapsible(
    const igraph_t *graph,
    const igraph_vector_int_t *r_nodes,
    igraph_vector_int_t *H_out
);

extern "C" igraph_error_t dag_get_ancestors(
    const igraph_t *graph,
    const igraph_vector_int_t *r_nodes,
    igraph_vector_int_t *ancestors_out
);

extern "C" igraph_error_t mcs_with_cliques(
    const igraph_t *graph,
    igraph_vector_int_t *alpha,
    igraph_vector_int_t *alpham1,
    igraph_vector_ptr_t *cliques
);

py::list get_minimal_collapsible_wrapper(py::object graph_obj, 
                                   const std::vector<int> &r_nodes) {
    // 获取 igraph_t 指针
    py::object graph_capsule = graph_obj.attr("__graph_as_capsule")();
    igraph_t *graph = static_cast<igraph_t*>(PyCapsule_GetPointer(graph_capsule.ptr(), nullptr));
    if (!graph) throw std::runtime_error("Invalid igraph capsule");

    // 将 r_nodes 转成 igraph_vector_int_t
    igraph_vector_int_t r_vec;
    igraph_vector_int_init(&r_vec, 0);
    for (int v : r_nodes) {
        igraph_vector_int_push_back(&r_vec, v);
    }

    // 调用 get_minimal_collapsible
    igraph_vector_int_t H_vec;
    igraph_vector_int_init(&H_vec, 0);

    igraph_error_t err = get_minimal_collapsible(graph, &r_vec, &H_vec);
    igraph_vector_int_destroy(&r_vec);

    if (err != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&H_vec);
        throw std::runtime_error("get_minimal_collapsible failed");
    }

    // 转换结果到 Python list
    py::list result;
    for (long i = 0; i < igraph_vector_int_size(&H_vec); ++i) {
        result.append(VECTOR(H_vec)[i]);
    }
    igraph_vector_int_destroy(&H_vec);

    return result;
}

py::list dag_get_minimal_collapsible_wrapper(py::object graph_obj,
                       const std::vector<int> &r_nodes) {
    py::object graph_capsule = graph_obj.attr("__graph_as_capsule")();
    igraph_t *graph = static_cast<igraph_t*>(PyCapsule_GetPointer(graph_capsule.ptr(), nullptr));
    if (!graph) throw std::runtime_error("Invalid igraph capsule");

    igraph_vector_int_t r_vec;
    igraph_vector_int_init(&r_vec, 0);
    for (int v : r_nodes) {
        igraph_vector_int_push_back(&r_vec, v);
    }

    igraph_vector_int_t H_vec;
    igraph_vector_int_init(&H_vec, 0);

    igraph_error_t err = dag_get_minimal_collapsible(graph, &r_vec, &H_vec);
    igraph_vector_int_destroy(&r_vec);

    if (err != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&H_vec);
        throw std::runtime_error("dag_get_minimal_collapsible failed, igraph error=" + std::to_string((int)err));
    }

    py::list result;
    for (long i = 0; i < igraph_vector_int_size(&H_vec); ++i) {
        result.append(VECTOR(H_vec)[i]);
    }
    igraph_vector_int_destroy(&H_vec);
    return result;
}

py::list dag_get_ancestors_wrapper(py::object graph_obj,
                                   const std::vector<int> &r_nodes) {
    py::object graph_capsule = graph_obj.attr("__graph_as_capsule")();
    igraph_t *graph = static_cast<igraph_t*>(PyCapsule_GetPointer(graph_capsule.ptr(), nullptr));
    if (!graph) throw std::runtime_error("Invalid igraph capsule");

    igraph_vector_int_t r_vec;
    igraph_vector_int_init(&r_vec, 0);
    for (int v : r_nodes) {
        igraph_vector_int_push_back(&r_vec, v);
    }

    igraph_vector_int_t anc_vec;
    igraph_vector_int_init(&anc_vec, 0);

    igraph_error_t err = dag_get_ancestors(graph, &r_vec, &anc_vec);
    igraph_vector_int_destroy(&r_vec);

    if (err != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&anc_vec);
        throw std::runtime_error("dag_get_ancestors failed, igraph error=" + std::to_string((int)err));
    }

    py::list result;
    for (long i = 0; i < igraph_vector_int_size(&anc_vec); ++i) {
        result.append(VECTOR(anc_vec)[i]);
    }
    igraph_vector_int_destroy(&anc_vec);
    return result;
}

py::list mcs_with_cliques_wrapper(py::object graph_obj,
                                         bool return_cliques = false) {
    py::object graph_capsule = graph_obj.attr("__graph_as_capsule")();
    igraph_t *graph = static_cast<igraph_t*>(PyCapsule_GetPointer(graph_capsule.ptr(), nullptr));
    if (!graph) throw std::runtime_error("Invalid igraph capsule");

    igraph_vector_int_t alpha, order;
    igraph_vector_int_init(&alpha, 0);
    igraph_vector_int_init(&order, 0);

    igraph_vector_ptr_t cliques;
    igraph_vector_ptr_init(&cliques, 0);

    igraph_error_t err = mcs_with_cliques(graph, &alpha, &order, 
                                                  return_cliques ? &cliques : NULL);

    if (err != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&alpha);
        igraph_vector_int_destroy(&order);
        if (return_cliques) {
            for (long i = 0; i < igraph_vector_ptr_size(&cliques); ++i) {
                igraph_vector_int_destroy((igraph_vector_int_t *)VECTOR(cliques)[i]);
                free(VECTOR(cliques)[i]);
            }
            igraph_vector_ptr_destroy(&cliques);
        }
        throw std::runtime_error("mcs_with_cliques failed");
    }

    // Convert alpha to Python list
    py::list py_alpha;
    for (long i = 0; i < igraph_vector_int_size(&alpha); ++i) {
        py_alpha.append(VECTOR(alpha)[i]);
    }

    // Convert order to Python list
    py::list py_order;
    for (long i = 0; i < igraph_vector_int_size(&order); ++i) {
        py_order.append(VECTOR(order)[i]);
    }

    // Convert cliques to nested Python list
    py::list py_cliques;
    if (return_cliques) {
        for (long i = 0; i < igraph_vector_ptr_size(&cliques); ++i) {
            igraph_vector_int_t *clique = (igraph_vector_int_t *)VECTOR(cliques)[i];
            py::list py_clique;
            for (long j = 0; j < igraph_vector_int_size(clique); ++j) {
                py_clique.append(VECTOR(*clique)[j]);
            }
            py_cliques.append(py_clique);
            igraph_vector_int_destroy(clique);
            free(clique);
        }
    }
    igraph_vector_ptr_destroy(&cliques);

    igraph_vector_int_destroy(&alpha);
    igraph_vector_int_destroy(&order);

    // Return tuple: (alpha, order, cliques)
    return py::make_tuple(py_alpha, py_order, py_cliques).cast<py::list>();
}

// 辅助函数，将 igraph_vector_ptr_t 转换为 Python list（嵌套 list）
py::list igraph_vector_ptr_to_pylist(igraph_vector_ptr_t *vec_ptr) {
    py::list result;
    igraph_integer_t n = igraph_vector_ptr_size(vec_ptr);
    for (igraph_integer_t i = 0; i < n; i++) {
        igraph_vector_int_t *v = (igraph_vector_int_t *)VECTOR(*vec_ptr)[i];
        py::list inner_list;
        for (int j = 0; j < igraph_vector_int_size(v); j++) {
            inner_list.append(VECTOR(*v)[j]);
        }
        result.append(inner_list);
    }
    return result;
}


// 递归分解包装，直接接收 Python igraph.Graph 对象，内部调用 __graph_as_capsule 获取指针
py::list decompose_atoms_wrapper(py::object graph_obj) {
    // 通过 __graph_as_capsule 获取底层指针
    py::object graph_capsule = graph_obj.attr("__graph_as_capsule")();
    igraph_t *g = static_cast<igraph_t *>(PyCapsule_GetPointer(graph_capsule.ptr(), nullptr));
    if (!g) throw std::runtime_error("Invalid igraph capsule");

    // 1. 初始化两个输出容器和 tree
    igraph_vector_ptr_t atoms;
    igraph_vector_ptr_t separators;
    igraph_t tree_out;
    igraph_vector_ptr_init(&atoms, 0);
    igraph_vector_ptr_init(&separators, 0);
    igraph_empty(&tree_out, 0, 1);

    igraph_error_t err = IGRAPH_SUCCESS;
    
    // 2. 调用 decompose_atoms 函数
    try {
        err = decompose_atoms(g, &atoms, &separators, &tree_out);
    } catch (...) {
        err = IGRAPH_EINTERNAL;
    }

    if (err != IGRAPH_SUCCESS) {
        // 3. 错误发生时，清理容器
        for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&atoms); i++) {
            igraph_vector_int_destroy((igraph_vector_int_t *)VECTOR(atoms)[i]);
            free(VECTOR(atoms)[i]);
        }
        igraph_vector_ptr_destroy(&atoms);
        
        for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&separators); i++) {
            igraph_vector_int_destroy((igraph_vector_int_t *)VECTOR(separators)[i]);
            free(VECTOR(separators)[i]);
        }
        igraph_vector_ptr_destroy(&separators);
        igraph_destroy(&tree_out);

        throw std::runtime_error("decompose_atoms failed with code " + std::to_string(err));
    }

    // 4. 将 atoms 和 separators 转换为 Python 列表
    py::list atoms_pylist = igraph_vector_ptr_to_pylist(&atoms);
    py::list separators_pylist = igraph_vector_ptr_to_pylist(&separators);

    // 5. 将 tree_out 转换为邻接表 (edge list)
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    igraph_get_edgelist(&tree_out, &edges, 0);
    py::list tree_edges;
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&edges); i += 2) {
        py::list edge;
        edge.append(VECTOR(edges)[i]);
        edge.append(VECTOR(edges)[i+1]);
        tree_edges.append(edge);
    }
    igraph_vector_int_destroy(&edges);
    
    // 6. 清理内存
    for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&atoms); i++) {
        igraph_vector_int_destroy((igraph_vector_int_t *)VECTOR(atoms)[i]);
        free(VECTOR(atoms)[i]);
    }
    igraph_vector_ptr_destroy(&atoms);
    
    for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&separators); i++) {
        igraph_vector_int_destroy((igraph_vector_int_t *)VECTOR(separators)[i]);
        free(VECTOR(separators)[i]);
    }
    igraph_vector_ptr_destroy(&separators);
    igraph_destroy(&tree_out);

    // 7. 返回 [atoms_list, separators_list, tree_edges]
    py::list result;
    result.append(atoms_pylist);
    result.append(separators_pylist);
    result.append(tree_edges);
    
    return result;
}




py::list SAHR_wrapper(py::object graph_obj, const std::vector<int> &r_nodes) {
    // 获取 igraph_t 指针
    py::object graph_capsule = graph_obj.attr("__graph_as_capsule")();
    igraph_t *graph = static_cast<igraph_t*>(PyCapsule_GetPointer(graph_capsule.ptr(), nullptr));
    if (!graph) throw std::runtime_error("Invalid igraph capsule");

    // 转换 r_nodes 为 igraph_integer_t 数组
    std::vector<igraph_integer_t> r_vec(r_nodes.begin(), r_nodes.end());
    igraph_integer_t *r_ptr = r_vec.data();
    igraph_integer_t r_size = static_cast<igraph_integer_t>(r_vec.size());

    // 输出参数
    igraph_integer_t *local2global = nullptr;
    igraph_integer_t result_size = 0;

    igraph_error_t err = SAHR(graph, r_ptr, r_size, &local2global, &result_size);
    if (err != IGRAPH_SUCCESS) {
        if (local2global) free(local2global);
        throw std::runtime_error("SAHR failed with code " + std::to_string(err));
    }

    // 转换结果为 Python list
    py::list result;
    for (igraph_integer_t i = 0; i < result_size; i++) {
        result.append(local2global[i]);
    }

    free(local2global);
    return result;
}



PYBIND11_MODULE(decom_h, m) {


    m.def("SAHR", &SAHR_wrapper,
    py::arg("graph"),
    py::arg("r_nodes"),
    "Run SAHR algorithm and return remaining nodes' global indices as a list.");

    m.doc() = "Example igraph C extension";

    m.def("decompose_atoms", &decompose_atoms_wrapper,
      py::arg("graph"),
      "Perform recursive atom decomposition of the graph using Atom Expansion.\n"
      "Returns: A list containing [atoms, separators, tree].\n"
      "  - atoms: List of lists, where each inner list is an atom (a vertex set).\n"
      "  - separators: List of lists, where each inner list is a clique minimal separator (a vertex set).\n"
      "  - tree: List of edges representing the atom tree.");

    m.def("get_minimal_collapsible", &get_minimal_collapsible_wrapper,
          py::arg("graph"),
          py::arg("r_nodes"),
          "Compute minimal collapsible set from given nodes.");
    

    m.def("mcs_with_cliques", &mcs_with_cliques_wrapper,
          py::arg("graph"),
          py::arg("return_cliques") = false,
          "Run MCS algorithm and optionally return the maximal cliques found.");


    m.def("dag_get_minimal_collapsible", &dag_get_minimal_collapsible_wrapper,
            py::arg("graph"),
            py::arg("r_nodes"),
            "Compute minimal collapsible set from given nodes in a directed graph.");

    m.def("dag_get_ancestors", &dag_get_ancestors_wrapper,
            py::arg("graph"),
            py::arg("r_nodes"),
            "Get all ancestors of the given root nodes in a directed graph.");

    m.def("components_forbidden", &components_forbidden_wrapper,
        py::arg("graph"),
        py::arg("forbidden_vertices") = std::vector<int>{},
        "Calculate connected components excluding forbidden vertices, "
        "and return each component with its boundary forbidden nodes.");

    m.def("close_separator", &close_separator_wrapper,
        py::arg("graph"),
        py::arg("vertex"),
        py::arg("forbidden_vertices") = std::vector<int>{},
        "Calculate forbidden boundary reachable from vertex without traversing forbidden vertices.");

    m.def("decompose_atoms", &decompose_atoms_wrapper,
        py::arg("graph"),
        "Perform atom decomposition on undirected graph using MCS ordering. Returns [atoms, separators].");

}

