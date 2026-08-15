#ifndef DECOM_H_H
#define DECOM_H_H


#include <igraph.h>

#ifdef __cplusplus
extern "C" {
#endif

igraph_error_t components_forbidden(
    const igraph_t *graph,
    igraph_vector_ptr_t *components,
    igraph_vector_ptr_t *boundaries,
    const igraph_vector_int_t *forbidden_vertices
);

// 修改后的 igraph_subcomponent_forbidden 函数声明
igraph_error_t close_separator(
    const igraph_t *graph,
    igraph_integer_t vertex,
    const igraph_vector_int_t *forbidden_vertices,
    igraph_vector_int_t *bound_b  // 用于存储遇到的禁忌节点
);

typedef struct {
    igraph_vector_int_t S;
    igraph_integer_t u_p_index;
} pair_t;

typedef struct {
    igraph_vector_int_t *V_nodes;
    igraph_vector_ptr_t *W;
} task_t;

igraph_error_t decompose_atoms(
    const igraph_t *graph,
    igraph_vector_ptr_t *atoms_out,
    igraph_vector_ptr_t *separators_out,
    igraph_t *tree_out
);

igraph_error_t get_minimal_collapsible(
    const igraph_t *graph,
    const igraph_vector_int_t *nodes,
    igraph_vector_int_t *H_out
);

igraph_error_t dag_get_minimal_collapsible(
    const igraph_t *graph,
    const igraph_vector_int_t *r_nodes,
    igraph_vector_int_t *H_out
);

igraph_error_t dag_get_ancestors(
    const igraph_t *graph,
    const igraph_vector_int_t *r_nodes,
    igraph_vector_int_t *ancestors_out
);

igraph_error_t mcs_with_cliques(
    const igraph_t *graph,
    igraph_vector_int_t *alpha,
    igraph_vector_int_t *alpham1,
    igraph_vector_ptr_t *cliques
);

igraph_error_t SAHR(
    igraph_t *g, 
    igraph_integer_t *r, 
    igraph_integer_t r_size, 
    igraph_integer_t **local2global, 
    igraph_integer_t *result_size
);




#ifdef __cplusplus
}
#endif

#endif // decom_h_H
