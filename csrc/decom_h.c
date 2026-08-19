//#define _CRTDBG_MAP_ALLOC //用于内存泄露检测
//#include <crtdbg.h>

#include <igraph.h>
#include <stdlib.h>
#include <stdint.h>  // for SIZE_MAX
#include <string.h>  // for memset
#include <stdbool.h>
#include <omp.h>
#include "uthash.h"
#include <stdio.h>
#include <time.h>
#include "decom_h.h"

static void debug_print_vector_int(const char *name, const igraph_vector_int_t *vec) {
    igraph_integer_t size = igraph_vector_int_size(vec);
    fprintf(stderr, "%s(size=%ld)=[", name, (long) size);
    for (igraph_integer_t i = 0; i < size; ++i) {
        fprintf(stderr, "%ld", (long) VECTOR(*vec)[i]);
        if (i + 1 < size) {
            fprintf(stderr, ", ");
        }
    }
    fprintf(stderr, "]\n");
}

static void debug_print_vector_ptr_of_int(const char *name, const igraph_vector_ptr_t *vecptr) {
    igraph_integer_t size = igraph_vector_ptr_size(vecptr);
    fprintf(stderr, "%s(count=%ld)\n", name, (long) size);
    for (igraph_integer_t i = 0; i < size; ++i) {
        igraph_vector_int_t *vec = (igraph_vector_int_t *) VECTOR(*vecptr)[i];
        if (!vec) {
            fprintf(stderr, "  [%ld]=NULL\n", (long) i);
            continue;
        }
        fprintf(stderr, "  [%ld] ", (long) i);
        debug_print_vector_int("", vec);
    }
}

static void debug_print_graph_edges(const char *name, const igraph_t *graph) {
    igraph_integer_t vcount = igraph_vcount(graph);
    igraph_integer_t ecount = igraph_ecount(graph);
    fprintf(stderr, "%s(vcount=%ld, ecount=%ld)\n", name, (long) vcount, (long) ecount);
    for (igraph_integer_t eid = 0; eid < ecount; ++eid) {
        igraph_integer_t from = 0;
        igraph_integer_t to = 0;
        if (igraph_edge(graph, eid, &from, &to) == IGRAPH_SUCCESS) {
            fprintf(stderr, "  (%ld -> %ld)\n", (long) from, (long) to);
        } else {
            fprintf(stderr, "  (%ld -> ?)\n", (long) eid);
        }
    }
}

/**
 * @brief 计算图中带有禁忌节点约束的弱连通分量及其边界节点集合。
 * 
 * 本函数在给定图中查找所有的弱连通分量，遍历过程中跳过指定的禁忌节点。
 * 每个连通分量由不包含禁忌节点的顶点集合构成，同时返回该连通分量与禁忌节点相邻的边界节点集合。
 * 禁忌节点既不会作为分量成员，也不会参与遍历扩展，但会被记录为相应连通分量的边界节点。
 * 
 * @param graph 输入的图结构，类型为 igraph_t*，要求图已初始化且合法。
 * @param components 已初始化的 igraph_vector_ptr_t* 指针，用于输出各连通分量，每个元素为指向 igraph_vector_int_t 的指针，存储该分量的节点集合。
 * @param boundaries 已初始化的 igraph_vector_ptr_t* 指针，用于输出对应连通分量的边界节点集合，每个元素为指向 igraph_vector_int_t 的指针。
 * @param forbidden_vertices 禁忌节点集合，类型为 igraph_vector_int_t*，指定的节点不会被访问或包含于分量。
 * 
 * @return 返回状态码，IGRAPH_SUCCESS 表示成功，其他值表示对应错误（如内存不足、无效顶点等）。
 * 
 * @note
 * - 禁忌节点在遍历过程中被标记以避免重复访问。
 * - 函数内部通过广度优先搜索(BFS)实现连通分量检测。
 * - 调用者需确保传入的 components 和 boundaries 已正确初始化，且调用后负责释放其中分配的内存。
 */


igraph_error_t components_forbidden(
    const igraph_t *graph,
    igraph_vector_ptr_t *components,
    igraph_vector_ptr_t *boundaries,
    const igraph_vector_int_t *forbidden_vertices)
{
    igraph_integer_t n = igraph_vcount(graph);

    // 访问标记数组：0未访问，1已访问, 2是禁忌节点
    char *visited = (char*)calloc(n, sizeof(char));
    if (!visited) {
        return IGRAPH_ENOMEM;
    }

    // 标记禁忌节点为已访问，防止遍历
    if (forbidden_vertices) {
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(forbidden_vertices); i++) {
            igraph_integer_t v = VECTOR(*forbidden_vertices)[i];
            visited[v] = 2;
        }
    }

     // 边界节点辅助标记，避免重复加入
    char *boundary_marked = (char*)calloc(n, sizeof(char));
    if (!boundary_marked) {
        free(visited);
        return IGRAPH_ENOMEM;
    }

    // 用于 BFS 的队列（动态数组）
    igraph_integer_t *queue = (igraph_integer_t*)malloc(n * sizeof(igraph_integer_t));
     if (!queue) {
        free(visited);
        free(boundary_marked);
        return IGRAPH_ENOMEM;
    }

    igraph_vector_int_t neighbors;
    igraph_vector_int_t parents;
    igraph_vector_int_init(&neighbors, 0);
    igraph_vector_int_init(&parents, 0);

    for (igraph_integer_t v = 0; v < n; v++) {
        if (visited[v]) continue;  // 已访问或禁忌跳过

        // 新连通分量和边界初始化
        igraph_vector_int_t *comp = (igraph_vector_int_t*)malloc(sizeof(igraph_vector_int_t));
        igraph_vector_int_t *bound = (igraph_vector_int_t*)malloc(sizeof(igraph_vector_int_t));
        if (!comp || !bound) {
            free(visited);
            free(boundary_marked);
            free(queue);
            igraph_vector_int_destroy(&neighbors);
            if (comp) free(comp);
            if (bound) free(bound);
            return IGRAPH_ENOMEM;
        }
        
        igraph_vector_int_init(comp, 0);
        igraph_vector_int_init(bound, 0);

        // BFS 初始化
        igraph_integer_t head = 0, tail = 0;
        queue[tail++] = v;
        visited[v] = 1;

        while (head < tail) {
            igraph_integer_t cur = queue[head++];

            igraph_vector_int_push_back(comp, cur);

            IGRAPH_CHECK(igraph_neighbors(graph, &neighbors, cur, IGRAPH_ALL));
            igraph_integer_t nei_count = igraph_vector_int_size(&neighbors);

            for (igraph_integer_t i = 0; i < nei_count; i++) {
                igraph_integer_t w = VECTOR(neighbors)[i];
                if (visited[w] == 2) {
                    // 禁忌节点，加入边界（避免重复）
                    if (!boundary_marked[w]) {
                        igraph_vector_int_push_back(bound, w);
                        boundary_marked[w] = 1;
                    }
                    continue; // 不加入队列
                }

                if (!visited[w]) {
                    visited[w] = 1;
                    queue[tail++] = w;
                }
            }
        }

        igraph_vector_ptr_push_back(components, comp);
        igraph_vector_ptr_push_back(boundaries, bound);

        // 清空边界辅助标记，为下一个连通分量做准备
        igraph_integer_t bsize = igraph_vector_int_size(bound);
        for (igraph_integer_t i = 0; i < bsize; i++) {
            boundary_marked[VECTOR(*bound)[i]] = 0;
        }

        
    }

    free(visited);
    free(boundary_marked);
    free(queue);
    igraph_vector_int_destroy(&neighbors);

    return IGRAPH_SUCCESS;
}

// 工具函数：unique + sort, 去重
igraph_integer_t vector_int_unique(igraph_vector_int_t *v) {
    if (igraph_vector_int_size(v) == 0) return IGRAPH_SUCCESS;
    igraph_vector_int_sort(v);
    igraph_integer_t  write_pos = 1;
    for (igraph_integer_t  i = 1; i < igraph_vector_int_size(v); i++) {
        if (VECTOR(*v)[i] != VECTOR(*v)[i - 1]) {
            VECTOR(*v)[write_pos++] = VECTOR(*v)[i];
        }
    }
    igraph_vector_int_resize(v, write_pos);
    return IGRAPH_SUCCESS;
}

// 返回 1 如果 S ⊆ U (按元素比较)，否则 0
static igraph_bool_t vector_int_is_subset(const igraph_vector_int_t *S, const igraph_vector_int_t *U) {
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(S); i++) {
        igraph_integer_t sv = VECTOR(*S)[i];
        igraph_bool_t found = 0;
        for (igraph_integer_t j = 0; j < igraph_vector_int_size(U); j++) {
            if (VECTOR(*U)[j] == sv) { found = 1; break; }
        }
        if (!found) return 0;
    }
    return 1;
}

/**
 * @brief 计算从指定起点出发、忽略禁忌节点的单点连通分量及其与禁忌节点相邻的边界节点。
 * 
 * 函数以指定顶点 `vertex` 作为起点，执行广度优先搜索（BFS）遍历图中除禁忌节点外的可达节点，
 * 构成该起点所在的连通分量。同时收集所有与该连通分量邻接但被标记为禁忌的节点，作为边界节点返回。
 * 
 * @param graph 输入图，类型为 `const igraph_t*`，要求已初始化且合法。
 * @param vertex 起始节点编号，类型为 `igraph_integer_t`，必须在图节点范围内。
 * @param forbidden_vertices 禁忌节点集合，类型为 `const igraph_vector_int_t*`，遍历时忽略这些节点。
 * @param bound_b 输出向量，类型为 `igraph_vector_int_t*`，返回遇到的禁忌节点集合（即边界节点）。
 * 
 * @return 返回状态码，`IGRAPH_SUCCESS` 表示执行成功，其他值代表错误状态（如起点无效、内存不足等）。
 * 
 * @note
 * - 禁忌节点在遍历过程中被标记为不可访问，且不包含于连通分量。
 * - 边界节点即为直接邻接于连通分量但为禁忌节点的集合。
 * - 函数内部使用队列实现 BFS，效率较高，适合中大型图遍历。
 * - 调用者需确保输入参数有效，且负责管理 `bound_b` 的内存。
 */


igraph_error_t close_separator(
    const igraph_t *graph,
    igraph_integer_t vertex,  //节点b,算法会序
    const igraph_vector_int_t *forbidden_vertices,
    igraph_vector_int_t *bound_b  // 新增的参数，用于存储遇到的禁忌节点
) {
    igraph_integer_t n = igraph_vcount(graph);
    if (vertex < 0 || vertex >= n) {
        return IGRAPH_EINVVID;
    }

    // 标记数组：0未访问，1已访问或禁忌
    char *visited = calloc(n, sizeof(char));
    if (!visited) {
        return IGRAPH_ENOMEM;
    }

    // 标记禁忌节点为已访问，防止遍历
    if (forbidden_vertices) {
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(forbidden_vertices); i++) {
            igraph_integer_t v = VECTOR(*forbidden_vertices)[i];
            if (v >= 0 && v < n) {
                visited[v] = 2;
            }
        }
    }

    igraph_vector_int_t res;
    igraph_vector_int_init(&res, 0);
    igraph_vector_int_clear(bound_b);

    // 用动态数组做队列
    igraph_integer_t *queue = malloc(n * sizeof(igraph_integer_t));
    if (!queue) {
        free(visited);
        return IGRAPH_ENOMEM;
    }
    igraph_integer_t head = 0, tail = 0;

    // 初始化队列
    queue[tail++] = vertex;
    visited[vertex] = 1;

    igraph_vector_int_t neighbors;
    igraph_vector_int_t parents;
    igraph_vector_int_init(&neighbors, 0);
    igraph_vector_int_init(&parents, 0);

    while (head < tail) {
        igraph_integer_t cur = queue[head++];

        igraph_vector_int_push_back(&res, cur);

        IGRAPH_CHECK(igraph_neighbors(graph, &neighbors, cur, IGRAPH_ALL));
        igraph_integer_t nei_count = igraph_vector_int_size(&neighbors);
        for (igraph_integer_t i = 0; i < nei_count; i++) {
            igraph_integer_t w = VECTOR(neighbors)[i];

            // 如果是禁忌节点 (visited[w] == 2)，将其加入 bound_b，并标记为已处理
            if (visited[w] == 2) {
                igraph_vector_int_push_back(bound_b, w);  // 将禁忌节点 w 加入 bound_b
                visited[w] = 1;  // 标记禁忌节点为已处理
                continue;  // 跳过继续处理该节点，直接进入下一个邻居
            }

            if (!visited[w]) {
                visited[w] = 1;
                queue[tail++] = w;
            }
        }
    }

    igraph_vector_int_destroy(&neighbors);
    igraph_vector_int_destroy(&res);
    free(queue);
    free(visited);

    return IGRAPH_SUCCESS;
}


/**
 * get_minimal_collapsible - 基于节点集 r_nodes 对图 graph 进行 CMSA（Clique Minimal Separator Algorithm）凸包扩展，
 *              并输出扩展后的节点集 H_out。
 *
 * 输入参数：
 *   - graph: 指向已初始化的 igraph_t 图结构的指针，表示待处理的图。
 *   - r_nodes: 指向 igraph_vector_int_t 类型的向量，包含初始的节点集合（即起始的禁忌节点集合）。
 *
 * 输出参数：
 *   - H_out: 指向 igraph_vector_int_t 类型的向量，用于返回扩展后的节点集合结果。
 *
 * 功能说明：
 *   本函数实现对输入图的递归处理，通过计算禁忌节点集 H 中节点的弱连通分量及其边界，
 *   在边界节点中寻找非邻接节点对，并通过调用 close_separator 函数扩展禁忌节点集 H。
 *   该过程重复进行，直到不能再找到新的非邻接节点对，停止迭代。
 *
 * 具体步骤：
 *   1. 初始化节点集合 H 为 r_nodes 的拷贝。
 *   2. 在迭代中，调用 components_forbidden 函数计算带有禁忌节点的连通分量及边界节点。
 *   3. 遍历每个连通分量的边界节点集合，寻找非邻接节点对。
 *   4. 对于每对非邻接节点，调用 close_separator 函数计算最小分割点集合，并将结果加入节点集合 H。
 *   5. 当所有连通分量的边界节点均无非邻接节点对时，迭代终止。
 *
 * 返回值：
 *   - IGRAPH_SUCCESS 表示函数执行成功。
 *   - 其他返回码表示执行过程中发生错误。
 *
 * 注意事项：
 *   - 本函数在迭代过程中会频繁分配和释放内存，可能影响性能，建议结合预分配和向量重用以优化。
 *   - 输入图应已初始化且有效，r_nodes 中的节点应均在图顶点范围内。
 *   - 输出向量 H_out 会被清空并重新初始化。
 */


igraph_error_t get_minimal_collapsible(
    const igraph_t *graph,
    const igraph_vector_int_t *r_nodes,
    igraph_vector_int_t *H_out
) {
    igraph_vector_int_t H;                 // 当前节点集合 H
    igraph_vector_ptr_t components;       // 连通分量列表
    igraph_vector_ptr_t boundaries;       // 边界节点列表
    igraph_bool_t s = 1;                   // 是否继续迭代标志
    igraph_error_t ret = IGRAPH_SUCCESS;

    igraph_integer_t n = igraph_vcount(graph);

    // 初始化 H，预分配空间
    //igraph_vector_int_init(&H,0);
    //igraph_vector_int_reserve(&H, n);
    //igraph_vector_int_copy(&H, r_nodes);

    igraph_vector_int_init_copy(&H, r_nodes);

    

    igraph_vector_ptr_init(&components, 0);
    igraph_vector_ptr_init(&boundaries, 0);

    // 预分配循环内常用向量，单线程分配
    igraph_vector_int_t sub_nodes;
    igraph_vector_int_t neighbors_a, neighbors_b;
    igraph_vector_int_t sep_a_local, sep_b_local;

    igraph_vector_int_init(&sub_nodes,0);
    igraph_vector_int_reserve(&sub_nodes, n);

    igraph_vector_int_init(&neighbors_a,0);
    igraph_vector_int_reserve(&neighbors_a, n);

    igraph_vector_int_init(&neighbors_b,0);
    igraph_vector_int_reserve(&neighbors_b, n);

    igraph_vector_int_init(&sep_a_local,0);
    igraph_vector_int_reserve(&sep_a_local, n);

    igraph_vector_int_init(&sep_b_local,0);
    igraph_vector_int_reserve(&sep_b_local, n);

    igraph_vector_int_t map, invmap;
    igraph_vector_int_init(&map, 0);
    igraph_vector_int_init(&invmap, 0);

    while (s) {
        s = 0;

        // 清理上次结果
        for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&components); i++) {
            igraph_vector_int_t *comp = (igraph_vector_int_t*) VECTOR(components)[i];
            igraph_vector_int_destroy(comp);
            free(comp);
        }
        for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&boundaries); i++) {
            igraph_vector_int_t *bound = (igraph_vector_int_t*) VECTOR(boundaries)[i];
            igraph_vector_int_destroy(bound);
            free(bound);
        }

        igraph_vector_ptr_clear(&components);
        igraph_vector_ptr_clear(&boundaries);

        components_forbidden(graph, &components, &boundaries, &H);

        igraph_integer_t ncomp = igraph_vector_ptr_size(&components);

        // 开启并行循环，每个线程维护局部H_local，循环结束后合并
        #pragma omp parallel
        {
            igraph_vector_int_t H_local;
            igraph_vector_int_init(&H_local,0);
            igraph_vector_int_reserve(&H_local, n);

            #pragma omp for schedule(dynamic)
            for (igraph_integer_t cidx = 0; cidx < ncomp; cidx++) {
                igraph_vector_int_t *component = (igraph_vector_int_t*) VECTOR(components)[cidx];
                igraph_vector_int_t *boundary = (igraph_vector_int_t*) VECTOR(boundaries)[cidx];

                if (igraph_vector_int_size(boundary) < 2) {
                    continue;
                }

                igraph_vector_int_sort(boundary);
                igraph_integer_t bsize = igraph_vector_int_size(boundary);

                igraph_bool_t found = 0;

                for (igraph_integer_t j = 0; j < bsize && !found; j++) {
                    igraph_integer_t hj = VECTOR(*boundary)[j];
                    for (igraph_integer_t k = j + 1; k < bsize && !found; k++) {
                        igraph_integer_t hk = VECTOR(*boundary)[k];

                        igraph_bool_t are_adj = 0;
                        igraph_are_adjacent(graph, hj, hk, &are_adj);
                        if (!are_adj) {
                            found = 1;
                            s = 1; // 标记继续迭代
                            // 清空并复用sub_nodes
                            igraph_vector_int_clear(&sub_nodes);
                            for (igraph_integer_t ci = 0; ci < igraph_vector_int_size(component); ci++) {
                                igraph_vector_int_push_back(&sub_nodes, VECTOR(*component)[ci]);
                            }
                            igraph_vector_int_push_back(&sub_nodes, hj);
                            igraph_vector_int_push_back(&sub_nodes, hk);
                            //igraph_vector_int_sort(&sub_nodes);

                            igraph_t subgraph; 
                            igraph_vector_int_init(&map, 0);
                            igraph_vector_int_init(&invmap, 0);

                            // 创建子图
                            igraph_vs_t vs;
                            igraph_vs_vector(&vs, &sub_nodes);
                            IGRAPH_CHECK(igraph_induced_subgraph_map(graph, &subgraph, vs,
                                                                        IGRAPH_SUBGRAPH_AUTO, &map, &invmap));
                            igraph_vs_destroy(&vs);
                            // 通过 map 获取局部编号
                            igraph_integer_t a = VECTOR(map)[hj] - 1;
                            igraph_integer_t b = VECTOR(map)[hk] - 1;

                            // 固定采用 CMSA 路线：调用两次 close_separator
                            igraph_vector_int_clear(&neighbors_a);
                            igraph_vector_int_clear(&neighbors_b);
                            igraph_neighbors(&subgraph, &neighbors_a, a, IGRAPH_ALL);
                            igraph_neighbors(&subgraph, &neighbors_b, b, IGRAPH_ALL);


                            igraph_vector_int_clear(&sep_a_local);
                            igraph_vector_int_clear(&sep_b_local);

                            close_separator(&subgraph, b, &neighbors_a, &sep_a_local);
                            close_separator(&subgraph, a, &neighbors_b, &sep_b_local);


                            for (igraph_integer_t si = 0; si < igraph_vector_int_size(&sep_a_local); si++) {
                                igraph_integer_t gv = VECTOR(invmap)[VECTOR(sep_a_local)[si]];
                                igraph_vector_int_push_back(&H_local, gv);
                            }
                            for (igraph_integer_t si = 0; si < igraph_vector_int_size(&sep_b_local); si++) {
                                igraph_integer_t gv = VECTOR(invmap)[VECTOR(sep_b_local)[si]];
                                igraph_vector_int_push_back(&H_local, gv);
                            }
                            igraph_destroy(&subgraph);
                            igraph_vector_int_destroy(&map);
                            igraph_vector_int_destroy(&invmap);
                            
                        }
                    }
                }
            }
            // 并行区域结束前合并线程局部结果到全局H，使用临界区保护
            #pragma omp critical
            {
                for (igraph_integer_t i = 0; i < igraph_vector_int_size(&H_local); i++) {
                    igraph_vector_int_push_back(&H, VECTOR(H_local)[i]);
                }
            }

            igraph_vector_int_destroy(&H_local);
        } // omp parallel end

        // 归一化H，去重
        vector_int_unique(&H);
    }

    igraph_vector_int_clear(H_out);
    igraph_vector_int_init_copy(H_out, &H);

    // 释放组件和边界内存
    for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&components); i++) {
        igraph_vector_int_t *comp = (igraph_vector_int_t*) VECTOR(components)[i];
        igraph_vector_int_destroy(comp);
        free(comp);
    }
    for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&boundaries); i++) {
        igraph_vector_int_t *bound = (igraph_vector_int_t*) VECTOR(boundaries)[i];
        igraph_vector_int_destroy(bound);
        free(bound);
    }

    igraph_vector_ptr_destroy(&components);
    igraph_vector_ptr_destroy(&boundaries);

    igraph_vector_int_destroy(&H);
    igraph_vector_int_destroy(&map);
    igraph_vector_int_destroy(&invmap);

    igraph_vector_int_destroy(&sub_nodes);
    igraph_vector_int_destroy(&neighbors_a);
    igraph_vector_int_destroy(&neighbors_b);
    igraph_vector_int_destroy(&sep_a_local);
    igraph_vector_int_destroy(&sep_b_local);

    return ret;
}


/*
 * decompose_atoms - 对无向图进行原子分解，并构建原子分解树。
 *
 * 本函数首先对输入无向图 g 使用最大基数搜索（MCS）计算节点排序，
 * 然后基于该排序迭代进行递归分解。每次从当前子图中选取最小 alpha 值节点 v，
 * 对其闭包 N[G'][v] 进行 atom 扩展，生成一个 atom U 并将其加入 atoms 列表。
 * 对于当前任务队列中的每个pair (S, U_p)，若 S ⊆ U，则在 tree_out 中
 * 添加以 U 为节点的新边，建立 atom tree 的父子关系。
 *
 * 接着，计算 G' - U 的连通分量，每个分量 M 对应的分离集 S_M = N_{G'}(M)，
 * 将 S_M 加入 separators，同时构建该分量的新任务，继续递归分解。
 *
 * 参数:
 *   - g: 输入无向图，必须是无向图，否则返回 IGRAPH_EINVAL。
 *   - atoms: 输出原子集合向量指针，每个元素为一个原子顶点集合（igraph_vector_int_t *）。
 *   - separators: 输出分离子集向量指针，每个元素为一个分离集顶点集合。
 *   - tree_out: 输出 atom tree 图，节点对应 atoms 中的 atom 索引，边表示包含关系。
 *
 * 返回:
 *   - IGRAPH_SUCCESS 表示成功。
 *   - 其他错误码表示内部计算或内存分配失败。
 */

igraph_error_t decompose_atoms(
    const igraph_t *g,
    igraph_vector_ptr_t *atoms, // 存储 Atoms (A)
    igraph_vector_ptr_t *separators, // 存储 Separators (C)
    igraph_t *tree_out // 存储 atom tree
) {

        // Step 1: MCS
    igraph_vector_int_t mcs_alpha, mcs_order;
    igraph_vector_int_init(&mcs_alpha, 0);
    igraph_vector_int_init(&mcs_order, 0);
    igraph_maximum_cardinality_search(g, &mcs_alpha, &mcs_order);


    igraph_vector_int_destroy(&mcs_alpha);

    igraph_vector_int_t mcs_index;
    igraph_vector_int_init(&mcs_index, igraph_vcount(g));
    for (igraph_integer_t i = 0; i < igraph_vcount(g); i++) {
        VECTOR(mcs_index)[VECTOR(mcs_order)[i]] = i;
    }
    igraph_vector_int_destroy(&mcs_order);


    igraph_integer_t n = igraph_vcount(g);

    igraph_vector_ptr_t A_list;
    igraph_vector_ptr_init(&A_list, 0);

    igraph_vector_int_t *Vfull = (igraph_vector_int_t *)malloc(sizeof(igraph_vector_int_t));
    igraph_vector_int_init_range(Vfull, 0, n); 

    /* initial empty W for the whole graph */
    igraph_vector_ptr_t *initial_W = (igraph_vector_ptr_t*)malloc(sizeof(igraph_vector_ptr_t));
    igraph_vector_ptr_init(initial_W, 0);

    /* wrap into a task_t and push onto A_list */
    task_t *initial_task = (task_t*)malloc(sizeof(task_t));
    initial_task->V_nodes = Vfull;
    initial_task->W = initial_W;
    igraph_vector_ptr_push_back(&A_list, initial_task);


    igraph_vector_int_t map, invmap;
    igraph_vector_int_init(&map, 0);
    igraph_vector_int_init(&invmap, 0);

    while (!igraph_vector_ptr_empty(&A_list)) {   
        task_t *task = (task_t*)VECTOR(A_list)[igraph_vector_ptr_size(&A_list) - 1];
        igraph_vector_ptr_pop_back(&A_list);

        igraph_vector_int_t *A_global = task->V_nodes;
        igraph_vector_ptr_t *W_current = task->W;


        igraph_t sub_g; 
        igraph_vector_int_clear(&map); 
        igraph_vector_int_clear(&invmap);

        igraph_vs_t vs;
        igraph_vs_vector(&vs, A_global); 
        IGRAPH_CHECK(igraph_induced_subgraph_map(g, &sub_g, vs, IGRAPH_SUBGRAPH_AUTO, &map, &invmap));
        igraph_vs_destroy(&vs);


        igraph_integer_t min_v_global = VECTOR(*A_global)[0];
        for (igraph_integer_t j = 1; j < igraph_vector_int_size(A_global); j++) {
            igraph_integer_t cand = VECTOR(*A_global)[j];
            if (VECTOR(mcs_index)[cand] < VECTOR(mcs_index)[min_v_global]) {
                min_v_global = cand; 
            }
        }
        igraph_integer_t min_v_local = VECTOR(map)[min_v_global]-1;


        igraph_vector_int_t N_v;
        igraph_vector_int_init(&N_v, 0);
        igraph_neighbors(&sub_g, &N_v, min_v_local, IGRAPH_ALL);
        igraph_vector_int_push_back(&N_v, min_v_local);


        igraph_bool_t is_clique;
        igraph_vs_vector(&vs, &N_v);
        igraph_is_clique(&sub_g, vs, 0, &is_clique);
        igraph_vs_destroy(&vs);


        igraph_vector_int_t H;
        igraph_vector_int_init(&H, 0);

        if (is_clique) {
            igraph_vector_int_append(&H, &N_v);
        } else {
            get_minimal_collapsible(&sub_g, &N_v, &H);
        }
        igraph_vector_int_destroy(&N_v);


        igraph_vector_int_t *H_global = malloc(sizeof(igraph_vector_int_t));
        igraph_vector_int_init(H_global, igraph_vector_int_size(&H));
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(&H); i++) {
            VECTOR(*H_global)[i] = VECTOR(invmap)[VECTOR(H)[i]];
        }
        igraph_vector_ptr_push_back(atoms, H_global);


        /* 新产生的 atom 在 atom 集合中的索引（用于在 tree 中连接） */
        igraph_integer_t u_index = igraph_vector_ptr_size(atoms) - 1;
        /* 为 atom tree 增加对应顶点 */
        IGRAPH_CHECK(igraph_add_vertices(tree_out, 1, 0));

        /* 将 W_current 中 S ⊆ U 的 pairs 连接到当前新 atom（u_index），并从 W_current 中移除 */
        if (W_current) {
            igraph_integer_t wsize = igraph_vector_ptr_size(W_current);
            for (igraph_integer_t wi = 0; wi < wsize; wi++) {
                pair_t *pp = (pair_t*) VECTOR(*W_current)[wi];
                if (!pp) continue;
                if (vector_int_is_subset(&pp->S, H_global)) {
                    /* 添加边 u -> pp->u_p_index */
                    IGRAPH_CHECK(igraph_add_edge(tree_out, u_index, pp->u_p_index));
                    /* 释放该 pair 结构并从 W_current 中删除 */
                    igraph_vector_ptr_remove(W_current, wi);
                    igraph_vector_int_destroy(&pp->S);
                    free(pp);
                    wi--; wsize--;
                }
            }
        }


        igraph_vector_ptr_t comps, bounds;
        igraph_vector_ptr_init(&comps, 0);
        igraph_vector_ptr_init(&bounds, 0);

        components_forbidden(&sub_g, &comps, &bounds, &H);
        igraph_vector_int_destroy(&H);


        igraph_integer_t ncomps = igraph_vector_ptr_size(&comps);
        for (igraph_integer_t i = 0; i < ncomps; i++) {
            igraph_vector_int_t *comp = (igraph_vector_int_t *)VECTOR(comps)[i];
            igraph_vector_int_t *bound = (igraph_vector_int_t *)VECTOR(bounds)[i];

            // ************************************************
            // ** 核心修改：将 N_{G'}(M) 映射回全局索引并存储到 separators **
            // ************************************************
            igraph_vector_int_t *C_global = (igraph_vector_int_t *)malloc(sizeof(igraph_vector_int_t));
            if (!C_global) { /* 错误处理 */ }
            IGRAPH_CHECK(igraph_vector_int_init(C_global, igraph_vector_int_size(bound)));
            for (igraph_integer_t j = 0; j < igraph_vector_int_size(bound); j++) {
                VECTOR(*C_global)[j] = VECTOR(invmap)[VECTOR(*bound)[j]];
            }
            IGRAPH_CHECK(igraph_vector_ptr_push_back(separators, C_global)); // 存储 Separator C


            igraph_vector_int_t A;
            igraph_vector_int_init(&A, 0);
            igraph_vector_int_append(&A, comp);
            igraph_vector_int_append(&A, bound);

            igraph_vector_int_t *A_global_ptr = (igraph_vector_int_t *)malloc(sizeof(igraph_vector_int_t));
            igraph_vector_int_init(A_global_ptr, igraph_vector_int_size(&A));
            for (igraph_integer_t j = 0; j < igraph_vector_int_size(&A); j++) {
                VECTOR(*A_global_ptr)[j] = VECTOR(invmap)[VECTOR(A)[j]];
            }


            /* 为本组件 M 构建对应的 W_M 并创建 task 推入 A_list */
            igraph_vector_ptr_t *W_M = (igraph_vector_ptr_t*)malloc(sizeof(igraph_vector_ptr_t));
            igraph_vector_ptr_init(W_M, 0);

            /* 将 (S_M, u_index) 放入 W_M，S_M 使用 C_global（已映射为全局索引） */
            pair_t *p0 = (pair_t*)malloc(sizeof(pair_t));
            igraph_vector_int_init(&p0->S, igraph_vector_int_size(C_global));
            for (igraph_integer_t jj = 0; jj < igraph_vector_int_size(C_global); jj++) VECTOR(p0->S)[jj] = VECTOR(*C_global)[jj];
            p0->u_p_index = u_index;
            igraph_vector_ptr_push_back(W_M, p0);

            /* 将 W_current 中 S ⊆ (M ∪ S_M) 的 pairs 复制到 W_M（按伪代码分发） */
            if (W_current) {
                igraph_integer_t wsize2 = igraph_vector_ptr_size(W_current);
                for (igraph_integer_t wi2 = 0; wi2 < wsize2; wi2++) {
                    pair_t *pp2 = (pair_t*) VECTOR(*W_current)[wi2];
                    if (!pp2) continue;
                    if (vector_int_is_subset(&pp2->S, A_global_ptr)) {
                        pair_t *pp_copy = (pair_t*)malloc(sizeof(pair_t));
                        igraph_vector_int_init(&pp_copy->S, igraph_vector_int_size(&pp2->S));
                        for (igraph_integer_t kk = 0; kk < igraph_vector_int_size(&pp2->S); kk++) {
                            VECTOR(pp_copy->S)[kk] = VECTOR(pp2->S)[kk];
                        }
                        pp_copy->u_p_index = pp2->u_p_index;
                        igraph_vector_ptr_push_back(W_M, pp_copy);
                    }
                }
            }

            /* 将 W_M 包装为任务并推入 A_list */
            task_t *task_M = (task_t*)malloc(sizeof(task_t));
            task_M->V_nodes = A_global_ptr;
            task_M->W = W_M;
            igraph_vector_ptr_push_back(&A_list, task_M);
            igraph_vector_int_destroy(&A);                     
        }

        for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&comps); i++) {
            igraph_vector_int_destroy((igraph_vector_int_t *)VECTOR(comps)[i]);
            free(VECTOR(comps)[i]);
            igraph_vector_int_destroy((igraph_vector_int_t *)VECTOR(bounds)[i]);
            free(VECTOR(bounds)[i]);
        }
        igraph_vector_ptr_destroy(&comps);
        igraph_vector_ptr_destroy(&bounds);
        igraph_destroy(&sub_g);
        // 清理当前 A_global
        igraph_vector_int_destroy(A_global);
        free(A_global);

        /* 清理当前任务的 W_current */
        if (W_current) {
            for (igraph_integer_t wi3 = 0; wi3 < igraph_vector_ptr_size(W_current); wi3++) {
                pair_t *pp3 = (pair_t*) VECTOR(*W_current)[wi3];
                if (!pp3) continue;
                igraph_vector_int_destroy(&pp3->S);
                free(pp3);
            }
            igraph_vector_ptr_destroy(W_current);
            free(W_current);
        }

        /* 释放当前任务结构 */
        free(task);

    }
    // 循环结束，A_list 已空，map / invmap / mcs_index 清理
    igraph_vector_int_destroy(&mcs_index);
    igraph_vector_int_destroy(&map);
    igraph_vector_int_destroy(&invmap);

    return IGRAPH_SUCCESS;
}






igraph_error_t SAHR(igraph_t *g, igraph_integer_t *r, igraph_integer_t r_size, igraph_integer_t **local2global, igraph_integer_t *result_size) {

    igraph_t g_copy;
    IGRAPH_CHECK(igraph_copy(&g_copy, g));

    igraph_integer_t n = igraph_vcount(&g_copy);
    *local2global = (igraph_integer_t *) malloc(n * sizeof(int));
    if (!*local2global) { igraph_destroy(&g_copy); return IGRAPH_ENOMEM; }
    for (igraph_integer_t i = 0; i < n; i++) (*local2global)[i] = i;

    igraph_integer_t local2global_size = n;  // 维护当前 local2global 的有效长度

    // 初始化 M：不在 r 中的节点
    igraph_integer_t *M = (igraph_integer_t *) malloc(n * sizeof(int));
    if (!M) { free(*local2global); igraph_destroy(&g_copy); return IGRAPH_ENOMEM; }
    igraph_integer_t M_size = 0;
    for (igraph_integer_t i = 0; i < n; i++) {
        igraph_integer_t found = 0;
        for (igraph_integer_t j = 0; j < r_size; j++) if (i == r[j]) { found = 1; break; }
        if (!found) M[M_size++] = i;
    }

    igraph_vector_int_t N;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&N, 0);
    igraph_vector_int_t neighbors_vec;
    igraph_vector_int_init(&neighbors_vec, 0);
    igraph_vs_t vs;

    igraph_integer_t changed = 1;
    while (changed) {
        changed = 0;

        igraph_integer_t N_size = igraph_vector_int_size(&N);

        // 处理上一轮删除节点邻集
        for (igraph_integer_t idx = 0; idx < N_size; idx++) {
            igraph_integer_t node = VECTOR(N)[idx];

            if (node >= igraph_vcount(&g_copy)) continue;

            IGRAPH_CHECK(igraph_neighbors(&g_copy, &neighbors_vec, node, IGRAPH_ALL));

            igraph_t subgraph;
            igraph_vs_vector(&vs, &neighbors_vec);
            igraph_error_t ret = igraph_induced_subgraph(&g_copy, &subgraph, vs, IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH);
            if (ret != IGRAPH_SUCCESS) continue;

            igraph_integer_t n_sub = igraph_vcount(&subgraph);
            igraph_integer_t m_sub = igraph_ecount(&subgraph);
            igraph_integer_t complete_edges = n_sub * (n_sub - 1) / 2;

            if (m_sub == complete_edges) {
                igraph_vs_t vs_delete;
                igraph_vs_1(&vs_delete, node);
                IGRAPH_CHECK(igraph_delete_vertices(&g_copy, vs_delete));
                igraph_vs_destroy(&vs_delete);

                // 更新 M
                igraph_integer_t new_M_size = 0;
                for (igraph_integer_t j = 0; j < M_size; j++) {
                    if (M[j] == node) continue;
                    M[new_M_size++] = (M[j] > node) ? M[j] - 1 : M[j];
                }
                M_size = new_M_size;

                // 更新 local2global
                for (igraph_integer_t j = node; j < local2global_size - 1; j++)
                    (*local2global)[j] = (*local2global)[j + 1];
                local2global_size--;
                *result_size = local2global_size;

                // 更新 N
                igraph_vector_int_t new_N;
                igraph_vector_int_init(&new_N, 0);
                for (igraph_integer_t j = 0; j < igraph_vector_int_size(&neighbors_vec); j++) {
                    igraph_integer_t neighbor = VECTOR(neighbors_vec)[j];
                    if (neighbor > node) neighbor--;
                    for (igraph_integer_t k = 0; k < M_size; k++) {
                        if (M[k] == neighbor) { igraph_vector_int_push_back(&new_N, neighbor); break; }
                    }
                }
                igraph_vector_int_init(&N, igraph_vector_int_size(&new_N)); // 初始化 N 大小
                for (long i = 0; i < igraph_vector_int_size(&new_N); i++) {
                    VECTOR(N)[i] = VECTOR(new_N)[i];
                }
                igraph_vector_int_destroy(&new_N);

                igraph_destroy(&subgraph);
                changed = 1;
                break;
            }
            igraph_destroy(&subgraph);
        }
        if (changed) continue;

        // 遍历剩余 M
        for (igraph_integer_t i = 0; i < M_size; i++) {
            igraph_integer_t node = M[i];

            igraph_integer_t in_N = 0;
            for (igraph_integer_t k = 0; k < igraph_vector_int_size(&N); k++) {
                if (VECTOR(N)[k] == node) {
                    in_N = 1;
                    break;
                }
            }
            if (in_N) continue;

            if (node >= igraph_vcount(&g_copy)) continue;

            IGRAPH_CHECK(igraph_neighbors(&g_copy, &neighbors_vec, node, IGRAPH_ALL));
            igraph_t subgraph;
            igraph_vs_vector(&vs, &neighbors_vec);
            IGRAPH_CHECK(igraph_induced_subgraph(&g_copy, &subgraph, vs, IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH));

            igraph_integer_t n_sub = igraph_vcount(&subgraph);
            igraph_integer_t m_sub = igraph_ecount(&subgraph);
            igraph_integer_t complete_edges = n_sub * (n_sub - 1) / 2;

            if (m_sub == complete_edges) {
                igraph_vs_t vs_delete;
                igraph_vs_1(&vs_delete, node);
                IGRAPH_CHECK(igraph_delete_vertices(&g_copy, vs_delete));
                igraph_vs_destroy(&vs_delete);

                igraph_integer_t new_M_size = 0;
                for (igraph_integer_t j = 0; j < M_size; j++) {
                    if (M[j] == node) continue;
                    M[new_M_size++] = (M[j] > node) ? M[j] - 1 : M[j];
                }
                M_size = new_M_size;

                for (igraph_integer_t j = node; j < local2global_size - 1; j++)
                    (*local2global)[j] = (*local2global)[j + 1];
                local2global_size--;
                *result_size = local2global_size;

                igraph_vector_int_t new_N;
                igraph_vector_int_init(&new_N, 0);
                for (igraph_integer_t j = 0; j < igraph_vector_int_size(&neighbors_vec); j++) {
                    igraph_integer_t neighbor = VECTOR(neighbors_vec)[j];
                    if (neighbor > node) neighbor--;
                    for (igraph_integer_t k = 0; k < M_size; k++) {
                        if (M[k] == neighbor) { igraph_vector_int_push_back(&new_N, neighbor); break; }
                    }
                }
                igraph_vector_int_init(&N, igraph_vector_int_size(&new_N)); // 初始化 N 大小
                for (long i = 0; i < igraph_vector_int_size(&new_N); i++) {
                    VECTOR(N)[i] = VECTOR(new_N)[i];
                }
                igraph_vector_int_destroy(&new_N);

                igraph_destroy(&subgraph);
                changed = 1;
                break;
            }
            igraph_destroy(&subgraph);
        }
    }

    *result_size = local2global_size;

    free(M);
    igraph_vector_int_destroy(&neighbors_vec);
    igraph_vector_int_destroy(&N);
    igraph_destroy(&g_copy);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}





/*
 * 实现最大基数搜索（MCS）算法，返回节点的 MCS 序列和对应的编号。
 * BLAIR, J. R. and PEYTON, B. (1993). An introduction to chordal graphs and clique trees. In Graph theory
and sparse matrix computation 1–29. Springer.
 * 参数:
 *  - graph: 输入图的常量指针
 *  - alpha: 输出参数，存储每个节点在 MCS 序列中的编号
 *  - alpham1: 输出参数，存储 MCS 序列中每个位置对应的节点编号
 * 返回:
 *  - igraph_error_t 类型，标示函数执行状态
 */
igraph_error_t mcs_with_cliques(const igraph_t *graph,
                                      igraph_vector_int_t *alpha,
                                      igraph_vector_int_t *alpham1,
                                      igraph_vector_ptr_t *cliques) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t size;
    igraph_vector_int_t head, next, prev; /* doubly linked list with head */
    igraph_integer_t i, j, v, x, k, len, w, ws, nw, pw;
    igraph_adjlist_t adjlist;
    igraph_vector_int_t *neis;
    igraph_integer_t pre_card = 0, new_card = 0;
    igraph_vector_int_t current_clique, L_set;
    igraph_vector_int_t *clique_copy;
    igraph_bool_t cliques_enabled;

    cliques_enabled = (cliques != NULL);

    if (cliques_enabled) {
        igraph_vector_ptr_clear(cliques);
        IGRAPH_CHECK(igraph_vector_int_init(&L_set, 0));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &L_set);
        IGRAPH_CHECK(igraph_vector_int_init(&current_clique, 0));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &current_clique);
    }

    if (no_of_nodes == 0) {
        igraph_vector_int_clear(alpha);
        if (alpham1) {
            igraph_vector_int_clear(alpham1);
        }
        if (cliques_enabled) {
            IGRAPH_FINALLY_CLEAN(2);
            igraph_vector_int_destroy(&current_clique);
            igraph_vector_int_destroy(&L_set);
        }
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&size, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&head, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&next, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&prev, no_of_nodes);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    IGRAPH_CHECK(igraph_vector_int_resize(alpha, no_of_nodes));
    if (alpham1) {
        IGRAPH_CHECK(igraph_vector_int_resize(alpham1, no_of_nodes));
    }

    /***********************************************/
    /* for i in [0,n-1] -> set(i) := emptyset rof; */
    /***********************************************/

    /* nothing to do, 'head' contains all zeros */

    /*********************************************************/
    /* for v in vertices -> size(v):=0; add v to set(0) rof; */
    /*********************************************************/

    VECTOR(head)[0] = 1;
    for (v = 0; v < no_of_nodes; v++) {
        VECTOR(next)[v] = v + 2;
        VECTOR(prev)[v] = v;
    }
    VECTOR(next)[no_of_nodes - 1] = 0;
    /* size is already all zero */

    /***************/
    /* i:=n; j:=0; */
    /***************/

    i = no_of_nodes; j = 0;
    pre_card = 0;

    /**************/
    /* do i>=1 -> */
    /**************/

    while (i >= 1) {
        /********************************/
        /* v :=  delete any from set(j) */
        /********************************/

        v = VECTOR(head)[j] - 1;
        x = VECTOR(next)[v];
        VECTOR(head)[j] = x;
        if (x != 0) {
            VECTOR(prev)[x - 1] = 0;
        }

        /*************************************************/
        /* alpha(v) := i; alpham1(i) := v; size(v) := -1 */
        /*************************************************/

        VECTOR(*alpha)[v] = i - 1;
        if (alpham1) {
            VECTOR(*alpham1)[i - 1] = v;
        }
        VECTOR(size)[v] = -1;

        /* Calculate new_card = |ne(v) ∩ L_{i+1}| and manage cliques */
        if (cliques_enabled) {
            new_card = 0;
            neis = igraph_adjlist_get(&adjlist, v);
            len = igraph_vector_int_size(neis);
            for (k = 0; k < len; k++) {
                w = VECTOR(*neis)[k];
                if (VECTOR(size)[w] == -1) {
                    /* w is in L_set (already processed) */
                    new_card++;
                }
            }
            
            /* if new_card <= pre_card, start a new clique */
            if (new_card <= pre_card) {
                /* Save the previous clique if it's not empty */
                if (igraph_vector_int_size(&current_clique) > 0) {
                    /* Remove duplicates before saving */
                    igraph_vector_int_sort(&current_clique);
                    igraph_integer_t write_idx = 0;
                    for (int ciq_i = 0; ciq_i < igraph_vector_int_size(&current_clique); ciq_i++) {
                        if (ciq_i == 0 || VECTOR(current_clique)[ciq_i] != VECTOR(current_clique)[ciq_i - 1]) {
                            VECTOR(current_clique)[write_idx++] = VECTOR(current_clique)[ciq_i];
                        }
                    }
                    igraph_vector_int_resize(&current_clique, write_idx);
                    
                    clique_copy = (igraph_vector_int_t *)malloc(sizeof(igraph_vector_int_t));
                    if (!clique_copy) {
                        IGRAPH_ERROR("Memory allocation failed", IGRAPH_ENOMEM);
                    }
                    igraph_vector_int_init_copy(clique_copy, &current_clique);
                    igraph_vector_ptr_push_back(cliques, clique_copy);
                }
                /* Start new clique: clear and add neighbors in L_set */
                igraph_vector_int_clear(&current_clique);
            }
            
            /* Add neighbors in L_set to current clique */
            neis = igraph_adjlist_get(&adjlist, v);
            len = igraph_vector_int_size(neis);
            for (k = 0; k < len; k++) {
                w = VECTOR(*neis)[k];
                if (VECTOR(size)[w] == -1) {
                    /* w is in L_set */
                    igraph_vector_int_push_back(&current_clique, w);
                }
            }
            
            /* Add v itself to the clique */
            igraph_vector_int_push_back(&current_clique, v);
            
            /* Update pre_card */
            pre_card = new_card;
            
            /* Add v to L_set */
            igraph_vector_int_push_back(&L_set, v);
        }

        /********************************************/
        /* for {v,w} in E such that size(w) >= 0 -> */
        /********************************************/

        neis = igraph_adjlist_get(&adjlist, v);
        len = igraph_vector_int_size(neis);
        for (k = 0; k < len; k++) {
            w = VECTOR(*neis)[k];
            ws = VECTOR(size)[w];
            if (ws >= 0) {

                /******************************/
                /* delete w from set(size(w)) */
                /******************************/

                nw = VECTOR(next)[w];
                pw = VECTOR(prev)[w];
                if (nw != 0) {
                    VECTOR(prev)[nw - 1] = pw;
                }
                if (pw != 0) {
                    VECTOR(next)[pw - 1] = nw;
                } else {
                    VECTOR(head)[ws] = nw;
                }

                /******************************/
                /* size(w) := size(w)+1       */
                /******************************/

                VECTOR(size)[w] += 1;

                /******************************/
                /* add w to set(size(w))      */
                /******************************/

                ws = VECTOR(size)[w];
                nw = VECTOR(head)[ws];
                VECTOR(next)[w] = nw;
                VECTOR(prev)[w] = 0;
                if (nw != 0) {
                    VECTOR(prev)[nw - 1] = w + 1;
                }
                VECTOR(head)[ws] = w + 1;

            }
        }

        /***********************/
        /* i := i-1; j := j+1; */
        /***********************/

        i -= 1;
        j += 1;

        /*********************************************/
        /* do j>=0 and set(j)=emptyset -> j:=j-1; od */
        /*********************************************/

        if (j < no_of_nodes) {
            while (j >= 0 && VECTOR(head)[j] == 0) {
                j--;
            }
        }
    }

    /* Save the last clique if cliques output is enabled */
    if (cliques_enabled && igraph_vector_int_size(&current_clique) > 0) {
        /* Remove duplicates before saving final clique */
        igraph_vector_int_sort(&current_clique);
        igraph_integer_t write_idx = 0;
        for (int ciq_i = 0; ciq_i < igraph_vector_int_size(&current_clique); ciq_i++) {
            if (ciq_i == 0 || VECTOR(current_clique)[ciq_i] != VECTOR(current_clique)[ciq_i - 1]) {
                VECTOR(current_clique)[write_idx++] = VECTOR(current_clique)[ciq_i];
            }
        }
        igraph_vector_int_resize(&current_clique, write_idx);
        
        clique_copy = (igraph_vector_int_t *)malloc(sizeof(igraph_vector_int_t));
        if (!clique_copy) {
            IGRAPH_ERROR("Memory allocation failed", IGRAPH_ENOMEM);
        }
        igraph_vector_int_init_copy(clique_copy, &current_clique);
        igraph_vector_ptr_push_back(cliques, clique_copy);
    }

    igraph_adjlist_destroy(&adjlist);
    if (cliques_enabled) {
        igraph_vector_int_destroy(&current_clique);
        igraph_vector_int_destroy(&L_set);
        IGRAPH_FINALLY_CLEAN(2);
    }
    igraph_vector_int_destroy(&prev);
    igraph_vector_int_destroy(&next);
    igraph_vector_int_destroy(&head);
    igraph_vector_int_destroy(&size);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}











/*
 * 下面是有向图凸包找寻算法
 * 下面是一些辅助数据结构和函数，用于递归图分解过程中的工作空间管理和节点对标记。
 * 这些结构和函数主要用于存储和操作在分解过程中需要频繁访问的节点状态、掩码和队列等信息，
 * 以提高算法的效率和可读性。
 */



 // 精简版：获取马尔可夫边界（不包含自身）
/*
 * 计算 DAG 中节点 target 的 Markov blanket：
 *
 *   mb(target)
 *     = parents(target)
 *       ∪ children(target)
 *       ∪ parents(children(target))
 *       \ {target}
 *
 * blanket 必须已经初始化。
 * 函数会先清空 blanket。
 */
static igraph_error_t get_markov_blanket(
    const igraph_t *graph,
    igraph_integer_t target,
    igraph_vector_int_t *blanket
) {
    if (!graph || !blanket) {
        return IGRAPH_EINVAL;
    }

    const igraph_integer_t n = igraph_vcount(graph);

    if (target < 0 || target >= n) {
        return IGRAPH_EINVVID;
    }

    igraph_error_t ret = IGRAPH_SUCCESS;
    igraph_vector_int_t parents;
    igraph_vector_int_t children;

    ret = igraph_vector_int_init(&parents, 0);
    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }

    ret = igraph_vector_int_init(&children, 0);
    if (ret != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&parents);
        return ret;
    }

    igraph_vector_int_clear(blanket);

    /*
     * 1. target 的父节点
     */
    ret = igraph_neighbors(
        graph,
        &parents,
        target,
        IGRAPH_IN
    );
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    ret = igraph_vector_int_append(blanket, &parents);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    /*
     * 2. target 的子节点
     */
    ret = igraph_neighbors(
        graph,
        &children,
        target,
        IGRAPH_OUT
    );
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    ret = igraph_vector_int_append(blanket, &children);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    /*
     * 3. target 每个孩子的父节点
     *
     * 这些节点包括：
     *   - target 自己；
     *   - target 的配偶节点。
     *
     * target 最后会被删除。
     */
    for (igraph_integer_t i = 0;
         i < igraph_vector_int_size(&children);
         ++i) {

        const igraph_integer_t child =
            VECTOR(children)[i];

        ret = igraph_neighbors(
            graph,
            &parents,
            child,
            IGRAPH_IN
        );
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup;
        }

        ret = igraph_vector_int_append(
            blanket,
            &parents
        );
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup;
        }
    }

    /*
     * 4. 排序、去重，并删除 target 自身
     */
    igraph_vector_int_sort(blanket);

    {
        igraph_integer_t write_pos = 0;
        const igraph_integer_t size =
            igraph_vector_int_size(blanket);

        for (igraph_integer_t read_pos = 0;
             read_pos < size;
             ++read_pos) {

            const igraph_integer_t node =
                VECTOR(*blanket)[read_pos];

            /* Markov blanket 不包含目标节点自身 */
            if (node == target) {
                continue;
            }

            /* 去除重复节点 */
            if (write_pos == 0 ||
                VECTOR(*blanket)[write_pos - 1] != node) {

                VECTOR(*blanket)[write_pos++] = node;
            }
        }

        ret = igraph_vector_int_resize(
            blanket,
            write_pos
        );
    }

cleanup:
    igraph_vector_int_destroy(&parents);
    igraph_vector_int_destroy(&children);

    return ret;
}


typedef struct dag_moral_cache_t {
    const igraph_t *graph;
    igraph_integer_t n;
    igraph_adjlist_t parents;
    igraph_adjlist_t children;
    igraph_vector_ptr_t blankets;
    unsigned char *blanket_ready;
    igraph_bool_t parents_initialized;
    igraph_bool_t children_initialized;
    igraph_bool_t blankets_initialized;
} dag_moral_cache_t;


static void dag_moral_cache_destroy(dag_moral_cache_t *cache) {
    if (!cache) {
        return;
    }

    if (cache->blankets_initialized) {
        for (igraph_integer_t i = 0; i < igraph_vector_ptr_size(&cache->blankets); ++i) {
            igraph_vector_int_t *blanket =
                (igraph_vector_int_t *) VECTOR(cache->blankets)[i];
            if (blanket) {
                igraph_vector_int_destroy(blanket);
                free(blanket);
            }
        }
        igraph_vector_ptr_destroy(&cache->blankets);
    }

    if (cache->children_initialized) {
        igraph_adjlist_destroy(&cache->children);
    }
    if (cache->parents_initialized) {
        igraph_adjlist_destroy(&cache->parents);
    }

    free(cache->blanket_ready);
    memset(cache, 0, sizeof(*cache));
}


static igraph_error_t dag_moral_cache_init(
    dag_moral_cache_t *cache,
    const igraph_t *graph
) {
    if (!cache || !graph) {
        return IGRAPH_EINVAL;
    }

    memset(cache, 0, sizeof(*cache));
    cache->graph = graph;
    cache->n = igraph_vcount(graph);

    igraph_error_t ret = igraph_adjlist_init(
        graph,
        &cache->parents,
        IGRAPH_IN,
        IGRAPH_NO_LOOPS,
        IGRAPH_NO_MULTIPLE
    );
    if (ret != IGRAPH_SUCCESS) {
        goto fail;
    }
    cache->parents_initialized = 1;

    ret = igraph_adjlist_init(
        graph,
        &cache->children,
        IGRAPH_OUT,
        IGRAPH_NO_LOOPS,
        IGRAPH_NO_MULTIPLE
    );
    if (ret != IGRAPH_SUCCESS) {
        goto fail;
    }
    cache->children_initialized = 1;

    ret = igraph_vector_ptr_init(&cache->blankets, cache->n);
    if (ret != IGRAPH_SUCCESS) {
        goto fail;
    }
    cache->blankets_initialized = 1;
    for (igraph_integer_t i = 0; i < cache->n; ++i) {
        VECTOR(cache->blankets)[i] = NULL;
    }

    cache->blanket_ready = (unsigned char *) calloc(
        (size_t) (cache->n > 0 ? cache->n : 1),
        sizeof(unsigned char)
    );
    if (!cache->blanket_ready) {
        ret = IGRAPH_ENOMEM;
        goto fail;
    }

    return IGRAPH_SUCCESS;

fail:
    dag_moral_cache_destroy(cache);
    return ret;
}


static igraph_error_t dag_moral_cache_get_blanket(
    dag_moral_cache_t *cache,
    igraph_integer_t target,
    const igraph_vector_int_t **blanket_out
) {
    if (!cache || !blanket_out || target < 0 || target >= cache->n) {
        return IGRAPH_EINVVID;
    }

    if (!cache->blanket_ready[target]) {
        igraph_vector_int_t *blanket =
            (igraph_vector_int_t *) VECTOR(cache->blankets)[target];

        if (!blanket) {
            blanket = (igraph_vector_int_t *) malloc(sizeof(igraph_vector_int_t));
            if (!blanket) {
                return IGRAPH_ENOMEM;
            }

            igraph_error_t ret = igraph_vector_int_init(blanket, 0);
            if (ret != IGRAPH_SUCCESS) {
                free(blanket);
                return ret;
            }
            VECTOR(cache->blankets)[target] = blanket;
        } else {
            igraph_vector_int_clear(blanket);
        }

        igraph_vector_int_t *parents =
            igraph_adjlist_get(&cache->parents, target);
        igraph_vector_int_t *children =
            igraph_adjlist_get(&cache->children, target);

        igraph_error_t ret = igraph_vector_int_append(blanket, parents);
        if (ret != IGRAPH_SUCCESS) {
            return ret;
        }
        ret = igraph_vector_int_append(blanket, children);
        if (ret != IGRAPH_SUCCESS) {
            return ret;
        }

        for (igraph_integer_t i = 0; i < igraph_vector_int_size(children); ++i) {
            const igraph_integer_t child = VECTOR(*children)[i];
            igraph_vector_int_t *child_parents =
                igraph_adjlist_get(&cache->parents, child);

            ret = igraph_vector_int_append(blanket, child_parents);
            if (ret != IGRAPH_SUCCESS) {
                return ret;
            }
        }

        igraph_vector_int_sort(blanket);

        igraph_integer_t write_pos = 0;
        const igraph_integer_t size = igraph_vector_int_size(blanket);
        for (igraph_integer_t read_pos = 0; read_pos < size; ++read_pos) {
            const igraph_integer_t node = VECTOR(*blanket)[read_pos];
            if (node == target) {
                continue;
            }
            if (write_pos == 0 || node != VECTOR(*blanket)[write_pos - 1]) {
                VECTOR(*blanket)[write_pos++] = node;
            }
        }

        ret = igraph_vector_int_resize(blanket, write_pos);
        if (ret != IGRAPH_SUCCESS) {
            return ret;
        }

        cache->blanket_ready[target] = 1;
    }

    *blanket_out =
        (const igraph_vector_int_t *) VECTOR(cache->blankets)[target];
    return IGRAPH_SUCCESS;
}


static igraph_error_t dag_collect_ancestor_nodes_cached(
    dag_moral_cache_t *cache,
    const igraph_vector_int_t *seeds,
    igraph_vector_int_t *ancestors_out,
    igraph_integer_t *stack,
    unsigned char *seen
) {
    if (!cache || !seeds || !ancestors_out || !stack || !seen) {
        return IGRAPH_EINVAL;
    }

    const igraph_integer_t n = cache->n;
    memset(seen, 0, (size_t) n);
    igraph_vector_int_clear(ancestors_out);

    igraph_integer_t top = 0;
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(seeds); ++i) {
        const igraph_integer_t seed = VECTOR(*seeds)[i];
        if (seed < 0 || seed >= n) {
            return IGRAPH_EINVVID;
        }
        if (!seen[seed]) {
            seen[seed] = 1;
            stack[top++] = seed;
            igraph_error_t ret = igraph_vector_int_push_back(ancestors_out, seed);
            if (ret != IGRAPH_SUCCESS) {
                return ret;
            }
        }
    }

    while (top > 0) {
        const igraph_integer_t cur = stack[--top];
        igraph_vector_int_t *parents =
            igraph_adjlist_get(&cache->parents, cur);

        for (igraph_integer_t i = 0; i < igraph_vector_int_size(parents); ++i) {
            const igraph_integer_t parent = VECTOR(*parents)[i];
            if (!seen[parent]) {
                seen[parent] = 1;
                stack[top++] = parent;
                igraph_error_t ret =
                    igraph_vector_int_push_back(ancestors_out, parent);
                if (ret != IGRAPH_SUCCESS) {
                    return ret;
                }
            }
        }
    }

    igraph_vector_int_sort(ancestors_out);
    return IGRAPH_SUCCESS;
}



static igraph_error_t dag_ancestor_moral_subgraph_with_mapping(
    const igraph_t *graph,
    const igraph_vector_int_t *r_nodes,
    igraph_t *ug_graph,
    igraph_vector_int_t *global_to_local,
    igraph_vector_int_t *local_to_global,
    igraph_vector_int_t *r_local_out
) {
    if (!graph || !r_nodes || !ug_graph || !global_to_local ||
        !local_to_global || !r_local_out) {
        return IGRAPH_EINVAL;
    }

    if (!igraph_is_directed(graph)) {
        return IGRAPH_EINVAL;
    }

    igraph_error_t ret = IGRAPH_SUCCESS;
    const igraph_integer_t n = igraph_vcount(graph);

    if (igraph_vector_int_size(r_nodes) == 0) {
        igraph_vector_int_clear(global_to_local);
        igraph_vector_int_clear(local_to_global);
        igraph_vector_int_clear(r_local_out);
        return igraph_empty(ug_graph, 0, IGRAPH_UNDIRECTED);
    }

    igraph_vector_int_t ancestors;
    ret = igraph_vector_int_init(&ancestors, 0);
    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }

    igraph_integer_t *ancestor_stack = (igraph_integer_t *) malloc(
        (size_t) (n > 0 ? n : 1) * sizeof(igraph_integer_t)
    );
    unsigned char *ancestor_seen = (unsigned char *) calloc(
        (size_t) (n > 0 ? n : 1),
        sizeof(unsigned char)
    );

    dag_moral_cache_t cache;
    igraph_bool_t cache_initialized = 0;
    igraph_vs_t vs;
    igraph_bool_t vs_initialized = 0;
    igraph_t induced_graph;
    igraph_bool_t induced_initialized = 0;
    igraph_adjlist_t parents;
    igraph_adjlist_t children;
    igraph_bool_t parents_initialized = 0;
    igraph_bool_t children_initialized = 0;

    if (!ancestor_stack || !ancestor_seen) {
        ret = IGRAPH_ENOMEM;
        goto cleanup;
    }

    ret = dag_moral_cache_init(&cache, graph);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    cache_initialized = 1;

    ret = dag_collect_ancestor_nodes_cached(
        &cache,
        r_nodes,
        &ancestors,
        ancestor_stack,
        ancestor_seen
    );
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    igraph_vector_int_sort(&ancestors);
    vector_int_unique(&ancestors);

    igraph_vector_int_clear(global_to_local);
    igraph_vector_int_clear(local_to_global);
    igraph_vector_int_clear(r_local_out);

    ret = igraph_vector_int_resize(global_to_local, n);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    for (igraph_integer_t i = 0; i < n; ++i) {
        VECTOR(*global_to_local)[i] = 0;
    }

    ret = igraph_vector_int_append(local_to_global, &ancestors);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&ancestors); ++i) {
        const igraph_integer_t global_node = VECTOR(ancestors)[i];
        if (global_node < 0 || global_node >= n) {
            ret = IGRAPH_EINVVID;
            goto cleanup;
        }
        VECTOR(*global_to_local)[global_node] = i + 1;
    }

    ret = igraph_vs_vector(&vs, &ancestors);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    vs_initialized = 1;

    ret = igraph_induced_subgraph_map(
        graph,
        &induced_graph,
        vs,
        IGRAPH_SUBGRAPH_AUTO,
        global_to_local,
        local_to_global
    );
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    induced_initialized = 1;

    igraph_vs_destroy(&vs);
    vs_initialized = 0;

    ret = igraph_empty(ug_graph, igraph_vcount(&induced_graph), IGRAPH_UNDIRECTED);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    /* Collect all edges (skeleton + co-parent) into a single vector for batch addition */
    igraph_vector_int_t all_edges;
    ret = igraph_vector_int_init(&all_edges, 0);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    igraph_vector_int_t edge_list;
    ret = igraph_vector_int_init(&edge_list, 0);
    if (ret != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&all_edges);
        goto cleanup;
    }

    ret = igraph_get_edgelist(&induced_graph, &edge_list, 0);
    if (ret != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&edge_list);
        igraph_vector_int_destroy(&all_edges);
        goto cleanup;
    }

    /* Add skeleton edges (undirected version of directed edges) */
    const igraph_integer_t ecount = igraph_ecount(&induced_graph);
    for (igraph_integer_t ei = 0; ei < ecount; ++ei) {
        const igraph_integer_t u = VECTOR(edge_list)[2 * ei];
        const igraph_integer_t v = VECTOR(edge_list)[2 * ei + 1];
        if (u != v) {
            ret = igraph_vector_int_push_back(&all_edges, u);
            if (ret != IGRAPH_SUCCESS) {
                igraph_vector_int_destroy(&edge_list);
                igraph_vector_int_destroy(&all_edges);
                goto cleanup;
            }
            ret = igraph_vector_int_push_back(&all_edges, v);
            if (ret != IGRAPH_SUCCESS) {
                igraph_vector_int_destroy(&edge_list);
                igraph_vector_int_destroy(&all_edges);
                goto cleanup;
            }
        }
    }
    igraph_vector_int_destroy(&edge_list);

    ret = igraph_adjlist_init(
        &induced_graph,
        &parents,
        IGRAPH_IN,
        IGRAPH_NO_LOOPS,
        IGRAPH_NO_MULTIPLE
    );
    if (ret != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&all_edges);
        goto cleanup;
    }
    parents_initialized = 1;

    /* Add co-parent edges */
    for (igraph_integer_t child = 0; child < igraph_vcount(&induced_graph); ++child) {
        igraph_vector_int_t *parent_list = igraph_adjlist_get(&parents, child);
        const igraph_integer_t parent_count = igraph_vector_int_size(parent_list);
        for (igraph_integer_t i = 0; i < parent_count; ++i) {
            for (igraph_integer_t j = i + 1; j < parent_count; ++j) {
                const igraph_integer_t a = VECTOR(*parent_list)[i];
                const igraph_integer_t b = VECTOR(*parent_list)[j];
                ret = igraph_vector_int_push_back(&all_edges, a);
                if (ret != IGRAPH_SUCCESS) {
                    igraph_vector_int_destroy(&all_edges);
                    goto cleanup;
                }
                ret = igraph_vector_int_push_back(&all_edges, b);
                if (ret != IGRAPH_SUCCESS) {
                    igraph_vector_int_destroy(&all_edges);
                    goto cleanup;
                }
            }
        }
    }

    /* Add all edges at once with automatic deduplication */
    ret = igraph_add_edges(ug_graph, &all_edges, NULL);
    if (ret != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&all_edges);
        goto cleanup;
    }
    igraph_vector_int_destroy(&all_edges);

    for (igraph_integer_t i = 0; i < igraph_vector_int_size(r_nodes); ++i) {
        const igraph_integer_t global_node = VECTOR(*r_nodes)[i];
        const igraph_integer_t mapped = VECTOR(*global_to_local)[global_node];
        if (mapped <= 0) {
            ret = IGRAPH_EINVAL;
            goto cleanup;
        }
        ret = igraph_vector_int_push_back(r_local_out, mapped - 1);
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup;
        }
    }

    vector_int_unique(r_local_out);

cleanup:
    if (vs_initialized) {
        igraph_vs_destroy(&vs);
    }
    if (induced_initialized) {
        igraph_destroy(&induced_graph);
    }
    if (parents_initialized) {
        igraph_adjlist_destroy(&parents);
    }
    if (cache_initialized) {
        dag_moral_cache_destroy(&cache);
    }
    free(ancestor_stack);
    free(ancestor_seen);
    igraph_vector_int_destroy(&ancestors);
    return ret;
}


static igraph_error_t dag_components_forbidden_cached(
    dag_moral_cache_t *cache,
    igraph_vector_ptr_t *components,
    igraph_vector_ptr_t *boundaries,
    const igraph_vector_int_t *forbidden_vertices
) {
    if (!cache || !boundaries) {
        return IGRAPH_EINVAL;
    }

    const igraph_integer_t n = cache->n;
    igraph_error_t ret = IGRAPH_SUCCESS;

    unsigned char *state = (unsigned char *) calloc(
        (size_t) (n > 0 ? n : 1),
        sizeof(unsigned char)
    );
    igraph_integer_t *queue = (igraph_integer_t *) malloc(
        (size_t) (n > 0 ? n : 1) * sizeof(igraph_integer_t)
    );

    if (!state || !queue) {
        free(state);
        free(queue);
        return IGRAPH_ENOMEM;
    }

    if (forbidden_vertices) {
        for (igraph_integer_t i = 0;
             i < igraph_vector_int_size(forbidden_vertices);
             ++i) {
            const igraph_integer_t node = VECTOR(*forbidden_vertices)[i];
            if (node < 0 || node >= n) {
                ret = IGRAPH_EINVVID;
                goto cleanup;
            }
            state[node] = 2;
        }
    }

    for (igraph_integer_t start = 0; start < n; ++start) {
        if (state[start] != 0) {
            continue;
        }

        igraph_vector_int_t *comp = NULL;
        igraph_vector_int_t *bound = NULL;
        igraph_bool_t comp_initialized = 0;
        igraph_bool_t bound_initialized = 0;
        igraph_bool_t comp_stored = 0;

        if (components) {
            comp = (igraph_vector_int_t *) malloc(sizeof(igraph_vector_int_t));
            if (!comp) {
                ret = IGRAPH_ENOMEM;
                goto component_failed;
            }

            ret = igraph_vector_int_init(comp, 0);
            if (ret != IGRAPH_SUCCESS) {
                goto component_failed;
            }
            comp_initialized = 1;
        }

        bound = (igraph_vector_int_t *) malloc(sizeof(igraph_vector_int_t));
        if (!bound) {
            ret = IGRAPH_ENOMEM;
            goto component_failed;
        }

        ret = igraph_vector_int_init(bound, 0);
        if (ret != IGRAPH_SUCCESS) {
            goto component_failed;
        }
        bound_initialized = 1;

        igraph_integer_t head = 0;
        igraph_integer_t tail = 0;
        queue[tail++] = start;
        state[start] = 1;

        while (head < tail) {
            const igraph_integer_t cur = queue[head++];
            const igraph_vector_int_t *moral_neighbors = NULL;

            if (comp) {
                ret = igraph_vector_int_push_back(comp, cur);
                if (ret != IGRAPH_SUCCESS) {
                    goto component_failed;
                }
            }

            ret = dag_moral_cache_get_blanket(cache, cur, &moral_neighbors);
            if (ret != IGRAPH_SUCCESS) {
                goto component_failed;
            }

            for (igraph_integer_t i = 0;
                 i < igraph_vector_int_size(moral_neighbors);
                 ++i) {
                const igraph_integer_t nb = VECTOR(*moral_neighbors)[i];

                if (state[nb] == 2) {
                    ret = igraph_vector_int_push_back(bound, nb);
                    if (ret != IGRAPH_SUCCESS) {
                        goto component_failed;
                    }
                    state[nb] = 3;
                } else if (state[nb] == 0) {
                    state[nb] = 1;
                    queue[tail++] = nb;
                }
            }
        }

        for (igraph_integer_t i = 0; i < igraph_vector_int_size(bound); ++i) {
            state[VECTOR(*bound)[i]] = 2;
        }

        if (components) {
            ret = igraph_vector_ptr_push_back(components, comp);
            if (ret != IGRAPH_SUCCESS) {
                goto component_failed;
            }
            comp_stored = 1;
        }

        ret = igraph_vector_ptr_push_back(boundaries, bound);
        if (ret != IGRAPH_SUCCESS) {
            if (comp_stored) {
                (void) igraph_vector_ptr_pop_back(components);
                comp_stored = 0;
            }
            goto component_failed;
        }

        continue;

component_failed:
        if (comp_initialized) {
            igraph_vector_int_destroy(comp);
        }
        free(comp);

        if (bound_initialized) {
            igraph_vector_int_destroy(bound);
        }
        free(bound);

        goto cleanup;
    }

cleanup:
    free(queue);
    free(state);
    return ret;
}


igraph_error_t dag_components_forbidden(
    const igraph_t *graph,
    igraph_vector_ptr_t *components,
    igraph_vector_ptr_t *boundaries,
    const igraph_vector_int_t *forbidden_vertices
) {
    if (!graph || !boundaries) {
        return IGRAPH_EINVAL;
    }

    if (!igraph_is_directed(graph)) {
        return IGRAPH_EINVAL;
    }

    {
        dag_moral_cache_t cache;
        igraph_error_t ret = dag_moral_cache_init(&cache, graph);
        if (ret != IGRAPH_SUCCESS) {
            return ret;
        }

        ret = dag_components_forbidden_cached(
            &cache,
            components,
            boundaries,
            forbidden_vertices
        );

        dag_moral_cache_destroy(&cache);
        return ret;
    }

    const igraph_integer_t n = igraph_vcount(graph);
    igraph_error_t ret = IGRAPH_SUCCESS;

    /*
     * components 可以为 NULL。
     * dag_get_minimal_collapsible 只需要各分量的边界，
     * 此时不再为分量节点集合分配和写入向量。
     */

    /*
     * state[v]：
     *   0 = 非禁忌，未访问；
     *   1 = 非禁忌，已加入某个分量；
     *   2 = 禁忌，尚未加入当前分量的边界；
     *   3 = 禁忌，已经加入当前分量的边界。
     */
    unsigned char *state =
        (unsigned char *) calloc(
            (size_t) n,
            sizeof(unsigned char)
        );

    igraph_integer_t *queue =
        (igraph_integer_t *) malloc(
            (size_t) n * sizeof(igraph_integer_t)
        );

    if (n > 0 && (!state || !queue)) {
        free(state);
        free(queue);
        return IGRAPH_ENOMEM;
    }

    /* 当前节点在道义图中的邻居，即其 Markov blanket。 */
    igraph_vector_int_t moral_neighbors;

    ret = igraph_vector_int_init(&moral_neighbors, 0);
    if (ret != IGRAPH_SUCCESS) {
        free(state);
        free(queue);
        return ret;
    }

    /* 标记禁忌节点。 */
    if (forbidden_vertices) {
        for (igraph_integer_t i = 0;
             i < igraph_vector_int_size(forbidden_vertices);
             ++i) {

            const igraph_integer_t node =
                VECTOR(*forbidden_vertices)[i];

            if (node < 0 || node >= n) {
                ret = IGRAPH_EINVVID;
                goto cleanup;
            }

            state[node] = 2;
        }
    }

    /* 在隐式道义图 graph^m - forbidden_vertices 上寻找所有连通分量。 */
    for (igraph_integer_t start = 0; start < n; ++start) {
        if (state[start] != 0) {
            continue;
        }

        igraph_vector_int_t *comp = NULL;
        igraph_vector_int_t *bound = NULL;
        igraph_bool_t comp_initialized = 0;
        igraph_bool_t bound_initialized = 0;
        igraph_bool_t comp_stored = 0;

        if (components) {
            comp = (igraph_vector_int_t *) malloc(
                sizeof(igraph_vector_int_t)
            );
            if (!comp) {
                ret = IGRAPH_ENOMEM;
                goto component_failed;
            }

            ret = igraph_vector_int_init(comp, 0);
            if (ret != IGRAPH_SUCCESS) {
                goto component_failed;
            }
            comp_initialized = 1;
        }

        bound = (igraph_vector_int_t *) malloc(
            sizeof(igraph_vector_int_t)
        );
        if (!bound) {
            ret = IGRAPH_ENOMEM;
            goto component_failed;
        }

        ret = igraph_vector_int_init(bound, 0);
        if (ret != IGRAPH_SUCCESS) {
            goto component_failed;
        }
        bound_initialized = 1;

        igraph_integer_t head = 0;
        igraph_integer_t tail = 0;

        queue[tail++] = start;
        state[start] = 1;

        while (head < tail) {
            const igraph_integer_t cur = queue[head++];

            if (comp) {
                ret = igraph_vector_int_push_back(comp, cur);
                if (ret != IGRAPH_SUCCESS) {
                    goto component_failed;
                }
            }

            ret = get_markov_blanket(
                graph,
                cur,
                &moral_neighbors
            );
            if (ret != IGRAPH_SUCCESS) {
                goto component_failed;
            }

            for (igraph_integer_t i = 0;
                 i < igraph_vector_int_size(&moral_neighbors);
                 ++i) {

                const igraph_integer_t nb =
                    VECTOR(moral_neighbors)[i];

                if (state[nb] == 2) {
                    ret = igraph_vector_int_push_back(bound, nb);
                    if (ret != IGRAPH_SUCCESS) {
                        goto component_failed;
                    }

                    /* 同一禁忌节点在当前边界中只记录一次。 */
                    state[nb] = 3;
                } else if (state[nb] == 0) {
                    state[nb] = 1;
                    queue[tail++] = nb;
                }
            }
        }

        /* 状态 3 只对当前分量有效，下一个分量仍可共享同一边界节点。 */
        for (igraph_integer_t i = 0;
             i < igraph_vector_int_size(bound);
             ++i) {
            state[VECTOR(*bound)[i]] = 2;
        }

        if (components) {
            ret = igraph_vector_ptr_push_back(components, comp);
            if (ret != IGRAPH_SUCCESS) {
                goto component_failed;
            }
            comp_stored = 1;
        }

        ret = igraph_vector_ptr_push_back(boundaries, bound);
        if (ret != IGRAPH_SUCCESS) {
            if (comp_stored) {
                (void) igraph_vector_ptr_pop_back(components);
                comp_stored = 0;
            }
            goto component_failed;
        }

        /* 所有权已交给输出容器。 */
        continue;

component_failed:
        if (comp_initialized) {
            igraph_vector_int_destroy(comp);
        }
        free(comp);

        if (bound_initialized) {
            igraph_vector_int_destroy(bound);
        }
        free(bound);

        goto cleanup;
    }

cleanup:
    igraph_vector_int_destroy(&moral_neighbors);
    free(state);
    free(queue);

    return ret;
}


/*
 * 检查无序节点对 {a, b} 是否第一次出现。
 *
 * 返回：
 *   1：该节点对此前未出现，并已标记；
 *   0：该节点对此前已经出现，或 a == b。
 *
 * pair_seen == NULL 表示未启用缓存，此时不同节点对总是返回 1。
 *
 * 前提：
 *   0 <= a, b < n。
 */

static inline igraph_bool_t dag_pair_mark_new(
    unsigned char *pair_seen,
    igraph_integer_t n,
    igraph_integer_t a,
    igraph_integer_t b
) {
    if (a == b) {
        return 0;
    }

    /* 未启用缓存时，每次都处理该节点对 */
    if (!pair_seen || n < 2) {
        return 1;
    }

    /* 将 (a,b) 规范化为 a < b，使其表示无序节点对 */
    if (a > b) {
        igraph_integer_t tmp = a;
        a = b;
        b = tmp;
    }

    const size_t n_sz = (size_t) n;
    const size_t a_sz = (size_t) a;
    const size_t b_sz = (size_t) b;

    /*
     * 将上三角中的无序节点对映射到：
     * 0, 1, ..., n(n-1)/2 - 1。
     *
     * 行首偏移为 a(2n-a-1)/2。先除以 2 再相乘，
     * 避免中间乘积在 size_t 中溢出。
     */
    size_t row_factor_a = a_sz;
    size_t row_factor_b = 2u * n_sz - a_sz - 1u;

    if ((row_factor_a & 1u) == 0u) {
        row_factor_a /= 2u;
    } else {
        row_factor_b /= 2u;
    }

    const size_t index =
        row_factor_a * row_factor_b
        + b_sz - a_sz - 1u;

    unsigned char *byte = &pair_seen[index >> 3];
    const unsigned char bit =
        (unsigned char) (1u << (index & 7u));

    if ((*byte & bit) != 0) {
        return 0;
    }

    *byte |= bit;
    return 1;
}


/*
 * 为 an_graph 的全部无序节点对分配一位缓存。
 *
 * 共需 C(n,2) 位，即约 n(n-1)/16 字节。
 * 本函数先约分再相乘，避免计算 n(n-1)/2 时发生 size_t 溢出。
 */
static igraph_error_t dag_pair_cache_init(
    igraph_integer_t n,
    unsigned char **pair_seen_out
) {
    if (!pair_seen_out || n < 0) {
        return IGRAPH_EINVAL;
    }

    *pair_seen_out = NULL;

    if (n < 2) {
        return IGRAPH_SUCCESS;
    }

    size_t factor_a = (size_t) n;
    size_t factor_b = (size_t) (n - 1);

    /* n 与 n-1 中必有一个为偶数，先除以 2 再相乘。 */
    if ((factor_a & 1u) == 0u) {
        factor_a /= 2u;
    } else {
        factor_b /= 2u;
    }

    if (factor_b != 0u && factor_a > SIZE_MAX / factor_b) {
        return IGRAPH_ENOMEM;
    }

    const size_t pair_count = factor_a * factor_b;
    const size_t byte_count =
        pair_count / 8u + (pair_count % 8u != 0u);

    unsigned char *pair_seen =
        (unsigned char *) calloc(byte_count, sizeof(unsigned char));

    if (!pair_seen) {
        return IGRAPH_ENOMEM;
    }

    *pair_seen_out = pair_seen;
    return IGRAPH_SUCCESS;
}



static igraph_error_t dag_cmds_cached(
    dag_moral_cache_t *cache,
    igraph_integer_t u,
    igraph_integer_t v,
    igraph_vector_int_t *out_nodes
) {
    if (!cache || !out_nodes) {
        return IGRAPH_EINVAL;
    }

    const igraph_integer_t n = cache->n;
    if (u < 0 || u >= n || v < 0 || v >= n) {
        return IGRAPH_EINVVID;
    }

    igraph_error_t ret = IGRAPH_SUCCESS;
    unsigned char *state = (unsigned char *) calloc(
        (size_t) (n > 0 ? n : 1),
        sizeof(unsigned char)
    );
    igraph_integer_t *queue = (igraph_integer_t *) malloc(
        (size_t) (n > 0 ? n : 1) * sizeof(igraph_integer_t)
    );

    if (!state || !queue) {
        free(state);
        free(queue);
        return IGRAPH_ENOMEM;
    }

    igraph_vector_int_clear(out_nodes);

    const igraph_vector_int_t *forbidden_vertices = NULL;
    ret = dag_moral_cache_get_blanket(cache, u, &forbidden_vertices);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    for (igraph_integer_t i = 0;
         i < igraph_vector_int_size(forbidden_vertices);
         ++i) {
        const igraph_integer_t node = VECTOR(*forbidden_vertices)[i];
        state[node] = 2;
    }

    if (state[v] >= 2) {
        ret = IGRAPH_EINVAL;
        goto cleanup;
    }

    igraph_integer_t head = 0;
    igraph_integer_t tail = 0;
    queue[tail++] = v;
    state[v] = 1;

    while (head < tail) {
        const igraph_integer_t cur = queue[head++];
        const igraph_vector_int_t *moral_neighbors = NULL;

        ret = dag_moral_cache_get_blanket(cache, cur, &moral_neighbors);
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup;
        }

        for (igraph_integer_t i = 0;
             i < igraph_vector_int_size(moral_neighbors);
             ++i) {
            const igraph_integer_t nb = VECTOR(*moral_neighbors)[i];

            if (state[nb] >= 2) {
                if (state[nb] == 2) {
                    ret = igraph_vector_int_push_back(out_nodes, nb);
                    if (ret != IGRAPH_SUCCESS) {
                        goto cleanup;
                    }
                    state[nb] = 3;
                }
                continue;
            }

            if (state[nb] == 0) {
                state[nb] = 1;
                queue[tail++] = nb;
            }
        }
    }

    igraph_vector_int_sort(out_nodes);

cleanup:
    free(queue);
    free(state);
    return ret;
}


static igraph_error_t dag_cmds(
    const igraph_t *uv_an_graph,
    igraph_integer_t u,
    igraph_integer_t v,
    igraph_vector_int_t *out_nodes
) {
    if (!uv_an_graph || !out_nodes) {
        return IGRAPH_EINVAL;
    }

    const igraph_integer_t n =
        igraph_vcount(uv_an_graph);

    if (u < 0 || u >= n ||
        v < 0 || v >= n) {
        return IGRAPH_EINVVID;
    }

    {
        dag_moral_cache_t cache;
        igraph_error_t ret = dag_moral_cache_init(&cache, uv_an_graph);
        if (ret != IGRAPH_SUCCESS) {
            return ret;
        }

        ret = dag_cmds_cached(&cache, u, v, out_nodes);
        dag_moral_cache_destroy(&cache);
        return ret;
    }

    igraph_error_t ret = IGRAPH_SUCCESS;

    /*
     * state：
     *   0：非禁忌，未访问
     *   1：非禁忌，已进入 BFS
     *   2：禁忌，尚未输出
     *   3：禁忌，已经输出
     */
    unsigned char *state =
        (unsigned char *) calloc(
            (size_t) n,
            sizeof(unsigned char)
        );

    igraph_integer_t *queue =
        (igraph_integer_t *) malloc(
            (size_t) n *
            sizeof(igraph_integer_t)
        );

    if (!state || !queue) {
        free(state);
        free(queue);
        return IGRAPH_ENOMEM;
    }

    igraph_vector_int_t forbidden_vertices;
    igraph_vector_int_t moral_neighbors;

    ret = igraph_vector_int_init(
        &forbidden_vertices,
        0
    );
    if (ret != IGRAPH_SUCCESS) {
        free(state);
        free(queue);
        return ret;
    }

    ret = igraph_vector_int_init(
        &moral_neighbors,
        0
    );
    if (ret != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(
            &forbidden_vertices
        );
        free(state);
        free(queue);
        return ret;
    }

    igraph_vector_int_clear(out_nodes);

    /*
     * Step 1：
     * mb(u) 作为禁忌节点集合。
     */
    ret = get_markov_blanket(
        uv_an_graph,
        u,
        &forbidden_vertices
    );
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    for (igraph_integer_t i = 0;
         i < igraph_vector_int_size(&forbidden_vertices);
         ++i) {

        const igraph_integer_t node =
            VECTOR(forbidden_vertices)[i];

        state[node] = 2;
    }

    /*
     * u、v 非邻接时，v 不应属于 mb(u)。
     */
    if (state[v] >= 2) {
        ret = IGRAPH_EINVAL;
        goto cleanup;
    }

    /*
     * Step 2：
     * 从 v 开始在隐式道义图上进行 BFS。
     */
    igraph_integer_t head = 0;
    igraph_integer_t tail = 0;

    queue[tail++] = v;
    state[v] = 1;

    while (head < tail) {
        const igraph_integer_t cur =
            queue[head++];

        /*
         * cur 的 Markov blanket
         * 就是其道义邻居集合。
         */
        ret = get_markov_blanket(
            uv_an_graph,
            cur,
            &moral_neighbors
        );
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup;
        }

        for (igraph_integer_t i = 0;
             i < igraph_vector_int_size(&moral_neighbors);
             ++i) {

            const igraph_integer_t nb =
                VECTOR(moral_neighbors)[i];

            /*
             * 遇到禁忌节点：
             * 加入边界输出，但不继续穿过。
             */
            if (state[nb] >= 2) {
                if (state[nb] == 2) {
                    ret =
                        igraph_vector_int_push_back(
                            out_nodes,
                            nb
                        );

                    if (ret != IGRAPH_SUCCESS) {
                        goto cleanup;
                    }

                    state[nb] = 3;
                }

                continue;
            }

            /*
             * 非禁忌节点加入当前分量。
             */
            if (state[nb] == 0) {
                state[nb] = 1;
                queue[tail++] = nb;
            }
        }
    }

    igraph_vector_int_sort(out_nodes);

cleanup:
    igraph_vector_int_destroy(
        &forbidden_vertices
    );
    igraph_vector_int_destroy(
        &moral_neighbors
    );

    free(state);
    free(queue);

    return ret;
}


static void dag_clear_int_vector_ptr(igraph_vector_ptr_t *vectors) {
    if (!vectors) {
        return;
    }

    for (igraph_integer_t i = 0;
         i < igraph_vector_ptr_size(vectors);
         ++i) {

        igraph_vector_int_t *vec =
            (igraph_vector_int_t *) VECTOR(*vectors)[i];

        if (vec) {
            igraph_vector_int_destroy(vec);
            free(vec);
        }
    }

    igraph_vector_ptr_clear(vectors);
}


/*
 * 直接使用 igraph_subcomponent(..., IGRAPH_IN) 收集 seeds 的全部祖先。
 *
 * 本实现完全基于节点向量，不再构造 seed_mask、allowed_mask、out_mask
 * 或额外的 DAG workspace。
 *
 * ancestors_out 与 tmp_ancestors 必须已经初始化。
 */
static igraph_error_t dag_collect_ancestor_nodes(
    const igraph_t *graph,
    const igraph_vector_int_t *seeds,
    igraph_vector_int_t *ancestors_out,
    igraph_vector_int_t *tmp_ancestors
) {
    if (!graph || !seeds || !ancestors_out || !tmp_ancestors) {
        return IGRAPH_EINVAL;
    }

    const igraph_integer_t n = igraph_vcount(graph);
    igraph_error_t ret = IGRAPH_SUCCESS;

    igraph_vector_int_clear(ancestors_out);

    for (igraph_integer_t i = 0;
         i < igraph_vector_int_size(seeds);
         ++i) {

        const igraph_integer_t seed = VECTOR(*seeds)[i];

        if (seed < 0 || seed >= n) {
            return IGRAPH_EINVVID;
        }

        igraph_vector_int_clear(tmp_ancestors);

        ret = igraph_subcomponent(
            graph,
            tmp_ancestors,
            seed,
            IGRAPH_IN
        );
        if (ret != IGRAPH_SUCCESS) {
            return ret;
        }

        ret = igraph_vector_int_append(
            ancestors_out,
            tmp_ancestors
        );
        if (ret != IGRAPH_SUCCESS) {
            return ret;
        }
    }

    vector_int_unique(ancestors_out);

    return IGRAPH_SUCCESS;
}


igraph_error_t dag_get_ancestors(
    const igraph_t *graph,
    const igraph_vector_int_t *r_nodes,
    igraph_vector_int_t *ancestors_out
) {
    if (!graph || !r_nodes || !ancestors_out) {
        return IGRAPH_EINVAL;
    }

    igraph_vector_int_t tmp_ancestors;
    igraph_error_t ret =
        igraph_vector_int_init(&tmp_ancestors, 0);

    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }

    ret = dag_collect_ancestor_nodes(
        graph,
        r_nodes,
        ancestors_out,
        &tmp_ancestors
    );

    igraph_vector_int_destroy(&tmp_ancestors);
    return ret;
}


/*
 * 对外层祖先图 an_graph 中的一对非邻接节点 aa、bb 执行一次 CMDS 扩展。
 *
 * 编号层次：
 *   原图编号
 *       <- outer invmap ->
 *   an_graph 编号
 *       <- uv_invmap ->
 *   uv_an_graph 编号
 *
 * 本函数只负责第二层映射：把 dag_cmds 返回的 uv_an_graph 局部编号
 * 通过 uv_invmap 映射回 an_graph 编号，并加入 h_nodes。
 */
static igraph_error_t dag_expand_nonadjacent_pair(
    const igraph_t *an_graph,
    dag_moral_cache_t *an_cache,
    igraph_integer_t aa,
    igraph_integer_t bb,
    igraph_vector_int_t *h_nodes,
    igraph_vector_int_t *uv_seeds,
    igraph_vector_int_t *uv_anc_nodes,
    igraph_vector_int_t *tmp_ancestors,
    igraph_vector_int_t *s_a,
    igraph_vector_int_t *s_b,
    igraph_vector_int_t *uv_map,
    igraph_vector_int_t *uv_invmap,
    igraph_integer_t *ancestor_stack,
    unsigned char *ancestor_seen,
    igraph_bool_t *added_new_node
) {
    if (!an_graph || !an_cache || !h_nodes || !uv_seeds ||
        !uv_anc_nodes || !tmp_ancestors ||
        !s_a || !s_b || !uv_map || !uv_invmap ||
        !ancestor_stack || !ancestor_seen || !added_new_node) {
        return IGRAPH_EINVAL;
    }

    (void) tmp_ancestors;

    const igraph_integer_t an = igraph_vcount(an_graph);
    if (aa < 0 || aa >= an || bb < 0 || bb >= an || aa == bb) {
        return IGRAPH_EINVVID;
    }

    igraph_error_t ret = IGRAPH_SUCCESS;
    igraph_t uv_an_graph;
    igraph_bool_t uv_an_graph_initialized = 0;
    dag_moral_cache_t uv_cache;
    igraph_bool_t uv_cache_initialized = 0;
    igraph_vs_t vs;
    igraph_bool_t vs_initialized = 0;

    const igraph_integer_t old_h_size =
        igraph_vector_int_size(h_nodes);

    *added_new_node = 0;

    /* 1. 构造种子 {aa, bb}。 */
    igraph_vector_int_clear(uv_seeds);

    ret = igraph_vector_int_push_back(uv_seeds, aa);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    ret = igraph_vector_int_push_back(uv_seeds, bb);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    /* 2. 在外层祖先图中收集 An({aa, bb})。 */
    ret = dag_collect_ancestor_nodes_cached(
        an_cache,
        uv_seeds,
        uv_anc_nodes,
        ancestor_stack,
        ancestor_seen
    );
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    /* 3. 构造 uv 祖先诱导子图，并保留双向编号映射。 */
    igraph_vector_int_clear(uv_map);
    igraph_vector_int_clear(uv_invmap);

    ret = igraph_vs_vector(&vs, uv_anc_nodes);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    vs_initialized = 1;

    ret = igraph_induced_subgraph_map(
        an_graph,
        &uv_an_graph,
        vs,
        IGRAPH_SUBGRAPH_AUTO,
        uv_map,
        uv_invmap
    );
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    uv_an_graph_initialized = 1;

    igraph_vs_destroy(&vs);
    vs_initialized = 0;

    ret = dag_moral_cache_init(&uv_cache, &uv_an_graph);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    uv_cache_initialized = 1;

    /*
     * map[old] == new + 1；若 old 不在子图中则为 0。
     * aa、bb 是祖先集合的种子，因此理论上一定在 uv_an_graph 中。
     */
    const igraph_integer_t a_uv = VECTOR(*uv_map)[aa] - 1;
    const igraph_integer_t b_uv = VECTOR(*uv_map)[bb] - 1;
    const igraph_integer_t uv_n = igraph_vcount(&uv_an_graph);

    if (a_uv < 0 || a_uv >= uv_n ||
        b_uv < 0 || b_uv >= uv_n) {
        ret = IGRAPH_EINVAL;
        goto cleanup;
    }

    /* 4. 在同一个 uv 祖先图上分别计算 S_a 和 S_b。 */
    igraph_vector_int_clear(s_a);
    igraph_vector_int_clear(s_b);

    ret = dag_cmds_cached(
        &uv_cache,
        a_uv,
        b_uv,
        s_a
    );
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    ret = dag_cmds_cached(
        &uv_cache,
        b_uv,
        a_uv,
        s_b
    );
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    /* 5. uv 局部编号 -> an_graph 编号，并加入当前 H。 */
    for (igraph_integer_t i = 0;
         i < igraph_vector_int_size(s_a);
         ++i) {

        const igraph_integer_t uv_node = VECTOR(*s_a)[i];
        if (uv_node < 0 || uv_node >= uv_n) {
            ret = IGRAPH_EINVVID;
            goto cleanup;
        }

        ret = igraph_vector_int_push_back(
            h_nodes,
            VECTOR(*uv_invmap)[uv_node]
        );
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup;
        }
    }

    for (igraph_integer_t i = 0;
         i < igraph_vector_int_size(s_b);
         ++i) {

        const igraph_integer_t uv_node = VECTOR(*s_b)[i];
        if (uv_node < 0 || uv_node >= uv_n) {
            ret = IGRAPH_EINVVID;
            goto cleanup;
        }

        ret = igraph_vector_int_push_back(
            h_nodes,
            VECTOR(*uv_invmap)[uv_node]
        );
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup;
        }
    }

    /*
     * 不再使用 in_h mask。
     * 通过“追加 -> 排序去重 -> 比较长度”判断 H 是否真正扩大。
     */
    vector_int_unique(h_nodes);

    *added_new_node =
        igraph_vector_int_size(h_nodes) > old_h_size;

cleanup:
    if (vs_initialized) {
        igraph_vs_destroy(&vs);
    }
    if (uv_cache_initialized) {
        dag_moral_cache_destroy(&uv_cache);
    }
    if (uv_an_graph_initialized) {
        igraph_destroy(&uv_an_graph);
    }

    return ret;
}


igraph_error_t dag_get_minimal_collapsible(
    const igraph_t *graph,
    const igraph_vector_int_t *r_nodes,
    igraph_vector_int_t *H_out
) {
    if (!graph || !r_nodes || !H_out) {
        return IGRAPH_EINVAL;
    }

    if (!igraph_is_directed(graph)) {
        return IGRAPH_EINVAL;
    }

    igraph_t ug_graph;
    igraph_bool_t ug_graph_initialized = 0;
    igraph_vector_int_t preprocessed_r;
    igraph_bool_t preprocessed_r_initialized = 0;
    igraph_vector_int_t preprocessed_H;
    igraph_bool_t preprocessed_H_initialized = 0;
    igraph_vector_int_t global_to_local;
    igraph_vector_int_t local_to_global;
    igraph_vector_int_t r_local;
    const igraph_vector_int_t *r_nodes_work = r_nodes;

    igraph_error_t ret = igraph_vector_int_init(&preprocessed_r, 0);
    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }
    preprocessed_r_initialized = 1;

    igraph_vector_int_init(&preprocessed_H, 0);
    igraph_vector_int_init(&global_to_local, 0);
    igraph_vector_int_init(&local_to_global, 0);
    igraph_vector_int_init(&r_local, 0);
    preprocessed_H_initialized = 1;

    if (igraph_vector_int_size(r_nodes) > 0) {
        ret = dag_ancestor_moral_subgraph_with_mapping(
            graph,
            r_nodes,
            &ug_graph,
            &global_to_local,
            &local_to_global,
            &r_local
        );
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup_preprocess;
        }
        ug_graph_initialized = 1;

        ret = get_minimal_collapsible(&ug_graph, &r_local, &preprocessed_H);
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup_preprocess;
        }

        igraph_vector_int_clear(&preprocessed_r);
        for (igraph_integer_t i = 0; i < igraph_vector_int_size(&preprocessed_H); ++i) {
            const igraph_integer_t local_node = VECTOR(preprocessed_H)[i];
            if (local_node < 0 || local_node >= igraph_vector_int_size(&local_to_global)) {
                ret = IGRAPH_EINVVID;
                goto cleanup_preprocess;
            }
            ret = igraph_vector_int_push_back(&preprocessed_r, VECTOR(local_to_global)[local_node]);
            if (ret != IGRAPH_SUCCESS) {
                goto cleanup_preprocess;
            }
        }

        vector_int_unique(&preprocessed_r);
        r_nodes_work = &preprocessed_r;
    }

cleanup_preprocess:
    if (preprocessed_H_initialized) {
        igraph_vector_int_destroy(&preprocessed_H);
    }
    if (ret != IGRAPH_SUCCESS) {
        if (preprocessed_r_initialized) {
            igraph_vector_int_destroy(&preprocessed_r);
        }
        if (ug_graph_initialized) {
            igraph_destroy(&ug_graph);
        }
        igraph_vector_int_destroy(&global_to_local);
        igraph_vector_int_destroy(&local_to_global);
        igraph_vector_int_destroy(&r_local);
        return ret;
    }

    igraph_bool_t is_dag = 0;
    ret = igraph_is_dag(graph, &is_dag);
    if (ret != IGRAPH_SUCCESS) {
        return ret;
    }
    if (!is_dag) {
        return IGRAPH_EINVAL;
    }

    ret = IGRAPH_SUCCESS;
    const igraph_integer_t n = igraph_vcount(graph);
    unsigned char *pair_seen = NULL;
    unsigned char *ancestor_seen = NULL;
    igraph_integer_t *ancestor_stack = NULL;

    igraph_vector_int_t anc_nodes;
    igraph_vector_int_t tmp_ancestors;
    igraph_vector_int_t h_nodes;
    igraph_vector_int_t outer_map;
    igraph_vector_int_t outer_invmap;

    igraph_vector_int_t uv_seeds;
    igraph_vector_int_t uv_anc_nodes;
    igraph_vector_int_t s_a;
    igraph_vector_int_t s_b;
    igraph_vector_int_t uv_map;
    igraph_vector_int_t uv_invmap;

    igraph_vector_ptr_t boundaries;

    igraph_bool_t anc_nodes_initialized = 0;
    igraph_bool_t tmp_ancestors_initialized = 0;
    igraph_bool_t h_nodes_initialized = 0;
    igraph_bool_t outer_map_initialized = 0;
    igraph_bool_t outer_invmap_initialized = 0;
    igraph_bool_t uv_seeds_initialized = 0;
    igraph_bool_t uv_anc_nodes_initialized = 0;
    igraph_bool_t s_a_initialized = 0;
    igraph_bool_t s_b_initialized = 0;
    igraph_bool_t uv_map_initialized = 0;
    igraph_bool_t uv_invmap_initialized = 0;
    igraph_bool_t boundaries_initialized = 0;

    igraph_t an_graph;
    igraph_bool_t an_graph_initialized = 0;
    dag_moral_cache_t graph_cache;
    igraph_bool_t graph_cache_initialized = 0;
    dag_moral_cache_t an_cache;
    igraph_bool_t an_cache_initialized = 0;
    igraph_vs_t vs;
    igraph_bool_t vs_initialized = 0;

    ret = igraph_vector_int_init(&anc_nodes, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    anc_nodes_initialized = 1;

    ret = igraph_vector_int_init(&tmp_ancestors, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    tmp_ancestors_initialized = 1;

    ret = igraph_vector_int_init(&h_nodes, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    h_nodes_initialized = 1;

    ret = igraph_vector_int_init(&outer_map, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    outer_map_initialized = 1;

    ret = igraph_vector_int_init(&outer_invmap, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    outer_invmap_initialized = 1;

    ret = igraph_vector_int_init(&uv_seeds, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    uv_seeds_initialized = 1;

    ret = igraph_vector_int_init(&uv_anc_nodes, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    uv_anc_nodes_initialized = 1;

    ret = igraph_vector_int_init(&s_a, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    s_a_initialized = 1;

    ret = igraph_vector_int_init(&s_b, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    s_b_initialized = 1;

    ret = igraph_vector_int_init(&uv_map, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    uv_map_initialized = 1;

    ret = igraph_vector_int_init(&uv_invmap, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    uv_invmap_initialized = 1;

    ret = igraph_vector_ptr_init(&boundaries, 0);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    boundaries_initialized = 1;

    /* 空初始集的闭包仍为空。 */
    if (igraph_vector_int_size(r_nodes_work) == 0) {
        igraph_vector_int_clear(H_out);
        goto cleanup;
    }

    ancestor_stack = (igraph_integer_t *) malloc(
        (size_t) (n > 0 ? n : 1) * sizeof(igraph_integer_t)
    );
    ancestor_seen = (unsigned char *) calloc(
        (size_t) (n > 0 ? n : 1),
        sizeof(unsigned char)
    );
    if (!ancestor_stack || !ancestor_seen) {
        ret = IGRAPH_ENOMEM;
        goto cleanup;
    }

    ret = dag_moral_cache_init(&graph_cache, graph);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    graph_cache_initialized = 1;

    /*
     * 第一层：构造原图中的 An(R)。
     * 后续 h_nodes 与 boundaries 全部使用 an_graph 的局部编号。
     */
    ret = dag_collect_ancestor_nodes_cached(
        &graph_cache,
        r_nodes_work,
        &anc_nodes,
        ancestor_stack,
        ancestor_seen
    );
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    ret = igraph_vs_vector(&vs, &anc_nodes);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    vs_initialized = 1;

    ret = igraph_induced_subgraph_map(
        graph,
        &an_graph,
        vs,
        IGRAPH_SUBGRAPH_AUTO,
        &outer_map,
        &outer_invmap
    );
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    an_graph_initialized = 1;

    igraph_vs_destroy(&vs);
    vs_initialized = 0;

    const igraph_integer_t an = igraph_vcount(&an_graph);

    ret = dag_moral_cache_init(&an_cache, &an_graph);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    an_cache_initialized = 1;

    /*
     * pair_seen 的编号空间固定为 an_graph 的局部编号。
     * 缓存在整个闭包扩张过程中保持，不随迭代轮次清空。
     */
    ret = dag_pair_cache_init(an, &pair_seen);
    if (ret != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    /* 预留常用向量容量，避免循环内反复扩容。 */
    ret = igraph_vector_int_reserve(&h_nodes, an);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    ret = igraph_vector_int_reserve(&uv_anc_nodes, an);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    ret = igraph_vector_int_reserve(&tmp_ancestors, an);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    ret = igraph_vector_int_reserve(&s_a, an);
    if (ret != IGRAPH_SUCCESS) goto cleanup;
    ret = igraph_vector_int_reserve(&s_b, an);
    if (ret != IGRAPH_SUCCESS) goto cleanup;

    /* 原图 R 编号 -> an_graph 局部编号。 */
    for (igraph_integer_t i = 0;
         i < igraph_vector_int_size(r_nodes_work);
         ++i) {

        const igraph_integer_t original_node =
            VECTOR(*r_nodes_work)[i];

        if (original_node < 0 || original_node >= n) {
            ret = IGRAPH_EINVVID;
            goto cleanup;
        }

        const igraph_integer_t mapped =
            VECTOR(outer_map)[original_node];

        if (mapped <= 0) {
            /* R 中的节点必须属于 An(R)。 */
            ret = IGRAPH_EINVAL;
            goto cleanup;
        }

        ret = igraph_vector_int_push_back(
            &h_nodes,
            mapped - 1
        );
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup;
        }
    }

    vector_int_unique(&h_nodes);
    
    

    for (;;) {
        igraph_bool_t changed = 0;

        /* 清理并复用上一轮的边界容器。 */
        dag_clear_int_vector_ptr(&boundaries);

        ret = dag_components_forbidden_cached(
            &an_cache,
            NULL,
            &boundaries,
            &h_nodes
        );
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup;
        }

        const igraph_integer_t group_count =
            igraph_vector_ptr_size(&boundaries);

        for (igraph_integer_t cidx = 0;
             cidx < group_count && !changed;
             ++cidx) {

            igraph_vector_int_t *boundary =
                (igraph_vector_int_t *) VECTOR(boundaries)[cidx];

            if (!boundary ||
                igraph_vector_int_size(boundary) < 2) {
                continue;
            }

            igraph_vector_int_sort(boundary);
            const igraph_integer_t bsize =
                igraph_vector_int_size(boundary);

            for (igraph_integer_t i = 0;
                 i < bsize && !changed;
                 ++i) {

                const igraph_integer_t aa =
                    VECTOR(*boundary)[i];

                for (igraph_integer_t j = i + 1;
                     j < bsize && !changed;
                     ++j) {

                    const igraph_integer_t bb =
                        VECTOR(*boundary)[j];

                    /*
                     * 同一无序节点对可能同时出现在多个分量边界中，
                     * 也可能在 H 扩张后的后续轮次再次出现。
                     *
                     * An({aa,bb}) 与 dag_cmds 的结果只依赖固定的
                     * an_graph 和无序节点对 {aa,bb}，因此全程只需处理一次。
                     * 把标记放在邻接检查之前，也可避免重复检查直接邻接对。
                     */
                    if (!dag_pair_mark_new(
                            pair_seen,
                            an,
                            aa,
                            bb)) {
                        continue;
                    }

                    /*
                     * 这里只排除 aa -> bb 或 bb -> aa 的直接有向边。
                     * igraph_are_adjacent() 在有向图中按方向查询，
                     * 因此必须检查两个方向。
                     *
                     * 不能用外层 an_graph 的 Markov blanket 判断：
                     * aa、bb 可能因共同子节点而在 An(R)^m 中相邻，
                     * 但该共同子节点通常不属于 An({aa,bb})，因此不应跳过。
                     */
                    igraph_bool_t adjacent_ab = 0;
                    igraph_bool_t adjacent_ba = 0;

                    ret = igraph_are_adjacent(
                        &an_graph,
                        aa,
                        bb,
                        &adjacent_ab
                    );
                    if (ret != IGRAPH_SUCCESS) {
                        goto cleanup;
                    }

                    if (!adjacent_ab) {
                        ret = igraph_are_adjacent(
                            &an_graph,
                            bb,
                            aa,
                            &adjacent_ba
                        );
                        if (ret != IGRAPH_SUCCESS) {
                            goto cleanup;
                        }
                    }

                    if (adjacent_ab || adjacent_ba) {
                        continue;
                    }

                    /*
                     * 第二层：构造 An({aa,bb})，调用新版 dag_cmds，
                     * 再通过 uv_invmap 把结果映射回 an_graph 后加入 H。
                     */
                    ret = dag_expand_nonadjacent_pair(
                        &an_graph,
                        &an_cache,
                        aa,
                        bb,
                        &h_nodes,
                        &uv_seeds,
                        &uv_anc_nodes,
                        &tmp_ancestors,
                        &s_a,
                        &s_b,
                        &uv_map,
                        &uv_invmap,
                        ancestor_stack,
                        ancestor_seen,
                        &changed
                    );
                    if (ret != IGRAPH_SUCCESS) {
                        goto cleanup;
                    }
                }
            }
        }

        if (!changed) {
            break;
        }
    }

    /* 最后一层：an_graph 局部编号 -> 原图编号。 */
    igraph_vector_int_sort(&h_nodes);
    igraph_vector_int_clear(H_out);

    for (igraph_integer_t i = 0;
         i < igraph_vector_int_size(&h_nodes);
         ++i) {

        const igraph_integer_t local_node =
            VECTOR(h_nodes)[i];

        if (local_node < 0 || local_node >= an) {
            ret = IGRAPH_EINVVID;
            goto cleanup;
        }

        ret = igraph_vector_int_push_back(
            H_out,
            VECTOR(outer_invmap)[local_node]
        );
        if (ret != IGRAPH_SUCCESS) {
            goto cleanup;
        }
    }

cleanup:
    if (vs_initialized) {
        igraph_vs_destroy(&vs);
    }

    if (boundaries_initialized) {
        dag_clear_int_vector_ptr(&boundaries);
        igraph_vector_ptr_destroy(&boundaries);
    }

    if (an_cache_initialized) {
        dag_moral_cache_destroy(&an_cache);
    }

    if (an_graph_initialized) {
        igraph_destroy(&an_graph);
    }

    if (graph_cache_initialized) {
        dag_moral_cache_destroy(&graph_cache);
    }

    free(pair_seen);
    free(ancestor_seen);
    free(ancestor_stack);

    if (uv_invmap_initialized) igraph_vector_int_destroy(&uv_invmap);
    if (uv_map_initialized) igraph_vector_int_destroy(&uv_map);
    if (s_b_initialized) igraph_vector_int_destroy(&s_b);
    if (s_a_initialized) igraph_vector_int_destroy(&s_a);
    if (uv_anc_nodes_initialized) igraph_vector_int_destroy(&uv_anc_nodes);
    if (uv_seeds_initialized) igraph_vector_int_destroy(&uv_seeds);
    if (outer_invmap_initialized) igraph_vector_int_destroy(&outer_invmap);
    if (outer_map_initialized) igraph_vector_int_destroy(&outer_map);
    if (h_nodes_initialized) igraph_vector_int_destroy(&h_nodes);
    if (tmp_ancestors_initialized) igraph_vector_int_destroy(&tmp_ancestors);
    if (anc_nodes_initialized) igraph_vector_int_destroy(&anc_nodes);
    if (preprocessed_r_initialized) {
        igraph_vector_int_destroy(&preprocessed_r);
    }
    if (ug_graph_initialized) {
        igraph_destroy(&ug_graph);
    }

    return ret;
}
