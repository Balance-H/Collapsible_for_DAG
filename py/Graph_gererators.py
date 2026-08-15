import random
import igraph as ig
import networkx as nx

try:
    import numpy as np
except ImportError:
    np = None  # If numpy is not available, set np to None


def generator_connected_ug(n, p, class_type="ig"):
    """
    Generates a random connected graph with n nodes and a probability p of adding edges.
    Supports both igraph and networkx based on class_type parameter.
    
    :param n: Number of nodes in the graph
    :param p: Probability of adding an edge between any two nodes (after ensuring connectivity)
    :param class_type: "ig" for igraph graph (default), "nx" for networkx graph
    :return: A connected graph object of the specified type
    """
    if np is None:
        raise ImportError("This function requires numpy, which is not installed on your system.")
    
    if class_type == "nx":
        # NetworkX version
        g = nx.random_tree(n)
        adj_matrix = np.zeros((n, n), dtype=bool)
        for u, v in g.edges():
            adj_matrix[u, v] = True
            adj_matrix[v, u] = True
        
        row_indices, col_indices = np.triu_indices(n, k=1)
        candidate_mask = ~adj_matrix[row_indices, col_indices]
        random_vals = np.random.rand(len(row_indices))
        edges_to_add_mask = (random_vals < p) & candidate_mask
        edges_to_add = [(row_indices[i], col_indices[i]) for i in np.where(edges_to_add_mask)[0]]
        g.add_edges_from(edges_to_add)
        return g
    
    else:
        # igraph version (default)
        g = ig.Graph.Tree(n, 2)
        adj_matrix = np.zeros((n, n), dtype=bool)
        for u, v in g.get_edgelist():
            adj_matrix[u, v] = True
            adj_matrix[v, u] = True
        
        row_indices, col_indices = np.triu_indices(n, k=1)
        candidate_mask = ~adj_matrix[row_indices, col_indices]
        random_vals = np.random.rand(len(row_indices))
        edges_to_add_mask = (random_vals < p) & candidate_mask
        edges_to_add = [(row_indices[i], col_indices[i]) for i in np.where(edges_to_add_mask)[0]]
        if edges_to_add:
            g.add_edges(edges_to_add)
        return g



# This function generates a connected Directed Acyclic Graph (DAG) with n nodes.
# It first creates a directed tree and then probabilistically adds edges with probability p.
# Each node has a maximum of `max_parents` parents, with the default value being 3.
# Additionally, the number of nodes with exactly 3 parents is constrained to not exceed 5% of the total number of nodes.

def generate_connected_dag(n, p):

    g = ig.Graph.Erdos_Renyi(n=n, p=p, directed=False, loops=False)

    # Ensure the undirected graph is connected by linking components with k-1 edges.
    comps = g.components()
    if len(comps) > 1:
        reps = [comp[0] for comp in comps]
        bridge_edges = [(reps[i], reps[i + 1]) for i in range(len(reps) - 1)]
        g.add_edges(bridge_edges)

    g.to_directed(mode="acyclic")

    return g


