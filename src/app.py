#!/usr/bin/env python3
"""
app.py – Graph Isomorphism and Automorphism Counting

Usage:
    python main.py <file_or_folder> [GI|Aut|GIAut]

Detects task from filename suffix (GI, Aut, GIAut). Defaults to GIAut.
Handles .gr (single graph) and .grl (graph list) files.
"""

import sys
import os
import time
from .graph import Graph
from .graph_io import load_graph

sys.setrecursionlimit(10000)


# ============================================================================
# Convert Graph object to flat adjacency list for fast processing
# ============================================================================

def graph_to_adjacency_list(graph):
    """Convert a Graph object to adjacency list (list of list of int)."""
    vertices = graph.vertices
    vertex_to_index = {vertex: index for index, vertex in enumerate(vertices)}
    neighbours = [[] for _ in range(len(vertices))]
    for vertex in vertices:
        index = vertex_to_index[vertex]
        for neighbour in vertex.neighbours:
            neighbours[index].append(vertex_to_index[neighbour])
    return neighbours


# ============================================================================
# Colour refinement
# ============================================================================

def color_refine(neighbours, colors):
    """
    Run colour refinement until stable colouring is reached.
    neighbours: adjacency list (list of list of int)
    colors: initial colouring (list of int)
    Returns the stable colouring as a list.
    """
    num_vertices = len(neighbours)
    current_colors = list(colors)

    while True:
        # Build a signature for each vertex: (own color, sorted neighbor colors)
        signatures = [None] * num_vertices
        for vertex in range(num_vertices):
            neighbor_list = neighbours[vertex]
            if neighbor_list:
                neighbor_colors = tuple(sorted(current_colors[neighbour] for neighbour in neighbor_list))
                signatures[vertex] = (current_colors[vertex], neighbor_colors)
            else:
                signatures[vertex] = (current_colors[vertex], ())

        # Sort vertices by signature and assign new color numbers
        sorted_vertices = sorted(range(num_vertices), key=lambda v: signatures[v])
        new_colors = [0] * num_vertices
        color_id = 0
        new_colors[sorted_vertices[0]] = 0
        for i in range(1, num_vertices):
            if signatures[sorted_vertices[i]] != signatures[sorted_vertices[i - 1]]:
                color_id += 1
            new_colors[sorted_vertices[i]] = color_id

        # If nothing changed, the colouring is stable
        if new_colors == current_colors:
            return current_colors
        current_colors = new_colors


# ============================================================================
# Merge two graphs into one adjacency list
# ============================================================================

def merge_graphs(graph_g, graph_h):
    """
    Merge two adjacency lists into one combined graph.
    Vertices 0..n_g-1 belong to G, vertices n_g..n_g+n_h-1 belong to H.
    Returns (merged_neighbours, num_vertices_g).
    """
    num_vertices_g = len(graph_g)
    merged_neighbours = [list(neighbor_list) for neighbor_list in graph_g]
    for vertex in range(len(graph_h)):
        # Shift H's vertex indices by num_vertices_g
        merged_neighbours.append([num_vertices_g + neighbour for neighbour in graph_h[vertex]])
    return merged_neighbours, num_vertices_g


# ============================================================================
# Branching (individualisation-refinement) — Algorithm 2 from the lectures
# ============================================================================

def count_isomorphisms(neighbours, num_vertices_g, colors, only_one=False):
    """
    Count the number of isomorphisms between G (vertices 0..num_vertices_g-1)
    and H (vertices num_vertices_g..total-1) in the merged graph.

    neighbours: merged adjacency list
    num_vertices_g: number of vertices in G
    colors: initial colouring
    only_one: if True, return as soon as we find at least 1 (for GI decision)
    """
    total_vertices = len(neighbours)

    # Step 1: Run colour refinement
    stable_colors = color_refine(neighbours, colors)

    # Step 2: Count colors in G and H separately to check balance
    colors_in_g = {}
    colors_in_h = {}
    for vertex in range(num_vertices_g):
        color = stable_colors[vertex]
        colors_in_g[color] = colors_in_g.get(color, 0) + 1
    for vertex in range(num_vertices_g, total_vertices):
        color = stable_colors[vertex]
        colors_in_h[color] = colors_in_h.get(color, 0) + 1

    # If unbalanced: no isomorphism possible
    if colors_in_g != colors_in_h:
        return 0

    # If every color class has size 1: we found a bijection (one isomorphism)
    if all(count == 1 for count in colors_in_g.values()):
        return 1

    # Step 3: Choose the smallest non-trivial color class (size >= 2 in G)
    best_color = None
    best_size = None
    for color, count in colors_in_g.items():
        if count >= 2 and (best_size is None or count < best_size):
            best_size = count
            best_color = color

    # Pick the first vertex x in G with that color
    vertex_x = None
    for vertex in range(num_vertices_g):
        if stable_colors[vertex] == best_color:
            vertex_x = vertex
            break

    # Collect all vertices y in H with that same color
    candidates_y = [vertex for vertex in range(num_vertices_g, total_vertices)
                    if stable_colors[vertex] == best_color]

    # Step 4: Branch — try mapping x to each candidate y
    new_unique_color = max(stable_colors) + 1
    num_isomorphisms = 0

    for vertex_y in candidates_y:
        # Give x and y the same new unique color (individualisation)
        branched_colors = list(stable_colors)
        branched_colors[vertex_x] = new_unique_color
        branched_colors[vertex_y] = new_unique_color

        # Recurse with the new colouring
        num_isomorphisms += count_isomorphisms(neighbours, num_vertices_g,
                                               branched_colors, only_one)
        # Early exit for GI decision: we only need to find one
        if only_one and num_isomorphisms > 0:
            return num_isomorphisms

    return num_isomorphisms


# ============================================================================
# Public solving functions
# ============================================================================

def are_isomorphic(graph_g, graph_h):
    """Check if two graphs are isomorphic (True/False)."""
    if len(graph_g) != len(graph_h):
        return False
    if len(graph_g) == 0:
        return True
    merged_neighbours, num_vertices_g = merge_graphs(graph_g, graph_h)
    initial_colors = [0] * len(merged_neighbours)
    return count_isomorphisms(merged_neighbours, num_vertices_g,
                              initial_colors, only_one=True) > 0


def count_automorphisms(graph):
    """Count |Aut(G)| by counting isomorphisms from G to itself."""
    if len(graph) == 0:
        return 1
    merged_neighbours, num_vertices_g = merge_graphs(graph, graph)
    initial_colors = [0] * len(merged_neighbours)
    return count_isomorphisms(merged_neighbours, num_vertices_g,
                              initial_colors, only_one=False)


def solve_gi(adjacency_lists):
    """
    Partition graphs into isomorphism classes.
    Returns list of groups (each group is a sorted list of graph indices).
    """
    num_graphs = len(adjacency_lists)
    visited = [False] * num_graphs
    groups = []
    for i in range(num_graphs):
        if visited[i]:
            continue
        group = [i]
        visited[i] = True
        for j in range(i + 1, num_graphs):
            if visited[j]:
                continue
            if are_isomorphic(adjacency_lists[i], adjacency_lists[j]):
                group.append(j)
                visited[j] = True
        groups.append(group)
    return groups


def solve_aut(adjacency_lists):
    """Compute |Aut(G)| for each graph. Returns list of (index, count)."""
    return [(i, count_automorphisms(graph))
            for i, graph in enumerate(adjacency_lists)]


def solve_gi_aut(adjacency_lists):
    """
    Solve both GI and #Aut together.
    Only computes automorphisms once per isomorphism class.
    Returns list of (group_indices, aut_count).
    """
    groups = solve_gi(adjacency_lists)
    result = []
    for group in groups:
        aut_count = count_automorphisms(adjacency_lists[group[0]])
        result.append((group, aut_count))
    return result


# ============================================================================
# File handling & main
# ============================================================================

def detect_task(filename):
    """Detect the task (GI, Aut, GIAut) from the filename."""
    base = os.path.splitext(filename)[0]
    if base.endswith("GIAut"):
        return "GIAut"
    if base.endswith("Aut"):
        return "Aut"
    if base.endswith("GI"):
        return "GI"
    return "GIAut"


def solve_file(filepath, task_override=None):
    """Load and solve a single .gr or .grl file."""
    filename = os.path.basename(filepath)
    extension = os.path.splitext(filepath)[1].lower()
    task = task_override if task_override else detect_task(filename)

    # Load graphs from file on variable graphs
    if extension == ".gr":
        with open(filepath, 'r') as file:
            graph = load_graph(file)
        graphs = [graph]
    elif extension == ".grl":
        with open(filepath, 'r') as file:
            graphs = load_graph(file, read_list=True)
    else:
        print(f"Skipping unknown file type: {filepath}")
        return

    # Convert Graph objects to adjacency lists for fast processing
    adjacency_lists = [graph_to_adjacency_list(graph) for graph in graphs]

    print(f"{filename}")
    start_time = time.time()

    # Route to solving method
    if task == "GI":
        groups = solve_gi(adjacency_lists)
        print("Sets of isomorphic graphs:")
        for group in groups:
            print(group)
    elif task == "Aut":
        automorphisms = solve_aut(adjacency_lists)
        print("Graphs with automorphisms:")
        for index, aut_count in automorphisms:
            print(f"{index}: {aut_count}")
    else:
        results = solve_gi_aut(adjacency_lists)
        print("Sets of isomorphic graphs and automorphisms:")
        for group, aut_count in results:
            print(f"{group}: {aut_count}")

    elapsed = time.time() - start_time
    print(f"(solved in {elapsed:.3f}s)\n")


def main():
    if len(sys.argv) < 2:
        print("Usage: python main.py <file_or_folder> [GI|Aut|GIAut]")
        print("  Task is auto-detected from filename, or you can override with 2nd argument.")
        sys.exit(1)

    # Parse command-line arguments
    path = sys.argv[1]
    task_override = sys.argv[2] if len(sys.argv) >= 3 else None

    # Check the directory given
    if os.path.isdir(path):
        # Get files 
        files = sorted(
            filename for filename in os.listdir(path)
            if os.path.splitext(filename)[1].lower() in (".gr", ".grl")
        )
        if not files:
            print(f"No .gr or .grl files found in {path}")
            sys.exit(1)
        for filename in files:
            # Call solver on each file
            solve_file(os.path.join(path, filename), task_override)
    elif os.path.isfile(path):
        # Call solver on the single file
        solve_file(path, task_override)
    else:
        print(f"Path not found: {path}")
        sys.exit(1)


