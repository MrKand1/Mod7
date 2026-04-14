def parse_grl(filepath):
    with open(filepath, 'r') as fh:
        raw = fh.readlines()
    chunks = [[]]
    for line in raw:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        if stripped[0] == '-':
            chunks.append([])
            continue
        chunks[-1].append(stripped)
    graphs = []
    for chunk in chunks:
        if not chunk:
            continue
        num_vertices = int(chunk[0])
        neighbours = [[] for _ in range(num_vertices)]
        for edge_line in chunk[1:]:
            if ':' in edge_line:
                edge_line = edge_line[:edge_line.index(':')]
            parts = edge_line.split(',')
            u, v = int(parts[0]), int(parts[1])
            neighbours[u].append(v)
            neighbours[v].append(u)
        graphs.append(neighbours)
    return graphs

def get_part_sizes(graph_idx, colors, offsets, graph_sizes):
    base = offsets[graph_idx]
    size = graph_sizes[graph_idx]
    if size == 0:
        return ()
    counts = {}
    for i in range(size):
        color = colors[base+i]
        if color not in counts:
            counts[color] = 0
        counts[color] += 1
    return tuple(sorted(counts.values()))

def get_color_counts(graph_idx, colors, offsets, graph_sizes):
    base = offsets[graph_idx]
    size = graph_sizes[graph_idx]
    if size == 0:
        return ()
    histogram = {}
    for i in range(size):
        color = colors[base+i]
        if not color in histogram:
            histogram[color] = 0
        histogram[color] += 1
    return tuple(sorted(histogram.items()))

def basic_colorref(path_to_grl):
    graphs = parse_grl(path_to_grl)
    num_graphs = len(graphs)
    if num_graphs == 0:
        return []
    graph_sizes = [len(g) for g in graphs]
    total_vertices = sum(graph_sizes)
    if total_vertices == 0:
        return [(list(range(num_graphs)), [], 0, True)]
    offsets = []
    offset = 0
    for size in graph_sizes:
        offsets.append(offset)
        offset += size
    neighbours = [None] * total_vertices
    for gr_index in range(num_graphs):
        base = offsets[gr_index]
        for v in range(graph_sizes[gr_index]):
            neighbours[base+v] = [base+u for u in graphs[gr_index][v]]
        
    colors = [0] * total_vertices
    prev_sizes = [get_part_sizes(gr_index, colors, offsets, graph_sizes) for gr_index in range(num_graphs)]
    last_change = [0] * num_graphs
    iteration = 0

    while True:
        signatures = [None] * total_vertices
        for v in range(total_vertices):
            neighbor_list = neighbours[v]
            if neighbor_list:
                neighbor_colors = sorted(colors[u] for u in neighbor_list)
                signatures[v] = (colors[v], tuple(neighbor_colors))
            else:
                signatures[v] = (colors[v], ())

        order = sorted(range(total_vertices), key=lambda v: signatures[v])
        new_colors = [0] * total_vertices
        color_id = 0
        new_colors[order[0]] = 0
        for i in range(1, total_vertices):
            if signatures[order[i]] != signatures[order[i-1]]:
                color_id += 1
            new_colors[order[i]] = color_id
        changed = False
        for v in range(total_vertices):
            if new_colors[v] != colors[v]:
                changed = True
                break
        if not changed:
            break
        colors = new_colors
        iteration += 1
        for gr_index in range(num_graphs):
            part_sizes = get_part_sizes(gr_index, colors, offsets, graph_sizes)
            if part_sizes != prev_sizes[gr_index]:
                last_change[gr_index] = iteration
                prev_sizes[gr_index] = part_sizes
    groups = {}
    for gr_index in range(num_graphs):
        key = get_color_counts(gr_index, colors, offsets, graph_sizes)
        if key not in groups:
            groups[key] = []
        groups[key].append(gr_index)
    result = []
    for histogram, gr_list in groups.items():
        gr_indices = sorted(gr_list)
        sorted_sizes = sorted(count for _, count in histogram)
        iterations = last_change[gr_indices[0]]
        discrete = all(count == 1 for _, count in histogram)
        result.append((gr_indices, sorted_sizes, iterations, discrete))
    result.sort(key=lambda t: t[0])
    return result