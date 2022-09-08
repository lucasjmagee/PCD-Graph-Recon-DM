import numpy as np
from sklearn.metrics.pairwise import pairwise_distances
from math import inf, exp
import os
import csv
import sys


def __sorted_insert(sorted_list, element):
    if len(sorted_list) == 0:
        return [element]
    for i in range(len(sorted_list)):
        ith = sorted_list[i]
        if ith < element:
            continue
        if ith == element:
            print('error dup element')
            sys.exit()
        return sorted_list[:i] + [element] + sorted_list[i:]
    #print('need last statement')
    #sys.exit()
    return sorted_list + [element]


def __gaussian_kernel(dist, bandwidth):
    return exp(-(dist ** 2) / (2 * bandwidth ** 2))


def density_estimation(input_filename, output_filename, k=15):
    with open(input_filename, 'r') as input_file:
        content = input_file.readlines()
        input_file.close()

    # content = content[1:]
    cells = [c.strip().split(' ') for c in content]
    # cells = [c for c in cells if c[4] in VALID_LABELS]
    print('number of cells:', len(cells))
    print('dims:', len(cells[0]))


    cells_array = np.asarray(cells)
    print('computing distances')
    # sklearn_dm = pairwise_distances(cells)
    sklearn_dm = pairwise_distances(cells, metric='l1')

    print('computing bandwidth')
    val_sum = 0
    count = 1
    for row in sklearn_dm:
        # print('working on', count)
        count += 1
        vals = sorted(row)[:k + 1]
        # print(vals)
        val_sum += sum(vals)

    bandwidth = val_sum / (len(cells) * k)
    # print(BANDWIDTH)

    print('calculating function values')
    function_vals = []
    for i in range(len(cells)):
        #print('func: ', i, '/', len(cells))
        val = 0
        distances = sklearn_dm[i]
        for dist in distances:
            val += __gaussian_kernel(dist, bandwidth)
        function_vals.append(val)

    with open(output_filename, 'w') as output_file:
        for f in function_vals:
            output_file.write(str(f) + '\n')
        output_file.close()


def compute_rips_complex_edges(input_filename, output_filename, thresh, metric='euclidean'):
    with open(input_filename, 'r') as input_file:
        content = input_file.readlines()
        input_file.close()
    cells = [c.strip().split(' ') for c in content]
    cells = [[float(n) for n in c] for c in cells]

    sklearn_dm = pairwise_distances(cells, metric=metric)

    edges = []
    for i in range(len(cells)):
        for j in range(i + 1, len(cells)):
            edges.append((i, j, sklearn_dm[i, j]))

    cnt = 0
    with open(output_filename, 'w') as output_file:
        for e in edges:
            if e[2] > thresh:
                continue
            cnt += 1
            output_file.write(str(e[0]) + ' ' + str(e[1]) + ' ' + str(e[2]) + '\n')
        output_file.close()


def build_baseline_complex(rips_edge_filename, density_filename, output_filename, threshold=inf):
    densities = []
    with open(density_filename, 'r') as density_file:
        reader = csv.reader(density_file, delimiter=' ')
        for row in reader:
            densities.append(float(row[0]))
        density_file.close()

    # content = content[1:]
    densities = [(i, densities[i]) for i in range(len(densities))]

    # print(cells[0])
    # print(cells[1])
    # sys.exit()

    print('reading in edges')
    weighted_edges = []
    with open(rips_edge_filename, 'r') as weighted_file:
        reader = csv.reader(weighted_file, delimiter=' ')
        for row in reader:
            weighted_edges.append([float(row[2]), int(row[0]), int(row[1])])
        weighted_file.close()
    print('read in', len(weighted_edges), 'weighted edges')

    possible_edges = weighted_edges

    print('possible edges:', len(possible_edges))

    possible_edges.sort()

    valid_edges = []

    tris = []
    filtration = []
    adjacent_map = {}

    for i in range(len(densities)):
        adjacent_map[i] = []

    edge_count = 1
    for e in possible_edges:
        if edge_count % 10000 == 0:
            print('edges:', edge_count, e[0])
        # print('working on edge', edge_count)
        edge_count += 1
        if e[0] > threshold:
            break
        v_e = (e[1], e[2])
        valid_edges.append(v_e)
        filtration.append(v_e)
        vi = e[1]
        vj = e[2]

        i_adj = adjacent_map[vi]
        j_adj = adjacent_map[vj]

        for n in i_adj:
            if n not in j_adj:
                continue
            tri_verts = [vi, vj, n]
            tri_verts.sort()
            t = (tri_verts[0], tri_verts[1], tri_verts[2])
            tris.append(t)
            filtration.append(t)
        adjacent_map[vi] = __sorted_insert(i_adj, vj)
        adjacent_map[vj] = __sorted_insert(j_adj, vi)

    print(len(valid_edges), '/', len(possible_edges), 'are valid')

    nV = len(densities)
    nE = len(valid_edges)
    nT = len(tris)
    nV = np.asarray(nV)
    nE = np.asarray(nE)
    nT = np.asarray(nT)

    # cells = np.asarray(cells)
    densities = np.asarray(densities)
    valid_edges = np.asarray(valid_edges)
    tris = np.asarray(tris)

    with open(output_filename, 'wb') as f:
        nV.astype(np.int32).tofile(f)
        densities.astype('d').tofile(f)
        nE.astype(np.int32).tofile(f)
        valid_edges.astype(np.int32).tofile(f)
        nT.astype(np.int32).tofile(f)
        tris.astype(np.int32).tofile(f)
        f.close()


def compute_persistence_baseline(input_filename, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    command = "./dmpcd/code/persistence_baseline/spt_cpp " + input_filename + ' ' + output_dir + ' 0 1'
    os.system(command)


def compute_graph_reconstruction_baseline(density_filename, edge_persistence_filename, persistence_threshold, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    print('checking...')
    command = "./dmpcd/code/baseline_morse_src/a.out " + density_filename + ' ' + edge_persistence_filename + ' ' + str(persistence_threshold) + ' ' + output_dir
    os.system(command)
