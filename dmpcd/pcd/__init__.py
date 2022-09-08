import numpy as np
from sklearn.metrics.pairwise import pairwise_distances
from math import inf
import os
import csv
import sys


def __choose_first_in_perm(distance_matrix):
    return 0


def __greedy_perm(distance_matrix, n_cells):
    perm = []
    intersection_radii = []
    available = list(range(n_cells))
    min_to_perm_dict = {}

    while len(perm) < n_cells:
        # print('perm len:', len(perm))
        if len(perm) == 0:
            first_element = __choose_first_in_perm(distance_matrix)
            perm.append(first_element)
            available.remove(first_element)
            intersection_radii.append(inf)
            for a in available:
                min_to_perm_dict[a] = distance_matrix[a][first_element]
            continue
        max_dist_to_prefix = -1
        max_index = -1
        for a in available:
            dist_to_perm = min_to_perm_dict[a]
            if dist_to_perm > max_dist_to_prefix:
                max_dist_to_prefix = dist_to_perm
                max_index = a
        perm.append(max_index)
        available.remove(max_index)
        intersection_radii.append(max_dist_to_prefix)
        for a in available:
            min_to_perm_dict[a] = min(distance_matrix[a][max_index], min_to_perm_dict[a])
    return perm, intersection_radii

'''
def __perturbation(p, a, intersection_radius_dict):
    intersection_radius = intersection_radius_dict[p]
    if a <= intersection_radius / epsilon:
        return 0
    elif a < intersection_radius / (epsilon * (1 - epsilon)):
        return a - intersection_radius / epsilon
    else:
        return a * epsilon


def __perturbed_distance(p, q, a, distance_matrix, intersection_radii_dict):
    return distance_matrix[p, q] + __perturbation(p, a, intersection_radii_dict) \
           + __perturbation(q, a, intersection_radii_dict)
'''

# sparse weighted rips filtration functions
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


def build_sparse_weighted_rips_filtration(feature_filename, output_dir, k=15, metric='euclidean', epsilon=.99, cutoff=inf):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    with open(feature_filename, 'r') as feature_file:
        content = feature_file.readlines()
        feature_file.close()

    verts = [c.strip().split(' ') for c in content]
    verts = [[float(n) for n in c] for c in verts]

    print('number of point clouds:', len(verts))
    print('dims:', len(verts[0]))

    vert_array = np.asarray(verts)
    sklearn_dm = pairwise_distances(verts, metric=metric)

    permutation, intersection_radii = __greedy_perm(sklearn_dm, len(verts))

    perm_index_dict = {}
    og_intersection_dict = {}
    for i in range(len(verts)):
        index = permutation.index(i)
        perm_index_dict[i] = index
        og_intersection_dict[i] = intersection_radii[index]

    edges = []
    for p in range(len(verts)):
        for q in range(p + 1, len(verts)):
            edges.append([p, q])

    filtration_net_alphas = []
    for e in edges:
        p = e[0]
        q = e[1]
        min_net_val = min(og_intersection_dict[p], og_intersection_dict[q])
        filtration_net_alpha = min_net_val / (epsilon * (1 - epsilon))
        filtration_net_alphas.append(filtration_net_alpha)

    print('calculating alpha in which edges fit distance requirement')
    total_not_in_filt = 0
    distance_alphas = []
    for i in range(len(edges)):
        # print('working on', i,'/',len(edges))
        ith_edge = edges[i]
        p = ith_edge[0]
        q = ith_edge[1]

        # print('\n\np', p)
        # print('q', q)

        metric_distance = sklearn_dm[p][q]
        p_lambda = og_intersection_dict[p]
        q_lambda = og_intersection_dict[q]
        max_to_stay_in_net = filtration_net_alphas[i]

        p_lower = p_lambda / epsilon
        q_lower = q_lambda / epsilon

        p_upper = p_lambda / (epsilon * (1 - epsilon))
        q_upper = q_lambda / (epsilon * (1 - epsilon))

        alpha = metric_distance / 2.0
        if alpha <= p_lower and alpha <= q_lower and alpha <= max_to_stay_in_net:
            # print('1 alpha', alpha, metric_distance, p_lower, q_lower, max_to_stay_in_net)
            distance_alphas.append(alpha)
            if alpha < 0:
                print('error', 1)
            continue

        if p_lower <= q_lower:
            if p_lower >= max_to_stay_in_net:
                total_not_in_filt += 1
                distance_alphas.append(inf)
                continue

            alpha = metric_distance - p_lambda / epsilon
            if alpha < p_upper and alpha < q_lower and alpha <= max_to_stay_in_net:
                assert (alpha >= p_lower)
                # print('2 alpha', alpha, metric_distance, p_lower, p_upper, q_lower, max_to_stay_in_net)
                distance_alphas.append(alpha)
                if alpha < 0:
                    print('error', 2)
                    sys.exit()
                continue
            if p_upper <= q_lower:

                alpha = metric_distance / (2 - epsilon)
                if alpha <= q_lower and alpha <= max_to_stay_in_net:
                    assert (alpha >= p_upper)
                    distance_alphas.append(alpha)
                    if alpha < 0:
                        print('error', 3)
                        sys.exit()
                    continue
                alpha = (metric_distance - q_lambda / epsilon) / (1 - epsilon)
                # alpha = (-EPSILON * metric_distance + 1) / (EPSILON ** 2 - EPSILON)
                # print('num:', (-EPSILON * metric_distance + 1))
                # print('den:', EPSILON ** 2 - EPSILON)
                '''
                if alpha < 0:
                    #print(4, 'less than zero')
                    alpha = q_lower
                '''
                if alpha < q_upper and alpha <= max_to_stay_in_net and q_lower != inf:
                    if q_lower >= max_to_stay_in_net:
                        total_not_in_filt += 1
                        distance_alphas.append(inf)
                        continue
                    assert (alpha > q_lower)
                    '''
                    if alpha < q_lower:
                        #print('4 switch!', i)
                        alpha = q_lower
                    '''
                    distance_alphas.append(alpha)
                    # print('4 alpha', alpha, metric_distance, p_upper, q_lower, q_upper, max_to_stay_in_net, q_lambda)
                    if alpha < 0:
                        '''
                        print('md:', metric_distance)
                        print('md - qlambda / E:', metric_distance - q_lambda / EPSILON)
                        print('alpha:', (metric_distance - q_lambda / EPSILON) / (1 - EPSILON))
                        print('error', 4, p, q, alpha, metric_distance, q_lambda, EPSILON, 1 - EPSILON, p_upper, q_lower)
                        '''
                        print('error', 4)
                        sys.exit()
                    continue
                alpha = metric_distance / (2 - 2 * epsilon)
                if alpha <= max_to_stay_in_net and q_lower != inf:
                    assert (alpha >= q_upper)
                    distance_alphas.append(alpha)
                    if alpha < 0:
                        print('error', 5)
                        sys.exit()
                    continue
            else:  # q_lower < p_upper
                if q_lower >= max_to_stay_in_net:
                    total_not_in_filt += 1
                    distance_alphas.append(inf)
                    continue
                if metric_distance <= p_lambda / epsilon + q_lambda / epsilon and max_to_stay_in_net > q_lower and max_to_stay_in_net > p_upper:
                    distance_alphas.append(q_lower)
                    if alpha < 0:
                        print('error', 6)
                        sys.exit()
                    continue
                if p_upper <= q_upper:
                    alpha = (metric_distance - q_lambda / epsilon) / (1 - epsilon)
                    if alpha < q_upper and alpha <= max_to_stay_in_net:
                        assert (alpha >= p_upper)
                        distance_alphas.append(alpha)
                        if alpha < 0:
                            print('error', 7)
                            sys.exit()
                        continue
                    alpha = metric_distance / (2 - 2 * epsilon)
                    if alpha <= max_to_stay_in_net:
                        assert (alpha >= q_upper)
                        distance_alphas.append(alpha)
                        if alpha < 0:
                            print('error', 8)
                            sys.exit()
                        continue
                else:  # q_upper < p_upper
                    alpha = (metric_distance - p_lambda / epsilon) / (1 - epsilon)
                    if alpha < p_upper and alpha <= max_to_stay_in_net:
                        assert (alpha >= q_upper)
                        distance_alphas.append(alpha)
                        if alpha < 0:
                            print('error', 9)
                            sys.exit()
                        continue
                    alpha = metric_distance / (2 - 2 * epsilon)
                    if alpha <= max_to_stay_in_net:
                        assert (alpha >= q_upper)
                        distance_alphas.append(alpha)
                        if alpha < 0:
                            print('error', 10)
                            sys.exit()
                        continue

        else:  # q_lower < p_lower
            if q_lower >= max_to_stay_in_net:
                total_not_in_filt += 1
                distance_alphas.append(inf)
                continue

            alpha = metric_distance - q_lambda / epsilon
            if alpha < q_upper and alpha < p_lower and alpha <= max_to_stay_in_net:
                assert (alpha >= q_lower)
                distance_alphas.append(alpha)
                if alpha < 0:
                    print('error', 11)
                    sys.exit()
                continue
            if p_lower <= q_upper:
                if p_lower >= max_to_stay_in_net:
                    total_not_in_filt += 1
                    distance_alphas.append(inf)
                    continue
                if metric_distance <= p_lambda / epsilon + q_lambda / epsilon and max_to_stay_in_net > q_lower and max_to_stay_in_net > p_upper:
                    distance_alphas.append(p_lower)
                    if alpha < 0:
                        print('error', 12)
                        sys.exit()
                    continue
                if p_upper <= q_upper:
                    alpha = (metric_distance - q_lambda / epsilon) / (1 - epsilon)
                    if alpha < q_upper and alpha <= max_to_stay_in_net:
                        assert (alpha >= p_upper)
                        distance_alphas.append(alpha)
                        if alpha < 0:
                            print('error', 13)
                            sys.exit()
                        continue
                    alpha = metric_distance / (2 - 2 * epsilon)
                    if alpha <= max_to_stay_in_net:
                        assert (alpha >= q_upper)
                        distance_alphas.append(alpha)
                        if alpha < 0:
                            print('error', 14)
                            sys.exit()
                        continue
                else:  # q_upper < p_upper
                    alpha = (metric_distance - p_lambda / epsilon) / (1 - epsilon)
                    if alpha < p_upper and alpha <= max_to_stay_in_net:
                        assert (alpha >= q_upper)
                        distance_alphas.append(alpha)
                        if alpha < 0:
                            print('error', 15)
                            sys.exit()
                        continue
                    alpha = metric_distance / (2 - 2 * epsilon)
                    if alpha <= max_to_stay_in_net:
                        assert (alpha >= p_upper)
                        distance_alphas.append(alpha)
                        if alpha < 0:
                            print('error', 16)
                            sys.exit()
                        continue
            else:  # q_upper < p_lower
                alpha = metric_distance / (2 - epsilon)
                # print('alpha before:', alpha)
                if alpha <= p_lower and alpha <= max_to_stay_in_net:
                    assert (alpha >= q_upper)
                    distance_alphas.append(alpha)
                    if alpha < 0:
                        print('error', 17)
                        sys.exit()
                    continue
                alpha = (metric_distance - p_lambda / epsilon) / (1 - epsilon)
                '''
                if alpha < 0:
                    #print(18, 'less than zero')
                    alpha = p_lower
                '''
                if p_lower >= max_to_stay_in_net:
                    total_not_in_filt += 1
                    distance_alphas.append(inf)
                    continue
                if alpha < p_upper and alpha <= max_to_stay_in_net and p_lower != inf:
                    assert (alpha >= p_lower)
                    distance_alphas.append(alpha)
                    # print('18 alpha', alpha, metric_distance, p_lower, q_upper, max_to_stay_in_net)
                    if alpha < 0:
                        '''
                        print('error p_lambda:', p_lambda)
                        print()
                        print('error', 18, p, q, alpha, metric_distance, p_lambda, EPSILON, 1 - EPSILON, q_upper, p_lower)
                        '''
                        sys.exit()
                    continue
                alpha = metric_distance / (2 - 2 * epsilon)
                if alpha <= max_to_stay_in_net and p_lower != inf:
                    assert (alpha >= p_upper)
                    distance_alphas.append(alpha)
                    if alpha < 0:
                        print('error', 19)
                        sys.exit()
                    continue

        # if you get here the edge will not be added
        # print('edge', i, 'will not be in filtration')
        total_not_in_filt += 1
        distance_alphas.append(inf)

    sparse_edges = []
    for i in range(len(edges)):
        if distance_alphas[i] == inf:
            continue
        sparse_edges.append([edges[i][0], edges[i][1], distance_alphas[i],
                             sklearn_dm[edges[i][0], edges[i][1]]])

    sparse_verts = set()
    for e in sparse_edges:
        sparse_verts.add(e[0])
        sparse_verts.add(e[1])

    sparse_verts = list(sparse_verts)
    sparse_verts.sort()

    print(len(sparse_verts), '/', len(verts), 'in filtration')

    print('calculating weights...')
    weights = {}
    for v in sparse_verts:
        # print('working on weights for cell', c)
        distances = sklearn_dm[v]
        distances.sort()
        weight = 0
        for i in range(1, k + 1):
            weight += distances[i] ** 2
        weight = (weight / k) ** .5
        weights[v] = weight

    del sklearn_dm

    smaller = 0
    larger = 0
    same = 0

    for e in sparse_edges:
        i = e[0]
        j = e[1]
        distance = e[3]
        w_i = weights[i]
        w_j = weights[j]
        # print('weight 1:', w_i)
        # print('weight 2: ', w_j)

        # alpha = Symbol('alpha')
        # answer = max(solve(distance - ((alpha ** 2 - w_i ** 2) ** .5 + (alpha ** 2 - w_j ** 2) ** .5), alpha))

        '''
        assert(abs(answer - ((distance ** 4) + (w_i ** 4) + (w_j ** 4) + (2 * (distance ** 2) * (w_i ** 2))
                 + (2 * (distance ** 2) * (w_j ** 2)) - (2 * (w_i ** 2) * (w_j ** 2))) ** .5
                / (2 * distance)) < EPS)
        '''

        answer = ((distance ** 4) + (w_i ** 4) + (w_j ** 4) + (2 * (distance ** 2) * (w_i ** 2)) + (
                    2 * (distance ** 2) * (w_j ** 2)) - (2 * (w_i ** 2) * (w_j ** 2))) ** .5 / (2 * distance)

        if not abs(distance - ((answer ** 2 - w_i ** 2) ** .5 + (answer ** 2 - w_j ** 2) ** .5)) < .0001:
            answer = max(w_i, w_j)
        if not abs(distance - ((answer ** 2 - w_i ** 2) ** .5 + (answer ** 2 - w_j ** 2) ** .5)) < .0001:
            print(distance, ((answer ** 2 - w_i ** 2) ** .5 + (answer ** 2 - w_j ** 2) ** .5), answer, w_i, w_j)

        # sys.exit()

        # assert(abs(distance - ((answer ** 2 - w_i ** 2) ** .5 + (answer ** 2 - w_j ** 2) ** .5)) < EPS)
        if answer < 0:
            print('had to negate, investigate')
            answer = -answer

        if distance / 2 < answer:
            larger += 1
        elif distance / 2 > answer:
            smaller += 1
        else:
            same += 1

        answer = max(answer, e[2], w_i, w_j)
        e.append(answer)

    sparse_edges.sort(key=lambda x: x[4])

    tris = []
    filtration = []
    adjacent_map = {}

    for v in sparse_verts:
        adjacent_map[v] = []

    valid_edges = []

    edge_count = 1
    for e in sparse_edges:
        # print('edges:', edge_count, e)
        # print('working on edge', edge_count)
        edge_count += 1
        if e[4] > cutoff:
            break

        '''
        if e[3] < e[4]:
            larger += 1
        elif e[3] > e[4]:
            smaller += 1
        else:
            same += 1
        '''

        v_e = (e[0], e[1], e[4])
        valid_edges.append(v_e)
        filtration.append(v_e)
        vi = e[0]
        vj = e[1]

        i_adj = adjacent_map[vi]
        j_adj = adjacent_map[vj]

        for n in i_adj:
            if n not in j_adj:
                continue
            tri_verts = [vi, vj, n]
            tri_verts.sort()
            t = (tri_verts[0], tri_verts[1], tri_verts[2], e[4])
            tris.append(t)
            filtration.append(t)
        adjacent_map[vi] = __sorted_insert(i_adj, vj)
        adjacent_map[vj] = __sorted_insert(j_adj, vi)

    print(len(valid_edges), '/', len(sparse_edges), 'are below cutoff')

    nV = len(weights)
    nE = len(valid_edges)
    nT = len(tris)
    nV = np.asarray(nV)
    nE = np.asarray(nE)
    nT = np.asarray(nT)

    valid_edges = np.asarray(valid_edges)
    tris = np.asarray(tris)

    weight_vert_pair = [[weights[i], i] for i in range(len(verts)) if i in sparse_verts]
    weight_vert_pair.sort()

    cell_mapping = {}
    for i in range(len(weight_vert_pair)):
        cell_mapping[weight_vert_pair[i][1]] = i

    print('writing weight information for later visualization')
    weights_filename = os.path.join(output_dir, 'weights.txt')
    with open(weights_filename, 'w') as weights_file:
        for i in range(len(weights)):
            if i not in sparse_verts:
                weights_file.write(str(-1) + '\n')
                continue
            weights_file.write(str(weights[i]) + '\n')
        weights_file.close()

    cell_ind = 0
    filtration_ind = 0
    counter = 0

    output_filename = os.path.join(output_dir, 'sparse_weighted_rips_filtration.txt')
    with open(output_filename, 'w') as output_file:
        while counter < len(weight_vert_pair) + len(filtration):
            current_cell = weight_vert_pair[cell_ind]
            # print(current_cell)
            # sys.exit()
            current_simp = filtration[filtration_ind]
            if current_cell[0] <= current_simp[len(current_simp) - 1]:
                output_file.write('0 ' + str(current_cell[1]) + ' ' + str(current_cell[0]) + '\n')
                # final_filtration.append([current_cell[1], current_cell[0]])
                cell_ind += 1
                # print(cell_ind)
                if cell_ind == len(weight_vert_pair):
                    for i in range(filtration_ind, len(filtration)):
                        ith_simp = filtration[i]
                        # print('ith simp', ith_simp)
                        if len(ith_simp) == 3:
                            # print('edge')
                            output_file.write('1 ')
                        else:
                            # print('tri')
                            assert (len(ith_simp) == 4)
                            output_file.write('2 ')
                        for j in range(len(ith_simp) - 1):
                            output_file.write(str(cell_mapping[ith_simp[j]]) + ' ')
                        output_file.write(str(ith_simp[len(ith_simp) - 1]) + '\n')
                    break
                counter += 1
            else:
                if len(current_simp) == 3:
                    output_file.write('1 ')
                else:
                    assert (len(current_simp) == 4)
                    output_file.write('2 ')

                for i in range(len(current_simp) - 1):
                    output_file.write(str(cell_mapping[current_simp[i]]) + ' ')
                output_file.write(str(current_simp[len(current_simp) - 1]) + '\n')
                # final_filtration.append(current_simp)
                filtration_ind += 1

                if filtration_ind == len(filtration):
                    break
                    print('ERROR', cell_ind)
                    sys.exit()

                counter += 1
        output_file.close()


def reorder_weights(input_filename, output_filename):
    weights = []
    with open(input_filename, 'r') as input_file:
        reader = csv.reader(input_file, delimiter=' ')
        for row in reader:
            weights.append(float(row[0]))
        input_file.close()
    weights.sort()
    with open(output_filename, 'w') as output_file:
        for w in weights:
            output_file.write(str(w) + '\n')
        output_file.close()


def compute_persistence_swr(input_filename, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    command = "./dmpcd/code/persistence_swr/spt_cpp " + input_filename + ' ' + output_dir + ' 0 1'
    os.system(command)


def reorder_verts_by_weight(weights_filename, verts_filename, output_filename):
    weights = []
    with open(weights_filename, 'r') as input_file:
        reader = csv.reader(input_file, delimiter=' ')
        for row in reader:
            weights.append(float(row[0]))
        input_file.close()

    with open(verts_filename, 'r') as input_file:
        content = input_file.readlines()
        input_file.close()
    cells = content
    pairs = [(weights[i], cells[i]) for i in range(len(cells))]
    pairs.sort()

    with open(output_filename, 'w') as output_file:
        for p in pairs:
            output_file.write(p[1])
        output_file.close()


def reorder_verts_and_annos_by_weight(weights_filename, verts_filename, anno_filename, output_vert_filename, output_anno_filename):
    weights = []
    with open(weights_filename, 'r') as input_file:
        reader = csv.reader(input_file, delimiter=' ')
        for row in reader:
            weights.append(float(row[0]))
        input_file.close()

    with open(verts_filename, 'r') as input_file:
        content = input_file.readlines()
        input_file.close()
    cells = content

    annos = []
    with open(anno_filename, 'r') as input_anno:
        reader = csv.reader(input_anno, delimiter=' ')
        for row in reader:
            annos.append(row[0])
        input_anno.close()

    pairs = [(weights[i], cells[i], annos[i]) for i in range(len(cells))]
    pairs.sort()

    with open(output_vert_filename, 'w') as output_file:
        for p in pairs:
            output_file.write(p[1])
        output_file.close()

    with open(output_anno_filename, 'w') as output_file:
        for p in pairs:
            output_file.write(p[2])
        output_file.close()


def compute_graph_reconstruction(sorted_weights_filename, edge_persistence_filename, persistence_threshold, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    command = "./dmpcd/code/phat_morse_src/a.out " + sorted_weights_filename + ' ' + edge_persistence_filename + ' ' + str(persistence_threshold) + ' ' + output_dir
    os.system(command)
