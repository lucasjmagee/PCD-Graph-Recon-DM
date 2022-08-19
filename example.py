import dmpcd as dm
import os


def pcd_test():
    input_filename = "data/1-circle/features.txt"
    output_dir = "results/1-circle-pcd/"
    k = 15
    metric = 'euclidean'
    epsiilon = .99

    persistence_threshold = .25

    dm.build_sparse_weighted_rips_filtration(input_filename, output_dir, k, metric, epsiilon)
    filtration_filename = os.path.join(output_dir, 'sparse_weighted_rips_filtration.txt')
    weights_filename = os.path.join(output_dir, 'weights.txt')

    dm.compute_persistence_swr(filtration_filename, output_dir)
    edge_filename = os.path.join(output_dir, "edge_for_morse_only.txt")

    sorted_weights_filename = os.path.join(output_dir, "sorted-weights.txt")
    dm.reorder_weights(weights_filename, sorted_weights_filename)

    morse_dir = os.path.join(output_dir, str(persistence_threshold) + '/')
    dm.compute_graph_reconstruction(sorted_weights_filename, edge_filename, persistence_threshold, morse_dir)

    result_edge_filename = os.path.join(morse_dir, 'dimo_edge.txt')
    sorted_feature_filename = os.path.join(output_dir, 'sorted-feature.txt')
    dm.reorder_verts_by_weight(weights_filename, input_filename, sorted_feature_filename)

    dm.visualize_results_2d(sorted_feature_filename, result_edge_filename)


def baseline_test():
    input_filename = "data/1-circle/features.txt"
    output_dir = "results/1-circle-baseline/"
    k = 15
    rips_alpha = .25
    persistence_threshold = 998

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    density_filename = os.path.join(output_dir, 'density.txt')
    dm.density_estimation(input_filename, density_filename)

    rips_edge_filename = os.path.join(output_dir, 'rips-edge.txt')
    dm.compute_rips_complex_edges(input_filename, rips_edge_filename, rips_alpha, metric='euclidean')

    complex_filename = os.path.join(output_dir, 'complex.bin')
    dm.build_baseline_complex(rips_edge_filename, density_filename, complex_filename)

    dm.compute_persistence_baseline(complex_filename, output_dir)
    edge_persistence_filename = os.path.join(output_dir, 'edge_for_morse_only.txt')

    morse_dir = os.path.join(output_dir, str(persistence_threshold) + '/')
    if not os.path.exists(morse_dir):
        os.mkdir(morse_dir)
    dm.compute_graph_reconstruction_baseline(density_filename, edge_persistence_filename, persistence_threshold, morse_dir)
    morse_edge_filename = os.path.join(morse_dir, 'dimo_edge.txt')

    dm.visualize_results_2d(input_filename, morse_edge_filename)


if __name__ == '__main__':
    pcd_test()
    baseline_test()
