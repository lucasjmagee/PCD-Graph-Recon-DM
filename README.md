# Point Cloud Discrete Morse Graph Reconstruction

dmpcd is a python package meant for executing the discrete Morse graph reconstruction algorithm on Point Cloud Datasets - designed with full mouse brain imaging data in mind.  The package includes functions which allow the user to perform the baseline DM approach and the more recently developed PCD DM approach.  For the basline approach, the package allows users to compute the Rips filtration, estimate density at each point, compute persistence of lower star filtration, and perform discrete Morse graph reconstruction.  For the PCD DM approach, users can compute the sparse weighted Rips filtration, compute persistence of the filtration, and perform discrete Morse graph reconstruction.  Data can be labeled or unlabeled, with visualization tools included for both cases.

* [Installation Intructions](#installation-instructions)
  * [System Requirements](#system-requirements)
  * [Required Python Libraries](#required-python-libraries)
  * [Compiling Code](#compiling-code)
* [dmpcd.pcd Functions](#dmpcd.pcd-functions)
* [dmpcd.baseline Functions](#dmpcd.baseline-functions)
* [Other dmpcd Functions](#other-dmpcd-functions)
* [Separate Programs](#separate-programs)
* [Example Use of Pipeline](#example-use-of-pipeline) 
  * [PCD](#pcd)
  * [Baseline](#baseline)

## Installation Instructions
### System Requirements
- Python 3.8.8 (or newer)
- g++ 9.4.0 (or newer)

### Required Python Libraries
- scitkit-learn - pip install scikit-learn (https://scikit-learn.org/)

### Compiling Code

Sparse Weighted Rips Filtration Persistence Module (PCD)
    
    > cd DiMo3d/code/persistence_swr/
    > g++ DiMoSC.cpp -I./phat/include -std=c++11 -o spt_cpp

Lower Star Filtration Persistence Module (baseline)

    > cd DiMo3d/code/spt_cpp/
    > g++ DiMoSC.cpp -I./phat/include -std=c++11 -o spt_cpp
    
Discrete Morse Graph Reconstruction (PCD)
    > cd DiMo3d/code/phat_morse_src/
    > g++ ComputeGraphReconstruction.cpp
    
Discrete Morse Graph Reconstruction (Baseline)
    > cd DiMo3d/code/baseline_morse_src/
    > g++ ComputeGraphReconstruction.cpp


## dmpcd.pcd Functions

### dmpcd.pcd.build_sparse_weighted_rips_filtration(feature_filename, output_dir, k=15, metric='euclidean', epsilon=.99, cutoff=inf)

#### Description
Divide the input domain into overlapping sub-rectangular prisms.

#### Input
- feature_filename - file where PCD is stored
- output_dir - directory where filtration file and weight file will be written to
- k - k value used for weight estimation (inverse density)
- metric - metric used for sparse weighted Rips filtration, default is euclidean
- epsilon - sparsification parameter
- cutoff - radius cut off for filtration

#### Output

Output dir is made containing a file with the filtration that will be passed Sparse Weighted Rips Filtration Persistence Module, and a weights file containing the weights of each point.

#### Example

    >import dmpcd as dm
    >import dmpcd.pcd as pcd

    >feature_filename = ???data/1-circle/features.txt???
    >output_dir = ???results/1-circle/???
    >k=15
    >metric='euclidean'
    >epsilon=.99
    >cutoff=inf
    >dm.pcd.build_sparse_weighted_rips_filtration(feature_filename, output_dir, k, metric, epsilon, cutoff)
    
### dmpcd.pcd.compute_persistence_swr(input_filename, output_dir)

#### Description
Compute persistence of sparse weighted Rips filtrations

#### Input
- input_filename - file where filtration is written
- output_dir - where persistence values of edges in filtration are written

#### Output

File containing persistence values and dimension of each edge in input filtration.  This is later used for DM graph reconstruction.

#### Example

    >import dmpcd as dm
    >import dmpcd.pcd as pcd

    >feature_filename = ???data/1-circle/features.txt???
    >output_dir = ???results/1-circle/???
    >k=15
    >metric='euclidean'
    >epsilon=.99
    >cutoff=inf
    >dm.pcd.build_sparse_weighted_rips_filtration(feature_filename, output_dir, k, metric, epsilon, cutoff)
    >dm.pcd.compute_persistence_swr(os.path.join(output_dir, "sparse_weighted_rips_filtration.txt"), output_dir)

### dmpcd.pcd.reorder_weights(input_filename, output_filename)

#### Description
order weights in ascending order - this is needed for input into PCD graph reconstruction.

#### Input
- input_filename - file where weights are written - rows match rows of original input_file for build_sparse_weighted_rips_filtration
- output_filename - file where sorted weights will be written

#### Output

File containing persistence values and dimension of each edge in input filtration.  This is later used for DM graph reconstruction.

#### Example

    >import dmpcd as dm
    >import dmpcd.pcd as pcd

    >feature_filename = ???data/1-circle/features.txt???
    >output_dir = ???results/1-circle/???
    >k=15
    >metric='euclidean'
    >epsilon=.99
    >cutoff=inf
    >dm.pcd.build_sparse_weighted_rips_filtration(feature_filename, output_dir, k, metric, epsilon, cutoff)
    >weights_filename = os.path.join(output_dir, 'weights.txt')
    >sorted_weights_filename = os.path.join(output_dir, "sorted-weights.txt")
    >dm.pcd.reorder_weights(weights_filename, sorted_weights_filename)
    
### dmpcd.pcd.compute_graph_reconstruction(sorted_weights_filename, edge_persistence_filename, persistence_threshold, output_dir)

#### Description
Compute DM graph reconstruction

#### Input
- sorted_weights_filename - list of vert weights in sorted order
- edge_persistence_filename - persistence values and type for each edge in domain
- persistence threshold - threshold value used by algorithm to simplify output
- output_dir - directory graph (.txt) is written to

#### Output

Graph reconstruction of PCD

#### Example

    >import dmpcd as dm
    >import dmpcd.pcd as pcd

    >input_filename = "data/1-circle/features.txt"
    >output_dir = "results/1-circle-pcd/"
    >k = 15
    >metric = 'euclidean'
    >epsiilon = .99
    >persistence_threshold = .25

    >dm.pcd.build_sparse_weighted_rips_filtration(input_filename, output_dir, k, metric, epsiilon)
    >filtration_filename = os.path.join(output_dir, 'sparse_weighted_rips_filtration.txt')
    >weights_filename = os.path.join(output_dir, 'weights.txt')

    >dm.pcd.compute_persistence_swr(filtration_filename, output_dir)
    >edge_filename = os.path.join(output_dir, "edge_for_morse_only.txt")

    >sorted_weights_filename = os.path.join(output_dir, "sorted-weights.txt")
    >dm.pcd.reorder_weights(weights_filename, sorted_weights_filename)

    >morse_dir = os.path.join(output_dir, str(persistence_threshold) + '/')
    >dm.pcd.compute_graph_reconstruction(sorted_weights_filename, edge_filename, persistence_threshold, morse_dir)

### dmpcd.pcd.reorder_verts_by_weight(weights_filename, verts_filename, output_filename)

#### Description
reorder original data by weights computed in dmpcd.build_sparse_weighted_rips_filtration

#### Input
- weights_filename - weight file outputted by dmpcd.build_sparse_weighted_rips_filtration
- verts_filename - original PCD dataset file for dmpcd.build_sparse_weighted_rips_filtration
- output_filename - PCD dataset with points listed in ascending weight order

#### Output

reordered PCD dataset based on weight values

#### Example

    >import dmpcd as dm
    >import dmpcd.pcd as pcd

    >input_filename = "data/1-circle/features.txt"
    >output_dir = "results/1-circle-pcd/"
    >k = 15
    >metric = 'euclidean'
    >epsiilon = .99
    >persistence_threshold = .25

    >dm.pcd.build_sparse_weighted_rips_filtration(input_filename, output_dir, k, metric, epsiilon)
    >filtration_filename = os.path.join(output_dir, 'sparse_weighted_rips_filtration.txt')
    >weights_filename = os.path.join(output_dir, 'weights.txt')

    >dm.pcd.compute_persistence_swr(filtration_filename, output_dir)
    >edge_filename = os.path.join(output_dir, "edge_for_morse_only.txt")

    >sorted_weights_filename = os.path.join(output_dir, "sorted-weights.txt")
    >dm.pcd.reorder_weights(weights_filename, sorted_weights_filename)

    >morse_dir = os.path.join(output_dir, str(persistence_threshold) + '/')
    >dm.pcd.compute_graph_reconstruction(sorted_weights_filename, edge_filename, persistence_threshold, morse_dir)

    >result_edge_filename = os.path.join(morse_dir, 'dimo_edge.txt')
    >sorted_feature_filename = os.path.join(output_dir, 'sorted-feature.txt')
    >dm.pcd.reorder_verts_by_weight(weights_filename, input_filename, sorted_feature_filename)
    
### dmpcd.pcd.reorder_verts_and_annos_by_weight(weights_filename, verts_filename, anno_filename, output_vert_filename, output_anno_filename)

#### Description
reorder original data (and annotations) by weights computed in dmpcd.build_sparse_weighted_rips_filtration

#### Input
- weights_filename - weight file outputted by dmpcd.build_sparse_weighted_rips_filtration
- verts_filename - original PCD dataset file for dmpcd.build_sparse_weighted_rips_filtration
- anno_filename - annotations of individual points in PCD
- output_vert_filename - PCD dataset with points listed in ascending weight order
- output_anno_filename - PCD dataset annotations with labels listed in ascending weight order


#### Output

reordered PCD dataset based on weight values

#### Example

    >import dmpcd as dm
    >import dmpcd.pcd as pcd

    >input_filename = "data/1-circle/features.txt"
    >output_dir = "results/1-circle-pcd/"
    >k = 15
    >metric = 'euclidean'
    >epsiilon = .99
    >persistence_threshold = .25

    >dm.pcd.build_sparse_weighted_rips_filtration(input_filename, output_dir, k, metric, epsiilon)
    >filtration_filename = os.path.join(output_dir, 'sparse_weighted_rips_filtration.txt')
    >weights_filename = os.path.join(output_dir, 'weights.txt')

    >dm.pcd.compute_persistence_swr(filtration_filename, output_dir)
    >edge_filename = os.path.join(output_dir, "edge_for_morse_only.txt")

    >sorted_weights_filename = os.path.join(output_dir, "sorted-weights.txt")
    >dm.pcd.reorder_weights(weights_filename, sorted_weights_filename)

    >morse_dir = os.path.join(output_dir, str(persistence_threshold) + '/')
    >dm.pcd.compute_graph_reconstruction(sorted_weights_filename, edge_filename, persistence_threshold, morse_dir)

    >result_edge_filename = os.path.join(morse_dir, 'dimo_edge.txt')
    >sorted_feature_filename = os.path.join(output_dir, 'sorted-feature.txt')
    >sorted_anno_filename = os.path.join(output_dir, 'sorted-anno.txt')
    >dm.pcd.reorder_verts_and_annos_by_weight(weights_filename, input_filename, "/path/to/anno.txt", sorted_feature_filename, output_anno_filename)

 
## dmpcd.baseline Functions

### dmpcd.baseline.density_estimation(input_filename, output_filename, k=15)

#### Description
Guassian kernel density estimation for each point in domain

#### Input
- input_filename - verts in PCD
- output_filename - file containing density estimation fo each point in PCD
- k - density factors in k-nearest neighbors

#### Output

a file containing the density estimation of each point in PCD

#### Example

    >import dmpcd as dm
    >import dmpcd.baseline as baseline
    
    >input_filename = "data/1-circle/features.txt"
    >output_dir = "results/1-circle-baseline/"
    >k = 15
    
    >density_filename = os.path.join(output_dir, 'density.txt', k)
    >dm.baseline.density_estimation(input_filename, density_filename)

### dmpcd.baseline.compute_rips_complex_edges(input_filename, output_filename, thresh, metric='euclidean')

#### Description
compute edges in Rips complex for given threshold

#### Input
- input_filename - verts in PCD
- output_filename - file containing edges in RIPS complex
- thresh - Radius parameter for Rips complex
- metric - metric used to compute Rips complex

#### Output

a file containing the edges in Rips complex

#### Example

    >import dmpcd as dm
    >import dmpcd.baseline as baseline
    
    >input_filename = "data/1-circle/features.txt"
    >output_dir = "results/1-circle-baseline/"
    >rips_alpha = .25
    
    >rips_edge_filename = os.path.join(output_dir, 'rips-edge.txt')
    >dm.baseline.compute_rips_complex_edges(input_filename, rips_edge_filename, rips_alpha, metric='euclidean')

### dmpcd.baseline.build_baseline_complex(rips_edge_filename, density_filename, output_filename, threshold=inf):

#### Description
build Rips complex with density values (.bin)

#### Input
- rips_edge_filename - file containing edges in Rips complex
- density_filename - file containing density of each point in PCD
- output filename - Complex (.bin) that can be used as input for persistence computation
- threshold - Rips cutoff, if rips_edge_filename was generated with cutoff this does not need to be used again

#### Output

a file containing the edges in Rips complex

#### Example

    >import dmpcd as dm
    >import dmpcd.baseline as baseline
    
    >input_filename = "data/1-circle/features.txt"
    >output_dir = "results/1-circle-baseline/"
    >k = 15
    >rips_alpha = .25
    >persistence_threshold = 998

    >if not os.path.exists(output_dir):
    >    os.mkdir(output_dir)

    >density_filename = os.path.join(output_dir, 'density.txt')
    >dm.baseline.density_estimation(input_filename, density_filename, k)

    >rips_edge_filename = os.path.join(output_dir, 'rips-edge.txt')
    >dm.baseline.compute_rips_complex_edges(input_filename, rips_edge_filename, rips_alpha, metric='euclidean')

    >complex_filename = os.path.join(output_dir, 'complex.bin')
    >dm.baseline.build_baseline_complex(rips_edge_filename, density_filename, complex_filename)


### dmpcd.baseline.compute_persistence_baseline(input_filename, output_dir)

#### Description
build Rips complex with density values (.bin)

#### Input
- input_filename - file where filtration is written
- output_dir - where persistence values of edges in filtration are written

#### Output

a file containing the persistence ifno of the edges in domain

#### Example

    >import dmpcd as dm
    >import dmpcd.baseline as baseline
    
    >input_filename = "data/1-circle/features.txt"
    >output_dir = "results/1-circle-baseline/"
    >k = 15
    >rips_alpha = .25
    >persistence_threshold = 998

    >if not os.path.exists(output_dir):
    >    os.mkdir(output_dir)

    >density_filename = os.path.join(output_dir, 'density.txt')
    >dm.baseline.density_estimation(input_filename, density_filename, k)

    >rips_edge_filename = os.path.join(output_dir, 'rips-edge.txt')
    >dm.baseline.compute_rips_complex_edges(input_filename, rips_edge_filename, rips_alpha, metric='euclidean')

    >complex_filename = os.path.join(output_dir, 'complex.bin')
    >dm.baseline.build_baseline_complex(rips_edge_filename, density_filename, complex_filename)

    >dm.baseline.compute_persistence_baseline(complex_filename, output_dir)
    >edge_persistence_filename = os.path.join(output_dir, 'edge_for_morse_only.txt')
    

### dmpcd.baseline.compute_graph_reconstruction_baseline(density_filename, edge_persistence_filename, persistence_threshold, output_dir)

#### Description
compute DM graph reconstruction

#### Input
- density_filename - file containing density of each point in PCD 
- edge_persistence_filename - file containing persistence values of edges in filtration 
- persistence_threshold - threshold used by algorithm to simplify output graph
- output_dir - directory output graph (.txt) is stored

#### Output

a file containing the persistence ifno of the edges in domain

#### Example

    >import dmpcd as dm
    >import dmpcd.baseline as baseline
    
    >input_filename = "data/1-circle/features.txt"
    >output_dir = "results/1-circle-baseline/"
    >k = 15
    >rips_alpha = .25
    >persistence_threshold = 998

    >if not os.path.exists(output_dir):
    >    os.mkdir(output_dir)

    >density_filename = os.path.join(output_dir, 'density.txt')
    >dm.baseline.density_estimation(input_filename, density_filename, k)

    >rips_edge_filename = os.path.join(output_dir, 'rips-edge.txt')
    >dm.baseline.compute_rips_complex_edges(input_filename, rips_edge_filename, rips_alpha, metric='euclidean')

    >complex_filename = os.path.join(output_dir, 'complex.bin')
    >dm.baseline.build_baseline_complex(rips_edge_filename, density_filename, complex_filename)

    >dm.baseline.compute_persistence_baseline(complex_filename, output_dir)
    >edge_persistence_filename = os.path.join(output_dir, 'edge_for_morse_only.txt')

    >morse_dir = os.path.join(output_dir, str(persistence_threshold) + '/')
    >if not os.path.exists(morse_dir):
    >    os.mkdir(morse_dir)
    >dm.baseline.compute_graph_reconstruction_baseline(density_filename, edge_persistence_filename, persistence_threshold, morse_dir)
    >morse_edge_filename = os.path.join(morse_dir, 'dimo_edge.txt')
    
## Other dmpcd Functions
    
### dmpcd.visualize_results_2d(vert_filename, edge_filename)

#### Description
visualize DM graph reconstruction on 2D embedding of points

#### Input
- vert_filename - verts in PCD (needs to be sorted by weight)
- edge_filename - edges of DM graph

#### Output

image of morse graph on top of 2D embedding of the dataset

#### Example

    >import dmpcd as dm
    >import dmpcd.pcd as pcd
    >import dmpcd.baseline as baseline

    >input_filename = "data/1-circle/features.txt"
    >output_dir = "results/1-circle-pcd/"
    >k = 15
    >metric = 'euclidean'
    >epsiilon = .99
    >persistence_threshold = .25

    >dm.pcd.build_sparse_weighted_rips_filtration(input_filename, output_dir, k, metric, epsiilon)
    >filtration_filename = os.path.join(output_dir, 'sparse_weighted_rips_filtration.txt')
    >weights_filename = os.path.join(output_dir, 'weights.txt')

    >dm.pcd.compute_persistence_swr(filtration_filename, output_dir)
    >edge_filename = os.path.join(output_dir, "edge_for_morse_only.txt")

    >sorted_weights_filename = os.path.join(output_dir, "sorted-weights.txt")
    >dm.pcd.reorder_weights(weights_filename, sorted_weights_filename)

    >morse_dir = os.path.join(output_dir, str(persistence_threshold) + '/')
    >dm.pcd.compute_graph_reconstruction(sorted_weights_filename, edge_filename, persistence_threshold, morse_dir)

    >result_edge_filename = os.path.join(morse_dir, 'dimo_edge.txt')
    >sorted_feature_filename = os.path.join(output_dir, 'sorted-feature.txt')
    >dm.pcd.reorder_verts_by_weight(weights_filename, input_filename, sorted_feature_filename)
    
    >dm.visualize_results_2d(sorted_feature_filename, result_edge_filename)

## Separate Programs

### Sparse Weighted Rips Filtration Persistence Module (PCD) (./dmpcd/code/persistence_swr/spt_cpp

#### Description

Compute persistence diagram of sparse weighted Rips filtration.  This is a modified version of the code found at (https://github.com/wangjiayuan007/graph_recon_DM).

#### Python Function

dmpcd.compute_persistence_swr

#### Input
- input_filename - path to sparse weighted Rips filtration file
- output_dir - directory where persistence results will be written to

#### Output:

output_dir/edge_for_morse_only.txt - contains all persistence information needed to run graph reconstruction

### Lower Star Filtration Persistence Module (Baseline) (./dmpcd/code/persistence_baseline/spt_cpp)

#### Description

Compute persistence diagram lower star filtration of given simplicial complex. This is a modified version of the code found at (https://github.com/wangjiayuan007/graph_recon_DM).

#### Python Function

dmpcd.compute_persistence_baseline

#### Input
- input_filename - path to complex file (.bin)
- output_dir - directory where persistence results will be written to

#### Output:

output_dir/edge_for_morse_only.txt - contains all persistence information needed to run graph reconstruction


### Discrete Morse Graph Reconstruction Module (PCD) (./dmpcd/code/phat_morse_src/a.out)

#### Description

Compute discrete Morse graph reconstruction of PCD.

#### Python Function

dmpcd.compute_graph_reconstruction

#### Input
- weights_filename - a sorted list of weights.
- input_edge_filename - persistence information for each edge in domain
- persistence_threshold - threshold used by the algorithm to simplify output
- output_dir - directory where file containing edges that are part of graph is written to.

#### Output:

output_dir/dimo_edge.txt - DM graph edges

### Discrete Morse Graph Reconstruction Module (Baseline) (./dmpcd/code/baseline_morse_src/a.out)

#### Description

Compute discrete Morse graph reconstruction of PCD.

#### Python Function

dmpcd.compute_graph_reconstruction_baseline

#### Input
- density_filename - a list of densities for each point in the dataset
- input_edge_filename - persistence information for each edge in domain
- persistence_threshold - threshold used by the algorithm to simplify output
- output_dir - directory where file containing edges that are part of graph is written to.

#### Output:

output_dir/dimo_edge.txt - DM graph edges


## Example Use of Pipeline

### PCD

    >import dmpcd as dm

    >input_filename = "data/1-circle/features.txt"
    >output_dir = "results/1-circle-pcd/"
    >k = 15
    >metric = 'euclidean'
    >epsiilon = .99
    >persistence_threshold = .25

    >dm.build_sparse_weighted_rips_filtration(input_filename, output_dir, k, metric, epsiilon)
    >filtration_filename = os.path.join(output_dir, 'sparse_weighted_rips_filtration.txt')
    >weights_filename = os.path.join(output_dir, 'weights.txt')

    >dm.compute_persistence_swr(filtration_filename, output_dir)
    >edge_filename = os.path.join(output_dir, "edge_for_morse_only.txt")

    >sorted_weights_filename = os.path.join(output_dir, "sorted-weights.txt")
    >dm.reorder_weights(weights_filename, sorted_weights_filename)

    >morse_dir = os.path.join(output_dir, str(persistence_threshold) + '/')
    >dm.compute_graph_reconstruction(sorted_weights_filename, edge_filename, persistence_threshold, morse_dir)

    >result_edge_filename = os.path.join(morse_dir, 'dimo_edge.txt')
    >sorted_feature_filename = os.path.join(output_dir, 'sorted-feature.txt')
    >dm.reorder_verts_by_weight(weights_filename, input_filename, sorted_feature_filename)
    
    >dm.visualize_results_2d(sorted_feature_filename, result_edge_filename)

### Baseline

    >import dmpcd as dm
    
    >input_filename = "data/1-circle/features.txt"
    >output_dir = "results/1-circle-baseline/"
    >k = 15
    >rips_alpha = .25
    >persistence_threshold = 998

    >if not os.path.exists(output_dir):
    >    os.mkdir(output_dir)

    >density_filename = os.path.join(output_dir, 'density.txt')
    >dm.density_estimation(input_filename, density_filename, k)

    >rips_edge_filename = os.path.join(output_dir, 'rips-edge.txt')
    >dm.compute_rips_complex_edges(input_filename, rips_edge_filename, rips_alpha, metric='euclidean')

    >complex_filename = os.path.join(output_dir, 'complex.bin')
    >dm.build_baseline_complex(rips_edge_filename, density_filename, complex_filename)

    >dm.compute_persistence_baseline(complex_filename, output_dir)
    >edge_persistence_filename = os.path.join(output_dir, 'edge_for_morse_only.txt')

    >morse_dir = os.path.join(output_dir, str(persistence_threshold) + '/')
    >if not os.path.exists(morse_dir):
    >    os.mkdir(morse_dir)
    >dm.compute_graph_reconstruction_baseline(density_filename, edge_persistence_filename, persistence_threshold, morse_dir)
    >morse_edge_filename = os.path.join(morse_dir, 'dimo_edge.txt')

    >dm.visualize_results_2d(input_filename, morse_edge_filename)


