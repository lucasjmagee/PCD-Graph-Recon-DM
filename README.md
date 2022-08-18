# Point Cloud Discrete Morse Graph Reconstruction

dmpcd is a python package meant for executing the discrete Morse graph reconstruction algorithm on Point Cloud Datasets - designed with full mouse brain imaging data in mind.  The package includes functions which allow the user to perform the baseline DM approach and the more recently developed PCD DM approach.  For the basline approach, the package allows users to compute the Rips filtration, estimate density at each point, compute persistence of lower star filtration, and perform discrete Morse graph reconstruction.  For the PCD DM approach, users can compute the sparse weighted Rips filtration, compute persistence of the filtration, and perform discrete Morse graph reconstruction.  Data can be labeled or unlabeled, with visualization tools included for both cases.

* [Installation Intructions](#installation-instructions)
  * [System Requirements](#system-requirements)
  * [Required Python Libraries](#required-python-libraries)
  * [Compiling Code](#compiling-code)
* [dmpcd Functions](#dmpcd-functions)
* [Separate Programs](#separate-programs)
* [Example Use of Pipeline](#example-use-of-pipeline) 

## Installation Instructions
### System Requirements
- Python 3.8.8 (or newer)
- g++ 9.4.0 (or newer)

### Required Python Libraries
- scitkit-learn - pip install scikit-learn (https://scikit-learn.org/)

### Compiling Code

Sparse Weighted Rips Filtration Persistence Module
    
    > cd DiMo3d/code/graph_recon_RIPS_v2/spt_cpp/
    > g++ DiMoSC.cpp -I./phat/include -std=c++11 -o spt_cpp

Lower Star Filtration Persistence Module

    > cd DiMo3d/code/spt_cpp/
    > g++ DiMoSC.cpp -I./phat/include -std=c++11 -o spt_cpp
    
Discrete Morse Graph Reconstruction

## dmpcd Functions

### dmpcd.build_sparse_weighted_rips_filtration(feature_filename, output_dir, k=15, metric='euclidean', epsilon=.99, cutoff=inf)

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

    >feature_filename = “data/1-circle/features.txt”
    >output_dir = “results/1-circle/”
    >k=15
    >metric='euclidean'
    >epsilon=.99
    >cutoff=inf
    >dm.build_sparse_weighted_rips_filtration(feature_filename, output_dir, k, metric, epsilon, cutoff)
    
### dmpcd.compute_persistence_swr(input_filename, output_dir)

#### Description
Compute persistence of sparse weighted Rips filtrations

#### Input
- input_filename - file where filtration is written
- output_dir - where persistence values of edges in filtration are written

#### Output

File containing persistence values and dimension of each edge in input filtration.  This is later used for DM graph reconstruction.

#### Example

    >import dmpcd as dm

    >feature_filename = “data/1-circle/features.txt”
    >output_dir = “results/1-circle/”
    >k=15
    >metric='euclidean'
    >epsilon=.99
    >cutoff=inf
    >dm.build_sparse_weighted_rips_filtration(feature_filename, output_dir, k, metric, epsilon, cutoff)
    >dm.compute_persistence_swr(os.path.join(output_dir, "sparse_weighted_rips_filtration.txt"), output_dir)


