#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Example script, part of MDAnalysis
"""
Example: Comparing a trajectories from different methods
========================================================

Example implementation of Path Similarity Analysis that shows how to read in a
set of trajectories, compute (discrete) Fréchet distances, and plot a heat
map-dendrogram.

This example uses the apo AdK transition between its open and closed crystal
structures as a testbed system (see [Seyler2014]). Trajectories are generated
given the known structural endpoints. A selection of ten different sampling
methods (three transitions each), plus the path of linear interpolation, were
used to generate a total of 31 transitions closed to open transition paths. The
(discrete) Fréchet distances are computed (between each unique pair of
trajectories) and stored in a distance matrix. (See [Seyler2015] for further
applications of and information on PSA.

The distance matrix is stored in a data file `discrete_frechet.dat` and a numpy
file `discrete_frechet.npy`, and the heat map-dendrogram showing Ward
hierarchical clustering of the distance matrix is also written to
`psadata/plots/df_war_psa-short.pdf` (requires :mod:`matplotlib`).

[Seyler2014]    S.L. Seyler and O. Beckstein, Sampling large conformational
                transitions: adenylate kinase as a testing ground. Mol Simul 40
                (2014), 855–877, doi:10.1080/08927022.2014.919497

[Seyler2015]    S.L. Seyler, A. Kumar, M.F. Thorpe, and O. Beckstein, Path
                Similarity Analysis: a Method for Quantifying Macromolecular
                Pathways. `arXiv:1505.04807v1`_ [q-bio.QM], 2015.

.. SeeAlso:: :mod:`MDAnalysis.analysis.psa`

"""

from MDAnalysis import Universe
from MDAnalysis.analysis.psa import PSAnalysis
from pair_id import PairID

if __name__ == '__main__':

    print("Building collection of simulations...")
    method_names = ['DIMS', 'FRODA', 'GOdMD', 'MDdMD', 'rTMD-F', 'rTMD-S',      \
                    'ANMP', 'iENM', 'MAP', 'MENM-SD', 'MENM-SP',                \
                    'Morph', 'LinInt']
    labels = [] # Heat map labels
    simulations = [] # List of simulation topology/trajectory filename pairs
    universes = [] # List of MDAnalysis Universes representing simulations

    # Build list of simulations, each represented by a pair of filenames
    # ([topology filename], [trajectory filename]). Generate corresponding label
    # list.
    for method in method_names:
        # Note: DIMS uses the PSF topology format
        topname = 'top.psf' if 'DIMS' in method or 'TMD' in method else 'top.pdb'
        pathname = 'fitted_psa.dcd'
        method_dir = 'methods/{}'.format(method)
        if method is not 'LinInt':
            for run in xrange(1, 4): # 3 runs per method
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(method_dir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))
        else: # only one LinInt trajectory
            topology = '{}/{}'.format(method_dir, topname)
            trajectory = '{}/{}'.format(method_dir, pathname)
            labels.append(method)
            simulations.append((topology, trajectory))

    # Generate simulation list represented as Universes. Each item, sim, in
    # simulations is a topology/trajectory filename pair that is unpacked into
    # an argument list with the "splat" ("*") operator.
    for sim in simulations:
        universes.append(Universe(*sim))

    print("Initializing Path Similarity Analysis...")
    psa_hpa = PSAnalysis(universes, path_select='name CA', labels=labels)

    print("Generating Path objects from trajectories...")
    psa_hpa.generate_paths()

    print("Performing full Hausdorff pairs analysis for all pairs of paths...")
    psa_hpa.run_pairs_analysis(neighbors=True, hausdorff_pairs=True)

    #------------------------------------------------
    # Generate nearest neighbor distance plots
    #------------------------------------------------

    # Add three runs per method, except for LinInt (only has one)
    identifier = PairID()
    for name in method_names:
        run_ids = [1] if 'LinInt' in name else [1,2,3]
        identifier.add_sim(name, run_ids)

    # Get the PairID:
    #    The comparison between a pair of simulations is assigned a unique
    #    PairID. Given the order in which simulations are added to PSA, the
    #    comparison between a pair of simulations can be identified by
    #    (distance) matrix indices. The PairID is the index in the corresponding
    #    distance vector of a given pair of simulations.
    s1, s2, s3 = 'DIMS 1', 'DIMS 2', 'rTMD-F 3'
    pid1 = identifier.get_pair_id(s1, s2)
    pid2 = identifier.get_pair_id(s2, s3)

    print("Plotting nearest neighbors as a function of normalized progress (by" \
        + " frame for:")
    print("    1. comparison {:d}, {} to {}...".format(pid1, s1, s2))
    psa_hpa.plot_nearest_neighbors(filename='nn_dims1_dims2.pdf', idx=pid1,     \
                                   labels=(s1, s2))

    print("    2. comparison {:d}, {} to {}...".format(pid2, s2, s3))
    psa_hpa.plot_nearest_neighbors(filename='nn_dims2_tmds3.pdf', idx=pid2,     \
                                   labels=(s2, s3))
