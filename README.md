# LMM:Local Microenvironment Map
The goal of LMM is to cluster multiplexed images data based on local microenvironment organization of cells of interest, for example, a subset of cells that expressed one or some key markers and/or cells from a specific type.

You can get feature matrix by:

`> ExtractFeatures.R`

and for clustering in cell level:

`> ClusteringCells.R`

and for clustering in sample level:

`> ClusteringSamples.R`

An example with downstream analysis can be found in Codes and DownstreamAnalysis.R file. All returned pluts and data can be found in Result directory. The example is for Keren et al. paper:

Keren L, Bosse M, Marquez D, Angoshtari R, Jain S, Varma S, Yang SR, Kurian A, Van Valen D, West R, Bendall SC, Angelo M. A Structured Tumor-Immune Microenvironment in Triple Negative Breast Cancer Revealed by Multiplexed Ion Beam Imaging. Cell. 2018 Sep 6;174(6):1373-1387.e19. doi: 10.1016/j.cell.2018.08.039. PMID: 30193111; PMCID: PMC6132072.

