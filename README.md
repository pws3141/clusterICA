# clusterICA
clusterICA: clustering method for ICA

An implementation of Independent Component Analysis, using clustering and optimisation of random projections of the data. It combines ordering of random points on the projective space and clustering of ≈the best points using divisive k-means until some tolerance condition is met, with an optimisation step for each cluster. Each loading is found sequentially, and Householder rotations are used to make the current loading orthogonal to all previous ones.

Requires:
> jvcoords: uses whiten function within Independent Component Analysis
