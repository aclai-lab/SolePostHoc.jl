using Clustering

# make a random dataset with 1000 random 5-dimensional points
X = rand(5, 1000)

# cluster X into 20 clusters using K-means
R = kmeans(X, 20; maxiter = 200, display = :iter)

@assert nclusters(R) == 20 # verify the number of clusters

a = assignments(R) # get the assignments of points to clusters
c = counts(R) # get the cluster sizes
M = R.centers # get the cluster centers
