Hierarchical clustering by correlation

This script creates gene clusters using the AgglomerativeClustering function from the scikit-learn package,
based on the Pearson correlation between the gene phylogenetic profiles (the rows in the input NPP matrix).

Installation and running:
Use Python 3+ (Developed with Python 3.8.5)
Install package scikit-learn (Developed with scikit-learn 1.0.2)

To use the code for creating clusters from NPP files:

1. Download "clusterNPP.py" and "upgrade_clean_data.py" in home directory
2.  Create in home directory folder "NPP" and locate the NPP file(s) in this directory.
3. Create in home directory folder "Hierarchical_by_corr_clusters" .
4. "clusterNPP" saves clustering results in this folder using similar name
 (including string with the cutting range of the dendrogram) and with csv suffix.
5. In  "clusterNPP.py", verify that fname1 [,fname2,fname3] contain the appropriate NPP filename(s),
  located in the "NPP" folder
6. To create the Clusters csv file run "clusterNPP.py"

 
Pseudo code
Hierarchical_clustering_by_correlation
- Input: NPP matrix file
- Remove rows with 0 std from the NPP matrix
- for each threshold distance T in the range 0.05 - 1.05 (with steps of 0.05):
  - Create a model of "AgglomerativeClustering" for the distance T with linkage='complete', affinity='correlation'
    (and n_clusters=None)
  - Fit the model for the NPP matrix
  - Calculate the Silhouette Coefficient values
  - Keep only clusters with size between 5 and 30 genes
  - For each cluster, calculate the mean and the std of the correlations between all cluster gene pairs
  - Each cluster is saved in the clusters list only if

The result is a table of clusters with the following columns:
        cluster_label - The label of the cluster within all clusters produced from the same distance threshold
        distance_threshold  - distance threshold
        model - The parameters used for the AgglomerativeClustering function
        num_samples_in_cluster -number of genes in cluster
        listGenesInCluster - Cluster gene names
        silhouette_avg - Avergae of the Silhouette Coefficient values of all clusters in the same distance
        CrossCorrMean - Average of the NPP correlations between all cluster gene pairs
        CrossCorrStd - Standard deviation of the NPP correlations between all cluster gene pairs
        silhoutte_sampls_values - Silhouette Coefficient values of all cluster genes
        silhoutte_sampls_avg - The average of the Silhouette Coefficient values of all cluster genes


Silhouette explanation:
The Silhouette Coefficient is calculated using the mean intra-cluster distance (a) and the mean nearest-cluster distance (b) for each sample. The Silhouette Coefficient for a sample is (b - a) / max(a, b). To clarify, b is the distance between a sample and the nearest cluster that the sample is not a part of. Note that Silhouette Coefficient is only defined if number of labels is 2 <= n_labels <= n_samples - 1.
