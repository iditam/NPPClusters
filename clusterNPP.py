# -*- coding: utf-8 -*-

#  4 Use a subset of only 3 classes from CIFAR10 (32x32 RGB images)

# Explore the database using unsupervised tools
# Try to cluster into 3 different clusters
# There is ground truth, but don’t use it for the fit, only for the evaluation and data exploration! (e.g. you can use Rand index)
# No need to split into Train & Test

#!pip install scikit-learn-extra
import sys
import numpy as np
import pandas as pd
import os
from pathlib import Path
import pickle
import statistics
from scipy.stats import hypergeom,ttest_1samp

import matplotlib
from matplotlib.pyplot import plot
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_samples, silhouette_score

min_genes_per_cluster=5
max_genes_per_cluster=20

# from  upgrade_clean_data import clean_df

param_grid = {
    'n_clusters': range (2,30)
}

###################save & restore pkl #####################
def save_ws_pkl(fname, listofOBJ):
    with open(fname, 'wb') as f:
        pickle.dump(listofOBJ, f)

def restore_ws_pkl(fname):
    with open(fname, 'rb') as f:
        listofOBJ = pickle.load(f)
    return (listofOBJ)
###########################silhouette ###############################

def modellist_above_max_Silhouette(silhouette4num_clusters):
    idx=np.where(np.array(silhouette4num_clusters['silhouette_avg']==silhouette4num_clusters['silhouette_avg'].max()))[0][0]
    print('thr', silhouette4num_clusters.iloc[idx]['distance_threshold'])
    return list(silhouette4num_clusters.iloc[idx:-1]['model']) #idx:-1

def model4max_Silhouette(silhouette4num_clusters):
    idx=np.where(np.array(silhouette4num_clusters['silhouette_avg']==silhouette4num_clusters['silhouette_avg'].max()))[0][0]
    print('thr', silhouette4num_clusters.iloc[idx]['distance_threshold'])
    return silhouette4num_clusters.iloc[idx]['model']

def isClusteralreadyexsist(clusterslist, cluster2check):
    if len(clusterslist)>0:
        d1=clusterslist[clusterslist['num_samples_in_cluster'] == cluster2check.iloc[0]['num_samples_in_cluster']]
        if len(d1)>0:
            for i in range(len(d1)):
                if set(d1.iloc[i]['listGenesInCluster'])==set(cluster2check.iloc[0]['listGenesInCluster']):
                    return True
    return False

def estimate_n_clusters_using_silhouette(X,fprfix,dstpath):
    #X- dataframe
    #fprfix - calde name
    # dst path for estimatetion result
    # Description: estimate number of distance_thr in the case of hirracial clustring distance='correlation'
    #return model with the value estimated ditance_thr and num_clusters=nan
    X1=X.copy()
    forcepkl=1
    silhouettefname=fprfix +'_'+'silhouette4_num_clusters'

    affinity = 'correlation'
    linkage = 'complete' #single,complete
    min_distance=0.05 #0.01
    max_distance=1.10 #1.99 1.10
    distance_step=0.05 #0.01
    silhouette_fullpath=os.path.join (dstpath,silhouettefname+str(min_distance)+'-'+str(max_distance)+'.pkl')

    param_grid = [ {'linkage': ['single','complete','average']}]
    #silhouette4num_clusters= pd.DataFrame({'distance_threshold': [], 'n_clusters_': [],'model': [],'silhouette_avg': []})
    dfcluster=pd.DataFrame()
    if not os.path.isfile(silhouette_fullpath) or forcepkl:
        clusters_list=[]
        for distance_threshold in np.arange(min_distance, max_distance, distance_step): # best was 0.35
            print('distance_threshold', distance_threshold)
            clustring_model = AgglomerativeClustering(distance_threshold=distance_threshold, n_clusters=None,
                                                      linkage=linkage, affinity=affinity)
            clustring_model.fit(X)
            X1['label']=clustring_model.labels_
            #R=plot_dendrogram(clustring_model)
            #if clustering_model.n_clusters_ < last_num_clusters:
            if len(np.unique( clustring_model.labels_)) > 1:
                silhouette_avg = silhouette_score(X, clustring_model.labels_, metric='correlation')
                sample_silhouette_values = silhouette_samples(X, clustring_model.labels_,metric='correlation')
                X1['SampSil']=sample_silhouette_values
                #for i in range(clustring_model.n_clusters_):
                for lbl in np.unique(clustring_model.labels_):
                   print('lbl',lbl)
                   SingleCluster = X1[X1['label'] == lbl].sort_values(by='SampSil', ascending=False)
                   num_samples_in_cluster= SingleCluster.shape[0]
                   # skip clusters with more than max_genes_per_cluster
                   if num_samples_in_cluster >max_genes_per_cluster or  \
                       num_samples_in_cluster < min_genes_per_cluster: # skip on clusters with less than 2 items...
                       continue
                   listGenesInCluster=list(SingleCluster.index)
                   silhoutte_sampls_values=list(SingleCluster['SampSil'])
                   X2=SingleCluster.copy()
                   X2.drop(['label', 'SampSil'],axis=1,inplace=True)
                   XcrossCorr=np.corrcoef(X2)
                   uper_triangle=XcrossCorr[np.triu(m=XcrossCorr, k=1) > 0]
                   AvgCrossCorr=uper_triangle.mean()
                   corr_med = statistics.median(uper_triangle)
                   corr_std = 0
                   if len(XcrossCorr) - 1 > 1:
                       corr_std = statistics.stdev(uper_triangle)
                   cluster2check=pd.DataFrame.from_dict(
                                        {'cluster_label': lbl,
                                         'distance_threshold': distance_threshold,
                                         'model': clustring_model,
                                         'num_samples_in_cluster':num_samples_in_cluster,
                                         'listGenesInCluster':[listGenesInCluster],
                                         'silhouette_avg': silhouette_avg,
                                         'CrossCorrMean': AvgCrossCorr,
                                         'CrossCorrStd': corr_std,
                                         'silhoutte_sampls_values':[silhoutte_sampls_values],
                                         'silhoutte_sampls_avg':X1[X1['label']==lbl]['SampSil'].mean()})
                   if isClusteralreadyexsist(dfcluster, cluster2check):
                       continue
                   dfcluster=pd.concat([dfcluster,cluster2check ], ignore_index=True)
            else:
                silhouette_avg=-1 # single label , silhouette not possible
        save_ws_pkl(silhouette_fullpath, [dfcluster])
    else:
        dfcluster = restore_ws_pkl(silhouette_fullpath)[0]
    m = model4max_Silhouette(dfcluster)
    models_list=modellist_above_max_Silhouette(dfcluster)
    return m,models_list,dfcluster

########################### Hierarchical_clustring ###############################

def rm_rows_with_zero_std(X):
    rows2rm=[]
    for i in range(X.shape[0]):
        if np.std(X.iloc[i])<10**-8:
            rows2rm.append(X.index[i])
    rmlbls=list(set(list(X.index)) - set(rows2rm))
    if len(rows2rm)>0:
        X=X.loc[rmlbls]
    return  X

def Hierarchical_clustering_by_correlation(X,clade,ResultsPath):
    X=rm_rows_with_zero_std(X)
    clustring_model,models_list,dfcluster=estimate_n_clusters_using_silhouette(X, fprfix=clade, dstpath=ResultsPath)
    dfcluster.sort_values(by=['CrossCorrMean'], inplace=True, ascending=[False])
    return dfcluster ,ResultsPath

################################### main #############################
fname1='NPP_classic_no_medium_pg_no_low_score_genes.csv'
fname2='NPP_no_medium_pg_no_low_score_genes.csv'
fname3='Homo_sapiens_NPP_UniProt_062018_Filtered_7057_genes_fixed.csv'
ResultsPath='Hierarchical_by_corr_clusters'

for fname in [fname2]: #[fname3,fname2,fname1]:
    fname=os.path.join('NPP',fname)
    NPPfile=pd.read_csv(fname)
    NPPfile.rename(columns={"Unnamed: 0": "Gname"},inplace=True)
    list_of_clustring_alg=['Hierarchical_correlation'] # 'DBSCAN' 'Hierarchical_euclidean',
    ################## change nan values to newmin ##################
    NPPfile_new=NPPfile.copy()
    NPPfile=NPPfile.set_index('Gname')
    NPPfile=NPPfile.astype(float)
    #NPPfile=clean_df(X)
    ##########################################################
    if 'Hierarchical_correlation' in list_of_clustring_alg:
        clade_names=['Eukaryota']
        for clade in clade_names:
            #NPP_with_parasite = NPP_clade_only_with_parasite(NPPfile, [clade])
            if len(NPPfile)>0:
               #NPP_with_parasite=rm_low_score_genes(NPP_with_parasite,clade)
               dfcluster,dstpath=Hierarchical_clustering_by_correlation(NPPfile,clade,ResultsPath)
               distance_win = '_range_dist_' + str(round(dfcluster['distance_threshold'].min(), 2)) + '-' + str(
                   round(dfcluster['distance_threshold'].max(), 2))
               # add destination dir name for results
               orgfname=Path(os.path.normpath(fname).split(os.sep)[1]).stem
               dfcluster.to_csv(os.path.join(dstpath ,'Clusters_of_'+ orgfname+distance_win + '.csv'))
