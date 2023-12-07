# cluster embeddings using dbscan from sckit-learn
import sklearn
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.model_selection import train_test_split
import numpy as np

def cluster_embeddings_dbscan(embeddings, labels, eps=0.3, min_samples=2):
    """
    Cluster embeddings using DBSCAN from scikit-learn
    :param embeddings: embeddings to cluster
    :param labels: labels for each embedding
    :param eps: epsilon parameter for DBSCAN
    :param min_samples: min_samples parameter for DBSCAN
    :return: labels for each embedding
    """
    embeddings = embeddings.cpu().numpy()
    print("Clustering")
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(embeddings)
    print("Done clustering")
    predicted_labels = db.labels_
    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(predicted_labels)) - (1 if -1 in predicted_labels else 0)
    print("True labels", set(labels))
    print("Predicted labels", set(predicted_labels))
    homogeneity = metrics.homogeneity_score(predicted_labels, labels)
    completeness = metrics.completeness_score(predicted_labels, labels)
    
    return db.labels_, n_clusters_, homogeneity, completeness

def cluster_embeddings_kmeans(k, embeddings, labels):
    """
    Cluster embeddings using KMeans from scikit-learn
    :param embeddings: embeddings to cluster
    :param labels: labels for each embedding
    :param eps: epsilon parameter for DBSCAN
    :param min_samples: min_samples parameter for DBSCAN
    :return: labels for each embedding
    """
    embeddings = embeddings.cpu().numpy()
    print("Clustering")
    db = KMeans(n_clusters=k).fit(embeddings)
    print("Done clustering")
    predicted_labels = db.labels_
    # log number of samples in each cluster
    print("Number of samples in each cluster:")
    print(np.bincount(predicted_labels))
    
    # group true labels by predicted cluster
    true_labels_by_cluster = {}
    for i in range(len(predicted_labels)):
        if predicted_labels[i] not in true_labels_by_cluster:
            true_labels_by_cluster[predicted_labels[i]] = []
        true_labels_by_cluster[predicted_labels[i]].append(labels[i]) 
    
    # print the true labels in each cluster
    for cluster in true_labels_by_cluster:
        print("Cluster " + str(cluster) + ":")
        print(set(true_labels_by_cluster[cluster]))
        # print the number of times each label appears in the cluster
        for lbl in set(true_labels_by_cluster[cluster]):
            print(lbl + ": " + str(true_labels_by_cluster[cluster].count(lbl)))
        print("Done")
    
    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(predicted_labels)) - (1 if -1 in predicted_labels else 0)

    homogeneity = metrics.homogeneity_score(predicted_labels, labels)
    completeness = metrics.completeness_score(predicted_labels, labels)
    
    return db.labels_, n_clusters_, homogeneity, completeness
