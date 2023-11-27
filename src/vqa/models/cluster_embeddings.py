# cluster embeddings using dbscan from sckit-learn
import sklearn
from sklearn.cluster import DBSCAN
from sklearn import metrics


def cluster_embeddings(embeddings, labels, eps=0.5, min_samples=2):
    """
    Cluster embeddings using DBSCAN from scikit-learn
    :param embeddings: embeddings to cluster
    :param labels: labels for each embedding
    :param eps: epsilon parameter for DBSCAN
    :param min_samples: min_samples parameter for DBSCAN
    :return: labels for each embedding
    """
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(embeddings)
    # calculate number of clusters
    n_clusters_ = len(set(db.labels_)) - (1 if -1 in db.labels_ else 0) 

    # homogeneous, completeness
    homogeneity = metrics.homogeneity_score(labels, db.labels_)
    completeness = metrics.completeness_score(labels, db.labels_)

    return db.labels_, n_clusters_, homogeneity, completeness
