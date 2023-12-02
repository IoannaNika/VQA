# cluster embeddings using dbscan from sckit-learn
import sklearn
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.model_selection import train_test_split


def cluster_embeddings(embeddings, labels, eps=0.3, min_samples=2):
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
    # split embeddings to train and test
    embeddings_train, embeddings_test, labels_train, labels_test = train_test_split(embeddings, labels, test_size=0.2, random_state=42)
    # train DBSCAN
    print("Training DBSCAN...")
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(embeddings_train, labels_train)
    # predict labels for test embeddings
    print("Training done")
    print("Predicting labels...")
    db.labels_ = db.fit_predict(embeddings_test)
    print("Predicting done")
    # calculate number of clusters
    n_clusters_ = len(set(db.labels_)) - (1 if -1 in db.labels_ else 0) 

    # homogeneous, completeness
    homogeneity = metrics.homogeneity_score(labels_test, db.labels_)
    completeness = metrics.completeness_score(labels_test, db.labels_)
    print("Calculated metrics")
    return db.labels_, n_clusters_, homogeneity, completeness
