# cluster embeddings using dbscan from sckit-learn
import sklearn
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn import metrics
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

def visualize_with_tsne(genomic_region, embeddings, labels, true_labels, perplexity=2, n_iter=1000):
    """
    Visualize embeddings using t-SNE from scikit-learn
    :param embeddings: embeddings to visualize
    :param labels: labels for each embedding
    :param perplexity: perplexity parameter for t-SNE
    :param n_iter: number of iterations for t-SNE
    """
    embeddings = np.asarray(embeddings, dtype=np.float32)
    print("Visualizing")
    pca = sklearn.decomposition.PCA(n_components=2) 
    tsne =  sklearn.manifold.TSNE(n_components=2, perplexity=perplexity, n_iter=n_iter)
    embeddings_2d = tsne.fit_transform(embeddings)
    print("Done visualizing")
    # plot the embeddings
    plt.figure(figsize=(6, 5))


    true_labels_set =list(set(true_labels))
    # assign colors to each unique label 
    colors = []
    for label in true_labels:
        colors.append(true_labels_set.index(label))

    plt.scatter(embeddings_2d[:, 0], embeddings_2d[:, 1], c=labels, cmap='jet')
    plt.colorbar()
    # save the plot
    plt.savefig("tsne_" + genomic_region + ".png")

    plt.scatter(embeddings_2d[:, 0], embeddings_2d[:, 1], c=colors, cmap='jet')
    plt.colorbar()
    # save the plot
    plt.savefig("true_tsne_" + genomic_region + ".png")
    return 

def cluster_embeddings_dbscan(embeddings, labels, eps=1, min_samples=2, genomic_region = "", produce_plots=True, verbose=True):
    """
    Cluster embeddings using DBSCAN from scikit-learn
    :param embeddings: embeddings to cluster
    :param labels: labels for each embedding
    :param eps: epsilon parameter for DBSCAN
    :param min_samples: min_samples parameter for DBSCAN
    :return: labels for each embedding
    """
    #embeddings = embeddings.cpu().numpy()
    if verbose:
        print("Clustering")

    embeddings_filtered = []
    labels_filtered = []

    # get rid of embeddings for which their label appears only once
    labels_counter = Counter(labels)
    
    for embedding, label in zip(embeddings, labels):
        if labels_counter[label] < 2:
            continue
        embeddings_filtered.append(embedding)
        labels_filtered.append(label)

    if verbose:
        print("Embedding length: ", len(embeddings[0]))
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(embeddings_filtered)
    if verbose:
        print("Done clustering")
    predicted_labels = db.labels_
    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(predicted_labels)) - (1 if -1 in predicted_labels else 0)
    if verbose:
        print("True labels", set(labels_filtered))
        print("Predicted labels", set(predicted_labels))
    homogeneity = metrics.homogeneity_score(labels_filtered, predicted_labels)
    completeness = metrics.completeness_score(labels_filtered, predicted_labels)
    if verbose:
        print("Done clustering")
    predicted_labels = db.labels_

    # group true labels by predicted cluster
    true_labels_by_cluster = {}
    for i in range(len(predicted_labels)):
        if predicted_labels[i] not in true_labels_by_cluster:
            true_labels_by_cluster[predicted_labels[i]] = []
        true_labels_by_cluster[predicted_labels[i]].append(labels_filtered[i]) 
    
    if verbose:
        # print the true labels in each cluster
        for cluster in true_labels_by_cluster:
            print("Cluster " + str(cluster) + ":")
            print(set(true_labels_by_cluster[cluster]))
            # print the number of times each label appears in the cluster
            for lbl in set(true_labels_by_cluster[cluster]):
                print(lbl + ": " + str(true_labels_by_cluster[cluster].count(lbl)))
            print("Done")
        
    # genomic_region = "_".join(labels_filtered[0].split("_")[1:])
    # plot the embeddings
    if produce_plots == True:
        visualize_with_tsne(genomic_region, embeddings_filtered, predicted_labels, labels_filtered)

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
    # embeddings = embeddings.cpu().numpy()
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

    homogeneity = metrics.homogeneity_score(labels, predicted_labels)
    completeness = metrics.completeness_score(labels, predicted_labels)
    
    return db.labels_, n_clusters_, homogeneity, completeness
