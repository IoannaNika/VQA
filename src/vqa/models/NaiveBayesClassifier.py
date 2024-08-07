from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
import numpy as np
from sklearn.utils import shuffle
import json
from sklearn.metrics import accuracy_score
import pickle
import pandas as pd
import matplotlib.pyplot as plt

ref_set_path = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/tuples_pacbio_sars_cov_2_rev_compl/dataset/samples.tsv"

reference_set = pd.read_csv(ref_set_path, sep='\t', header=0)

pos = list(reference_set[reference_set["label"] == "positive"]["edit_distance"])
neg = list(reference_set[(reference_set["label"] == "negative") | (reference_set["label"] == "hard negative")]["edit_distance"])


fig = plt.figure()
ax = fig.add_subplot(111)
n1, bins1, _ = ax.hist(pos, bins=100, alpha=0.5, label='positive')
n2, bins2, _ = ax.hist(neg, bins=100, alpha=0.5, label='negative')
ax.set_xlabel('Edit distance')
ax.set_ylabel('Frequency')
ax.legend(loc='upper right')
# save the plot
fig.savefig('Experiments/lumc_baselines/edit_distance/ed_distribution.png')


x_pos = [np.asarray([i]) for i in pos]
y_pos = np.ones(len(pos)).tolist()

x_neg = [np.asarray([i]) for i in neg]
y_neg = np.zeros(len(neg)).tolist()

x = x_pos + x_neg
y = y_pos + y_neg

X = np.asarray(x)
y = np.asarray(y)
print(len(X))

X, y = shuffle(X, y, random_state=0)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

gnb = GaussianNB()
y_pred = gnb.fit(X_train, y_train).predict(X_test)

print("Accuracy: ", accuracy_score(y_test, y_pred))


# save the model
filename = 'Experiments/lumc_baselines/edit_distance/model.pkl'

pickle.dump(gnb, open(filename, 'wb'))
