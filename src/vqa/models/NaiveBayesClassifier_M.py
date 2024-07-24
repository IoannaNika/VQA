from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB, MultinomialNB
import numpy as np
from sklearn.utils import shuffle
import json
from sklearn.metrics import accuracy_score
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import f1_score


ref_set_path = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/tuples_pacbio_sars_cov_2_rev_compl/dataset/samples_mutations.tsv"

reference_set = pd.read_csv(ref_set_path, sep='\t', header=0)

pos = list(reference_set[reference_set["label"] == "positive"]["mutations"])
neg = list(reference_set[reference_set["label"] != "positive"]["mutations"])

# plot distribution of edit distances for the positive and negative samples
# plot them on the same plot
fig = plt.figure()
ax = fig.add_subplot(111)
n1, bins1, _ = ax.hist(pos, bins=100, alpha=0.5, label='positive')
n2, bins2, _ = ax.hist(neg, bins=100, alpha=0.5, label='negative')
ax.set_xlabel('Number of mutations')
ax.set_ylabel('Frequency')
ax.legend(loc='upper right')
# save the plot
fig.savefig('Experiments/lumc_baselines/mutations/mutations_distribution.png')

x_pos = [np.asarray([i]) for i in pos]
y_pos = np.ones(len(pos)).tolist()

x_neg = [np.asarray([i]) for i in neg]
y_neg = np.zeros(len(neg)).tolist()

x = x_pos + x_neg
y = y_pos + y_neg

X = np.asarray(x)
y = np.asarray(y)

X, y = shuffle(X, y, random_state=0)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

nb = GaussianNB()
nb.fit(X_train, y_train)
y_pred = nb.predict(X_test)
print("Accuracy: ", accuracy_score(y_test, y_pred))


# save the model
filename = 'Experiments/lumc_baselines/mutations/model.pkl'

pickle.dump(nb, open(filename, 'wb'))
