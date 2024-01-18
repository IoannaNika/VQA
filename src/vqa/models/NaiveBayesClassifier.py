from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
import numpy as np
from sklearn.utils import shuffle
import json
from sklearn.metrics import accuracy_score


f= open("pos")
pos = json.load(f)
x_pos = [np.asarray([i]) for i in pos]
y_pos = np.ones(len(pos)).tolist()

f= open("neg")
neg = json.load(f)
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