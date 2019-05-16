# Authors: Vlad Niculae, Mathieu Blondel
# License: BSD 3 clause
"""
=========================
Multilabel classification
=========================

This example simulates a multi-label document classification problem. The
dataset is generated randomly based on the following process:

    - pick the number of labels: n ~ Poisson(n_labels)
    - n times, choose a class c: c ~ Multinomial(theta)
    - pick the document length: k ~ Poisson(length)
    - k times, choose a word: w ~ Multinomial(theta_c)

In the above process, rejection sampling is used to make sure that n is more
than 2, and that the document length is never zero. Likewise, we reject classes
which have already been chosen.  The documents that are assigned to both
classes are plotted surrounded by two colored circles.

The classification is performed by projecting to the first two principal
components found by PCA and CCA for visualisation purposes, followed by using
the :class:`sklearn.multiclass.OneVsRestClassifier` metaclassifier using two
SVCs with linear kernels to learn a discriminative model for each class.
Note that PCA is used to perform an unsupervised dimensionality reduction,
while CCA is used to perform a supervised one.

Note: in the plot, "unlabeled samples" does not mean that we don't know the
labels (as in semi-supervised learning) but that the samples simply do *not*
have a label.
"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd
import seaborn as sn
import pickle

from sklearn.datasets import load_digits
from sklearn import model_selection
from sklearn.pipeline import Pipeline
from tempfile import mkdtemp
from shutil import rmtree
from sklearn.externals.joblib import Memory
#from sklearn.metrics import make_scorer
#from sklearn.metrics import accuracy_score
#from sklearn.metrics import confusion_matrix
from sklearn import metrics

#from sklearn.datasets import make_multilabel_classification
#from sklearn.multiclass import OneVsRestClassifier
#from sklearn.preprocessing import LabelBinarizer
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA, NMF
from sklearn.cross_decomposition import CCA
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis 
from sklearn import svm
from sklearn.externals import joblib


COLORS = np.array([
                  '#e6194b',
                  '#3cb44b',
                  '#ffe119',
                  '#0082c8',
                  '#f58231',
                  '#911eb4',
                  '#46f0f0',
                  '#f032e6',
                  '#d2f53c',
                  '#fabebe',
                  '#008080',
                  '#e6beff',
                  '#aa6e28',
                  '#fffac8',
                  '#800000',
                  '#aaffc3',
                  '#808000',
                  '#ffd8b1',
                  '#000080',
                  '#808080',
                  '#FFFFFF',
                  '#000000'
                 ])

#################################################
# read in data


set1 = pd.read_table('../result.rsem.TET.instance/set.TE.noovp.100K.txt', dtype={'patient': str} )
patient = set1["patient"]
X = set1.drop(columns=['patient'])
Xlabels = X.columns
X = X.values
#patient2label = pd.read_table('/DataDrives/dd4/xiaogang/patientToCancerMappings.tsv')
#patient2label[['TCGA','center', 'patient']] = patient2label['barcode'].str.split('-',expand=True)
#patient2label = patient2label.set_index(['patient'])
patient2label = patient.str.split('.',expand=True)
label = patient2label[1]
Y = pd.get_dummies(label)
Yidx = np.where(Y!=0)[1]
Ylabelidx = np.where(Y.drop_duplicates()!=0)[1]
Ylabels = dict(zip(Ylabelidx, Y.columns[Ylabelidx]))
print (Ylabels)
print (Yidx)
Y = Y.values
# dividing X, y into train and test data
X_train, X_test, y_train, y_test = model_selection.train_test_split(X, Yidx, random_state = 0)
#print (X_train)
#print (y_train)
print (X_train.shape)
print (y_train.shape)




################################################
# plot the decision surface for SVM classifiers 


def make_meshgrid(x, y, h=.02):
    """Create a mesh of points to plot in

    Parameters
    ----------
    x: data to base x-axis meshgrid on
    y: data to base y-axis meshgrid on
    h: stepsize for meshgrid, optional

    Returns
    -------
    xx, yy : ndarray
    """
    x_min, x_max = x.min() - 1, x.max() + 1
    y_min, y_max = y.min() - 1, y.max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                         np.arange(y_min, y_max, h))
    return xx, yy


def plot_contours(ax, clf, xx, yy, **params):
    """Plot the decision boundaries for a classifier.

    Parameters
    ----------
    ax: matplotlib axes object
    clf: a classifier
    xx: meshgrid ndarray
    yy: meshgrid ndarray
    params: dictionary of params to pass to contourf, optional
    """
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    out = ax.contourf(xx, yy, Z, **params)
    return out


#################################################
## plot LDA results
#def plot_scikit_lda(X, title):
#    ax = plt.subplot(111)
#    for label,marker,color in zip(
#        range(1,4),('^', 's', 'o'),('blue', 'red', 'green')):
#
#        plt.scatter(x=X[:,0][y == label],
#                    y=X[:,1][y == label] * -1, # flip the figure
#                    marker=marker,
#                    color=color,
#                    alpha=0.5,
#                    label=label_dict[label])
#
#    plt.xlabel('LD1')
#    plt.ylabel('LD2')
#
#    leg = plt.legend(loc='upper right', fancybox=True)
#    leg.get_frame().set_alpha(0.5)
#    plt.title(title)
#
#    # hide axis ticks
#    plt.tick_params(axis="both", which="both", bottom="off", top="off",  
#            labelbottom="on", left="off", right="off", labelleft="on")
#
#    # remove axis spines
#    ax.spines["top"].set_visible(False)  
#    ax.spines["right"].set_visible(False)
#    ax.spines["bottom"].set_visible(False)
#    ax.spines["left"].set_visible(False)    
#
#    plt.grid()
#    plt.tight_layout
#    plt.show()
#
#
#



####################################################
## dividing X, y into train and test data
#X_train, X_test, y_train, y_test = model_selection.train_test_split(X, Yidx, random_state = 0)
##print (X_train)
##print (y_train)
#print (X_train.shape)
#print (y_train.shape)
#
## decision tree
## training a DescisionTreeClassifier
#dtree_model = DecisionTreeClassifier(max_depth = 2).fit(X_train, y_train)
#dtree_predictions = dtree_model.predict(X_test)
# 
#dt_accuracy = dtree_model.score(X_test, y_test)
#print (dt_accuracy)
#
## creating a confusion matrix
##cm = metrics.confusion_matrix(y_test.argmax(axis=1), dtree_predictions.argmax(axis=1))
#cm = metrics.confusion_matrix(y_test, dtree_predictions)
#
## plotting the confusion matrix
#plt.figure(figsize = (10,7))
#sn.set(font_scale=1.4)#for label size
#sn.heatmap(cm, annot=True,annot_kws={"size": 16})# font size
#plt.savefig('dt.cm.png')
#plt.close()
#
#
#
###################################################
## dividing X, y into train and test data
#X_train, X_test, y_train, y_test = model_selection.train_test_split(X, Yidx, random_state = 0)
# 
## training a KNN classifier
#knn = KNeighborsClassifier(n_neighbors = 7).fit(X_train, y_train)
# 
## accuracy on X_test
#knn_accuracy = knn.score(X_test, y_test)
#print (knn_accuracy)
# 
## creating a confusion matrix
#knn_predictions = knn.predict(X_test) 
#cm = metrics.confusion_matrix(y_test, knn_predictions)
#
## plotting the confusion matrix
#plt.figure(figsize = (10,7))
#sn.set(font_scale=1.4)#for label size
#sn.heatmap(cm, annot=True,annot_kws={"size": 16})# font size
#plt.savefig('knn.cm.png')
#plt.close()
#
#
#
####################################################
## dividing X, y into train and test data
#X_train, X_test, y_train, y_test = model_selection.train_test_split(X, Yidx, random_state = 0)
# 
## training a Naive Bayes classifier
#gnb = GaussianNB().fit(X_train, y_train)
#gnb_predictions = gnb.predict(X_test)
# 
## accuracy on X_test
#nb_accuracy = gnb.score(X_test, y_test)
#print (nb_accuracy)
# 
## creating a confusion matrix
#cm = metrics.confusion_matrix(y_test, gnb_predictions)
#
## plotting the confusion matrix
#plt.figure(figsize = (10,7))
#sn.set(font_scale=1.4)#for label size
#sn.heatmap(cm, annot=True,annot_kws={"size": 16})# font size
#plt.savefig('nb.cm.png')
#plt.close()
#
#
#
#
###################################################
## dividing X, y into train and test data
#X_train, X_test, y_train, y_test = model_selection.train_test_split(X, Yidx, random_state = 0)
#
#sklearn_lda = LinearDiscriminantAnalysis(n_components=2)
#lda = sklearn_lda.fit(X_train, y_train)
#X_lda = lda.transform(X_train) #using the model to project X 
#
#Z = lda.transform(X_test) #using the model to project Z
#lda_predictions = lda.predict(X_test) #gives you the predicted label for each sample
#z_prob = lda.predict_proba(X_test) #the probability of each sample to belong to each class
#
## model accuracy for X_test  
#lda_accuracy = lda.score(X_test, y_test)
#print (lda_accuracy)
# 
## creating a confusion matrix
#cm = metrics.confusion_matrix(y_test, lda_predictions)
#
## plotting the confusion matrix
#plt.figure(figsize = (10,7))
#sn.set(font_scale=1.4)#for label size
#sn.heatmap(cm, annot=True,annot_kws={"size": 16})# font size
#plt.savefig('lda.cm.png')
#plt.close()
#
##plot_scikit_lda(X_lda, title='Default LDA via scikit-learn')
#
#
####################################################
## dividing X, y into train and test data
#
## training a linear SVM classifier
#svm_model_linear = svm.SVC(kernel = 'linear', C = 1).fit(X_train, y_train)
#svm_predictions = svm_model_linear.predict(X_test)
# 
## model accuracy for X_test  
#svm_accuracy = svm_model_linear.score(X_test, y_test)
#print (svm_accuracy)
# 
## creating a confusion matrix
#cm = metrics.confusion_matrix(y_test, svm_predictions)
#
## plotting the confusion matrix
#plt.figure(figsize = (10,7))
#sn.set(font_scale=1.4)#for label size
#sn.heatmap(cm, annot=True,annot_kws={"size": 16})# font size
#plt.savefig('svm.cm.png')
#plt.close()
#


################################################
## prepare configuration for cross validation test harness
#seed = 7
## prepare models
#models = []
#models.append(('LR', LogisticRegression()))
#models.append(('LDA', LinearDiscriminantAnalysis(n_components=2)))
#models.append(('KNN', KNeighborsClassifier()))
#models.append(('CART', DecisionTreeClassifier()))
#models.append(('NB', GaussianNB()))
#models.append(('SVM', svm.SVC(kernel = 'linear', C = 1)))
## evaluate each model in turn
#results = []
#names = []
#scoring = 'accuracy'
#for name, model in models:
#  kfold = model_selection.KFold(n_splits=10, random_state=seed)
#  cv_results = model_selection.cross_val_score(model, X, Yidx, cv=kfold, scoring=scoring)
#  results.append(cv_results)
#  names.append(name)
#  msg = "%s: %f (%f)" % (name, cv_results.mean(), cv_results.std())
#  print(msg)
## boxplot algorithm comparison
#fig = plt.figure()
#fig.suptitle('Algorithm Comparison')
#ax = fig.add_subplot(111)
#plt.boxplot(results)
#ax.set_xticklabels(names)
#plt.show()
#plt.savefig('cv.boxplot.png')
#plt.close()
#
#
#

























################################################
# we create an instance of SVM and fit out data. We do not scale our
# data since we want to plot the support vectors


cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, verbose=10)
cached_pipe = Pipeline([
    ('reduce_dim', PCA(iterated_power=7)),
    ('classify',  svm.SVC(kernel='linear'))],
    memory=memory)

N_FEATURES_OPTIONS = [10, 100, 300, 1000, 2000]
C_OPTIONS = [0.01, 0.1, 1, 10, 100]
gamma_range = np.logspace(-9, 3, 13)
CLASSIFIERS = [
          svm.SVC(kernel='linear'), 
          #svm.LinearSVC(),
          #svm.SVC(kernel='rbf', gamma=0.7),
          #svm.SVC(kernel='poly', degree=3)      
          ]
param_grid = [ 
    {   
        'reduce_dim__n_components': N_FEATURES_OPTIONS,
        'classify__C': C_OPTIONS
        #'classify__gamma': gamma_range
    }  
]
reducer_labels = ['PCA', 'NMF']
classifier_labels = ['SVClinear']
#classifier_labels = ['SVClinear', 'SVCrbf', 'SVCpoly']
#classifier_labels = ['SVClinear', 'LinearSVC', 'SVCrbf', 'SVCpoly']

grid = model_selection.GridSearchCV(cached_pipe, cv=5, n_jobs=1, param_grid=param_grid)
grid.fit(X_train, y_train)
joblib.dump(grid, 'grid.set3.pkl')


grid = joblib.load('grid.set3.pkl')
y_predictions = grid.predict(X_test)

report = metrics.classification_report( y_test, y_predictions )
print(report)

print("Best parameters set found on development set:")
print()
print(grid.best_params_)

print("With a Best Score of:")
print()
print(grid.best_score_)
results = grid.cv_results_

# Delete the temporary cache before exiting
rmtree(cachedir)



mean_scores = np.array(grid.cv_results_['mean_test_score'])
print (mean_scores.shape)
# scores are in the order of param_grid iteration, which is alphabetical
mean_scores = mean_scores.reshape(-1, 2, len(N_FEATURES_OPTIONS) )
print (mean_scores.shape)
# select score for best C
mean_scores = mean_scores.max(axis=0)
bar_offsets = (np.arange(len(N_FEATURES_OPTIONS)) *
               (len(reducer_labels) + 1) + .5)

plt.figure()
COLORS = 'bgrcmyk'
for i, (label, reducer_scores) in enumerate(zip(reducer_labels, mean_scores)):
    plt.bar(bar_offsets + i, reducer_scores, label=label, color=COLORS[i])

plt.title("Comparing feature reduction techniques")
plt.xlabel('Reduced number of features')
plt.xticks(bar_offsets + len(reducer_labels) / 2, N_FEATURES_OPTIONS)
plt.ylabel('tissue classification accuracy')
plt.ylim((0, 1))
plt.legend(loc='upper left')

plt.show()
plt.savefig('grid.reducer.png')
plt.close()




mean_scores = np.array(grid.cv_results_['mean_test_score'])
# scores are in the order of param_grid iteration, which is alphabetical
print (mean_scores.shape)
mean_scores = mean_scores.reshape(len(CLASSIFIERS), len(C_OPTIONS), -1 )
print (mean_scores.shape)
# select score for best C
mean_scores = mean_scores.max(axis=2)
bar_offsets = (np.arange(len(CLASSIFIERS)) *
               (len(C_OPTIONS) + 1) + .5)


plt.figure()
COLORS = 'bgrcmyk'
for i, (label, classifier_scores) in enumerate(zip(classifier_labels, mean_scores)):
    plt.bar(bar_offsets + i, classifier_scores, label=label, color=COLORS[i])

plt.title("Comparing classifiers")
plt.xlabel('SVM kernels')
plt.xticks(bar_offsets + len(C_OPTIONS) / 2, C_OPTIONS) 
plt.ylabel('tissue classification accuracy')
plt.ylim((0, 1))
plt.legend(loc='upper left')

plt.show()
plt.savefig('grid.classifier.png')
plt.close()



mean_scores = np.array(grid.cv_results_['mean_test_score'])
# scores are in the order of param_grid iteration, which is alphabetical
print (mean_scores.shape)
mean_scores = mean_scores.reshape(len(C_OPTIONS), len(N_FEATURES_OPTIONS) )
#mean_scores = mean_scores.reshape(len(C_OPTIONS), -1, len(N_FEATURES_OPTIONS) )
print (mean_scores.shape)
# select score for best C
#mean_scores = mean_scores.max(axis=0)
bar_offsets = (np.arange(len(C_OPTIONS)) *
               (len(N_FEATURES_OPTIONS) + 1) + .5)
#bar_offsets = (np.arange(len(gamma_range)) *
#               (len(N_FEATURES_OPTIONS) + 1) + .5)
print (bar_offsets)
print (mean_scores.shape)
print (mean_scores)


plt.figure()
COLORS = 'bgrcmyk'
for i, (label, classifier_scores) in enumerate(zip(C_OPTIONS, mean_scores)):
    plt.bar(bar_offsets + i, classifier_scores, label=label, color=COLORS[i])
#for i, (label, classifier_scores) in enumerate(zip(gamma_range, mean_scores)):
#    plt.bar(bar_offsets + i, classifier_scores, label=label, color=COLORS[i])

plt.title("Comparing parameters")
plt.xlabel('N Features')
plt.xticks(bar_offsets + len(N_FEATURES_OPTIONS) / 2, N_FEATURES_OPTIONS) 
plt.ylabel('test accuracy')
plt.ylim((0, 1))
plt.legend(loc='upper left')

plt.show()
plt.savefig('grid.parameters.linear.png')
plt.close()


###########################################################



#############################################################

def f_importances(coef, names, cl):
    imp = coef
    imp,names = zip(*sorted(zip(imp,names)))
    plt.figure()
    plt.barh(range(20), imp[:10]+imp[-10:], align='center')
    plt.yticks(range(20), names[:10]+names[-10:])
    plt.show()  
    plt.savefig('important_features.'+cl+'.pdf')
    plt.close()
    return names


features_names = Xlabels
clf = svm.LinearSVC(C=0.01)
clf.fit(X_train, y_train)
y_predictions = clf.predict(X_test)
report = metrics.classification_report( y_test, y_predictions )
print(report)

for cl in sorted(Ylabels.keys()):
  cl_names = f_importances(clf.coef_[cl,], features_names, Ylabels[cl])
  print(Ylabels[cl])
  print(cl_names[:50])
  print(cl_names[-50:])
  print()



#print (grid)
results = grid.cv_results_
#print (results)

scaler = grid.best_estimator_.named_steps['reduce_dim']
classifier = grid.best_estimator_.named_steps['classify']
print (classifier)

importances = classifier.feature_importances_
features = grid.named_steps['feat']
print (X.columns[features.transform(np.arange(len(X.columns)))])





##################################################################
C = 1.0  # SVM regularization parameter
models = (svm.SVC(kernel='linear', C=C),
        svm.LinearSVC(C=C),
        svm.SVC(kernel='rbf', gamma=0.7, C=C),
        svm.SVC(kernel='poly', degree=3, C=C))
models = (clf.fit(X_train, y_train) for clf in models)

# title for the plots
titles = ('SVC with linear kernel',
        'LinearSVC (linear kernel)',
        'SVC with RBF kernel',
        'SVC with polynomial (degree 3) kernel')

# Set-up 2x2 grid for plotting.
fig, sub = plt.subplots(2, 2)
plt.subplots_adjust(wspace=0.4, hspace=0.4)

X0, X1 = X_train[:, 0], X_train[:, 1]
xx, yy = make_meshgrid(X0, X1)

for clf, title, ax in zip(models, titles, sub.flatten()):
  plot_contours(ax, clf, xx, yy,
                cmap=plt.cm.coolwarm, alpha=0.8)
  ax.scatter(X0, X1, c=y, cmap=plt.cm.coolwarm, s=20, edgecolors='k')
  ax.set_xlim(xx.min(), xx.max())
  ax.set_ylim(yy.min(), yy.max())
  ax.set_xlabel('Sepal length')
  ax.set_ylabel('Sepal width')
  ax.set_xticks(())
  ax.set_yticks(())
  ax.set_title(title)

plt.show()
plt.savefig('svm.kernels.png')
plt.close()
