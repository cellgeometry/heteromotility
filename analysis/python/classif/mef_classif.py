'''
Classify WT vs. MycRas MEF based on Heteromotility features utilizing
round-robin training
'''
import sys, csv
import numpy as np
import pandas as pd
import itertools

from sklearn import preprocessing
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectFromModel, SelectPercentile, f_classif, chi2, VarianceThreshold
from sklearn.decomposition import PCA
from sklearn.linear_model import Lasso
from sklearn import cross_validation
from sklearn.model_selection import GridSearchCV

from multiprocessing import Pool
from functools import partial

data_path = '../../data/mef/mycwt_classification_plate_df.csv'
# Load data
df = pd.read_csv(data_path)
y = df['class'] # load labels into `y`; 1 = mycras, 2 = wt
e = df['plate'] # load exp numbers into a vector `e`
myc_e = e[e <= 2]
wt_e = e[e > 2]
X_raw = df.iloc[:,:-2] # load remaining features into `X`
X = preprocessing.scale(X_raw) # center and scale data, mean=0, unit var


estimator_svm = ('svm', SVC())
estimator_rf = ('rf100', RandomForestClassifier(n_estimators=100))
selector_none = ('None', SelectPercentile(f_classif, percentile=100))
selector_pca = ('PCA', PCA(copy=True, n_components=None))
selector_anova = ('ANOVA', SelectPercentile(f_classif, percentile=60))

clf_svm_highdim = Pipeline([selector_anova, estimator_svm])
clf_svm_lowdim = Pipeline([selector_pca, selector_anova, estimator_svm])
clf_rf_anova = Pipeline([selector_anova, estimator_rf])
clf_rf_none = Pipeline([selector_none, estimator_rf])
clf_rf_pca = Pipeline([selector_pca, estimator_rf])


def class_balance(X, y):
    '''
    Class balance a data set by random selection from majority class.

    Parameters
    ----------
    X : ndarray. feature data.
    y : ndarray. labels.

    Returns
    -------
    Xb : ndarray. class balanced feature data.
    y : ndarray. class balanced labels.
    '''
    # type case in case pd.DataFrames are provided
    X = np.array(X)
    y = np.array(y)

    classes = np.unique(y)
    nb_class = np.zeros(len(classes))
    for i in range(len(classes)):
        nb_class[i] = (y==classes[i]).sum()
    nb_min = np.min(nb_class).astype('int32')

    # create balanced data set and labels
    Xb = np.zeros([nb_min*len(classes), X.shape[1]])
    yb = np.zeros(nb_min*len(classes))

    for i in range(len(classes)):
        c = np.where(y==classes[i])[0]
        ridx = np.random.choice(c, size=nb_min, replace=False)
        Xb[i*nb_min:(i+1)*nb_min, :] = X[ridx,:]
        yb[i*nb_min:(i+1)*nb_min] = y[ridx]

    return Xb, yb

def train_cv(clf, X, y, n_iter=100):
    '''
    Train a classifier object on data `X` labels `y` with 5-fold CV

    Parameters
    ----------
    clf : sklearn classifier object with .fit() method
    X : ndarray. N x M feature data.
    y : ndarray. N x 1 labels.
    n_iter : integer. Number of cross-val splits.

    Returns
    -------
    score : mean of cross val scores.
    '''
    # generate cross-val splits
    cv = cross_validation.StratifiedShuffleSplit(y, test_size=0.2, n_iter=n_iter)
    scores = cross_validation.cross_val_score(clf, X, y, cv=cv)
    return np.mean(scores)

def get_test_set(train, e):
    '''
    Get the remaining experiments in `e` not in `train`

    Parameters
    ----------
    train : array-like listing experiment labels for training.
    e : array-like listing all experiment labels.
    '''

    in_train = np.in1d(e, train)
    test = e[~in_train]
    return test

def round_robin(clf, X, y, e0, e1, balance=True, n_iter=10, all_scores=False):
    '''
    Perform round-robin training, leaving one experiment out from each category
    for classification

    Parameters
    ----------
    clf : scikit-learn classifier object.
    X : ndarray. feature data in N x M format.
    y : ndarray. N x 1 labels.
    e0 : ndarray. experiment labels for category `0`
    e1 : ndarray. experiment labels for category `1`
    balance : boolean. class balance data before training.
        Set to false to use a class weighted classifier.
    n_iter : iterations for CV scoring.
    all_scores : boolean. return score arrays rather than means.
    '''

    e0_exps = np.unique(e0)
    e1_exps = np.unique(e1)

    e0_combs = list(itertools.combinations(e0_exps, len(e0_exps)-1))
    e1_combs = list(itertools.combinations(e1_exps, len(e1_exps)-1))

    train_sets = list(itertools.product(e0_combs, e1_combs))

    e = np.concatenate([e0, e1])

    cv_scores = []
    pred_scores = []
    # for all models export
    models = []
    y_tests = []
    y_preds = []

    for s in range(len(train_sets)):
        train0_exp, train1_exp = train_sets[s]
        test0_exp = get_test_set(train0_exp, e0_exps)
        test1_exp = get_test_set(train1_exp, e1_exps)

        X_train, y_train = [],[]
        X_test, y_test = [], []

        # Build up training and test arrays in order of experiment class
        for i in np.unique(e):

            if i in train0_exp:
                X_train.append(X[e==i,:])
                y_train.append(y[e==i])
            elif i in test0_exp:
                X_test.append(X[e==i,:])
                y_test.append(y[e==i])
            else:
                pass

            if i in train1_exp:
                X_train.append(X[e==i,:])
                y_train.append(y[e==i])
            elif i in test1_exp:
                X_test.append(X[e==i,:])
                y_test.append(y[e==i])
            else:
                pass

        # concat arrays
        X_train = np.concatenate(X_train)
        y_train = np.concatenate(y_train)
        X_test = np.concatenate(X_test)
        y_test = np.concatenate(y_test)

        # Class balance arrays
        if balance:
            X_train, y_train = class_balance(X_train, y_train)
            X_test, y_test = class_balance(X_test, y_test)
        else:
            clf.class_weight = 'balanced'

        # Get cross val score
        score = train_cv(clf, X_train, y_train, n_iter=n_iter)
        cv_scores.append(score)
        # Fit clf
        clf.fit(X_train, y_train)
        # evaluate clf
        y_pred = clf.predict(X_test)
        pred_scores.append( (y_test == y_pred).sum() / y_test.shape[0] )

        if all_scores:
            models.append(clf)
            y_preds.append(y_pred)
            y_tests.append(y_test)

    if all_scores:
        return cv_scores, pred_scores, models, y_preds, y_tests
    else:
        return np.mean(cv_scores), np.mean(pred_scores)

def rr_predict(clf, X, y, e0, e1, balance=True, all_scores=False):
    e0_exps = np.unique(e0)
    e1_exps = np.unique(e1)

    e0_combs = list(itertools.combinations(e0_exps, len(e0_exps)-1))
    e1_combs = list(itertools.combinations(e1_exps, len(e1_exps)-1))

    train_sets = list(itertools.product(e0_combs, e1_combs))

    e = np.concatenate([e0, e1])

    cv_scores = []
    pred_scores = []
    # for all models export
    models = []
    y_tests = []
    y_preds = []

    for s in range(len(train_sets)):
        train0_exp, train1_exp = train_sets[s]
        test0_exp = get_test_set(train0_exp, e0_exps)
        test1_exp = get_test_set(train1_exp, e1_exps)

        X_train, y_train = [],[]
        X_test, y_test = [], []


        # Build up training and test arrays in order of experiment class
        for i in np.unique(e):

            if i in train0_exp:
                X_train.append(X[e==i,:])
                y_train.append(y[e==i])
            elif i in test0_exp:
                X_test.append(X[e==i,:])
                y_test.append(y[e==i])
            else:
                pass

            if i in train1_exp:
                X_train.append(X[e==i,:])
                y_train.append(y[e==i])
            elif i in test1_exp:
                X_test.append(X[e==i,:])
                y_test.append(y[e==i])
            else:
                pass

        # concat arrays
        X_train = np.concatenate(X_train)
        y_train = np.concatenate(y_train)
        X_test = np.concatenate(X_test)
        y_test = np.concatenate(y_test)

        # Class balance arrays
        if balance:
            X_train, y_train = class_balance(X_train, y_train)
            X_test, y_test = class_balance(X_test, y_test)
        else:
            clf.class_weight = 'balanced'

        # evaluate clf
        y_pred = clf.predict(X_test)
        pred_scores.append( (y_test == y_pred).sum() / y_test.shape[0] )

        if all_scores:
            models.append(clf)
            y_preds.append(y_pred)
            y_tests.append(y_test)

    if all_scores:
        return pred_scores, models, y_preds, y_tests
    else:
        return np.mean(pred_scores)

# c = [clf_svm_lowdim, clf_svm_highdim, clf_rf_none, clf_rf_anova, clf_rf_pca]
# cv_scores = []
# pred_scores = []
# for clf in c:
#     a, b = round_robin(clf, X, y, myc_e, wt_e)
#     cv_scores.append(a)
#     pred_scores.append(b)
#


_prc = np.arange(1,10) * 10
_C = np.logspace(-4, 4, 20)
_kernel = ('linear', 'rbf', 'poly', 'sigmoid')

comb = list(itertools.product(_prc, _C, _kernel))

# for prc, C, kernel in comb:
#     s0 = selector_pca
#     s1 = ('ANOVA', SelectPercentile(f_classif, percentile=prc))
#     e = ('svm', SVC(C=C, kernel=kernel))
#     clf = Pipeline([s0, s1, e])
#     a, b = round_robin(clf, X, y, myc_e, wt_e)
#     cv_scores.append(a)
#     pred_scores.append(b)

def run_round_robin_params(p, X, y, e0, e1, verbose=True):
    '''
    Run round robin on a param set `p`
    `p` : tuple of percentile for ANOVA, C for SVM, SVM kernel
    '''
    prc, C, kernel = p
    s0 = selector_pca
    s1 = ('ANOVA', SelectPercentile(f_classif, percentile=prc))
    e = ('svm', SVC(C=C, kernel=kernel))
    #clf = Pipeline([s0, s1, e])
    clf = Pipeline([s1, e])
    a, b = round_robin(clf, X, y, e0, e1)
    if verbose:
        print('-----')
        print('Tested   : ', p)
        print('Accuracy : ', b)
        print('-----')
    return (a, b)

part_rrr = partial(run_round_robin_params, X=X, y=y, e0=myc_e, e1=wt_e)

mpool = Pool()
# scores = list(mpool.map(part_rrr, comb))
#
# scores = np.array(scores)
# best_fit = np.argmax(scores[:,1])
# best_params = comb[best_fit]

# Redo search with a narrow param set to try and further increase accuracy
_prc = np.linspace(50,70,25)
_C = np.linspace(0.15,0.9,50)
_kernel = (['rbf'])
comb_narrow = list(itertools.product(_prc, _C, _kernel))

scores_agg = []
for i in range(3):

    scores_narrow = list(mpool.map(part_rrr, comb_narrow))
    scores_agg.append(scores_narrow)

scores_narrow=pd.DataFrame(np.array(scores_narrow))

df = pd.DataFrame(comb_narrow)
df = pd.concat([df, scores_narrow], axis=1)
df.columns = ['ANOVA_prc', 'SVM_C', 'SVM_kernel', 'cv_score', 'round_robin']
df.to_csv('round_robin_scores.csv', index=False)

# Get round robin split accuracies
best_pidx = np.argmax(df.iloc[:,-1])
best_params = tuple(df.iloc[best_pidx,:-2])
prc, C, kernel = best_params

s1 = ('ANOVA', SelectPercentile(f_classif, percentile=prc))
e = ('svm', SVC(C=C, kernel=kernel))
clf = Pipeline([s1, e])

cv_scores, pred_scores, models, y_tests, y_preds = round_robin(clf, X, y, myc_e, wt_e, all_scores=True)

# save outputs
pickle.dump(y_tests, open('mycwt_classif_svm_round_robin_tests.pickle', 'wb'))
pickle.dump(y_tests, open('mycwt_classif_svm_round_robin_preds.pickle', 'wb'))
pickle.dump(models, open('mycwt_classif_svm_round_robin_fitted.pickle', 'wb'))

np.savetxt('mycwt_round_robin_scores.csv', pred_scores, delimiter=',')

# plot confusion matrix
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
cm = confusion_matrix(np.concatenate(y_tests), np.concatenate(y_preds))

cm = pd.DataFrame(cm)
cm.columns = ['MycRas', 'WT']
cm.index = ['MycRas', 'WT']

cmap = sns.cubehelix_palette(32, reverse=False, start=2.8, rot=.1)

g = sns.heatmap(cm, cbar_kws={'label':'Count'}, annot=True, fmt='d', cmap=cmap)
plt.xlabel('Predicted Class')
plt.ylabel('True Class')
plt.title('MEF Classification Confusion Matrix')

# violin plot of predicted scores

g = sns.violinplot(data=pred_scores)

def optimize_clf(clf, X, y):
    '''Optimize clf round-robin score with GridSearch'''

    parameters = {svm__kernel:('linear', 'rbf', 'poly', 'sigmoid'),
                    svm__C:np.logspace(-4, 4, 50),
                    svm__shrinking:(True, False)}

    gso = GridSearchCV(clf, parameters)

def min_features(selector, X, y):
    selector = selector.fit(X, y)
    min_score = np.percentile(selector.scores_, 100-selector.percentile)
    features_idx = selector.scores_ > min_score
    return features_idx


for prc in [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
    selector = SelectPercentile(f_classif, prc)
    features_idx = min_features(selector, X, y)
    features = X_raw.columns[features_idx]

    np.savetxt('mycwt_min_features_bool' + str(prc).zfill(3) + '.csv', features_idx, delimiter=',')
    features.to_series().to_csv('mycwt_min_features_names' + str(prc).zfill(3) + '.csv', index=False)

# show accuracy decreases with feature elimination
_prc = np.arange(1,10) * 10
_C = np.linspace(0.4, 0.6, 5)
_kernel = ['rbf']
comb_feature_red = list(itertools.product(_prc, _C, _kernel))

scores_feat_red = []
for p in comb_feature_red:
    s = part_rrr(p)
    scores_feat_red += [s]
    print(s)

comb_df = pd.DataFrame(comb_feature_red)
scores_df = pd.DataFrame(scores_feat_red)

df = scores_df
df.columns = ['cv', 'round_robin']
df['prc'] = comb_df.iloc[:,0]

# plot
import matplotlib.pyplot as plt
import seaborn as sns
sns.lmplot(x='prc', y='round_robin', data=df)
plt.title('MEF Classification Accuracy')
plt.xlabel('Percentage of Features Used')
plt.ylabel('Round Robin Accuracy')
plt.tight_layout()

from scipy.stats import linregress
f = linregress(np.array(df['prc']), np.array(df['round_robin']))
print(f.rvalue, f.pvalue)




# find feature importances by removing each feature

def predict_all_models(m_list, X, y, e0, e1):
    '''
    m_list : list of clf model objects.

    Returns
    -------
    mean_acc : mean accuracy across models.
    '''
    acc = []
    for i in range(len(m_list)):
        acc.append( rr_predict(m_list[i], X, y, e0, e1) )
    return np.mean(acc)


def leave_features_out(clf, X, y, e0, e1):
    '''
    Find the most important features based on a decrease in accuracy when
    features are removed.
    '''

    base_acc = np.zeros(10)
    for i in range(10):
        _, base_acc[i] = round_robin(clf, X, y, e0, e1)
    mean_base_acc = np.mean(base_acc)

    feature_accs = np.zeros(X.shape[1])
    for i in range(X.shape[1]):
        X_rm = np.delete(X, i, axis=1)

        accs = np.zeros(10)
        for j in range(10):
            _, accs[j] = round_robin(clf, X_rm, y, e0, e1)
        feature_accs[i] = np.mean(accs)

    diff_accs = feature_accs - np.mean(base_acc)
    return diff_accs

feats = min_features(s1[1], X, y)
diff_accs = leave_features_out(clf, X[:,feats], y, myc_e, wt_e)
idx = np.argsort(diff_accs*-1)

top10 = diff_accs[idx[-10:]]
top10_names = X_raw.columns[feats][idx[-10:]]
top10_df = pd.concat([pd.DataFrame(top10_names), pd.DataFrame(top10*-1)], 1)

top10_df.columns = ['Feature', 'Decrease in Accuracy']
top10_df = top10_df.sort_values('Decrease in Accuracy')
top10_df.to_csv('feature_removal_top10_names_accs.csv', index = False)
