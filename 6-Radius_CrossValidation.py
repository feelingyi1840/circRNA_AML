###### Radius positioning of circRNA in AML prediction -- Cross Validation ######

import pandas as pd
import random
import math
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import numpy as np
import seaborn as sns

from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import roc_curve, auc, precision_recall_curve, roc_auc_score, average_precision_score
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier
from sklearn.ensemble import GradientBoostingClassifier, AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import lightgbm as lgb
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore')


# Randomly select n circRNA from each bin
def select(data, bin_start, bin_end, n):
    bin_data = data[(data['distance'] > bin_start) & (data['distance'] <= bin_end)]
    if len(bin_data) > 0:
        sample_size = min(n, len(bin_data))
        selected_rows = bin_data.sample(sample_size)
        return selected_rows
    else:
        return pd.DataFrame()


# Perform final clustering using the optimal k values, and calculate silhouette scores
def calculate_silhouette_scores(selected_rows,n):
    select_data = selected_rows[['X', 'Y', 'Z']]
    best_k = 0
    best_score = -1
    max_k = min(n, 10)
    if max_k > 1:  # At least two clusters are needed
        for k in range(2, max_k + 1):
            try:
                kmeans = KMeans(n_clusters=k, random_state=42)
                labels = kmeans.fit_predict(select_data)
                score = silhouette_score(select_data, labels)

                if score > best_score:
                    best_k = k
                    best_score = score

            except Exception as e:
                print(f"k = {k} encountered an error: {e}")

    return best_score,best_k


def evaluate_model(models,model, X, y, cv):
    roc_aucs = []
    pr_aucs = []
    fprs, tprs, precisions, recalls = [], [], [], []
    for train_ix, test_ix in cv.split(X, y):
        X_train, X_test = X.iloc[train_ix], X.iloc[test_ix]
        y_train, y_test = y.iloc[train_ix], y.iloc[test_ix]
        fitted_model = models[model].fit(X_train, y_train)
        yhat = fitted_model.predict_proba(X_test)[:, 1]
        fpr, tpr, _ = roc_curve(y_test, yhat)
        precision, recall, _ = precision_recall_curve(y_test, yhat)
        fprs.append(fpr)
        tprs.append(tpr)
        precisions.append(precision)
        recalls.append(recall)
        roc_auc = roc_auc_score(y_test, yhat)
        pr_auc = average_precision_score(y_test, yhat)
        roc_aucs.append(roc_auc)
        pr_aucs.append(pr_auc)

    return np.mean(roc_aucs), np.mean(pr_aucs)


def model_auc(data):
    X = data.iloc[:, 1:]
    y = data.iloc[:, 0]

    models = {
        "Logistic Regression": LogisticRegression(max_iter=1000),
        "SVM": SVC(probability=True),
        "Random Forest": RandomForestClassifier(),
        "K-Nearest Neighbors": KNeighborsClassifier(),
        "XGBoost": XGBClassifier(),
        "Gradient Boosting Classifier": GradientBoostingClassifier(),
        "Decision Tree Classifier": DecisionTreeClassifier(),
        "Naive Bayes": GaussianNB(),
        "AdaBoost": AdaBoostClassifier(),
        "Multi-layer Perceptron": MLPClassifier(max_iter=1000)
    }

    cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=100, random_state=42)

    results = {}
    for name in models:
        roc_auc, pr_auc = evaluate_model(models,name, X, y, cv)
        tag1 = name + " ROC-AUC"
        results[tag1] = [roc_auc]
        tag2 = name + " PR-AUC"
        results[tag2] = [pr_auc]

    return results


# Load circ_3D_co.csv
circ_3d_data = pd.read_csv('all_bin_co_with_distance.csv')
circ_3d = circ_3d_data[['circRNA_ID', 'X', 'Y', 'Z', 'distance']]

# Load bin information
bins = pd.DataFrame({
    'Bin Start': [0.000000, 1.185943, 1.581257, 1.976571, 2.371886, 2.767200, 3.162514, 3.557828, 3.953143, 4.348457, 4.743771, 5.139085, 5.534400, 5.929714, 6.325028, 6.720343, 7.115657, 7.510971],
    'Bin End': [1.185943, 1.581257, 1.976571, 2.371886, 2.767200, 3.162514, 3.557828, 3.953143, 4.348457, 4.743771, 5.139085, 5.534400, 5.929714, 6.325028, 6.720343, 7.115657, 7.510971, 7.906285]
})

file_path = 'Panel-All_cv.csv'
data = pd.read_csv(file_path)


for n in [3,5,10, 20, 50]: #n = 3、5、10、20、50
    for bin_index, bin_range in bins.iterrows():
        output_file0 = f"result/bin_{bin_index}_n=_{n}_results_cv_"

        for i in range(101):
            selected_rows = select(circ_3d, bin_range['Bin Start'], bin_range['Bin End'], n)
            print(bin_range['Bin Start'], bin_range['Bin End'])
            if selected_rows.empty:
                continue

            score, _ = calculate_silhouette_scores(selected_rows, n)

            geometry = {
                'num': [i],
                'silhouette_score': [score]
            }
            geometry = pd.DataFrame(geometry)

            ids_to_select = selected_rows["circRNA_ID"].tolist()
            columns_to_keep = data.columns.isin(ids_to_select)
            B_filtered = data.loc[:, columns_to_keep]
            B_filtered.insert(0, 'cancer_State', data['cancer_State'])
            model_result = model_auc(B_filtered)
            output_auc = pd.DataFrame(model_result)

            results = pd.concat([geometry, output_auc], axis=1)
            if i == 0:
                output = results
            elif i % 10 == 0:
                output_file = output_file0 + str(i) + ".csv"
                output.to_csv(output_file, index=False, encoding='utf-8-sig')
                output = results
            else:
                output = pd.concat([output, results], ignore_index=True)