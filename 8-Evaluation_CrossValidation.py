###### Evaluation for different Panels -- Cross Validation ######

import pandas as pd
import numpy as np
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import (roc_auc_score, average_precision_score, accuracy_score, 
                             precision_score, recall_score, f1_score, cohen_kappa_score,
                             roc_curve, precision_recall_curve)
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from xgboost import XGBClassifier

models = {
    "Logistic Regression": LogisticRegression(max_iter=1000),
    "SVM": SVC(probability=True),
    "K-Nearest Neighbors": KNeighborsClassifier(),
    "Random Forest": RandomForestClassifier(),
    "XGBoost": XGBClassifier(),
    "Gradient Boosting Classifier": GradientBoostingClassifier(),
}

cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=100, random_state=42)

def evaluate_model(model, X, y):
    accuracies = []
    precisions = []
    recalls = []
    f1_scores = []
    roc_aucs = []
    pr_aucs = []
    kappa_scores = []
    
    for train_ix, test_ix in cv.split(X, y):
        X_train, X_test = X.iloc[train_ix], X.iloc[test_ix]
        y_train, y_test = y.iloc[train_ix], y.iloc[test_ix]
        
        fitted_model = model.fit(X_train, y_train)
        
        yhat_proba = fitted_model.predict_proba(X_test)[:, 1]
        yhat_class = fitted_model.predict(X_test)
        
        accuracy = accuracy_score(y_test, yhat_class)
        precision = precision_score(y_test, yhat_class, zero_division=0)
        recall = recall_score(y_test, yhat_class, zero_division=0)
        f1 = f1_score(y_test, yhat_class, zero_division=0)
        roc_auc = roc_auc_score(y_test, yhat_proba)
        pr_auc = average_precision_score(y_test, yhat_proba)
        kappa = cohen_kappa_score(y_test, yhat_class)
        
        accuracies.append(accuracy)
        precisions.append(precision)
        recalls.append(recall)
        f1_scores.append(f1)
        roc_aucs.append(roc_auc)
        pr_aucs.append(pr_auc)
        kappa_scores.append(kappa)

    return {
        # Mean
        'Accuracy_mean': np.mean(accuracies),
        'Precision_mean': np.mean(precisions),
        'Recall_mean': np.mean(recalls),
        'F1_score_mean': np.mean(f1_scores),
        'AUPR_mean': np.mean(pr_aucs),
        'AUC_mean': np.mean(roc_aucs),
        'Kappa_score_mean': np.mean(kappa_scores),
        
        # Var
        'Accuracy_var': np.var(accuracies),
        'Precision_var': np.var(precisions),
        'Recall_var': np.var(recalls),
        'F1_score_var': np.var(f1_scores),
        'AUPR_var': np.var(pr_aucs),
        'AUC_var': np.var(roc_aucs),
        'Kappa_score_var': np.var(kappa_scores),
        
        # Std
        'Accuracy_std': np.std(accuracies),
        'Precision_std': np.std(precisions),
        'Recall_std': np.std(recalls),
        'F1_score_std': np.std(f1_scores),
        'AUPR_std': np.std(pr_aucs),
        'AUC_std': np.std(roc_aucs),
        'Kappa_score_std': np.std(kappa_scores)
    }

#-------------------------------------------------------------------
file = 'Panel-TopN.csv'
data = pd.read_csv(file)
X = data.iloc[:, 1:]
X.columns = [col.replace('|', '_').replace(':', '_') for col in X.columns]
y = data.iloc[:, 0]
print(file)

results_list = []

for name, model in models.items():
    print(f"Evaluating {name}...")
    metrics = evaluate_model(model, X, y)
    metrics['Model'] = name
    results_list.append(metrics)

results_df = pd.DataFrame(results_list)

print("\nFinal Results:")
print(results_df)

results_df.to_csv('Panel-TopN_model_evaluation_results.csv')


#-------------------------------------------------------------------
file = 'Panel-Pathway_51.csv'
data = pd.read_csv(file)
X = data.iloc[:, 1:]
X.columns = [col.replace('|', '_').replace(':', '_') for col in X.columns]
y = data.iloc[:, 0]
print(file)

results_list = []

for name, model in models.items():
    print(f"Evaluating {name}...")
    metrics = evaluate_model(model, X, y)
    metrics['Model'] = name
    results_list.append(metrics)

results_df = pd.DataFrame(results_list)

print("\nFinal Results:")
print(results_df)

results_df.to_csv('Panel-Pathway_model_evaluation_results.csv')


#-------------------------------------------------------------------
file = 'Panel-3DG.csv'
data = pd.read_csv(file)
X = data.iloc[:, 1:]
X.columns = [col.replace('|', '_').replace(':', '_') for col in X.columns]
y = data.iloc[:, 0]
print(file)

results_list = []

for name, model in models.items():
    print(f"Evaluating {name}...")
    metrics = evaluate_model(model, X, y)
    metrics['Model'] = name
    results_list.append(metrics)

results_df = pd.DataFrame(results_list)

print("\nFinal Results:")
print(results_df)

results_df.to_csv('Panel-3DG-Cluster_model_evaluation_results.csv')


#-------------------------------------------------------------------
file = 'Panel-Radius5.csv'
data = pd.read_csv(file)
X = data.iloc[:, 1:]
X.columns = [col.replace('|', '_').replace(':', '_') for col in X.columns]
y = data.iloc[:, 0]
print(file)

results_list = []

for name, model in models.items():
    print(f"Evaluating {name}...")
    metrics = evaluate_model(model, X, y)
    metrics['Model'] = name
    results_list.append(metrics)

results_df = pd.DataFrame(results_list)

print("\nFinal Results:")
print(results_df)

results_df.to_csv('Panel-3DG-Radius5_model_evaluation_results.csv')


#-------------------------------------------------------------------
file = 'Panel-Radius1.csv'
data = pd.read_csv(file)
X = data.iloc[:, 1:]
X.columns = [col.replace('|', '_').replace(':', '_') for col in X.columns]
y = data.iloc[:, 0]
print(file)

results_list = []

for name, model in models.items():
    print(f"Evaluating {name}...")
    metrics = evaluate_model(model, X, y)
    metrics['Model'] = name
    results_list.append(metrics)

results_df = pd.DataFrame(results_list)

print("\nFinal Results:")
print(results_df)

results_df.to_csv('Panel-3DG-Radius1_model_evaluation_results.csv')