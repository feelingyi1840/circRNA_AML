###### Evaluation for different Panels -- Train Test ######

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.metrics import roc_curve, auc, precision_recall_curve, f1_score, average_precision_score,RocCurveDisplay, PrecisionRecallDisplay
import matplotlib.pyplot as plt
from imblearn.over_sampling import SMOTE
from sklearn.neural_network import MLPClassifier  

def gaussian(dist, sigma = 10.0):
    """ Input a distance and return it`s weight"""
    weight = np.exp(-dist**2/(2*sigma**2))
    return weight

import numpy as np
file_paths = {
    'Panel-TopN.csv': 'Panel-TopN',
    'Panel-Pathway_51.csv': 'Panel-Pathway',
    'Panel-3DG.csv': 'Panel-3DG-Cluster',
    'Panel-Radius5-test.csv': 'Panel-3DG-Radius5',
    'Panel-Radius1-test.csv': 'Panel-3DG-Radius1'
}

models = {
    'Logistic Regression': LogisticRegression(C=7),
    'Random Forest': RandomForestClassifier(n_estimators=1200, max_depth=None, max_features= None, min_samples_leaf = 1, min_samples_split = 2),
    'Support Vector Machine': SVC(C=2, kernel='rbf', probability=True),
    'K-Nearest Neighbors': KNeighborsClassifier(weights=gaussian),
    'Gradient Boosting': GradientBoostingClassifier(),
    'XGBoost': XGBClassifier(use_label_encoder=False, eval_metric='logloss'),
    'MLP': MLPClassifier(hidden_layer_sizes=(30, 20), max_iter=1000, random_state=42),
}

for model_name, model in models.items():
    y_true = {}
    y_scores = {}
    for file_path, file in file_paths.items():
        data = pd.read_csv(file_path)
        X = data.iloc[:, 1:]  # Feature
        y = data.iloc[:, 0]   # Label

        X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=94, test_size=26, random_state=42, stratify=y)
        smote = SMOTE(random_state=42)
        X_train, y_train = smote.fit_resample(X_train, y_train)

        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)[:,1]

        y_true[file] = y_test
        y_scores[file] = y_proba

    plt.figure(1, figsize=(8, 6))
    for label in y_true.keys():
        fpr, tpr, _ = roc_curve(y_true[label], y_scores[label])
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr,label=f'{label} AUC: {roc_auc:.3f}')

    plt.plot([0, 1], [0, 1], linestyle='--', color='gray')

    plt.title(f"{model_name}"+' ROC AUC Curve',fontsize=18)
    plt.xlabel('False Positive Rate',fontsize=16)
    plt.ylabel('True Positive Rate',fontsize=16)
    plt.legend(loc='lower right',fontsize=14)
    plt.grid()
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])

    plt.savefig(f'{model_name}_ROC.png', dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(2, figsize=(8, 6))
    for label in y_true.keys():
        precision, recall, _ = precision_recall_curve(y_true[label], y_scores[label])
        avg_precision = average_precision_score(y_true[label], y_scores[label])
        plt.plot(recall, precision, label=f'{label} AUPR: {avg_precision:.3f}')

    plt.title(f"{model_name}"+' AUPR Curve',fontsize=18)
    plt.xlabel('Recall',fontsize=16)
    plt.ylabel('Precision',fontsize=16)
    plt.legend(loc='lower left', fontsize=14)
    plt.grid()

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])

    plt.savefig(f'{model_name}_AUPR.png', dpi=300, bbox_inches='tight')
    plt.close()


### More Matrics ###
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.metrics import roc_curve, auc, precision_recall_curve, f1_score, average_precision_score,RocCurveDisplay, PrecisionRecallDisplay
import matplotlib.pyplot as plt
from imblearn.over_sampling import SMOTE
from sklearn.neural_network import MLPClassifier  


from sklearn.metrics import (roc_auc_score, average_precision_score, accuracy_score, 
                             precision_score, recall_score, f1_score, cohen_kappa_score,
                             roc_curve, precision_recall_curve)

file_paths = {
    'Panel-TopN.csv': 'Panel-TopN',
    'Panel-Pathway_51.csv': 'Panel-Pathway',
    'Panel-3DG.csv': 'Panel-3DG-Cluster',
    'Panel-Radius5-test.csv': 'Panel-3DG-Radius5',
    'Panel-Radius1-test.csv': 'Panel-3DG-Radius1'
}

models = {
    'Logistic Regression': LogisticRegression(C=7),
    'Random Forest': RandomForestClassifier(n_estimators=1200, max_depth=None, max_features= None, min_samples_leaf = 1, min_samples_split = 2),
    'Support Vector Machine': SVC(C=2, kernel='rbf', probability=True),
    'K-Nearest Neighbors': KNeighborsClassifier(weights=gaussian),
    'Gradient Boosting': GradientBoostingClassifier(),
    'XGBoost': XGBClassifier(eval_metric='logloss'),
    'MLP': MLPClassifier(hidden_layer_sizes=(30, 20), max_iter=1000, random_state=42),
}

for file_path, file in file_paths.items():
    print(f'Processing: {file}')
    data = pd.read_csv(file_path)
    X = data.iloc[:, 1:]  # Feature
    y = data.iloc[:, 0]   # Label

    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=94, test_size=26, random_state=42, stratify=y)
    smote = SMOTE(random_state=42)
    X_train, y_train = smote.fit_resample(X_train, y_train)

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    
    eva = None
    for model_name, model in models.items():
        model.fit(X_train, y_train)
        yhat_proba = model.predict_proba(X_test)[:, 1]
        yhat_class = model.predict(X_test)
        
        model_eva = {}
        model_eva['Model'] = model_name
        model_eva['Accuracy'] = accuracy_score(y_test, yhat_class)
        model_eva['Precision'] = precision_score(y_test, yhat_class, zero_division=0)
        model_eva['Recall'] = recall_score(y_test, yhat_class, zero_division=0)
        model_eva['F1_score'] = f1_score(y_test, yhat_class, zero_division=0)
        model_eva['AUC'] = roc_auc_score(y_test, yhat_proba)
        model_eva['AUPR'] = average_precision_score(y_test, yhat_proba)
        model_eva['Kappa_score'] = cohen_kappa_score(y_test, yhat_class)

        model_eva_df = pd.DataFrame([model_eva])
        if eva is None:
            eva = model_eva_df
        else:
            eva = pd.concat([eva,model_eva_df],ignore_index=True)
    print(eva)
    eva.to_csv(f'{file}_model_evaluation_results.csv')