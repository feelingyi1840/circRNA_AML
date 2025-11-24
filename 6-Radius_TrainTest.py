###### Radius positioning of circRNA in AML prediction -- Train Test ######

import pandas as pd
import random
import math
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, precision_recall_curve, roc_auc_score, average_precision_score
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from imblearn.over_sampling import SMOTE

# Suppress all warnings
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
def calculate_silhouette_scores(selected_rows, n):
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

    return best_score, best_k

def model_auc(data):
    X = data.iloc[:, 1:]  # Feature
    y = data.iloc[:, 0]   # Label

    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=94, test_size=26, random_state=42, stratify=y)
    smote = SMOTE(random_state=42)
    X_train, y_train = smote.fit_resample(X_train, y_train)

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    models = {
        'Logistic Regression': LogisticRegression(),
        'Random Forest': RandomForestClassifier(),
        'Support Vector Machine': SVC(probability=True),
        'K-Nearest Neighbors': KNeighborsClassifier(),
        'Gradient Boosting': GradientBoostingClassifier(),
        'XGBoost': XGBClassifier(use_label_encoder=False, eval_metric='logloss')
    }

    results = {}

    for model_name, model in models.items():
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        auc = roc_auc_score(y_test, y_pred)

        results[model_name] = [auc]

    return results


# Load circ_3D_co.csv
circ_3d_data = pd.read_csv('all_bin_co_with_distance.csv')
circ_3d = circ_3d_data[['circRNA_ID', 'X', 'Y', 'Z', 'distance']]

# Load bin information
bins = pd.DataFrame({
    'Bin Start': [0.000000, 1.185943, 1.581257, 1.976571, 2.371886, 2.767200, 3.162514, 3.557828, 3.953143, 4.348457, 4.743771, 5.139085, 5.534400, 5.929714, 6.325028, 6.720343, 7.115657, 7.510971],
    'Bin End': [1.185943, 1.581257, 1.976571, 2.371886, 2.767200, 3.162514, 3.557828, 3.953143, 4.348457, 4.743771, 5.139085, 5.534400, 5.929714, 6.325028, 6.720343, 7.115657, 7.510971, 7.906285]
})

file_path = 'Panel-All_train.csv'
data = pd.read_csv(file_path)

for n in [10, 20, 50]: #n = 3、5、10、20、50
    for bin_index, bin_range in bins.iterrows():
        output_file = f"bin_{bin_index}_n=_{n}_results.csv"

        for i in range(1001):
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
            else:
                output = pd.concat([output, results], ignore_index=True)

            if i == 1000:
                output.to_csv(output_file, index=False, encoding='utf-8-sig')