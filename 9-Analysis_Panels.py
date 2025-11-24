# Function: Analysis of different circRNA panels -- Part 1.
# Programmer: Zhangli Yuan, Yi Shi, 20251124.

#### Comparison of correlation coefficients ####
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

file_paths = {
    'control_3d_spearman.csv': 'Panel-TopN',
    'go_path_51.csv': 'Panel-Pathway',
    'circ_3d_spearman.csv': 'Panel-3DG-Cluster',
    'circ_radius_5_spearman.csv': 'Panel-3DG-Radius5',
}

data_com = pd.DataFrame()
for file_path, file in file_paths.items():
    data = pd.read_csv(file_path)
    data["Group"] = file
    data0 = data[['gene_correlations','Group']]
    data_com = pd.concat([data_com, data0], ignore_index=True)
 
plt.figure(figsize=(10, 6))

sns.violinplot(x='Group', y='gene_correlations', data=data_com, hue="Group")  
plt.xlabel(None)
plt.ylabel(None)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title('Spearman correlation coefficient between circRNA and label in each panel',fontsize=18)  
plt.show()


#### Heatmap ####
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

file_paths = {
    'Panel-TopN.csv': 'Panel-TopN',
    'Panel-Pathway_51.csv': 'Panel-Pathway',
    'Panel-3DG.csv': 'Panel-3DG-Cluster',
    'Panel-Radius5.csv': 'Panel-3DG-Radius5',
    'Panel-Radius1.csv': 'Panel-3DG-Radius1',
}

for file_path, file in file_paths.items():
    data = pd.read_csv(file_path)
    data = data.iloc[:,1:]
    new_df = data.corr()
    mask = np.triu(np.ones_like(new_df, dtype=bool))
    plt.figure(figsize=(25, 18))  
    sns.heatmap(new_df, 
                annot=True, fmt=".1f",
                cmap="coolwarm", center=0,   
                xticklabels=False, yticklabels=False)
    plt.title(file, fontsize=100)
    plt.figure(1)

    triu_indices = np.triu_indices_from(new_df, k=1)  
    abs_df = new_df.abs()
    average_correlation = (abs_df.values[triu_indices]).mean() 
    print(file,f'Average Correlation Coefficient (excluding self-correlation): {average_correlation:.3f}')  
