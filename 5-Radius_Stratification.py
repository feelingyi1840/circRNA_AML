###### Radius Stratification ######

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('all_bin_co.csv')

def calculate_distance(row):
    return np.sqrt(row['X']**2 + row['Y']**2 + row['Z']**2)

df['distance'] = df.apply(calculate_distance, axis=1)
df.to_csv('all_bin_co_with_distance.csv', index=False)

plt.figure(figsize=(10, 6))
plt.hist(df['distance'], bins=30, edgecolor='k')
plt.title('Distribution of Distances')
plt.xlabel('Distance')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()

max_distance = df['distance'].max()
bins = np.linspace(0, max_distance, 21)
merged_bins = np.concatenate(([bins[0]], bins[3:]))

# Calculate and display descriptive statistics
distance_stats = df['distance'].describe()
counts, bin_edges = np.histogram(df['distance'], bins=merged_bins)
print(merged_bins)

bin_counts_df = pd.DataFrame({
    'Bin Start': bin_edges[:-1],
    'Bin End': bin_edges[1:],
    'Count': counts
})
print(bin_counts_df)