import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
from scipy import stats

# Use glob to find all *_B_factors.csv files in the specified directory
b_factor_files = glob.glob(os.path.join('*_B_factors.csv'))

# Initialize a list to store the DataFrames
b_factor_dfs = []

# Loop through each file and read it into a DataFrame
for file in b_factor_files:
    df = pd.read_csv(file)
    # Optionally, add a new column to identify the source file
    b_factor_dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
combined_b_factor_df = pd.concat(b_factor_dfs, ignore_index=True)
combined_b_factor_df['chain'] = combined_b_factor_df['chain'].str.replace(r"[\[\]']+", '', regex=True).str.strip()
combined_b_factor_df['resn'] = combined_b_factor_df['resn'].str.replace(r"[\[\]']+", '', regex=True).str.strip()


# Display the first few rows of the combined DataFrame
print(combined_b_factor_df.head())
