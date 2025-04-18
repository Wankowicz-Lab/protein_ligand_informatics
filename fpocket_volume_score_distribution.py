import pandas as pd
import matplotlib.pyplot as plt
import os
import glob

# Define the directory containing the CSV files
directory = '/dors/wankowicz_lab/all_pdb/fpocket'  # Path to the directory with your CSV files
# Define the directory to save the histograms
save_directory = '/dors/wankowicz_lab/ellas'  # Replace with the directory where you want to save the plots

# Ensure the save directory exists
os.makedirs(save_directory, exist_ok=True)

# Get all the pocket CSV files ending with 'fpocket_pocket1.csv'
all_files = glob.glob(os.path.join(directory, '*_pocket1.csv'))  # This line fetches the specific files

li = []

# Iterate through the list of files and read them in
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    
    # Extract PDB ID from the filename (optional)
    df['PDB'] = filename.split('/')[-1].split('_')[0]
    
    li.append(df)

# Concatenate all the dataframes into a single one
pocket_info = pd.concat(li, axis=0, ignore_index=True)

# Reshape the dataframe: pivot the table so that Attributes become columns
pocket_info_pivot = pocket_info.pivot_table(index=['PDB', 'Pocket'], columns='Attribute', values='Value', aggfunc='first')

# Flatten the columns to avoid multi-level columns
pocket_info_pivot.columns = [col for col in pocket_info_pivot.columns]

# Reset index to make PDB and Pocket regular columns
pocket_info_pivot.reset_index(inplace=True)

# Show the first few rows to verify the transformation
print(pocket_info_pivot.head())

# Initialize lists to store 'Volume' and 'Score' from all files
all_volumes = []
all_scores = []

# Iterate over the reshaped dataframe to extract 'Volume' and 'Score' values
for index, row in pocket_info_pivot.iterrows():
    # Check if 'Volume' and 'Score' columns exist and extract values
    if 'Volume' in row and pd.notnull(row['Volume']):
        all_volumes.append(row['Volume'])
    if 'Score' in row and pd.notnull(row['Score']):
        all_scores.append(row['Score'])

# Convert lists to float for plotting
all_volumes = list(map(float, all_volumes))
all_scores = list(map(float, all_scores))

# Plot and save histogram for 'Volume'
if all_volumes:
    plt.figure(figsize=(10, 5))
    # Set dynamic x-axis limits based on min and max values
    volume_min, volume_max = min(all_volumes), max(all_volumes)
    plt.hist(all_volumes, bins=50, range=(volume_min, volume_max), color='skyblue', edgecolor='black', alpha=0.7)
    plt.title('Distribution of Volume')
    plt.xlabel('Volume')
    plt.ylabel('Frequency')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(save_directory, 'volume_histogram.png'))
    plt.close()  # Close the figure to avoid overlap
else:
    print("No 'Volume' data available for plotting.")

# Plot and save histogram for 'Score'
if all_scores:
    plt.figure(figsize=(10, 5))
    # Set dynamic x-axis limits based on min and max values
    score_min, score_max = min(all_scores), max(all_scores)
    plt.hist(all_scores, bins=50, range=(score_min, score_max), color='salmon', edgecolor='black', alpha=0.7)
    plt.title('Distribution of Score')
    plt.xlabel('Score')
    plt.ylabel('Frequency')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(save_directory, 'score_histogram.png'))
    plt.close()  # Close the figure to avoid overlap
else:
    print("No 'Score' data available for plotting.")
