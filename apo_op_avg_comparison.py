import pandas as pd
import matplotlib.pyplot as plt
import os

# File paths
all_path = "/.../all_op_combined.csv"
binding_path = "/...//binding_combined.csv"

# Load data
df_all = pd.read_csv(all_path)
df_binding = pd.read_csv(binding_path)

# Make sure s2calc is numeric
df_all['s2calc'] = pd.to_numeric(df_all['s2calc'], errors='coerce')
df_binding['s2calc'] = pd.to_numeric(df_binding['s2calc'], errors='coerce')

# Drop any NaNs just in case
df_all = df_all.dropna(subset=['s2calc'])
df_binding = df_binding.dropna(subset=['s2calc'])

# Calculate average s2calc
avg_all = df_all['s2calc'].mean()
avg_binding = df_binding['s2calc'].mean()
delta = avg_all - avg_binding

# Print summary
print(f"Average s2calc (All): {avg_all:.4f}")
print(f"Average s2calc (Binding): {avg_binding:.4f}")
print(f"Δs2calc (All - Binding): {delta:.4f}")

# Save summary
with open("/dors/wankowicz_lab/ellas/apo_op/close_resi_filtered/avg_s2calc_delta.txt", "w") as f:
    f.write(f"Average s2calc (All): {avg_all:.4f}\n")
    f.write(f"Average s2calc (Binding): {avg_binding:.4f}\n")
    f.write(f"Δs2calc (All - Binding): {delta:.4f}\n")

# Boxplot
plt.figure(figsize=(8, 5))
plt.boxplot([df_all['s2calc'], df_binding['s2calc']],
            labels=["All Residues", "Binding Site Residues"],
            patch_artist=True,
            boxprops=dict(facecolor='lightblue'),
            medianprops=dict(color='red'))
plt.title("Boxplot of s2calc")
plt.ylabel("s2calc")
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig("/.../s2calc_boxplot.png", dpi=300)
plt.show()

# Histogram of delta (optional, just to visualize the shift)
plt.figure(figsize=(8, 5))
plt.hist(df_all['s2calc'], bins=30, alpha=0.6, label='All', color='gray')
plt.hist(df_binding['s2calc'], bins=30, alpha=0.6, label='Binding', color='purple')
plt.title("Histogram of s2calc Values")
plt.xlabel("s2calc")
plt.ylabel("Frequency")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig("/.../s2calc_histogram.png", dpi=300)
plt.show()
