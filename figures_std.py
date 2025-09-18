import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


datasets = ["citeseq", "meth", "lawlor", "pbmc", "tcga"]
result_dir = "/result"
features_to_include = [1, 3, 5, 10, 15, 20]

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.5)

fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
axes = axes.flatten()


handles, labels = None, None

for i, dataset in enumerate(datasets):
    file_path = os.path.join(result_dir, f"{dataset}.csv")
    if not os.path.exists(file_path):
        print(f"Warning: file {file_path} not found, skipping.")
        continue

    df = pd.read_csv(file_path, index_col=0)
    filtered_df = df[df['nb_genes'].isin(features_to_include)]

    ax = axes[i]
    sns.pointplot(
        data=filtered_df,
        x="nb_genes",
        y="mcc",
        hue="method",
        palette="colorblind",
        markers=["o", "s", "D", "v", "^", "<", ">"],
        linestyles=["-", "--", "-.", ":", "-", "--", "-."],
        errorbar="sd",
        capsize=.1,
        ax=ax
    )

    ax.set_title(dataset.upper(), fontsize=14, weight="bold")
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)

    
    if handles is None and labels is None:
        handles, labels = ax.get_legend_handles_labels()
    ax.get_legend().remove()


fig.delaxes(axes[-1])

# Common X and Y labels
fig.text(0.5, 0.07, "Number of Top Features", ha='center', fontsize=16, weight="bold")
fig.text(0.02, 0.5, "Matthews Correlation Coefficient (MCC)", va='center', rotation='vertical',
         fontsize=16, weight="bold")

fig.legend(
    handles, labels,
    loc="lower right",
    bbox_to_anchor=(0.95, 0.05),
    ncol=2,  
    fontsize=12,
    title="Method",
    title_fontsize=13,
    frameon=True,
    facecolor='white',
    edgecolor='black',
    framealpha=1
)

plt.tight_layout(rect=[0.03, 0.08, 1, 0.95])  
plt.savefig("comparison_all_datasets.png", dpi=300)
print("Plot saved as 'comparison_all_datasets.png'")
