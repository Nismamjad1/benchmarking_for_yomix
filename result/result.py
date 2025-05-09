import os
import pandas as pd

def save_results(results, marker_method, comaprsion_mode, method, selected_dataset, output_dir="/home/nisma/new_yomix/yomix/output"):
    """
    Saves benchmarking results in a structured format, including MCC scores and feature selections.
    
    Parameters:
    - results: Dictionary containing benchmark comparisons, signature sizes, and MCC scores.
    - output_dir: Directory to save CSV files.

    Example Structure of `results`:
    {
        "T_BRCA_vs_T_BLCA": {
            "1_features": 0.72,
            "3_features": 0.78,
            "10_features": 0.85,
            "20_features": 0.88
        },
        "T_LUAD_vs_T_LUSC": {
            "1_features": 0.65,
            "3_features": 0.70,
            "10_features": 0.82,
            "20_features": 0.90
        }
    }
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Convert dictionary into DataFrame
    all_results = []

    for benchmark, scores in results.items():
        row = {"Benchmark": benchmark}
        print(benchmark, type(scores)) 
        row.update(scores)  # Add MCC scores for each signature size
        all_results.append(row)

    # Create DataFrame
    df_results = pd.DataFrame(all_results)

   
    #  Different file naming based on method
    if marker_method == "cosg":
        output_csv = os.path.join(output_dir, f"benchmark_mcc_scores_{selected_dataset}_{marker_method}.csv")
    else:  # If using scanpy methods
        output_csv = os.path.join(output_dir, f"benchmark_mcc_scores_{selected_dataset}_{marker_method}_{method}_{comaprsion_mode}.csv")

    # Save DataFrame to CSV
    df_results.to_csv(output_csv, index=False)

    print(f" Results saved successfully at: {output_csv}")

