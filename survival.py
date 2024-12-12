# -------------------------

genes_threshold = 0.001

survival_plot_format = 'png'
survival_plot_size = (3.5, 3.5)
survival_plot_dpi = 600
survival_plot_fontsize = 6
survival_plot_linewidth = 1
survival_plot_censorheight = 5
survival_plot_color = 'tab10'


# -------------------------

import sys
import os
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_logrank(durations, events, groups):
    unique_times = np.unique(durations[events == 1])
    unique_groups = np.unique(groups)
    n_groups = len(unique_groups)
    
    # Degrees of freedom, observed - expected, variance-covariance matrix
    df = n_groups - 1
    O_E = np.zeros(n_groups)
    V = np.zeros((n_groups, n_groups))
    
    for t in unique_times:
        at_risk = durations >= t
        events_at_t = (durations == t) & (events == 1)
        n_events = sum(events_at_t)
        n_at_risk = sum(at_risk)
        if n_events == 0 or n_at_risk == 0:
            continue
        
        for i, group in enumerate(unique_groups):
            mask = groups == group
            n_group_at_risk = sum(at_risk & mask)
            observed = sum(events_at_t & mask)
            expected = n_events * n_group_at_risk / n_at_risk
            O_E[i] += observed - expected
            
            factor = n_events * (n_at_risk - n_events) / (n_at_risk * n_at_risk * (n_at_risk - 1))
            for j, group_j in enumerate(unique_groups):
                mask_j = groups == group_j
                n_j_at_risk = sum(at_risk & mask_j)
                if i == j:
                    V[i, i] += n_group_at_risk * (n_at_risk - n_group_at_risk) * factor
                else:
                    V[i, j] += -n_group_at_risk * n_j_at_risk * factor
    
    # Chi-square statistic (df = n-1)
    O_E = O_E[:-1]
    V = V[:-1, :-1]
    try:
        chi2 = np.dot(O_E, np.linalg.solve(V, O_E))
        p_value = stats.chi2.sf(chi2, df)
    except np.linalg.LinAlgError:
        return 0, 1
    
    return chi2, p_value

def find_optimal_threshold(values, durations, events):
    if len(values) < 3:
        return np.median(values)
    
    percentiles = np.percentile(values, np.arange(20, 81, 1))
    min_p_value = float('inf')
    optimal_threshold = np.median(values)
    
    for threshold in percentiles:
        groups = (values > threshold).astype(int)
        if len(np.unique(groups)) < 2:
            continue
            
        _, p_value = calculate_logrank(durations, events, groups)
        if p_value < min_p_value:
            min_p_value = p_value
            optimal_threshold = threshold
            
    return optimal_threshold

def generate_survival_plot(durations, events, groups, group_labels, time_label, output_file):
    plt.style.use('ggplot')
    plt.figure(figsize=survival_plot_size)
    
    colors = plt.get_cmap(survival_plot_color)
    unique_groups = np.unique(groups)
    
    for i, group in enumerate(unique_groups):
        mask = (groups == group)
        group_durations = durations[mask]
        group_events = events[mask]
        
        # Sort by duration
        sort_idx = np.argsort(group_durations)
        duration_sorted = group_durations[sort_idx]
        event_sorted = group_events[sort_idx]
        
        # Calculate survival probability
        unique_times = np.unique(np.append(0, duration_sorted))
        survival = np.ones(len(unique_times))
        at_risk = len(duration_sorted)
        idx = 0
        
        censored_times = []
        censored_survivals = []
        
        for j, t in enumerate(unique_times[1:], 1):
            while idx < len(duration_sorted) and duration_sorted[idx] <= t:
                if event_sorted[idx]:
                    survival[j:] *= (at_risk - 1) / at_risk
                else:
                    censored_times.append(duration_sorted[idx])
                    censored_survivals.append(survival[j-1])
                at_risk -= 1
                idx += 1
        
        # Plot survival curve and censored points
        label = f"{group_labels[i]} (n={sum(mask)})"
        plt.step(unique_times, survival, where='post', label=label, linewidth=survival_plot_linewidth, color=colors(i))
        if censored_times:
            plt.plot(censored_times, censored_survivals, '|', 
                    color=colors(i),
                    markersize=survival_plot_censorheight, 
                    markeredgewidth=survival_plot_linewidth)
    
    plt.xlabel(time_label, fontsize=survival_plot_fontsize)
    plt.ylabel('Survival Probability', fontsize=survival_plot_fontsize)
    plt.xticks(fontsize=survival_plot_fontsize)
    plt.yticks(fontsize=survival_plot_fontsize)
    plt.legend(fontsize=survival_plot_fontsize)
    
    plt.savefig(output_file + '.' + survival_plot_format, dpi=survival_plot_dpi, bbox_inches='tight')
    plt.close()

def survival(summarize_file, time_column, event_column, group_columns):
    # Read data
    print(f"[Read Data] ...", end='\r')
    df = pd.read_csv(summarize_file, low_memory=False)
    start_gene_index = df.columns.get_loc('START_GENE')
    metadata = df.iloc[:, :start_gene_index + 1]
    gene_data = df.iloc[:, start_gene_index + 1:]
    print("[Read Data] Complete                 ")

    # Analyze variables
    print("[Analyze Variables] ...", end='\r')
    output_dir = os.path.dirname(summarize_file)
    results = []
    for column in group_columns:
        if column in metadata.columns:
            values = metadata[column]
        elif column in gene_data.columns:
            values = gene_data[column]
        elif column == '-all':
            continue
        else:
            print(f"[Warning] Column {column} not found")
            continue
            
        # Survival data
        durations = pd.to_numeric(df[time_column], errors='coerce')
        events = pd.to_numeric(df[event_column], errors='coerce')
        mask = ~(durations.isna() | events.isna() | values.isna())
        durations = durations[mask].reset_index(drop=True).values
        events = events[mask].reset_index(drop=True).values
        values = values[mask].reset_index(drop=True)
        
        # Numerical column
        if pd.api.types.is_numeric_dtype(values):
            values = values.astype(float)
            threshold = find_optimal_threshold(values, durations, events)
            groups = (values > threshold).astype(int)
            if len(np.unique(groups)) < 2:
                continue
            
            chi2, p_value = calculate_logrank(durations, events, groups)

            generate_survival_plot(
                durations,
                events,
                groups,
                [f'Low', f'High'],
                time_column,
                os.path.join(output_dir, f'survival_KM_{column}'),
            )
            results.append({
                'variable': column,
                'type': 'numerical',
                'threshold': threshold,
                'df': 1,
                'p_value': p_value,
            })
            print(f"[Analyze Variables] {column}: {p_value:.3e}                 ")

        # Categorical column
        else:
            groups = values
            unique_groups = np.unique(groups)
            if len(unique_groups) < 2:
                continue
            
            chi2, p_value = calculate_logrank(durations, events, groups)
            
            output_file = os.path.join(output_dir, f'survival_KM_{column}')
            generate_survival_plot(
                durations,
                events,
                groups,
                unique_groups,
                time_column,
                output_file,
            )
            results.append({
                'variable': column,
                'type': 'categorical',
                'threshold': "",
                'df': len(unique_groups)-1,
                'p_value': p_value
            })
            print(f"[Analyze Variables] {column}: {p_value:.3e}                 ")
    
    # Save variable results
    results_df = pd.DataFrame(results)
    results_file = os.path.join(output_dir, 'survival_variables.csv')
    results_df.to_csv(results_file, index=False)
    
    # Analyze genes
    if '-all' not in group_columns:
        return
    print("[Analyze Genes] ...", end='\r')
    genes = [gene for gene in gene_data.columns]
    gene_results = []
    total_genes = len(genes)
    
    for i, gene in enumerate(genes, 1):
        values = gene_data[gene]
        durations = pd.to_numeric(df[time_column], errors='coerce')
        events = pd.to_numeric(df[event_column], errors='coerce')
        mask = ~(durations.isna() | events.isna() | values.isna())
        durations = durations[mask].reset_index(drop=True).values
        events = events[mask].reset_index(drop=True).values
        values = values[mask].reset_index(drop=True).astype(float)
        threshold = np.percentile(values, 50)
        groups = (values > threshold).astype(int)
        if len(np.unique(groups)) < 2:
            continue
        
        chi2, p_value = calculate_logrank(durations, events, groups)
        
        gene_results.append({
            'gene': gene,
            'threshold': threshold,
            'df': 1,
            'p_value': p_value
        })
        print(f"[Analyze Genes] {i}/{total_genes}                 ", end='\r')
    
    # Save gene results
    gene_results_df = pd.DataFrame(gene_results)
    _, adjusted_pvalues, _, _ = multipletests(gene_results_df['p_value'], method='fdr_bh')
    gene_results_df['adjusted_pvalue'] = adjusted_pvalues
    gene_results_df = gene_results_df.sort_values('p_value')
    results_file = os.path.join(output_dir, 'survival_genes.csv')
    gene_results_df.to_csv(results_file, index=False)

    significant_genes = gene_results_df[gene_results_df['adjusted_pvalue'] < genes_threshold]
    print(f"[Analyze Genes] Significant genes: {len(significant_genes)}                 ")

if __name__ == "__main__":
    '''Command: python3 survival.py <cohort/summarize.csv> <time> <event> <group_column1> <group_column2> ... '''
    summarize_file = sys.argv[1]
    time_column = sys.argv[2]
    event_column = sys.argv[3]
    group_columns = sys.argv[4:]
    
    survival(summarize_file, time_column, event_column, group_columns)