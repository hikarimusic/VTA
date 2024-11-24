import sys
import os
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_logrank(durations, events, groups):
    """Calculate log rank test statistic and p-value."""
    unique_times = np.unique(durations[events == 1])
    unique_groups = np.unique(groups)
    n_groups = len(unique_groups)
    
    if n_groups < 2:
        return 0, 1
    
    # Degrees of freedom for chi-square test
    df = n_groups - 1
    
    # Calculate O-E and variance for each group
    O_E = np.zeros(n_groups)  # Observed minus Expected
    V = np.zeros((n_groups, n_groups))  # Variance-covariance matrix
    
    for t in unique_times:
        # For each time point
        at_risk = durations >= t
        events_at_t = (durations == t) & (events == 1)
        n_events = sum(events_at_t)
        n_at_risk = sum(at_risk)
        
        # Skip if no events or no one at risk
        if n_events == 0 or n_at_risk == 0:
            continue
        
        # Calculate expected values for each group
        for i, group in enumerate(unique_groups):
            mask = groups == group
            n_group_at_risk = sum(at_risk & mask)
            observed = sum(events_at_t & mask)
            expected = n_events * n_group_at_risk / n_at_risk
            
            O_E[i] += observed - expected
            
            # Calculate variance
            factor = n_events * (n_at_risk - n_events) / (n_at_risk * n_at_risk * (n_at_risk - 1))
            for j, group_j in enumerate(unique_groups):
                mask_j = groups == group_j
                n_j_at_risk = sum(at_risk & mask_j)
                
                if i == j:
                    V[i, i] += n_group_at_risk * (n_at_risk - n_group_at_risk) * factor
                else:
                    V[i, j] += -n_group_at_risk * n_j_at_risk * factor
    
    # Calculate chi-square statistic
    # Use only n-1 groups (remove last) to avoid singularity
    O_E = O_E[:-1]
    V = V[:-1, :-1]
    
    try:
        chi2 = np.dot(O_E, np.linalg.solve(V, O_E))
        p_value = stats.chi2.sf(chi2, df)
    except np.linalg.LinAlgError:
        # Handle singular matrix case
        return 0, 1
    
    return chi2, p_value

def find_optimal_threshold(values, durations, events):
    """Find threshold that gives minimum p-value in survival analysis."""
    if len(np.unique(values)) < 3:
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

def generate_survival_plot(durations, events, groups, group_labels, title, output_file, p_value):
    """Generate Kaplan-Meier survival plot."""
    plt.style.use('ggplot')
    plt.figure(figsize=(2.8, 2.8))
    
    # Get tab10 color palette
    colors = plt.get_cmap('tab10')
    unique_groups = np.unique(groups)
    
    for i, group in enumerate(unique_groups):
        mask = (groups == group)
        group_durations = durations[mask]
        group_events = events[mask]
        
        # Sort by duration
        sort_idx = np.argsort(group_durations)
        duration_sorted = group_durations[sort_idx]
        event_sorted = group_events[sort_idx]
        
        # Calculate survival curve
        unique_times = np.unique(np.append(0, duration_sorted))  # Add 0 as starting point
        survival = np.ones(len(unique_times))
        at_risk = len(duration_sorted)
        idx = 0
        
        # Store censored points
        censored_times = []
        censored_survivals = []
        
        for j, t in enumerate(unique_times[1:], 1):  # Skip the first time point (0)
            while idx < len(duration_sorted) and duration_sorted[idx] <= t:
                if event_sorted[idx]:
                    survival[j:] *= (at_risk - 1) / at_risk
                else:
                    censored_times.append(duration_sorted[idx])
                    censored_survivals.append(survival[j-1])  # Use previous survival probability
                at_risk -= 1
                idx += 1
        
        label = f"{group_labels[i]} (n={sum(mask)})"
        plt.step(unique_times, survival, where='post', label=label, linewidth=1,
                color=colors(i % 10))  # Cycle through tab10 colors
        
        # Plot censored points
        if censored_times:
            plt.plot(censored_times, censored_survivals, '|', 
                    color=colors(i % 10),  # Match line color
                    markersize=5, 
                    markeredgewidth=1)
    
    if p_value < 0.001:
        p_value_text = f'p = {p_value:.2e}'
    else:
        p_value_text = f'p = {p_value:.3f}'
    
    plt.text(0.98, 0.98, p_value_text,
             horizontalalignment='right',
             verticalalignment='top',
             transform=plt.gca().transAxes,
             fontsize=6)
    
    plt.xlabel('Time (months)', fontsize=6)
    plt.ylabel('Survival Probability', fontsize=6)
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    plt.title(title, fontsize=7)
    plt.legend(fontsize=6)
    
    plt.savefig(output_file + '.pdf', format='pdf', dpi=600, bbox_inches='tight')
    plt.savefig(output_file + '.png', format='png', dpi=600, bbox_inches='tight')
    plt.close()

def survival(summarize_file, time_column, event_column, group_columns):
    """Main function for survival analysis."""
    # Read data
    print(f"[Read Data] ...", end='\r')
    df = pd.read_csv(summarize_file, low_memory=False)
    
    # Extract survival information
    if time_column not in df.columns or event_column not in df.columns:
        print(f"Error: Missing {time_column} or {event_column} columns")
        return
    
    durations = df[time_column].values
    events = df[event_column].values
    
    # Verify event values are binary
    unique_events = np.unique(events)
    if not set(unique_events).issubset({0, 1}):
        print(f"Error: Event column should only contain binary values (0 or 1)")
        print(f"Found values: {unique_events}")
        return
    
    # Extract gene data start position
    start_gene_index = df.columns.get_loc('START_GENE')
    metadata = df.iloc[:, :start_gene_index + 1]
    gene_data = df.iloc[:, start_gene_index + 1:]
    print("[Read Data] Complete")
    
    # Create output directory
    output_dir = os.path.dirname(summarize_file)
    
    # Analyze specified columns/genes
    print("[Analyzing Variables]")
    results = []
    for column in group_columns:
        # Determine if column is in metadata or gene data
        if column in metadata.columns:
            values = metadata[column]
            source = 'metadata'
        elif column in gene_data.columns:
            values = gene_data[column]
            source = 'gene'
        else:
            print(f"[Warning] Column/Gene {column} not found")
            continue
        
        # Check if column is numerical or categorical
        if pd.api.types.is_numeric_dtype(values):
            values = values.values
            threshold = find_optimal_threshold(values, durations, events)
            groups = (values > threshold).astype(int)
            
            if len(np.unique(groups)) < 2:
                continue
            
            # Calculate log rank test
            chi2, p_value = calculate_logrank(durations, events, groups)
            
            output_file = os.path.join(output_dir, f'survival_{column}')
            generate_survival_plot(
                durations,
                events,
                groups,
                [f'Low (≤{threshold:.2f})', f'High (>{threshold:.2f})'],
                f'Survival Analysis by {column}',
                output_file,
                p_value
            )
            
            print(f"[Numerical] {column}: {p_value:.3e}")
            results.append({
                'variable': column,
                'type': f'numerical_{source}',
                'p_value': p_value,
                'threshold': threshold
            })
            
        else:
            groups = values.values
            unique_groups = np.unique(groups)
            
            if len(unique_groups) < 2:
                continue
            
            # Calculate log rank test
            chi2, p_value = calculate_logrank(durations, events, groups)
            
            output_file = os.path.join(output_dir, f'survival_{column}')
            generate_survival_plot(
                durations,
                events,
                groups,
                unique_groups,
                f'Survival Analysis by {column}',
                output_file,
                p_value
            )
            
            print(f"[Categorical] {column}: {p_value:.3e}")
            results.append({
                'variable': column,
                'type': f'categorical_{source}',
                'p_value': p_value
            })
    
    # Save variable results
    results_df = pd.DataFrame(results)
    if not results_df.empty:
        # Perform multiple testing correction if any genes were analyzed
        gene_results = results_df[results_df['type'].str.endswith('_gene')]
        if not gene_results.empty:
            _, adjusted_pvalues, _, _ = multipletests(gene_results['p_value'], method='fdr_bh')
            results_df.loc[gene_results.index, 'adjusted_p_value'] = adjusted_pvalues
        
        results_df = results_df.sort_values('p_value')
        results_file = os.path.join(output_dir, 'survival_variables.csv')
        results_df.to_csv(results_file, index=False)
    
    # Analyze remaining genes
    remaining_genes = [gene for gene in gene_data.columns if gene not in group_columns]
    if remaining_genes:
        print("\n[Analyzing Remaining Genes]")
        gene_results = []
        total_genes = len(remaining_genes)
        
        for i, gene in enumerate(remaining_genes, 1):
            values = gene_data[gene].values
            threshold = find_optimal_threshold(values, durations, events)
            groups = (values > threshold).astype(int)
            
            if len(np.unique(groups)) < 2:
                continue
            
            chi2, p_value = calculate_logrank(durations, events, groups)
            
            gene_results.append({
                'gene': gene,
                'p_value': p_value,
                'threshold': threshold
            })
            
            print(f"[Gene Analysis] {i}/{total_genes}", end='\r')
        
        # Create gene results DataFrame
        gene_results_df = pd.DataFrame(gene_results)
        if not gene_results_df.empty:
            # Multiple testing correction
            _, adjusted_pvalues, _, _ = multipletests(gene_results_df['p_value'], method='fdr_bh')
            gene_results_df['adjusted_p_value'] = adjusted_pvalues
            
            # Sort and save results
            gene_results_df = gene_results_df.sort_values('p_value')
            results_file = os.path.join(output_dir, 'survival_genes.csv')
            gene_results_df.to_csv(results_file, index=False)
            
            # Generate plots for significant genes
            significant_genes = gene_results_df[gene_results_df['adjusted_p_value'] < 0.05]
            print(f"\n[Significant Genes] {len(significant_genes)}")
            
            for _, row in significant_genes.iterrows():
                gene = row['gene']
                threshold = row['threshold']
                values = gene_data[gene].values
                groups = (values > threshold).astype(int)
                
                output_file = os.path.join(output_dir, f'survival_gene_{gene}')
                generate_survival_plot(
                    durations,
                    events,
                    groups,
                    [f'High (>{threshold:.2f})', f'Low (≤{threshold:.2f})'],
                    f'Survival Analysis by {gene}',
                    output_file,
                    row['adjusted_p_value']
                )

if __name__ == "__main__":
    '''Command: python3 survival.py <cohort/summarize.csv> <time> <event> <group_column1> <group_column2> ...'''
    if len(sys.argv) < 5:
        print("Usage: python3 survival.py <cohort/summarize.csv> <time> <event> <group_column1> <group_column2> ...")
        sys.exit(1)
    
    summarize_file = sys.argv[1]
    time_column = sys.argv[2]
    event_column = sys.argv[3]
    group_columns = sys.argv[4:]
    
    survival(summarize_file, time_column, event_column, group_columns)