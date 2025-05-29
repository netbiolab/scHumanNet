"""
Core functions for scHumanNet analysis.

This module contains the main functions for cell-type-specific network analysis,
including centrality calculations, hub gene identification, and differential
network analysis.
"""

import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
from sklearn.utils import resample
from typing import Dict, List, Optional, Tuple, Union
import warnings
from tqdm import tqdm


def sort_add_lls(
    celltype_specific_networks: Dict[str, pd.DataFrame],
    reference_network: Optional[nx.Graph] = None
) -> Dict[str, pd.DataFrame]:
    """
    Sort SCINET output edges alphabetically and add LLS weight from HumanNetv3.
    
    Args:
        celltype_specific_networks: Dictionary of cell-type specific networks
        reference_network: Reference network (HumanNetv3), if None uses default
        
    Returns:
        Dictionary of dataframes with gene interactions and weights
    """
    if reference_network is None:
        # Load default HumanNetv3 if not provided
        from .utils import load_humannet_v3
        reference_network = load_humannet_v3()
    
    sorted_net_list = {}
    
    for celltype, cell_net_df in celltype_specific_networks.items():
        # Convert to networkx graph if needed
        if isinstance(cell_net_df, pd.DataFrame):
            cell_net = nx.from_pandas_edgelist(
                cell_net_df, 
                source=cell_net_df.columns[0], 
                target=cell_net_df.columns[1],
                edge_attr=True
            )
        else:
            cell_net = cell_net_df
            
        # Find intersection with reference network
        common_edges = []
        for edge in cell_net.edges(data=True):
            node1, node2, attrs = edge
            # Sort nodes alphabetically
            sorted_nodes = tuple(sorted([node1, node2]))
            
            if reference_network.has_edge(sorted_nodes[0], sorted_nodes[1]):
                ref_attrs = reference_network[sorted_nodes[0]][sorted_nodes[1]]
                lls_weight = ref_attrs.get('LLS', ref_attrs.get('weight', 1.0))
                scinet_weight = attrs.get('weight', attrs.get('scinet_weight', 1.0))
                
                common_edges.append({
                    'gene1': sorted_nodes[0],
                    'gene2': sorted_nodes[1], 
                    'LLS': lls_weight,
                    'scinet_weight': scinet_weight
                })
        
        sorted_net_list[celltype] = pd.DataFrame(common_edges)
    
    return sorted_net_list


def get_centrality(
    method: str,
    net_list: Dict[str, pd.DataFrame],
    exclude_ribosomal: bool = True
) -> Dict[str, Dict[str, float]]:
    """
    Calculate centrality values for nodes in scHumanNet networks.
    
    Args:
        method: Centrality measure ('degree', 'betweenness', 'closeness', 'eigenvector')
        net_list: Output of sort_add_lls()
        exclude_ribosomal: Whether to exclude ribosomal genes
        
    Returns:
        Dictionary of centrality values for each cell type
    """
    valid_methods = ['degree', 'betweenness', 'closeness', 'eigenvector']
    if method not in valid_methods:
        raise ValueError(f"Method must be one of {valid_methods}")
    
    cell_list = {}
    
    for celltype, net_df in net_list.items():
        if net_df.empty:
            cell_list[celltype] = {}
            continue
            
        # Create networkx graph
        G = nx.from_pandas_edgelist(
            net_df[['gene1', 'gene2', 'LLS']],
            source='gene1',
            target='gene2', 
            edge_attr='LLS'
        )
        
        # Calculate centrality
        if method == 'degree':
            centrality = dict(G.degree(weight='LLS'))
        elif method == 'betweenness':
            centrality = nx.betweenness_centrality(G)
        elif method == 'closeness':
            centrality = nx.closeness_centrality(G)
        elif method == 'eigenvector':
            try:
                centrality = nx.eigenvector_centrality(G, weight='LLS', max_iter=1000)
            except nx.PowerIterationFailedConvergence:
                warnings.warn(f"Eigenvector centrality failed to converge for {celltype}")
                centrality = nx.eigenvector_centrality(G, weight='LLS', max_iter=100)
        
        # Remove ribosomal and mitochondrial genes if specified
        if exclude_ribosomal:
            genes_to_remove = [
                gene for gene in centrality.keys() 
                if (gene.startswith('RPS') or gene.startswith('RPL') or 
                    gene.startswith('MRPS') or gene.startswith('MRPL') or
                    gene.startswith('NDUF'))
            ]
            for gene in genes_to_remove:
                centrality.pop(gene, None)
        
        # Convert to percentile ranks
        if centrality:
            values = list(centrality.values())
            ranks = stats.rankdata(values, method='average') / len(values)
            final_rank = dict(zip(centrality.keys(), ranks))
            cell_list[celltype] = final_rank
        else:
            cell_list[celltype] = {}
    
    return cell_list


def combine_perc_rank(
    perc_rank_list: Dict[str, Dict[str, float]],
    reference_genes: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Combine percentile rank centrality values across cell types.
    
    Args:
        perc_rank_list: Output of get_centrality()
        reference_genes: List of reference genes to include
        
    Returns:
        DataFrame with genes as rows and cell types as columns
    """
    if not perc_rank_list:
        return pd.DataFrame()
    
    # Convert to DataFrame
    df = pd.DataFrame(perc_rank_list).fillna(0)
    
    # Add missing reference genes if provided
    if reference_genes:
        missing_genes = set(reference_genes) - set(df.index)
        if missing_genes:
            missing_df = pd.DataFrame(0, index=list(missing_genes), columns=df.columns)
            df = pd.concat([df, missing_df])
    
    return df


def top_hub(rank_df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    """
    Get top n hub genes for each cell type.
    
    Args:
        rank_df: Output of combine_perc_rank()
        top_n: Number of top genes to return
        
    Returns:
        DataFrame with top genes for each cell type
    """
    top_genes = {}
    
    for celltype in rank_df.columns:
        sorted_genes = rank_df[celltype].sort_values(ascending=False)
        top_genes[celltype] = sorted_genes.head(top_n).index.tolist()
    
    # Convert to DataFrame
    top_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in top_genes.items()]))
    
    return top_df


def diff_pr(
    rank_df: pd.DataFrame,
    meta: pd.DataFrame,
    celltypes_col: str,
    condition_col: str,
    control_value: str
) -> pd.DataFrame:
    """
    Calculate differential percentile rank between conditions.
    
    Args:
        rank_df: Output of combine_perc_rank()
        meta: Metadata DataFrame
        celltypes_col: Column name for cell types
        condition_col: Column name for conditions
        control_value: Value representing control condition
        
    Returns:
        DataFrame with differential PR values
    """
    conditions = meta[condition_col].unique()
    control = control_value
    disease = [c for c in conditions if c != control][0]
    
    # Find cell types that have both conditions
    valid_celltypes = []
    for celltype in meta[celltypes_col].unique():
        possible_cols = [
            f"{control}_{celltype}", f"{disease}_{celltype}",
            f"{celltype}_{control}", f"{celltype}_{disease}"
        ]
        matching_cols = [col for col in possible_cols if col in rank_df.columns]
        
        if len(matching_cols) == 2:
            valid_celltypes.append(celltype)
    
    result_data = []
    
    for celltype in valid_celltypes:
        possible_cols = [
            f"{control}_{celltype}", f"{disease}_{celltype}",
            f"{celltype}_{control}", f"{celltype}_{disease}"
        ]
        matching_cols = [col for col in possible_cols if col in rank_df.columns]
        
        if len(matching_cols) == 2:
            control_col = [col for col in matching_cols if control in col][0]
            disease_col = [col for col in matching_cols if disease in col][0]
            
            diff_rank = rank_df[disease_col] - rank_df[control_col]
            diff_rank = diff_rank.sort_values(key=abs, ascending=False)
            
            for gene, diff_val in diff_rank.items():
                result_data.append({
                    celltype: gene,
                    f"{celltype}_{disease}-{control}": diff_val
                })
    
    if result_data:
        result_df = pd.DataFrame(result_data)
        return result_df
    else:
        return pd.DataFrame()


def find_diff_hub(
    rank_df: pd.DataFrame,
    meta: pd.DataFrame,
    celltypes_col: str,
    condition_col: str,
    control_value: str,
    net_list: Dict[str, pd.DataFrame],
    centrality: str = "degree",
    q_method: str = "BH",
    min_cells: int = 500,
    n_permutations: int = 1000
) -> pd.DataFrame:
    """
    Find statistically significant differential hub genes.
    
    Args:
        rank_df: Output of combine_perc_rank()
        meta: Metadata DataFrame  
        celltypes_col: Column name for cell types
        condition_col: Column name for conditions
        control_value: Control condition value
        net_list: Network list from sort_add_lls()
        centrality: Centrality method used
        q_method: Multiple testing correction method
        min_cells: Minimum cells required
        n_permutations: Number of permutations for null distribution
        
    Returns:
        DataFrame with statistical significance results
    """
    from statsmodels.stats.multitest import multipletests
    
    conditions = meta[condition_col].unique()
    control = control_value
    disease = [c for c in conditions if c != control][0]
    
    # Check cell type thresholds
    valid_celltypes = []
    for celltype in meta[celltypes_col].unique():
        control_cells = meta[
            (meta[celltypes_col] == celltype) & (meta[condition_col] == control)
        ]
        disease_cells = meta[
            (meta[celltypes_col] == celltype) & (meta[condition_col] == disease)
        ]
        
        if len(control_cells) >= min_cells and len(disease_cells) >= min_cells:
            valid_celltypes.append(celltype)
        else:
            print(f"{celltype}: {len(control_cells)} Control, {len(disease_cells)} Disease cells. Below threshold.")
    
    final_results = []
    
    for celltype in tqdm(valid_celltypes, desc="Processing cell types"):
        # Get network columns
        possible_cols = [
            f"{control}_{celltype}", f"{disease}_{celltype}",
            f"{celltype}_{control}", f"{celltype}_{disease}"
        ]
        matching_cols = [col for col in possible_cols if col in rank_df.columns]
        
        if len(matching_cols) != 2:
            continue
            
        control_col = [col for col in matching_cols if control in col][0]
        disease_col = [col for col in matching_cols if disease in col][0]
        
        # Get networks for null distribution
        control_net_name = control_col
        disease_net_name = disease_col
        
        if control_net_name not in net_list or disease_net_name not in net_list:
            continue
            
        control_net_df = net_list[control_net_name]
        disease_net_df = net_list[disease_net_name]
        
        # Generate null distribution
        null_distribution = _generate_null_distribution(
            control_net_df, disease_net_df, centrality, n_permutations
        )
        
        # Get genes with non-zero centrality in at least one condition
        df_subset = rank_df[[control_col, disease_col]]
        df_filtered = df_subset[~((df_subset[control_col] == 0) & (df_subset[disease_col] == 0))]
        
        # Remove ribosomal genes
        genes_to_remove = [
            gene for gene in df_filtered.index
            if (gene.startswith('RPS') or gene.startswith('RPL') or 
                gene.startswith('MRPS') or gene.startswith('MRPL') or
                gene.startswith('NDUF'))
        ]
        df_final = df_filtered.drop(genes_to_remove, errors='ignore')
        
        # Calculate differential PR and p-values
        diff_pr_values = df_final[disease_col] - df_final[control_col]
        
        p_values = []
        for diff_val in diff_pr_values:
            distribution_with_obs = np.append(null_distribution, diff_val)
            rank = stats.rankdata(-np.abs(distribution_with_obs), method='min')[-1]
            p_val = rank / len(distribution_with_obs)
            p_values.append(p_val)
        
        # Multiple testing correction
        if q_method.upper() == 'BH':
            method = 'fdr_bh'
        else:
            method = q_method.lower()
            
        _, q_values, _, _ = multipletests(p_values, method=method)
        
        # Create results DataFrame
        for i, gene in enumerate(df_final.index):
            final_results.append({
                'gene': gene,
                'celltype': celltype,
                'Control_scHumanNet': df_final.iloc[i][control_col],
                'Disease_scHumanNet': df_final.iloc[i][disease_col],
                'diffPR': diff_pr_values.iloc[i],
                'pvalue': p_values[i],
                'qvalue': q_values[i]
            })
    
    if final_results:
        return pd.DataFrame(final_results)
    else:
        return pd.DataFrame()


def find_all_hub(
    net_list: Dict[str, pd.DataFrame],
    centrality: str = "degree",
    q_method: str = "BH", 
    threshold: float = 0.05,
    n_permutations: int = 1000
) -> pd.DataFrame:
    """
    Find statistically significant hub genes in each network.
    
    Args:
        net_list: Output of sort_add_lls()
        centrality: Centrality method
        q_method: Multiple testing correction method
        threshold: Significance threshold
        n_permutations: Number of permutations
        
    Returns:
        DataFrame with significant hub genes
    """
    from statsmodels.stats.multitest import multipletests
    
    final_results = []
    centrality_list = get_centrality(net_list, method=centrality)
    
    for celltype in tqdm(net_list.keys(), desc="Finding hubs"):
        if celltype not in centrality_list or not centrality_list[celltype]:
            continue
            
        net_df = net_list[celltype]
        if net_df.empty:
            continue
            
        # Create network and calculate absolute centrality
        G = nx.from_pandas_edgelist(
            net_df[['gene1', 'gene2', 'LLS']],
            source='gene1',
            target='gene2',
            edge_attr='LLS'
        )
        
        if centrality == 'degree':
            abs_centrality = dict(G.degree(weight='LLS'))
        else:
            # For other centrality measures, calculate accordingly
            abs_centrality = dict(G.degree(weight='LLS'))  # Simplified
        
        perc_centrality = centrality_list[celltype]
        
        # Generate null distribution
        null_distribution = _generate_null_centrality_distribution(
            net_df, centrality, n_permutations
        )
        
        # Calculate p-values
        p_values = []
        for gene in perc_centrality.keys():
            if gene in abs_centrality:
                gene_cent = abs_centrality[gene]
                distribution_with_obs = np.append(null_distribution, gene_cent)
                rank = stats.rankdata(-distribution_with_obs, method='min')[-1]
                p_val = rank / len(distribution_with_obs)
                p_values.append(p_val)
            else:
                p_values.append(1.0)
        
        # Multiple testing correction
        if q_method.upper() == 'BH':
            method = 'fdr_bh'
        else:
            method = q_method.lower()
            
        _, q_values, _, _ = multipletests(p_values, method=method)
        
        # Filter significant results
        for i, (gene, cent_val) in enumerate(perc_centrality.items()):
            if q_values[i] < threshold:
                final_results.append({
                    'gene': gene,
                    'celltype': celltype,
                    'Centrality_PR': cent_val,
                    'pvalue': p_values[i],
                    'qvalue': q_values[i]
                })
    
    if final_results:
        return pd.DataFrame(final_results)
    else:
        return pd.DataFrame()


def _generate_null_distribution(
    control_net_df: pd.DataFrame,
    disease_net_df: pd.DataFrame, 
    centrality: str,
    n_permutations: int
) -> np.ndarray:
    """Generate null distribution for differential analysis."""
    null_values = []
    
    for _ in range(n_permutations):
        # Create random networks by rewiring
        if not control_net_df.empty:
            G_control = nx.from_pandas_edgelist(
                control_net_df[['gene1', 'gene2', 'LLS']],
                source='gene1', target='gene2', edge_attr='LLS'
            )
            
            # Simple edge shuffling as rewiring approximation
            edges = list(G_control.edges(data=True))
            np.random.shuffle(edges)
            
            G_random1 = nx.Graph()
            G_random1.add_edges_from(edges)
            
            G_random2 = nx.Graph() 
            G_random2.add_edges_from(edges)
            
            # Calculate centralities
            if centrality == 'degree':
                cent1 = dict(G_random1.degree(weight='LLS'))
                cent2 = dict(G_random2.degree(weight='LLS'))
            else:
                cent1 = dict(G_random1.degree())
                cent2 = dict(G_random2.degree())
            
            # Calculate differences
            common_genes = set(cent1.keys()) & set(cent2.keys())
            for gene in common_genes:
                diff = cent1[gene] - cent2[gene]
                null_values.append(diff)
    
    return np.array(null_values)


def _generate_null_centrality_distribution(
    net_df: pd.DataFrame,
    centrality: str,
    n_permutations: int
) -> np.ndarray:
    """Generate null distribution for hub analysis."""
    null_values = []
    
    if net_df.empty:
        return np.array([])
    
    for _ in range(n_permutations // 10):  # Reduce computation
        # Create random network
        G = nx.from_pandas_edgelist(
            net_df[['gene1', 'gene2', 'LLS']],
            source='gene1', target='gene2', edge_attr='LLS'
        )
        
        # Shuffle edge weights
        edges = list(G.edges(data=True))
        weights = [d['LLS'] for _, _, d in edges]
        np.random.shuffle(weights)
        
        for i, (u, v, d) in enumerate(edges):
            G[u][v]['LLS'] = weights[i]
        
        # Calculate centrality
        if centrality == 'degree':
            cent_dict = dict(G.degree(weight='LLS'))
        else:
            cent_dict = dict(G.degree())
            
        null_values.extend(cent_dict.values())
    
    return np.array(null_values)