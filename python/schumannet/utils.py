"""
Utility functions for scHumanNet analysis.

This module contains helper functions for loading data, network operations,
and other utilities.
"""

import numpy as np
import pandas as pd
import networkx as nx
from typing import Dict, List, Optional, Tuple, Union
import pickle
import os


def load_humannet_v3(data_path: Optional[str] = None) -> nx.Graph:
    """
    Load HumanNetv3 reference network.
    
    Args:
        data_path: Path to HumanNetv3 data file
        
    Returns:
        NetworkX graph of HumanNetv3
    """
    if data_path is None:
        # Try to load from package data directory
        package_dir = os.path.dirname(__file__)
        data_path = os.path.join(package_dir, '..', '..', 'data', 'HNv3_XC.rda')
    
    # This would need to be implemented based on the actual data format
    # For now, return an empty graph as placeholder
    G = nx.Graph()
    
    # TODO: Implement actual loading of R data files
    # This might require rpy2 or conversion to pickle/json format
    print(f"Warning: HumanNetv3 loading not yet implemented. Using empty graph.")
    
    return G


def connectivity(
    gene_list: List[str],
    network: nx.Graph,
    method: str = "mean"
) -> float:
    """
    Calculate connectivity measure for a gene list in a network.
    
    Args:
        gene_list: List of genes to analyze
        network: NetworkX graph
        method: Method for calculating connectivity ('mean', 'sum', 'density')
        
    Returns:
        Connectivity value
    """
    if not gene_list or network.number_of_nodes() == 0:
        return 0.0
    
    # Get subgraph of genes present in network
    present_genes = [gene for gene in gene_list if gene in network.nodes()]
    
    if len(present_genes) < 2:
        return 0.0
    
    subgraph = network.subgraph(present_genes)
    
    if method == "mean":
        # Mean degree in subgraph
        degrees = [subgraph.degree(node) for node in subgraph.nodes()]
        return np.mean(degrees) if degrees else 0.0
    
    elif method == "sum":
        # Total edges in subgraph
        return subgraph.number_of_edges()
    
    elif method == "density":
        # Network density
        return nx.density(subgraph)
    
    else:
        raise ValueError(f"Unknown method: {method}")


def deconvolute_net(
    network_list: Dict[str, pd.DataFrame],
    gene_signatures: Dict[str, List[str]],
    method: str = "connectivity"
) -> pd.DataFrame:
    """
    Deconvolute networks using gene signatures.
    
    Args:
        network_list: Dictionary of networks
        gene_signatures: Dictionary of gene signatures
        method: Deconvolution method
        
    Returns:
        DataFrame with deconvolution scores
    """
    results = []
    
    for net_name, net_df in network_list.items():
        if net_df.empty:
            continue
            
        # Convert to networkx
        G = nx.from_pandas_edgelist(
            net_df[['gene1', 'gene2']],
            source='gene1',
            target='gene2'
        )
        
        for sig_name, gene_list in gene_signatures.items():
            if method == "connectivity":
                score = connectivity(gene_list, G, method="density")
            else:
                # Could implement other deconvolution methods
                score = connectivity(gene_list, G, method="mean")
            
            results.append({
                'network': net_name,
                'signature': sig_name,
                'score': score
            })
    
    return pd.DataFrame(results)


def filter_networks_by_size(
    network_list: Dict[str, pd.DataFrame],
    min_edges: int = 100,
    max_edges: Optional[int] = None
) -> Dict[str, pd.DataFrame]:
    """
    Filter networks by size criteria.
    
    Args:
        network_list: Dictionary of networks
        min_edges: Minimum number of edges
        max_edges: Maximum number of edges (optional)
        
    Returns:
        Filtered network dictionary
    """
    filtered_networks = {}
    
    for name, net_df in network_list.items():
        n_edges = len(net_df)
        
        if n_edges < min_edges:
            continue
            
        if max_edges is not None and n_edges > max_edges:
            continue
            
        filtered_networks[name] = net_df
    
    return filtered_networks


def calculate_overlap(
    gene_list1: List[str],
    gene_list2: List[str],
    method: str = "jaccard"
) -> float:
    """
    Calculate overlap between two gene lists.
    
    Args:
        gene_list1: First gene list
        gene_list2: Second gene list
        method: Overlap method ('jaccard', 'overlap_coefficient', 'dice')
        
    Returns:
        Overlap score
    """
    set1 = set(gene_list1)
    set2 = set(gene_list2)
    
    intersection = len(set1 & set2)
    
    if method == "jaccard":
        union = len(set1 | set2)
        return intersection / union if union > 0 else 0.0
    
    elif method == "overlap_coefficient":
        min_size = min(len(set1), len(set2))
        return intersection / min_size if min_size > 0 else 0.0
    
    elif method == "dice":
        return 2 * intersection / (len(set1) + len(set2)) if (len(set1) + len(set2)) > 0 else 0.0
    
    else:
        raise ValueError(f"Unknown method: {method}")


def network_stats(network: Union[nx.Graph, pd.DataFrame]) -> Dict[str, float]:
    """
    Calculate basic network statistics.
    
    Args:
        network: NetworkX graph or edge list DataFrame
        
    Returns:
        Dictionary of network statistics
    """
    if isinstance(network, pd.DataFrame):
        if network.empty:
            return {"nodes": 0, "edges": 0, "density": 0.0, "avg_degree": 0.0}
            
        G = nx.from_pandas_edgelist(
            network,
            source=network.columns[0],
            target=network.columns[1]
        )
    else:
        G = network
    
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    
    if n_nodes == 0:
        return {"nodes": 0, "edges": 0, "density": 0.0, "avg_degree": 0.0}
    
    density = nx.density(G)
    avg_degree = 2 * n_edges / n_nodes if n_nodes > 0 else 0.0
    
    stats_dict = {
        "nodes": n_nodes,
        "edges": n_edges,
        "density": density,
        "avg_degree": avg_degree
    }
    
    # Add connected components info
    if not G.is_directed():
        n_components = nx.number_connected_components(G)
        largest_cc_size = len(max(nx.connected_components(G), key=len)) if n_components > 0 else 0
        stats_dict.update({
            "connected_components": n_components,
            "largest_component_size": largest_cc_size
        })
    
    return stats_dict


def save_networks(
    network_list: Dict[str, pd.DataFrame],
    output_path: str,
    format: str = "csv"
) -> None:
    """
    Save networks to files.
    
    Args:
        network_list: Dictionary of networks
        output_path: Output directory path
        format: Output format ('csv', 'pickle', 'graphml')
    """
    os.makedirs(output_path, exist_ok=True)
    
    for name, net_df in network_list.items():
        if format == "csv":
            filepath = os.path.join(output_path, f"{name}.csv")
            net_df.to_csv(filepath, index=False)
        
        elif format == "pickle":
            filepath = os.path.join(output_path, f"{name}.pkl")
            net_df.to_pickle(filepath)
        
        elif format == "graphml":
            if not net_df.empty:
                G = nx.from_pandas_edgelist(
                    net_df,
                    source=net_df.columns[0],
                    target=net_df.columns[1],
                    edge_attr=True
                )
                filepath = os.path.join(output_path, f"{name}.graphml")
                nx.write_graphml(G, filepath)
        
        else:
            raise ValueError(f"Unknown format: {format}")


def load_networks(
    input_path: str,
    format: str = "csv",
    pattern: str = "*.csv"
) -> Dict[str, pd.DataFrame]:
    """
    Load networks from files.
    
    Args:
        input_path: Input directory path
        format: Input format ('csv', 'pickle', 'graphml')
        pattern: File pattern to match
        
    Returns:
        Dictionary of loaded networks
    """
    import glob
    
    if format == "csv":
        pattern = pattern or "*.csv"
        files = glob.glob(os.path.join(input_path, pattern))
        networks = {}
        
        for filepath in files:
            name = os.path.splitext(os.path.basename(filepath))[0]
            networks[name] = pd.read_csv(filepath)
        
        return networks
    
    elif format == "pickle":
        pattern = pattern or "*.pkl"
        files = glob.glob(os.path.join(input_path, pattern))
        networks = {}
        
        for filepath in files:
            name = os.path.splitext(os.path.basename(filepath))[0]
            networks[name] = pd.read_pickle(filepath)
        
        return networks
    
    else:
        raise ValueError(f"Format {format} not implemented yet")


def convert_r_data_to_pickle(r_data_path: str, output_path: str) -> None:
    """
    Convert R data files to pickle format.
    
    Args:
        r_data_path: Path to R data file (.rda)
        output_path: Output pickle file path
    """
    try:
        import rpy2.robjects as robjects
        from rpy2.robjects import pandas2ri
        
        pandas2ri.activate()
        
        # Load R data
        robjects.r['load'](r_data_path)
        
        # This would need to be customized based on the specific R objects
        # For now, just print available objects
        print("Available R objects:", list(robjects.r.ls()))
        
    except ImportError:
        print("rpy2 not available. Please install rpy2 to convert R data files.")
        print("Alternative: Convert R data to CSV format manually and use load_networks()")


def validate_network_format(network_df: pd.DataFrame) -> bool:
    """
    Validate network DataFrame format.
    
    Args:
        network_df: Network DataFrame to validate
        
    Returns:
        True if format is valid
    """
    if network_df.empty:
        return True
    
    # Check minimum columns
    if len(network_df.columns) < 2:
        return False
    
    # Check for required gene columns
    gene_cols = ['gene1', 'gene2']
    if not all(col in network_df.columns for col in gene_cols):
        # Try with first two columns
        if len(network_df.columns) >= 2:
            network_df.columns = list(network_df.columns[:2]) + list(network_df.columns[2:])
            return True
        return False
    
    # Check for weight columns
    weight_cols = ['LLS', 'weight', 'scinet_weight']
    has_weights = any(col in network_df.columns for col in weight_cols)
    
    return has_weights