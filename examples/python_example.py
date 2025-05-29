#!/usr/bin/env python3
"""
Example usage of scHumanNet Python package.

This script demonstrates how to use the scHumanNet Python package for
cell-type-specific network analysis.
"""

import pandas as pd
import numpy as np
from schumannet import (
    sort_add_lls,
    get_centrality,
    combine_perc_rank,
    find_all_hub,
    find_diff_hub,
    network_stats
)


def create_example_data():
    """Create example network and metadata for demonstration."""
    
    # Create example cell-type-specific networks
    np.random.seed(42)
    
    # Generate example gene names
    genes = [f"GENE{i:04d}" for i in range(1, 1001)]
    
    networks = {}
    
    # Create networks for different cell types and conditions
    for celltype in ['Neuron', 'Astrocyte', 'Microglia']:
        for condition in ['Control', 'Disease']:
            # Generate random network
            n_edges = np.random.randint(200, 500)
            
            network_data = []
            for _ in range(n_edges):
                gene1 = np.random.choice(genes)
                gene2 = np.random.choice(genes)
                if gene1 != gene2:
                    network_data.append({
                        'gene1': gene1,
                        'gene2': gene2,
                        'LLS': np.random.uniform(0.1, 1.0),
                        'scinet_weight': np.random.uniform(0.1, 1.0)
                    })
            
            network_name = f"{condition}_{celltype}"
            networks[network_name] = pd.DataFrame(network_data)
    
    # Create example metadata
    n_cells = 3000
    metadata = []
    
    for i in range(n_cells):
        celltype = np.random.choice(['Neuron', 'Astrocyte', 'Microglia'])
        condition = np.random.choice(['Control', 'Disease'])
        
        metadata.append({
            'cell_id': f"Cell_{i:04d}",
            'celltype': celltype,
            'condition': condition
        })
    
    metadata_df = pd.DataFrame(metadata)
    
    return networks, metadata_df


def main():
    """Run the example analysis."""
    
    print("scHumanNet Python Example")
    print("=" * 40)
    
    # Create example data
    print("Creating example data...")
    networks, metadata = create_example_data()
    
    print(f"Created {len(networks)} networks")
    print(f"Created metadata for {len(metadata)} cells")
    
    # Sort networks and add LLS weights
    print("\nSorting networks and adding LLS weights...")
    sorted_networks = sort_add_lls(networks)
    
    # Calculate network statistics
    print("\nCalculating network statistics...")
    for name, net_df in sorted_networks.items():
        stats = network_stats(net_df)
        print(f"{name}: {stats['nodes']} nodes, {stats['edges']} edges")
    
    # Calculate centralities
    print("\nCalculating degree centrality...")
    centralities = get_centrality('degree', sorted_networks)
    
    # Combine percentile ranks
    print("Combining percentile ranks...")
    rank_df = combine_perc_rank(centralities)
    print(f"Combined centrality matrix: {rank_df.shape}")
    
    # Find hub genes
    print("\nFinding significant hub genes...")
    hub_results = find_all_hub(
        sorted_networks,
        centrality='degree',
        threshold=0.05
    )
    
    print(f"Found {len(hub_results)} significant hub genes")
    if len(hub_results) > 0:
        print("Top 5 hub genes:")
        top_hubs = hub_results.nlargest(5, 'Centrality_PR')
        for _, row in top_hubs.iterrows():
            print(f"  {row['gene']} ({row['celltype']}): PR={row['Centrality_PR']:.3f}, q={row['qvalue']:.3e}")
    
    # Differential analysis
    print("\nRunning differential network analysis...")
    try:
        diff_results = find_diff_hub(
            rank_df=rank_df,
            meta=metadata,
            celltypes_col='celltype',
            condition_col='condition',
            control_value='Control',
            net_list=sorted_networks,
            min_cells=100  # Lower threshold for example
        )
        
        print(f"Found {len(diff_results)} differential results")
        if len(diff_results) > 0:
            # Show top differential genes
            top_diff = diff_results.nlargest(5, lambda x: abs(x), subset=['diffPR'])
            print("Top 5 differential genes:")
            for _, row in top_diff.iterrows():
                print(f"  {row['gene']} ({row['celltype']}): diffPR={row['diffPR']:.3f}, q={row['qvalue']:.3e}")
                
    except Exception as e:
        print(f"Differential analysis failed: {e}")
    
    print("\nExample completed successfully!")


if __name__ == "__main__":
    main()