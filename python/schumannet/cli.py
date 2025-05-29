"""
Command-line interface for scHumanNet.

This module provides a command-line interface for running scHumanNet analyses.
"""

import argparse
import sys
import os
import pandas as pd
import json
from typing import Dict, Any

from . import core, utils


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="scHumanNet: Cell-type-specific functional gene network analysis"
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Hub analysis command
    hub_parser = subparsers.add_parser('find-hubs', help='Find hub genes in networks')
    hub_parser.add_argument('--networks', required=True, help='Path to network files directory')
    hub_parser.add_argument('--method', default='degree', choices=['degree', 'betweenness', 'closeness', 'eigenvector'],
                           help='Centrality method')
    hub_parser.add_argument('--output', required=True, help='Output file path')
    hub_parser.add_argument('--threshold', type=float, default=0.05, help='Significance threshold')
    hub_parser.add_argument('--format', default='csv', choices=['csv', 'pickle'], help='Input format')
    
    # Differential analysis command
    diff_parser = subparsers.add_parser('diff-analysis', help='Differential network analysis')
    diff_parser.add_argument('--networks', required=True, help='Path to network files directory')
    diff_parser.add_argument('--metadata', required=True, help='Path to metadata file')
    diff_parser.add_argument('--celltype-col', required=True, help='Cell type column name')
    diff_parser.add_argument('--condition-col', required=True, help='Condition column name')
    diff_parser.add_argument('--control', required=True, help='Control condition value')
    diff_parser.add_argument('--output', required=True, help='Output file path')
    diff_parser.add_argument('--method', default='degree', help='Centrality method')
    diff_parser.add_argument('--min-cells', type=int, default=500, help='Minimum cells required')
    
    # Network stats command
    stats_parser = subparsers.add_parser('network-stats', help='Calculate network statistics')
    stats_parser.add_argument('--networks', required=True, help='Path to network files directory')
    stats_parser.add_argument('--output', required=True, help='Output file path')
    stats_parser.add_argument('--format', default='csv', help='Input format')
    
    # Convert data command
    convert_parser = subparsers.add_parser('convert', help='Convert R data to Python format')
    convert_parser.add_argument('--input', required=True, help='Input R data file')
    convert_parser.add_argument('--output', required=True, help='Output pickle file')
    
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        sys.exit(1)
    
    try:
        if args.command == 'find-hubs':
            run_hub_analysis(args)
        elif args.command == 'diff-analysis':
            run_diff_analysis(args)
        elif args.command == 'network-stats':
            run_network_stats(args)
        elif args.command == 'convert':
            run_convert(args)
        else:
            print(f"Unknown command: {args.command}")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


def run_hub_analysis(args):
    """Run hub gene analysis."""
    print(f"Loading networks from {args.networks}...")
    networks = utils.load_networks(args.networks, format=args.format)
    
    if not networks:
        raise ValueError("No networks found")
    
    print(f"Found {len(networks)} networks")
    print(f"Running hub analysis with {args.method} centrality...")
    
    # Sort and add LLS weights
    sorted_networks = core.sort_add_lls(networks)
    
    # Find hub genes
    results = core.find_all_hub(
        sorted_networks,
        centrality=args.method,
        threshold=args.threshold
    )
    
    print(f"Found {len(results)} significant hub genes")
    
    # Save results
    results.to_csv(args.output, index=False)
    print(f"Results saved to {args.output}")


def run_diff_analysis(args):
    """Run differential network analysis."""
    print(f"Loading networks from {args.networks}...")
    networks = utils.load_networks(args.networks, format='csv')
    
    print(f"Loading metadata from {args.metadata}...")
    metadata = pd.read_csv(args.metadata)
    
    if not networks:
        raise ValueError("No networks found")
    
    print(f"Found {len(networks)} networks")
    
    # Sort and add LLS weights
    sorted_networks = core.sort_add_lls(networks)
    
    # Calculate centralities
    print("Calculating centralities...")
    centralities = core.get_centrality(args.method, sorted_networks)
    
    # Combine percentile ranks
    rank_df = core.combine_perc_rank(centralities)
    
    # Run differential analysis
    print("Running differential analysis...")
    results = core.find_diff_hub(
        rank_df=rank_df,
        meta=metadata,
        celltypes_col=args.celltype_col,
        condition_col=args.condition_col,
        control_value=args.control,
        net_list=sorted_networks,
        centrality=args.method,
        min_cells=args.min_cells
    )
    
    print(f"Found {len(results)} differential results")
    
    # Save results
    results.to_csv(args.output, index=False)
    print(f"Results saved to {args.output}")


def run_network_stats(args):
    """Calculate and save network statistics."""
    print(f"Loading networks from {args.networks}...")
    networks = utils.load_networks(args.networks, format=args.format)
    
    if not networks:
        raise ValueError("No networks found")
    
    print(f"Calculating statistics for {len(networks)} networks...")
    
    stats_results = []
    for name, network in networks.items():
        stats = utils.network_stats(network)
        stats['network_name'] = name
        stats_results.append(stats)
    
    stats_df = pd.DataFrame(stats_results)
    
    # Save results
    stats_df.to_csv(args.output, index=False)
    print(f"Network statistics saved to {args.output}")


def run_convert(args):
    """Convert R data to Python format."""
    print(f"Converting {args.input} to {args.output}...")
    utils.convert_r_data_to_pickle(args.input, args.output)
    print("Conversion complete")


if __name__ == '__main__':
    main()