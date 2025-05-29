# scHumanNet API Reference

## Python API

### Core Functions

#### `sort_add_lls(celltype_specific_networks, reference_network=None)`
Sort SCINET output edges alphabetically and add LLS weight from HumanNetv3.

**Parameters:**
- `celltype_specific_networks` (Dict[str, pd.DataFrame]): Dictionary of cell-type specific networks
- `reference_network` (nx.Graph, optional): Reference network (HumanNetv3)

**Returns:**
- `Dict[str, pd.DataFrame]`: Dictionary of dataframes with gene interactions and weights

#### `get_centrality(method, net_list, exclude_ribosomal=True)`
Calculate centrality values for nodes in scHumanNet networks.

**Parameters:**
- `method` (str): Centrality measure ('degree', 'betweenness', 'closeness', 'eigenvector')
- `net_list` (Dict[str, pd.DataFrame]): Output of sort_add_lls()
- `exclude_ribosomal` (bool): Whether to exclude ribosomal genes

**Returns:**
- `Dict[str, Dict[str, float]]`: Dictionary of centrality values for each cell type

#### `combine_perc_rank(perc_rank_list, reference_genes=None)`
Combine percentile rank centrality values across cell types.

**Parameters:**
- `perc_rank_list` (Dict[str, Dict[str, float]]): Output of get_centrality()
- `reference_genes` (List[str], optional): List of reference genes to include

**Returns:**
- `pd.DataFrame`: DataFrame with genes as rows and cell types as columns

#### `find_all_hub(net_list, centrality='degree', q_method='BH', threshold=0.05)`
Find statistically significant hub genes in each network.

**Parameters:**
- `net_list` (Dict[str, pd.DataFrame]): Output of sort_add_lls()
- `centrality` (str): Centrality method
- `q_method` (str): Multiple testing correction method
- `threshold` (float): Significance threshold

**Returns:**
- `pd.DataFrame`: DataFrame with significant hub genes

#### `find_diff_hub(rank_df, meta, celltypes_col, condition_col, control_value, net_list, **kwargs)`
Find statistically significant differential hub genes.

**Parameters:**
- `rank_df` (pd.DataFrame): Output of combine_perc_rank()
- `meta` (pd.DataFrame): Metadata DataFrame
- `celltypes_col` (str): Column name for cell types
- `condition_col` (str): Column name for conditions
- `control_value` (str): Control condition value
- `net_list` (Dict[str, pd.DataFrame]): Network list from sort_add_lls()

**Returns:**
- `pd.DataFrame`: DataFrame with statistical significance results

### Utility Functions

#### `load_humannet_v3(data_path=None)`
Load HumanNetv3 reference network.

#### `connectivity(gene_list, network, method='mean')`
Calculate connectivity measure for a gene list in a network.

#### `network_stats(network)`
Calculate basic network statistics.

#### `save_networks(network_list, output_path, format='csv')`
Save networks to files.

#### `load_networks(input_path, format='csv')`
Load networks from files.

---

## R API

### Core Functions

#### `SortAddLLS(Celltype.specific.networks, reference.network)`
Sort SCINET output edges alphabetically and add LLS weight derived from HumanNetv3.

**Parameters:**
- `Celltype.specific.networks`: SCINET output from function run.SCINET.clusters()
- `reference.network`: reference network input, default is HumanNetv3

**Returns:**
- List of dataframe edgelist, with gene interaction and weights from SCINET and HumanNetv3

#### `GetCentrality(method, net.list)`
Get centrality values for each nodes of scHumanNet list.

**Parameters:**
- `method`: Centrality measure ('degree', 'betweenness', 'closeness', 'eigenvector')
- `net.list`: Output of SortAddLLS()

**Returns:**
- List of named vector, each value corresponding to node's centrality value

#### `CombinePercRank(perc.rank.list)`
Calculate percentile rank of each centrality measure.

**Parameters:**
- `perc.rank.list`: Output of GetCentrality()

**Returns:**
- Dataframe of celltypes and their centrality values

#### `FindAllHub(net.list, centrality='degree', q.method='BH', threshold=0.05)`
Find statistically significant hub genes in each scHumanNets.

**Parameters:**
- `net.list`: output of SortAddLLS
- `centrality`: centrality method
- `q.method`: multiple testing correction method
- `threshold`: significance threshold

**Returns:**
- Dataframe with Percentile Rank Centrality, gene, pvalue, qvalue and celltype

#### `FindDiffHub(rank.df.final, meta, celltypes, condition, control, net.list, **kwargs)`
Calculate statistical significance of diffPR values.

**Parameters:**
- `rank.df.final`: Output from CombinePercRank
- `meta`: metadata data.frame
- `celltypes`: column name for celltypes annotation
- `condition`: column name for disease and control annotation
- `control`: character string specifying control value
- `net.list`: Output of SortAddLLS

**Returns:**
- Dataframe with percentile rank of centrality, differences, and statistical significance

### Helper Functions

#### `TopHub(rank.df.final, top.n)`
Get top n genes in terms of centrality for each scHumanNet.

#### `DiffPR(rank.df.final, celltypes, condition, control, meta)`
Get difference of normalized centrality values.

---

## Command Line Interface

The Python package includes a command-line interface:

```bash
# Find hub genes
schumannet find-hubs --networks /path/to/networks --output results.csv

# Differential analysis  
schumannet diff-analysis --networks /path/to/networks --metadata meta.csv \
  --celltype-col celltype --condition-col condition --control Control \
  --output diff_results.csv

# Calculate network statistics
schumannet network-stats --networks /path/to/networks --output stats.csv

# Convert R data to Python format
schumannet convert --input data.rda --output data.pkl
```

---

## Data Formats

### Network Format
Networks should be provided as DataFrames/data.frames with columns:
- `gene1`: First gene in interaction
- `gene2`: Second gene in interaction  
- `LLS`: LLS weight from HumanNetv3
- `scinet_weight`: Weight from SCINET (optional)

### Metadata Format
Metadata should contain:
- Cell identifiers
- Cell type annotations
- Condition/group annotations
- Additional covariates (optional)