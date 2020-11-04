# CellPhoneDB for scanpy

Usage:

```
import cellphonedb as cphdb
outs=cphdb.statistical_analysis_scanpy(adata, adata.var_names, adata.obs_names, obs_key)
```

`obs_key` must be some key in adata.obs that groups the cells, such as 'celltype', 'louvain' or 'cluster'.

There is a jupyter notebook in here somewhere which shows that this works fine.
