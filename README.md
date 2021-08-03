# scTITANS
### Identifying key genes and cell subclusters for time-series single cell sequencing data

# Usage
For users with differnt needs, two modes are provided. In the 'full' mode, only a file containing gene expression matrix (file_data), the status of the matrix (indication), the needs of cell annotation (annotation),the name of root cluster and fdr level are required.

### `scTITANS_full`
In the 'full' mode, details about the parameters are shown as follows.
### Parameters
`data_file` file for gene expression matrix (row genes, column cells)
`cell_metadata` file for cell metadata (row cells);if not provided, provide "none" 
`gene_metadata` file for gene metadata (row genes,header must contain gene_short_name); if not provided, provide "none" 
`indication` status of the matrix (normalized, filtered, none); if raw data is used, provide "none"; if the gene matrix has been filter, provide "filtered"; if the gene matrix has been filtered and noralized, provide "normalized".
`annotation` if cell annotation is required, provide T; else, provide F
`sp` must be Human or Mouse. If anntation is "F", provide ""
`tissue` organ, such as Liver; if not sure, provide "none". If annotation is "F", provide ""
`root_type` the name of root cluster; if not sure, provide "none"
`fdr_level` fdr level

### Example

For users with a raw matrix obtained from human liver with CellRanger, the following command will be ok (no cell annotation). 
In this example, besieds the file for gene expression matrix, users chose no cell annotation and fdr level of 0.01 in identifying key genes and cell subclusters. Moreover, no root cell type is provided.

```
result = scTITANS_full(data_file,"none","none","none","F","","","none",0.01)
```

For users with a raw matrix obtained from human liver with CellRanger, the following command will be ok (cell annotation required). 
In this example, besieds the file for gene expression matrix, users chose cell annotation and a fdr leve of 0.01 in identifying key genes and cell subclusters. 
In this case, users must also provide the information for sp and organ.

```
result = scTITANS_full(data_file,"none","none","none","T","Human","Liver","none",0.01)
```

The above two command will both return two txt files named "SigGenes.fdr0.01.txt" and "SigClusters.fdr0.01.txt", respectively. If something wrong happens in identifying significant clusters, an error information will be displayed.

### `scTITANS_partial`
For users who have finished trajectory inference analysis, the 'partial' mode will much more suitable.
### Parameters
`data_file` file for gene expression matrix (row genes, column cells)
`cell_metadata` file for cell metadata (row cells) 
`gene_meta` file for gene metadata (row genes,header must contain gene_short_name) 
`root_type` the name of root cluster; if not sure, provide "none"
`fdr_level` fdr level

### Example
In this example, besides the file for gene expression matrix, users must also provide files for `cell_metadata` and `gene_metadata`. In the file for `cell_metada`, information for cell clusters and/or cell types are required.
This mode can work in cases where root cell type is provided or not. 

For a user with results from trajectory inference analysis and without information about root cell type, the following command will be ok. 
```
result = scTITANS_partial(data_file,cell_meta,gene_meta,"none",0.01)
```

The above command will both return two txt files named "SigGenes.fdr0.01.txt" and "SigClusters.fdr0.01.txt", respectively. If something wrong happens in identifying significant clusters, an error information will be displayed.

# Cite
Shao et al., Identify differential genes and cell subclusters from time-series scRNA-seq data using scTITANS, Computational and Structural Biotechnology Journal, Volume 19, 2021, Pages 4132-4141, [https://doi.org/10.1016/j.csbj.2021.07.016](https://doi.org/10.1016/j.csbj.2021.07.016)