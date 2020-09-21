#======== scTITANS_partial ========================
#Usage: For identifying key genes and cell subclusters for time-series single cell sequencing data 

#Parameters:
#(1) data_file: file for gene expression matrix (row genes, column cells)
#(2) cell_metadata: file for cell metadata (row cells) 
#(3) gene_meta: file for gene metadata (row genes,header must contain gene_short_name) 
#(4) root_type: the name of root cluster; if not sure, provide "none"
#(5) fdr_level: fdr level


#Example:
#For a user with results from trajectory inference analysis, the following command will be ok. 
# result = scTITANS_partial(data_file,cell_meta,gene_meta,"none",0.01)


scTITANS_partial <- function(data_file,cell_metadata,gene_metadata,root_type,fdr_level){
	if (!requireNamespace("BiocManager", quietly = TRUE))
    	install.packages("BiocManager")

	pcg = installed.packages()[,c('Package','Version','LibPath')]
	if (length(which(pcg == "monocle")) == 0){
		BiocManager::install("monocle")
	}
	
	if (length(which(pcg == "Biobase")) == 0){
		BiocManager::install("Biobase")
	}
	
	if (length(which(pcg == "edge")) == 0){
		BiocManager::install("edge")
	}
	
	if (length(which(pcg == "scran")) == 0){
		BiocManager::install("scran")
	}
	
	if (length(which(pcg == "scater")) == 0){
		BiocManager::install("scater")
	}
	
	if (length(which(pcg == "BiocGenerics")) == 0){
		BiocManager::install("BiocGenerics")	
	}

	if (length(which(pcg == "DelayedArray")) == 0){
		BiocManager::install("DelayedArray")
	}

	if (length(which(pcg == "DelayedMatrixStats")) == 0){
		BiocManager::install("DelayedMatrixStats")
	}

	if (length(which(pcg == "limma")) == 0){
		BiocManager::install("limma")	
	}

	if (length(which(pcg == "S4Vectors")) == 0){
		BiocManager::install("S4Vectors")	
	}

	if (length(which(pcg == "SingleCellExperiment")) == 0){
		BiocManager::install("SingleCellExperiment")	
	}

	if (length(which(pcg == "SummarizedExperiment")) == 0)(
		BiocManager::install("SummarizedExperiment")
	)

	if (length(which(pcg == "batchelor")) == 0){
		BiocManager::install("batchelor")
	}

	if (length(which(pcg == "Matrix.utils")) == 0){
		BiocManager::install("Matrix.utils")
	}

	if (length(which(pcg == "knitr")) == 0 ){
		install.packages("knitr")
	}

	if (length(which(pcg == "reshape2")) == 0){
		install.packages("reshape2")
	}
	
	if (length(which(pcg == "ggplot2")) ==0){
		install.packages("ggplot2")
	}
	
	if (length(which(pcg == "Seurat")) ==0){
		install.packages('Seurat')
	}
	
	if (length(which(pcg == "dplyr")) == 0){
		install.packages("dplyr")
	}
	
	if (length(which(pcg == "devtools")) == 0){
		install.packages("devtools")
	}
	
	if (length(which(pcg == "leidenbase")) == 0){
		devtools::install_github('cole-trapnell-lab/leidenbase')
	}
	
	if (length(which(pcg == "monocle3")) ==0){
		devtools::install_github('cole-trapnell-lab/monocle3')
	}

	library(monocle)
	library(Biobase)
	library(knitr)
	library(reshape2)
	library(ggplot2)
	library(edge)
	library(monocle3)
	library(scran)
	library(Seurat)
	library(dplyr)
	library(scater)

	if (length(which(pcg == "scCATCH")) == 0){
		install.packages(pkgs = 'devtools')
		devtools::install_github('ZJUFanLab/scCATCH')
	}
	
	library(scCATCH)

#====== read file and pre-processing ==================
	data = read.table(data_file,header=T,row.names=1,sep="\t")
	cell_meta = read.table(cell_metadata,header=T,row.names=1,sep="\t")
	gene_meta = read.table(gene_metadata,header=T,row.names=1,sep="\t")
	
	cds = new_cell_data_set(as.matrix(data),cell_metadata = cell_meta,gene_metadata=gene_meta)

	cds = preprocess_cds(cds,norm_method="none")
	cds = reduce_dimension(cds)
	cds = cluster_cells(cds)
	clusters = cds@clusters$UMAP$clusters
	#clusters = clusters(cds)
	clusterID = as.vector(clusters)
	names(clusterID) = names(clusters)
	pData(cds)$clusterID = clusterID

	ct = as.vector(pData(cds)$cell.type)
	names(ct) = clusterID

	#==== plot umap for cell clusters =======================
	cds = learn_graph(cds)

	rt = root_type
	if (rt!="none"){
		rt = tolower(rt)
		rt=gsub("[' ',',','.','*','-','--','_']","",rt)
		ctype = pData(cds)$cell.type
		celltype = c()
		for (jj in 1:length(ctype)){
			tmp = tolower(ctype[jj])
			tmp2=gsub("[' ',',','.','*','-','--','_']","",tmp)
			celltype = c(celltype,tmp2)
		}
		if (length(which(celltype==rt))>0){
			cds = order_cells(cds,root_cells=rownames(pData(cds))[which(celltype==rt)])

		}else{
			get_earliest_principal_node <- function(cds, clust_id="1"){
			cell_ids <- which(pData(cds)$clusterID == clust_id)
  
  			closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  			closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  			root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
			return(root_pr_nodes)
			}

			cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
		}
		
	}else{
			get_earliest_principal_node <- function(cds, clust_id="1"){
			cell_ids <- which(pData(cds)$clusterID == clust_id)
  
  			closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  			closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  			root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
			return(root_pr_nodes)
			}
			cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
	}
	
		#cds = order_cells(cds)
	#==== Identification of key genes ==============
	pseudotime <- pseudotime(cds)

	pseudo = pseudotime[which(pseudotime!="Inf")]
	index_sample = c()
	for (i in 1:length(pseudo)){
		tmp = which(colnames(data)==names(pseudo[i]))
		index_sample = c(index_sample,tmp)
	}

	data_pseudo = data[,index_sample]

	de_obj <- build_study(data = as.matrix(data_pseudo), tme = pseudo, sampling = "timecourse")
	full_model <- fullModel(de_obj)
	null_model <- nullModel(de_obj)
	full_matrix = fullMatrix(de_obj)
	null_matrix = nullMatrix(de_obj)
	#ef_obj <- fit_models(de_obj,stat.type = "lrt")
	de_lrt <- lrt(de_obj, nullDistn = "normal",pi0=1)
	sig_results <- qvalueObj(de_lrt)
	pvalues <- sig_results$pvalues
	qvalues <- sig_results$qvalues
	lfdr <- sig_results$lfdr
	pi0 <- sig_results$pi0
	index_sigGenes = which(qvalues < fdr_level)
	qvalue_sigGenes = qvalues[index_sigGenes]

	write.table(data_pseudo,sprintf("Data.pseudo.genes.fdr%s.txt",fdr_level),row.names=T,col.names=T,sep="\t")
	write.table(qvalues,sprintf("SigGenes.fdr%s.txt",fdr_level),row.names=T,col.names=F,sep="\t")

	#==== Fitting curves for top 20 genes ================
	r = sort(qvalue_sigGenes,decreasing=F,index.return=T)
	qvalue_sigGenes = qvalue_sigGenes[r$ix]
	write.table(qvalue_sigGenes,sprintf("SigGenes.fdr%s.txt",fdr_level),row.names=T,col.names=F,sep="\t")

	x <- names(qvalue_sigGenes)[1:20]
	qvalue_subset <- cds[x,colnames(data_pseudo)]
	#for (j in 1:20){
		#pdf(file=sprintf("Pseudotime.%s.pdf",x[j]))
		#plot_genes_in_pseudotime(qvalues_subset[j,],color_cells_by = "pseudotime")
	#}
	#dev.off()
	
	saveRDS(cds,file="result_scTITANS_full.rds")
	saveRDS(data,file="ExpressionData.rds")
	saveRDS(colData(cds),file="cell_metadata.rds")
	saveRDS(rowData(cds),file="gene_metadata.rds")
#=========== For key subclusters ===============
	r2 = sort(pseudo,decreasing=F,index.return=T)
	pseudo = pseudo[r2$ix]
	cellt = pData(cds)$cell.type
	index = c()
	for (i in 1:length(pseudo)){
		tmp = which(rownames(pData(cds))==names(pseudo)[i])
		index = c(index,tmp)
	}
	cellt = cellt[index]
	level_ct = as.vector(unique(cellt))

	rh = hist(pseudo,breaks=50)
	mids = rh$mids
	counts = rh$counts

	num = c()
	for (j in 1:length(counts)){
		if (j==1){
			ct_tmp = cellt[1:counts[1]]
		}else if (j==length(counts)){
			ct_tmp = cellt[(sum(counts[1:(length(counts)-1)])+1):sum(counts)]
		}else{
			ct_tmp = cellt[(sum(counts[1:j])+1):sum(counts[1:(j+1)])]
		}
		num_cell = c()
		for (k in 1:length(level_ct)){
			tmp = which(ct_tmp==level_ct[k])
			num_cell = c(num_cell,length(tmp))
		}
		num = cbind(num,num_cell)
	}

	rownames(num) = level_ct
	colnames(num) = mids

	names(mids) = mids
	de_obj = build_study(data=num,tme=mids,sampling="timecourse")

	full_model <- fullModel(de_obj)
	null_model <- nullModel(de_obj)
	full_matrix = fullMatrix(de_obj)
	null_matrix = nullMatrix(de_obj)
	#ef_obj <- fit_models(de_obj,stat.type = "lrt")
	de_lrt <- lrt(de_obj, nullDistn = "normal",pi0=1)
	sig_results <- qvalueObj(de_lrt)
	pvalues <- sig_results$pvalues
	qvalues <- sig_results$qvalues
	lfdr <- sig_results$lfdr
	pi0 <- sig_results$pi
	index_sigClusters = which(qvalues < fdr_level)
	qvalue_sigClusters = qvalues[index_sigClusters]

	write.table(qvalue_sigClusters,sprintf("SigClusters.fdr%s.txt",fdr_level),row.names=T,col.names=F,sep="\t")
}