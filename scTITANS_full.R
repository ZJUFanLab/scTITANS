#======== scTITANS ========================
#Usage: For identifying key genes and cell subclusters for time-series single cell sequencing data 

#Parameters:
#(1) data_file: file for gene expression matrix (row genes, column cells)
#(2) cell_metadata: file for cell metadata (row cells);if not provided, provide "none" 
#(3) gene_metadata: file for gene metadata (row genes,header must contain gene_short_name); if not provided, provide "none" 
#(4) indication: status of the matrix (normalized, filtered, none); if raw data is used, provide "none"
#(5) sp: must be Human or Mouse
#(6) tissue: organ, such as Liver; if not sure, provide "none"
#(7) root_type: the name of root cluster; if not sure, provide "none"
#(8) fdr_level: fdr level
#(9) annotation: if cell annotation is required, provide T; else, provide F

#Example:
#For a user with a raw matrix obtained from human liver with CellRanger, the following command will be ok (no cell annotation). 
# result = scTITANS_full(data_file,"none","none","none","Human","Liver","none",0.01,"F")


scTITANS_full <- function(data_file,cell_metadata,gene_metadata,indication,sp,organ,root_type,fdr_level,annotation){
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
	if (cell_metadata == "none"){
		cell_meta = matrix(colnames(data),ncol=1)
		rownames(cell_meta) = colnames(data)
		colnames(cell_meta) = "cells"
	}else{
		cell_meta = read.table(cell_file,header=T,row.names=1,sep="\t")
	}
	if (gene_metadata == "none"){
		gene_meta = matrix(rownames(data),ncol=1)
		rownames(gene_meta) = rownames(data)
		colnames(gene_meta) = "gene_short_name"
	}else{
		gene_meta = read.table(gene_file,header=T,row.names=1,sep="\t")
	}
	
	#========= check if the cell names in data and cell_meta are the same ==========================
	#====== match cells in data and cell_meta ==================
	index_data = c()
	index_meta = c()
	cells = c()
	for (jj in 1:ncol(data)){
		tmp = which(rownames(cell_meta) == colnames(data)[jj])
		if (length(tmp)!=0){
			cells = c(cells,colnames(data)[jj])
		}
	}

	for (jj in 1:length(cells)){
		t1 = which(colnames(data) == cells[jj])
		t2 = which(rownames(cell_meta) == cells[jj])
		index_data = c(index_data,t1)
		index_meta = c(index_meta,t2)
	}	

	data = data[,index_data]
	if (ncol(cell_meta) >1){
		cell_meta = cell_meta[index_meta,]
	}else{
		tmp = rownames(cell_meta)
		cell_meta = matrix(cell_meta[index_meta,],ncol=1)
		rownames(cell_meta) = tmp[index_meta]
		colnames(cell_meta) = "cell"
	}

	genes = c()
	for (ii in 1:nrow(gene_meta)){
		tmp = which(rownames(data) == rownames(gene_meta)[ii])
		if (length(tmp)!=0){
			genes = c(genes,rownames(gene_meta)[ii])
		}
	}

	indexg_data = c()
	indexg_gene = c()
	for (ii in 1:length(genes)){
		t1 = which(rownames(data)==genes[ii])
		t2 = which(rownames(gene_meta)==genes[ii])
		indexg_data = c(indexg_data,t1)
		indexg_gene = c(indexg_gene,t2)
	}
	data = data[indexg_data,]

	if (ncol(gene_meta) > 1){
		gene_meta = gene_meta[indexg_gene,]
	}else{
		tmp = rownames(gene_meta)
		gene_meta = matrix(gene_meta[indexg_gene,],ncol=1)
		rownames(gene_meta) = tmp[indexg_gene]
		colnames(gene_meta) = "gene_short_name"
	}

	data = as.matrix(data)
	meta_cell = c()	

	if (indication == "normalized"){
		data = data
	}

	if (indication == "filtered"){
		sce = SingleCellExperiment(assays=list(counts=data))
		clusters = quickCluster(sce)
		sce = computeSumFactors(sce,clusters=clusters)
		sce = logNormCounts(sce)
		data = logcounts(sce)
	}

	if (indication == "none"){
	#==== filter genes without any counts in any cells =========
		s1 = apply(data,1,sum)
		index1 = which(s1>0)
		data = data[index1,]
		if (ncol(gene_meta) > 1){
			gene_meta = gene_meta[index1,]
		}else{
			tmp = rownames(gene_meta)
			gene_meta = matrix(gene_meta[index1,],ncol=1)
			rownames(gene_meta) = tmp[index1]
			colnames(gene_meta) = "gene_short_name"
		}

	#====== a gene was defined detected if two or more transcripts were present in at least 10 cells =========
		s2 = c()
		for (i in 1:nrow(data)){
       			t1 = length(which(data[i,]>=2))
       			s2 = c(s2,t1)
		}
		index2 = which(s2>=10)
		data = data[index2,]
		if (ncol(gene_meta) > 1){
			gene_meta = gene_meta[index2,]
		}else{
			tmp = rownames(gene_meta)
			gene_meta = matrix(gene_meta[index2,],ncol=1)
			rownames(gene_meta) = tmp[index2]
			colnames(gene_meta) = "gene_short_name"
		}

	#=============cells outside the 5th and 95th percentile with respect to the number of genes detected and the number of UMIs were discarded =========
		ngene = c()
		nUMI = c()
		for (i in 1:ncol(data)){
	        t1 = length(which(data[,i]!=0))
	        ngene = c(ngene,t1)
		}
		r1 = quantile(ngene,c(0.05,0.95))
		index3 = which(ngene >= as.numeric(r1[1]) & ngene <= as.numeric(r1[2]))

		data = data[,index3]
		if (ncol(cell_meta) > 1){
			cell_meta = cell_meta[index3,]
		}else{
			tmp = rownames(cell_meta)
			cell_meta = matrix(cell_meta[index3,],ncol=1)
			rownames(cell_meta) = tmp[index3]
			colnames(cell_meta) = "cell"
		}

		nUMI = colSums(as.matrix(data))
		r2 = quantile(nUMI,c(0.05,0.95))
		index4 = which(nUMI >= as.numeric(r2[1]) & nUMI <= as.numeric(r2[2]))
		data = data[,index4]
		if (ncol(cell_meta) > 1){
			cell_meta = cell_meta[index4,]
		}else{
			tmp = rownames(cell_meta)
			cell_meta = matrix(cell_meta[index4,],ncol=1)
			rownames(cell_meta) = tmp[index4]
			colnames(cell_meta) = "cell"
		}
		
		#========cells with more than 10% of the their UMIs assigned to mitochondrial genes were filtered out =============
		proj = CreateSeuratObject(counts=data,project="proj",min.cells=3,min.features=200)
		proj[["percent.mt"]] = PercentageFeatureSet(proj,pattern="^MT-")
		index5 = which(proj[["percent.mt"]] <= 10)
		data = data[,index5]
		if (ncol(cell_meta) > 1){
			cell_meta = cell_meta[index5,]
		}else{
			tmp = rownames(cell_meta)
			cell_meta = matrix(cell_meta[index5,],ncol=1)
			rownames(cell_meta) = tmp[index5]
			colnames(cell_meta) = "cell"
		}

		sce = SingleCellExperiment(assays=list(counts=data))
		clusters = quickCluster(sce)
		sce = computeSumFactors(sce,clusters=clusters)
		sce = logNormCounts(sce)
		data = logcounts(sce)
		meta_cell = proj@meta.data
	}

	if (length(meta_cell)!=0){
		cell_meta = cbind(cell_meta,meta_cell)
	}

	cell_meta = data.frame(cell_meta)
	gene_meta = data.frame(gene_meta)

#======== Monocle3 =================
	#cds = new_cell_data_set(as(data,"sparseMatrix"),cell_metadata = cell_meta,gene_metadata=gene_meta)
	cds = new_cell_data_set(data,cell_metadata = cell_meta,gene_metadata=gene_meta)

	cds = preprocess_cds(cds,norm_method="none")
	cds = reduce_dimension(cds)
	cds = cluster_cells(cds)
	clusters = cds@clusters$UMAP$clusters
	#clusters = clusters(cds)
	clusterID = as.vector(clusters)
	names(clusterID) = names(clusters)
	pData(cds)$clusterID = clusterID

	#======== scCATCH for cell type annotation =============
	proj = CreateSeuratObject(counts=data,project="singlecell",min.cells=3,min.features=200)
	proj@meta.data$seurat_clusters = clusterID
	Idents(proj) = clusterID

	if (annotation == "T"){
		if (organ != "none"){
			clu_markers_all = findmarkergenes(proj,species=sp,match_CellMatch=TRUE,tissue=organ)
			clu_ann_all = scCATCH(clu_markers_all$clu_markers,species=sp,tissue=organ)
		}else{
			clu_markers_all = findmarkergenes(proj,species=sp,match_CellMatch=F)
			clu_ann_all = scCATCH(clu_markers_all$clu_markers,species=sp)
		}
	
		ct = clu_ann_all$cell_type
		names(ct) = clu_ann_all$cluster
		proj.anno = RenameIdents(proj,ct)
		proj@meta.data$cell_type = Idents(proj.anno)
		pData(cds)$cell.type = Idents(proj.anno)
	}else{
		ct = unique(clusterID)
		pData(cds)$cell.type = clusterID	
	}

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
	
    saveRDS(cds,file="result_scTITANS_full.rds")
	saveRDS(data,file="ExpressionData.rds")
	saveRDS(colData(cds),file="cell_metadata.rds")
	saveRDS(rowData(cds),file="gene_metadata.rds")
    #==== Identification of key genes ==============
	pseudotime <- pseudotime(cds)

	pseudo = pseudotime[which(pseudotime!="Inf")]
	index_sample = c()
	for (i in 1:length(pseudo)){
		tmp = which(colnames(data)==names(pseudo[i]))
		index_sample = c(index_sample,tmp)
	}

	data_pseudo = data[,index_sample]
	write.table(data_pseudo,"Expresson.with.pseudotime.txt",row.names=T,col.names=T,sep="\t")
	write.table(pseudo,"Pseudotime.txt",row.names=T,col.names=F,sep="\t")
	
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

	#==== Fitting curves for top 20 genes ================
	r = sort(qvalue_sigGenes,decreasing=F,index.return=T)
	qvalue_sigGenes = qvalue_sigGenes[r$ix]
	write.table(qvalue_sigGenes,sprintf("SigGenes.fdr%s.txt",fdr_level),row.names=T,col.names=F,sep="\t")

	x <- names(qvalue_sigGenes)[1:20]
	qvalue_subset <- cds[x,colnames(data_pseudo)]
	
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

	write.table(num,"NumberOfCells.in.eachcluster.along.pseudotime.txt",row.names=T,col.names=T,sep="\t")

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