rm(list=ls())
gc()

## package
{
  library(data.table)
  library(dplyr)
  library(stringr)
  library(openxlsx)
  library(clusterProfiler)
  library(tidyr)
  
  
  ##生信分析
  {
    #生存分析
    library(survival)
    library(survminer)
    library(phenoTest) #smooth HR
    ##GSEA
    library(enrichplot)
    library(GSVA)
    ##gsea
    library(fgsea)
    library(msigdbr)
    
    library(ggpubr) #计算boxplot显著性
    library(ggsignif)
    
    #计算OR值用的包
    {
      # #sscVis安装：rlang版本过低，安装过程中报错了；缺少impute包
      # remove.packages('rlang')
      # install.packages('rlang')
      # BiocManager::install("impute")
      # remotes::install_github("Japrin/sscVis")
      
      # library("sscVis")
      # library("data.table")
      library("grid")
      library("cowplot")
      library("ggrepel")
      library("readr")
      library("plyr")
      library("ggpubr")
      # library("ggplot2")
    }
    
    #GO simplify
    # library(devtools)
    # install_github("jokergoo/simplifyEnrichment")
    # library(simplifyEnrichment)
    

    
    ## 计算Gini 系数的R包
    library(DescTools)
  }
  
  ##单细胞分析
  {
    library(Seurat)
    
    # library(devtools)
    # install_github("immunogenomics/harmony")
    library(harmony)
    
    #去除双细胞
    # BiocManager::install('harmony',force = TRUE)
    # devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
    library(DoubletFinder)
    
    #RNA速率
    # library(devtools)
    # BiocManager::install("pcaMethods")
    # install_github("velocyto-team/velocyto.R")
    # library(velocyto.R)
    
    ## 分类树
    # BiocManager::install('clustree')
    # library(clustree)
    
    ##surat示例数据
    # devtools::install_github('satijalab/seurat-data')
    # library(SeuratData)
    
    #infercnv  
    ##安装记录：需要到下面这个网站下载JAGS软件，建立与JAGS库的连接，才能调用rjags包，才能调用infercnv
    # # https://sourceforge.net/projects/mcmc-jags/
    # BiocManager::install("rjags")
    # # install.packages('rjags')
    Sys.setenv(JAGS_HOME="C:\\Program Files\\JAGS\\JAGS-4.3.1")
    # install.packages("rjags","http://rforge.net/",type="source")
    library(rjags)
    # BiocManager::install('infercnv')
    library(infercnv)
    
    
    ##Seurat 数据转换成 h5ad数据格式
    # remotes::install_github('mojaveazure/seurat-disk')
    library(SeuratDisk)      # SaveH5Seurat
    

    ##读取hd5数据
    # BiocManager::install('hdf5r')
    library(hdf5r)
    library(Matrix)
    
    ##cellchat
    # devtools::install_github("sqjin/CellChat")
    # ComplexHeatmap
    # BiocNeighbors
    library(CellChat)
    
    # devtools::install_github("Japrin/STARTRAC")
    # 张泽民组的包 可以计算R o/e值
    library("Startrac")
    
    ##AUCell 单细胞计算基因集合得分
    # To support paralell execution:
    # BiocManager::install(c("doMC", "doRNG","doSNOW"))
    # BiocManager::install('AUCell')
    library(AUCell)
    
    #安装包 画单细胞密度图,只能画 基因的密度
    # BiocManager::install("Nebulosa")
    # devtools::install_github("powellgenomicslab/Nebulosa")
    #加载R包
    # library(Nebulosa)
    
    ## pseudobulk（muscat包）给count值加和
    # BiocManager::install("densvis")
    # #下载到 global文件夹本地安装了
    # BiocManager::install("scater")
    # #下载到 global文件夹本地安装了
    # install.packages('variancePartition')
    # #下载到 global文件夹本地安装了
    # install.packages('muscat')
    # library(muscat)
    
    ##张泽民团队R包，用来判断细胞注释是否需要细分
    # devtools::install_github("PaulingLiu/ROGUE")
    library(ROGUE)
  }
  
  ##单细胞VDJ分析
  {
    #TCR
    # BiocManager::install('scRepertoire')
    library(scRepertoire)
    
    ##TCR BCR 分析
    # BiocManager::install('alakazam')
    library(alakazam)
  }
  
  ##单细胞轨迹分析
  {
    #2和3 冲突，用二的时候需要关掉3
    ##monocle2
    # BiocManager::install("monocle")
    library(monocle)
    
    ##monocle3 
    # if (!requireNamespace("BiocManager", quietly = TRUE))
    #   install.packages("BiocManager")
    # BiocManager::install(version = "3.10")
    # BiocManager::install('terra')
    # BiocManager::install('SingleCellExperiment')
    # BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
    #                        'limma', 'S4Vectors', 'SingleCellExperiment',
    #                        'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
    # install.packages("devtools")
    # devtools::install_github('cole-trapnell-lab/leidenbase')
    # install.packages('igraph')
    # devtools::install_github('cole-trapnell-lab/monocle3')
    # library(monocle3)
    
    ##diffusion map-- destiny R包
    # BiocManager::install("destiny")
    library(SingleCellExperiment)
    library(destiny)
    # browseVignettes("destiny")
    
    ##CytoTRACE
    # ERROR: dependencies 'HiClimR', 'ccaPP', 'nnls', 'egg' are not available for package 'CytoTRACE'
    # devtools::install_local(file.path(cloud_path,'cytoTRACE',"CytoTRACE_0.3.3.tar.gz"))
    library(CytoTRACE) 
    
    #类似cytoTRACE 也是评估分化
    # https://github.com/aet21/SCENT
    # devtools::install_github("aet21/SCENT")
    # library(SCENT)
  }

  ##画图
  {
    library(ggplot2)
    library(ggthemes) #色盘
    library(scales) #查看色块
    
    library(patchwork) #拼图
    ##纵坐标截断画图
    # install.packages("gg.gap")
    # library(gg.gap)
    
    ##对半小提琴图
    # library(devtools)
    # install_github("JanCoUnchained/ggunchained")
    # library(ggunchained)
    
    #画对半小提琴图或boxplot
    # BiocManager::install('gghalves')
    library(gghalves)
    
    ##热图
    # library(pheatmap)
    library(ComplexHeatmap)
    
    ##饼图
    # install.packages("ggpie")
    library(ggpie)
    
    #画单细胞密度图用到的包
    # library(viridis)
    
    #桑吉图/冲击图 之前画不出来，需要升级R包
    # install.packages("ggalluvial")
    library(ggalluvial)
  }
  

}

## Paths
{
  creat_path <- function(path){
    if (!dir.exists(path)) {
      dir.create(path,recursive=T)
    }
  }
  
  #cloud_path
  cloud_path <- 'E:\\Hao_lab_work\\A_Cloud_File'
  
  global_data_path <- 'global_data'
  
  # function path
  function_path <- file.path("Script","Functions")
  
  # sc 原始数据路径
  sc_path <- file.path("SingleCell_Data")
  sc_rna_path <- file.path(sc_path,"scRNA")
  sc_rna_cellbender_path <- file.path(sc_path,"scRNA_cellbender")
  sc_rna_cellranger_path <- file.path(sc_path,"cellranger_res")
  sc_tcr_path <- file.path(sc_path,"scTCR")
  
  # 结果路径
  sc_res_path <- file.path('SingleCell_Result')
  cellbender_path <- file.path(sc_res_path,'cellbender_res')
  cellranger_path <- file.path(sc_res_path,'cellranger_res')
  
  path_obj <- file.path(cellbender_path,'Seurat_objects_afterN')
  
  final_path <- file.path(cellbender_path,'object_final')
  # fig1_path <- file.path(cellbender_path,'Fig1_main_res')
  # fig2_path <- file.path(cellbender_path,'Fig2_res')
  # fig3_path <- file.path(cellbender_path,'Fig3_res')
  # fig3_tcr_path <- file.path(cellbender_path,'Fig3_TCR_res')
  # fig4_path <- file.path(cellbender_path,'Fig4_res')
  # fig5_mayloid_path <- file.path(cellbender_path,'Fig5_Mayloid_res')
  
  # path_cd8 <- file.path(fig2_path,'CD8')
  # creat_path(path_cd8)
  # 
  # path_cd4 <- file.path(fig2_path,'CD4')
  # creat_path(path_cd4)
  
  afterN_path <- file.path(cellbender_path,'Analysis_afterN')
  creat_path(afterN_path)
  
  fig_T <- file.path(afterN_path,'Figure','Tcell')
  fig_all <- file.path(afterN_path,'Figure','All')
  fig_tcr <- file.path(afterN_path,'Figure','TCR')
  fig_myeloid_b <- file.path(afterN_path,'Figure','Myeloid')
  fig_B <- file.path(afterN_path,'Figure','Fig5 B')
  fig_6 <- file.path(afterN_path,'Figure','Fig6')
  fig_NK <- file.path(afterN_path,'Figure','NK cell')
  # creat_path(fig_NK)
  # creat_path(fig_all)
  # creat_path(fig_tcr)
}

## Functions
fucs <- list.files(function_path,full.names=T)
# 注：有的程序如果读不进来 用notepad改成utf-8编码
lapply(fucs,function(x) {source(x,encoding="utf-8");return("Yes")})

# local fucntion
{
  #
  groupMeans <- function(mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
      if (sparse) {
        Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
      }
      else {
        rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
      }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
  }
  scCluster <- function(obj,nfeature=2500,min.d=0.1,res=1,n.nei=30,TCRBCR=TRUE) {
    if (TCRBCR) {
      TCR.genes <- grep("^TR[AB][VJ]",rownames(obj),value = T)# keep GD TCR genes for GDT cells
      BCR.genes <- c(grep("^IG[KHL][VJC]",rownames(obj),value = T),
                     grep("^AC[0-9]",rownames(obj),value = T))# some RNA genes are also excluded.
    } else {TCR.genes <- c();BCR.genes <- c()}
    obj@assays$RNA@scale.data <- matrix() # in order to subset memory efficiently
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj,nfeatures = nfeature)
    var.genes <- VariableFeatures(obj)
    var.genes <- setdiff(var.genes,c(TCR.genes,BCR.genes))
    VariableFeatures(obj) <- var.genes
    obj <- ScaleData(obj,features = var.genes)
    obj <- RunPCA(obj, verbose = FALSE,features = VariableFeatures(obj),npcs=30)
    obj <- RunHarmony(obj, group.by.vars=c("orig.ident"),assay.use ="RNA")
    obj <- FindNeighbors(obj, dims=1:30,reduction = "harmony",k.param = 30)
    obj <- FindClusters(obj,resolution=res,random.seed=123,graph.name = 'RNA_snn')
    obj <- RunUMAP(obj,reduction = "harmony",seed.use = 123,dims=1:30,n.neighbors = n.nei,
                   umap.method='uwot',min.dist=min.d,spread=1)
    return(obj)
  }
  creat_jump <- function(object,files){
          object.meta <- data.table(cellid=colnames(object),object@meta.data)
          object.meta$umap1 <- object@reductions$umap@cell.embeddings[,1]
          object.meta$umap2 <- object@reductions$umap@cell.embeddings[,2]
          object.meta$seurat_clusters_use <- paste('cluster',object.meta$seurat_clusters)
          fwrite(object.meta,file = files,sep = '\t')
  }
  updata_meta <- function(object){
    
    object.meta <- data.table(cellid=colnames(object),object@meta.data)
    
    object.meta$CellType_2 <- all.metadata[match(object.meta$cellid,cellid),CellType_2]
    object.meta$CellType_n <- all.metadata[match(object.meta$cellid,cellid),CellType_n]
    object.meta$CellType_umap1 <- all.metadata[match(object.meta$cellid,cellid),CellType_umap1]
    object.meta$celltype_main <- all.metadata[match(object.meta$cellid,cellid),celltype_main]
    
    object.meta[grepl('yyh',orig.ident),sample:=str_replace(orig.ident,'yyh','_post')]
    object.meta[grepl('yyq',orig.ident),sample:=str_replace(orig.ident,'yyq','_pre')]
    ##treat
    object.meta[grep('yyq',orig.ident),treat:='Pre']
    object.meta[grep('yyh',orig.ident),treat:='Post']
    #response
    object.meta$response <- object.meta$orig.ident
    object.meta[grep('MA1',orig.ident),response:='NR']
    object.meta[grep('MA2',orig.ident),response:='NR']
    object.meta[grep('MA5',orig.ident),response:='NR']
    object.meta[grep('MA6',orig.ident),response:='NR']
    object.meta[grep('MA3',orig.ident),response:='R']
    object.meta[grep('MA4',orig.ident),response:='NR'] ##Updata
    object.meta[grep('MA7',orig.ident),response:='R']
    #combine label
    object.meta$combine_label <- object.meta$response
    object.meta[treat=='Pre'&response=='NR',combine_label:='pre_NR']
    object.meta[treat=='Pre'&response=='R',combine_label:='pre_R']
    object.meta[treat=='Post'&response=='NR',combine_label:='post_NR']
    object.meta[treat=='Post'&response=='R',combine_label:='post_R']
    object.meta$combine_label <- factor(object.meta$combine_label,levels = c('pre_R','post_R','pre_NR','post_NR'))
    object.meta$treat <- factor(object.meta$treat,levels = c('Pre','Post'))
    object.meta$response <- factor(object.meta$response,levels = c('R','NR'))
    
    ##sample order
    sample_order <- NULL
    for (i in 1:7) {
      tmp <- c(paste0('MA',i,'_pre'),paste0('MA',i,'_post'))
      sample_order <- c(sample_order,tmp)
    }
    object.meta$sample <- factor(object.meta$sample,levels = sample_order)
    
    object.meta$patient <- str_split_fixed(object.meta$sample,'_',2)[,1]
    
    object@meta.data <- dt_2_df(object.meta)
    return(object)
  }
  ColAssign <- function(Var,palettes="Classic 20", n = 20){
    require(ggthemes);require(RColorBrewer)
    pal <- tableau_color_pal(palette = palettes,direction = 1,type="regular")
    if (length(Var) > n) {
      palOut <- colorRampPalette(pal(n))(length(Var))
      names(palOut) <- Var
    } else if (length(Var) == n) {
      palOut <- pal(n)
      names(palOut) <- Var
    } else if (length(Var) < n) {
      palOut <- pal(n)
      palOut <- setdiff(palOut,c("#7f7f7f","#c7c7c7"))# remove grey colors
      #palOut <- sample(palOut)
      palOut <- c(palOut,c("#7f7f7f","#c7c7c7"))
      palOut <- palOut[1:length(Var)]
      names(palOut) <- Var
    }
    return(palOut)
  }
}

## color 样本 细胞类型 factor 顺序颜色固定
{
  sample_order <- c("MA1_pre","MA1_post","MA2_pre","MA2_post","MA4_pre","MA4_post",
                    "MA5_pre","MA5_post","MA6_pre","MA6_post","MA3_pre","MA3_post","MA7_pre","MA7_post")
  sample_order_ori <- c("MA1yyq","MA1yyh","MA2yyq","MA2yyh","MA4yyq","MA4yyh","MA5yyq","MA5yyh",
                        "MA6yyq","MA6yyh","MA3yyq","MA3yyh","MA7yyq","MA7yyh")
  
  color_sample <- c(
    "MA1_pre"="#1f77b4",
    "MA1_post"="#aec7e8",
    "MA2_pre"="#ff7f0e",
    "MA2_post"="#ffbb78",
    "MA3_pre"="#2ca02c",
    "MA3_post"="#98df8a",
    "MA4_pre"="#e377c2",
    "MA4_post"="#f7b6d2",
    "MA5_pre"="#9467bd",
    "MA5_post"="#c5b0d5",
    "MA6_pre"="#8c564b",
    "MA6_post"="#c49c94",
    "MA7_pre"="#d62728",
    "MA7_post"="#ff9896"
  )
  color_sample_new <- c(
    "MA1_pre"="#1f77b4",
    "MA1_post"="#aec7e8",
    "MA2_pre"="#ff7f0e",
    "MA2_post"="#ffbb78",
    "MA4_pre"="#e377c2",
    "MA4_post"="#f7b6d2",
    "MA5_pre"="#9467bd",
    "MA5_post"="#c5b0d5",
    "MA6_pre"="#8c564b",
    "MA6_post"="#c49c94",
    "MA3_pre"="#2ca02c",
    "MA3_post"="#98df8a",
    "MA7_pre"="#d62728",
    "MA7_post"="#ff9896"
  )
  
  color_patient <- c(
    "MA1"="#1f77b4",
    "MA2"="#ff7f0e",
    "MA3"="#2ca02c",
    "MA4"="#e377c2",
    "MA5"="#9467bd",
    "MA6"="#8c564b",
    "MA7"="#d62728"
  )
  
  color_combine <- c('pre_R'='#add995',
                     'post_R'='#00a233',
                     'pre_NR'='#ffa2a4',
                     'post_NR'='#ff0700')
  color_combine_fill <- c('pre_R'='#dbe9cf',
                     'post_R'='#b7cea2',
                     'pre_NR'='#ffd4d2',
                     'post_NR'='#ffa58e')
  
  color_combine <- c('pre_R'='#71afe5',
                     'post_R'='#106ebe',
                     'pre_NR'='#ffa2a4',
                     'post_NR'='#ff0700')
  color_combine_fill <- c('pre_R'='#a3c6e6',
                          'post_R'='#88ABCB',
                          'pre_NR'='#ffd4d2',
                          'post_NR'='#ffa58e')
  
  
  color_treat <- c('Pre'='#a5cde1',
                   'Post'='#2278b3')
  
  color_treat <- c('Pre'='#7AB774',
                   'Post'='#B0BBDE')
  
  color_pdsd_fill <- c('R'='#dbe9cf','NR'='#ffd3d3')
  color_pdsd <- c('R'='#b1d998','NR'='#ffa3a5')
  
  color_pdsd_fill <- c('R'='#a3c6e6','NR'='#ffd3d3')
  color_pdsd <- c('R'='#88ABCB','NR'='#ffa3a5') 
  
  #main
  {
    color_main <-c('T cells'="#1f77b4",
                   'NK cells'="#aec7e8",
                   'Myeloid'='#f68422',
                   # 'Mono/Macro cells'='#f68422',
                   # 'DC'='#007175',
                   # 'Neutrophils'='#8b7042',
                   'B cells'='#d5c8a0',
                   'Plasma cells'='#9d9d82',
                   'cancer cells'='#9668BC')
    
    color_main2 <-c('T cells'="#1f77b4",
                   'NK cells'="#aec7e8",
                   'B cells'='#d5c8a0',
                   'Plasma cells'='#9d9d82',
                   'Myeloid'='#f68422',
                   'cancer cells'='#9668BC')
    
    
    # color_all_celltype1 <- c('CD4_naive'='#c8e6c9',
    #                          'CD4_memory'='#9EBEE2',
    #                          'CD4_EFF'='#E85959',
    #                          'CD4_TFH'='#9E76C1',
    #                          'CD4_Treg'='#FF9332',
    #                          'CD8_naive_T'="#0a7b37",
    #                          'CD8_memory'="#5d9fcf",
    #                          'CD8_CTL'='#cc101c',
    #                          'CD8_ISG_T'='#f9be25',
    #                          'CD8_GD_T'='#C8969A',
    #                          'MAIT'='#FB9697',
    #                          'NKT'='#8D594F',
    #                          'NKT_AP-1'='#AF8DA0',
    #                          'NKT_ISG'='#AA7C73',
    #                          'Proliferative_T'='#1F77B4',
    #                          'NK cells'='#03635f',
    #                          "Monocyte"='#fbbb7c',
    #                          "Macrophage"='#f68422',
    #                          "Neutrophils"='#8b7042',
    #                          "cDC1"='#284852',
    #                          "cDC2"='#007175',
    #                          "pDC"='#3271ae',
    #                          'MDSC'='#b83570',
    #                          "Plasma cells"='#a6cee2',
    #                          "B cells"='#1f78b4',
    #                          "cancer cells"='#9668BC')
    
    
    # color_all_celltype1 <- c('CD4_naive'='#c8e6c9',
    #                          'CD4_memory'='#9EBEE2',
    #                          'CD4_EFF'='#f7dedc',
    #                          'CD4_TFH'='#9E76C1',
    #                          'CD4_Treg'='#FF9332',
    #                          'CD8_naive_T'="#0a7b37",
    #                          'CD8_memory'="#5d9fcf",
    #                          'CD8_CTL'='#cc101c',
    #                          'CD8_ISG_T'='#f9be25',
    #                          'CD8_GD_T'='#C8969A',
    #                          'MAIT'='#FB9697',
    #                          'NKT'='#5976ba',
    #                          'NKT_AP-1'='#6f94cd',
    #                          'NKT_ISG'='#88abda',
    #                          'Proliferative_T'='#547689',
    #                          'NK cells'='#b5c7e1',
    #                          "Monocyte"='#fbbb7c',
    #                          "Macrophage"='#f68422',
    #                          "Neutrophils"='#8b7042',
    #                          "cDC1"='#284852',
    #                          "cDC2"='#007175',
    #                          "pDC"='#3271ae',
    #                          'MDSC'='#b83570',
    #                          "Plasma cells"='#9d9d82',
    #                          "B cells"='#d5c8a0',
    #                          "cancer cells"='#9668BC')
    
    color_all_celltype1 <- c('CD4_Naive'='#c8e6c9',
                             'CD4_Mem'='#9EBEE2',
                             'CD4_GZMK'= '#BCBD22',
                             'CD4_EFF'='#ffa58d',
                             'CD4_Tfh'='#9E76C1',
                             'CD4_Treg'='#FF9332',
                             'CD8_Naive'="#0a7b37",
                             'CD8_Tcm'="#3BBEB2",
                             'CD8_GZMK'="#5D9FCF",
                             'CD8_ISG'='#f9be25',
                             'CD8_CTL'='#cc101c',
                             'GDT'='#C8969A',
                             'MAIT'='#FB9697',
                             'NKT'='#5976ba',
                             # 'NKT_CD16_low'='#98DF8A',
                             # 'NKT_ISG'='#ff6641',
                             'Prolif.T'='#547689',
                             'NK_CD16'="#1F78B4",
                             'NK_CD56'="#98D083",
                             'NK_ISG'="#FD9F35",
                             "Mono_CD14"='#e67f26',
                             "Mono_FCGR3A"='#87b7d8',
                             "Macro_C1QC"='#e078a8',
                             "Macro_NFKB"='#e1ab8c',
                             "Macro_CXCL10"='#f7cc76',
                             "Macro_LGALS3"='#86d4d8',
                             "Macro_MT1H"='#7bed9f',
                             "Neutrophils"='#8b7042',
                             "Neu_ISG"='#fdbada',
                             "cDC1"='#f7ce9c',
                             "cDC2"='#3478b3',
                             'cCD3_LAMP3'='#b83570',
                             "pDC"='#a7cfe7',
                             'Naive B'="#1f77b4",
                             'Memory B'="#aec7e8",
                             'Plasma'="#2ca02c",
                             'Plasmablast'="#98df8a"
                             )
    
    
    color_umap1 <- c('CD4_Naive'='#c8e6c9',
                     'CD4_Mem'='#9EBEE2',
                     'CD4_GZMK'= '#BCBD22',
                     'CD4_EFF'='#ffa58d',
                     'CD4_Tfh'='#9E76C1',
                     'CD4_Treg'='#FF9332',
                     'CD8_Naive'="#0a7b37",
                     'CD8_Tcm'="#3BBEB2",
                     'CD8_GZMK'="#5D9FCF",
                     'CD8_ISG'='#f9be25',
                     'CD8_CTL'='#cc101c',
                     'NKT'='#96CC86',
                     'GDT'='#C8969A',
                     'MAIT'='#FB9697',
                     'Prolif.T'='#547689',
                     'NK'='#B1C1D7',
                     "Monocyte"='#e67f26',
                     "Macrophage"='#e1ab8c',
                     "Neutrophil"='#8b7042',
                     "DC"='#6F9C5A',
                     "pDC"='#a7cfe7',
                     "Plasma"='#9d9d82',
                     "B cells"='#d5c8a0',
                     "cancer"='#9668BC')
    
    
    # color_main_2 <-c(
    #   'CD4 T cells'="#81d4fa",
    #   'CD8 T cells'="#1f77b4",
    #   'Proliferative_T'="#1a237e",
    #   'NKT cells'="#00bcd4",
    #   'NK cells'="#00695c",
    #   'Mono/Macro cells'='#ff8022',
    #   'DC'='#b05828',
    #   'Neutrophils'='#c4a730',
    #   'B cells'='#E04242',
    #   'Plasma cells'='#FC918F',
    #   'cancer cells'='#9668BC')
    
  }
  
  #T
  { 
    color_T <- c('CD4_Naive'='#c8e6c9',
                 'CD4_Mem'='#9EBEE2',
                 'CD4_GZMK'= '#BCBD22',
                 'CD4_EFF'='#ffa58d',
                 'CD4_Tfh'='#9E76C1',
                 'CD4_Treg'='#FF9332',
                 'CD8_Naive'="#0a7b37",
                 'CD8_Tcm'="#3BBEB2",
                 'CD8_GZMK'="#5D9FCF",
                 'CD8_ISG'='#f9be25',
                 'CD8_CTL'='#cc101c',
                 'NKT'='#96CC86',
                 'GDT'='#C8969A',
                 'MAIT'='#FB9697',
                 'Prolif.T'='#547689'
    )
    
    color_cd8_NKT <- c(
      'CD8_Naive'="#0a7b37",
      'CD8_Tcm'="#3BBEB2",
      'CD8_GZMK'="#5D9FCF",
      'CD8_ISG'='#f9be25',
      'CD8_CTL'='#cc101c',
      'NKT'='#96CC86'
    )
    
    color_cd8 <-c('CD8_Naive'="#0a7b37",
                   'CD8_Tcm'="#3BBEB2",
                   'CD8_GZMK'="#5D9FCF",
                   'CD8_ISG'='#f9be25',
                   'CD8_CTL'='#cc101c'
    )
    
    color_cd4 <-c('CD4_Naive'='#c8e6c9',
                   'CD4_Mem'='#9EBEE2',
                   'CD4_GZMK'= '#BCBD22',
                   'CD4_EFF'='#ffa58d',
                   'CD4_Tfh'='#9E76C1',
                   'CD4_Treg'='#FF9332')
    
  }
  

  color_Myeloid <- c("Mono_CD14"='#e67f26',
                     "Mono_FCGR3A"='#87b7d8',
                     "Macro_C1QC"='#e078a8',
                     "Macro_NFKB"='#e1ab8c',
                     "Macro_CXCL10"='#f7cc76',
                     "Macro_LGALS3"='#86d4d8',
                     "Macro_MT1H"='#9467BD',
                     "Neutrophils"='#8b7042',
                     "Neu_ISG"='#FF9896',
                     "cDC1"='#f7ce9c',
                     "cDC2"='#3478b3',
                     'cCD3_LAMP3'='#b83570',
                     "pDC"='#a7cfe7')
  
  color_macro <- c(
                     "Macro_C1QC"='#e078a8',
                     "Macro_NFKB"='#e1ab8c',
                     "Macro_CXCL10"='#f7cc76',
                     "Macro_LGALS3"='#86d4d8',
                     "Macro_MT1H"='#9467BD')

  color_B <- c('Naive B'="#1f77b4",
               'Memory B'="#aec7e8",
               # 'ISG Memory B'= "#ff7f0e",
               # 'Atypical Memory B'="#ffbb78",
               'Plasma'="#2ca02c",
               'Plasmablast'="#98df8a"
  )
  
  color_NK <- c('NK_CD16'="#1F78B4",
               'NK_CD56'="#98D083",
               'NK_ISG'="#FD9F35"
  )

}

## 1.1 data load
{
  # #all meta 2023年9月13日最新注释信息
  # load(file.path(path_obj,'all.metadata.filter.rda'))
  load('SingleCell_Result/cellbender_res/object_final/all.metadata_filter.rda')
  all.metadata
  
  #AUCell score
  load(file.path(afterN_path,'AUCell_score_all.rda'))
  
  # all
  object <- readRDS('SingleCell_Result/cellbender_res/object_final/object_final.rds')#105007 
  object <- updata_meta(object)
  # object$CellType_n <- factor(object$CellType_n,levels = names(color_all_celltype1))
  # object$CellType_umap1 <- factor(object$CellType_umap1,levels = names(color_umap1))
  object.meta <- data.table(cellid=colnames(object),object@meta.data)
  # saveRDS(object,file = 'SingleCell_Result/cellbender_res/object_final/object_final.rds')
  
  #CytoTRACE
  load(file.path(final_path,'CytoTRACE_T.rda'))
  
  # T
  object_T <- readRDS(file.path(final_path,'object_T.rds')) #44026
  object_T <- updata_meta(object_T)
  object_T$CytoTrace <- CytoTRACE_T[match(colnames(object_T),cellid),CytoTrace]
  object_T$TCR <- ifelse(is.na(object_T$barcode),"F",'T')
  object_T.meta <- data.table(cellid=colnames(object_T),object_T@meta.data)
  
  # CD8_NKT
  object_cd8_NKT <- readRDS(file = file.path(final_path,'object_cd8_NKT.rds'))
  object_cd8_NKT <- updata_meta(object_cd8_NKT)
  object_cd8_NKT$CytoTrace <- CytoTRACE_T[match(colnames(object_cd8_NKT),cellid),CytoTrace]
  object_cd8_NKT_meta <- data.table(cellid=colnames(object_cd8_NKT),object_cd8_NKT@meta.data)
  
  
  # cd4
  object_cd4T <- readRDS(file.path(final_path,'object_cd4T.rds')) #11777
  object_cd4T <- updata_meta(object_cd4T)
  object_cd4T$CytoTrace <- CytoTRACE_T[match(colnames(object_cd4T),cellid),CytoTrace]
  object_cd4T$CellType_n <- factor(object_cd4T$CellType_n,levels = c("CD4_Naive","CD4_Mem","CD4_GZMK","CD4_EFF","CD4_Tfh",'CD4_Treg'))
  object_cd4T_meta <- data.table(cellid=colnames(object_cd4T),object_cd4T@meta.data)
  
  # cd8
  object_cd8T <- readRDS(file.path(final_path,'object_cd8T.rds')) #15070
  object_cd8T <- updata_meta(object_cd8T)
  object_cd8T$CytoTrace <- CytoTRACE_T[match(colnames(object_cd8T),cellid),CytoTrace]
  object_cd8T$CellType_n <- factor(object_cd8T$CellType_n,levels = c("CD8_Naive","CD8_Tcm","CD8_GZMK","CD8_ISG","CD8_CTL"))
  object_cd8T_meta <- data.table(cellid=colnames(object_cd8T),object_cd8T@meta.data)
  
  # ## 计算cytoTRACE得分
  # {
  #   object_use <- object_T
  #   GeneCounts <- as.matrix(object_use@assays$RNA@counts)
  #   iOrd <- rowSums(GeneCounts>0)
  #   GeneCounts <- GeneCounts[iOrd>10,]#only keep genes expressing in more than 10 cell
  #   # ncores>1 报错，win系统不支持？
  #   object.CytoTRACE <- CytoTRACE(GeneCounts, enableFast = TRUE, ncores = 1, subsamplesize = 1000)
  #   
  #   CytoTRACE_T <- data.table(cellid=names(object.CytoTRACE$CytoTRACE),CytoTrace=object.CytoTRACE$CytoTRACE)
  #   save(CytoTRACE_T,file = file.path(path_obj,'CytoTRACE_T.rda'))
  # }
  
  # NK
  object_NK <- readRDS(file.path(final_path,'object_NK.rds'))
  object_NK <- updata_meta(object_NK)
  object_NK_meta <- data.table(cellid=colnames(object_NK),object_NK@meta.data)
  
  # B
  object_B <- readRDS(file.path(final_path,'object_B.rds'))
  object_B <- updata_meta(object_B)
  object_B$CellType_n <- factor(object_B$CellType_n,levels = names(color_B))
  object_B_meta <- data.table(cellid=colnames(object_B),object_B@meta.data)
  
  # 
  # Myeloid
  object_Myeloid <- readRDS(file.path(final_path,'object_Myeloid.rds'))
  object_Myeloid <- updata_meta(object_Myeloid)
  object_Myeloid$CellType_n <- factor(object_Myeloid$CellType_n,levels = names(color_Myeloid))
  object_Myeloid_meta <- data.table(cellid=colnames(object_Myeloid),object_Myeloid@meta.data)
  
  object_macro <- subset(object_Myeloid, subset=CellType_umap1=='Macrophage')
  object_macro_meta <- data.table(cellid=colnames(object_macro),object_macro@meta.data)
  # 
  # # cancer
  object_cancer <- readRDS('SingleCell_Result/cellbender_res/object_final/object_cancer.rds')
  object_cancer <- updata_meta(object_cancer)
  object_cancer_meta <- data.table(cellid=colnames(object_cancer),object_cancer@meta.data)
}

## 1.2 TCR load
{
  ##data
  TCR_combined <- readRDS(file = file.path(sc_res_path,'CombineTCR.rds'))
  #添加patient列，按照病人合并TCR数据，一个患者用药前后克隆扩增放一起统计
  for (i in names(TCR_combined)) {
    # i=names(TCR_combined)[1]
    print(i)
    data <- TCR_combined[[`i`]]
    data$patient <- str_sub(i,1,3)
    data$sample <- str_replace(data$sample,'yyh','_post')
    data$sample <- str_replace(data$sample,'yyq','_pre')
    
    TCR_combined[[`i`]] <- data
  }
  
  #添加治疗信息
  for (i in names(TCR_combined)) {
    # i=names(TCR_combined)[1]
    print(i)
    data <- TCR_combined[[`i`]]
    if (grepl('yyh',i)) {
      data$treat <- 'Post'
    }
    if (grepl('yyq',i)) {
      data$treat <- 'Pre'
    }
    TCR_combined[[`i`]] <- data
  }
  
  TCR_data <- do.call(rbind,TCR_combined)%>% data.table()
  # %>% select(CTgene,patient,treat,sample,barcode) 
  table(TCR_data$patient,TCR_data$treat)
  # TCR_data[,sample:=str_replace(sample,'yyh','_post')]
  # TCR_data[,sample:=str_replace(sample,'yyq','_pre')]
  # TCR_data$sample %>% table()
  TCR_data[,clone_id:=paste0(patient,'_',CTgene)]
  
  
  
  #T cell
  
  object_T_tcr   <- combineExpression(TCR_combined, object_T,group.by = 'patient',proportion = F,cloneCall = "gene+nt")
  object_T_tcr_meta <- data.table(cellid=colnames(object_T_tcr),object_T_tcr@meta.data)
  object_T_tcr_meta[!is.na(CTgene),clone_id:=paste0(patient,'_',CTgene)]
  object_T_tcr$clone_id <- object_T_tcr_meta$clone_id
  
  #
  clone_table <- object_T_tcr_meta[!is.na(clone_id),.(Freq.Patient=.N),by=.(clone_id)] %>% arrange(-Freq.Patient)
  object_T_tcr_meta <- left_join(object_T_tcr_meta,clone_table)
  object_T_tcr$Freq.Patient <- object_T_tcr_meta$Freq.Patient
  
  
  
  # object_T_tcr_meta <- data.table(cellid=colnames(object_T_tcr),object_T_tcr@meta.data)
  # object_T_tcr_meta$Clone_size <- 'None'
  # object_T_tcr_meta[Frequency<=3,Clone_size:='1-3']
  # object_T_tcr_meta[Frequency>=4&Frequency<=10,Clone_size:='4-10']
  # object_T_tcr_meta[Frequency>10&Frequency<=50,Clone_size:='11-50']
  # object_T_tcr_meta[Frequency>50&Frequency<=100,Clone_size:='51-100']
  # object_T_tcr_meta[Frequency>100,Clone_size:='>100']
  # object_T_tcr$Clone_size <- object_T_tcr_meta$Clone_size
  
  ##重新计算TCR 频率,不需要了，CD8，CD4 单独合并TCR数据就行了
  # {
  #   clone_table <- object_T_tcr_meta[,.(TCR_Frequency_new=.N),by=.(patient,CTgene)] %>% select(CTgene,TCR_Frequency_new)
  #   clone_table[is.na(CTgene),TCR_Frequency_new:=0]
  #   object_T_tcr_meta <- left_join(object_T_tcr_meta,clone_table,multiple = "all")
  #   object_T_tcr$TCR_Frequency_new <- object_T_tcr_meta$TCR_Frequency_new
  # }
  
  
  #CD8 Tcell
  # object_cd8T_tcr <- object_cd8T
  object_cd8_tcr_meta <- object_T_tcr_meta[CellType_2%in%c('CD8+ cells','NKT')&!is.na(clone_id),]
  
  #CD4 Tcell
  # object_cd4T_tcr <- object_cd4T
  object_cd4_tcr_meta <- object_T_tcr_meta[CellType_2=='CD4+ cells'&!is.na(clone_id),]

}

## 1.3 功能marker基因
{
  #樊嘉 2021 Cell 肝癌 marker
  sig_fj_hcc_Proliferation <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8","HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA")
  sig_fj_hcc_Hepatic <- c("ABCB11","ACOT12","ACOX2","ADH1A","AFP","ALB","ALDH1B1","ALDH1L1","ALDH4A1","ALDH6A1","ALDH7A1","APOA1","APOA2","APOB","APOC1","APOC2","APOC4","APOE","APOF","APOH","APOM","AQP8","ARG1","AVPR1A","B4GALNT1","BBOX1","CES1","CPB2","CPM","CPN1","CPS1","CPT2","CYP39A1","DPYS","ENPP2","ENPP3","ENPP5","EPHX1","FMO1","GPLD1","KMO","MAT1A","MIOX","MOCOS","NUDT18","PAPSS2","PEMT","PGLYRP2","PON1","PPARA","PPP1R3B","RARRES2","SERPINA1","SERPINA10","SERPINA6","SERPINC1","SERPIND1","SERPINF1","SERPINF2","SERPING1","SLC10A1","SLC19A2","SLC22A18","SLC25A1","SLC25A10","SLC25A13","SLC25A15","SLC25A39","SLC25A4","SLC27A2","SLC2A2","SLC31A1","SLC35F5","SLC36A2","SLC37A4","SLC38A4","SLC39A4","SLC39A5","SLC6A13","SLC9A3R1","STARD5","SULT1A1","THRB","TMEM30B","TTR")
  sig_fj_hcc_Immune_surveillance <- c("HLA-A","HLA-B","HLA-C","MICA","MICB")
  sig_fj_hcc_Immune_escape <- c("CD47","ADAM10","HLA-G","CD274","FASLG","CCL5","TGFB1","IL10","PTGER4")
  
  #STING
  sig_STING <- fread(file = file.path(fig_T,'STING.csv'))
  sig_STING <- sig_STING$SYMBOL
  
  #ISG score 38gene nature_medician
  sig_NM_ISG <- c("ADAR","DDX60","HERC6","IRF7","OASL","PSME2","STAT2",
                  "TRIM25","BST2","DHX58","IFI35","ISG15","OGFR","RSAD2",
                  "TDRD7","UBE2L6","CASP1","EIF2AK2","IFIH1","ISG20",
                  "PARP12","RTP4","TRAFD1","USP18","CMPK2","EPSTI1","IFIT2",
                  "MX1","PARP14","SAMD9L","TRIM14","CXCL10","GBP4","IFIT3","NMI","PNPT1","SP110","TRIM21")
  
  #上面的ISG相关基因有些不属于ISG，筛选
  sig_ISG_select <- c('IFIT1','IFIT2','IFIT3','IFIT5','IFI35','ISG15',
                      'IRF1','IRF3','IRF7','IRF9','OASL','ADAR','EIF2AK2','EPSTI1','GBP4','MX1','OGFR','SP110',
                      'STAT1','STAT2')
  
  #增值marker  张泽民 pancaner T
  sig_proliferation <- c("ZWINT","E2F1","FEN1","FOXM1","H2AFZ","HMGB2",
                         "MCM2","MCM3","MCM4","MCM5","MCM6","MKI67",   
                         "MYBL2","PCNA","PLK1","CCND1","AURKA","BUB1",
                         "TOP2A","TYMS","DEK","CCNB1","CCNE1")
  
  ##NC胃癌腹水 Supplementary Data 1 里的marker
  {
    ##T cell related
    proliferative_gene <- c("MCM5","MCM2","MCM4","MCM6","MKI67")
    naive_gene <- c("CCR7","LEF1","SELL","TCF7")
    cytotoxic_gene <- c("GNLY","IFNG","NKG7","PRF1","GZMA","GZMB","GZMH","GZMK",
                        "IL2","IL17A","IL17F","LAMTOR3","KLRK1","KLRD1","CTSW",
                        "CST7","CX3CR1","FGFBP2","S1PR5","FCGR3A","PLAC8","C1orf21",
                        "TGFBR3","PLEK","FGR","KLRF1","SPON2","CD300A","S1PR1","RIPOR2",
                        "EFHD2","STK38","C1orf162","SORL1","EMP3","ARL4C","BIN2",
                        "CCND3","FCRL6","SAMD3","TRDC","TYROBP","KLRG1","KLRB1")
    Treg_gene <- c("FOXP3","IL2RA","IKZF2","ENTPD1","IL10RA","RTKN2","CDC25B","NT5E",
                   "CD5","CTLA4","IZUMO1R","TNFRSF18","ITGAE","LAG3","TGFB1","LRRC32",
                   "TNFRSF4","SELL","STAT5A","STAT5B","LGALS1","IL10","IL12A","EBI3")
    inhibitory_gene <- c("HAVCR2","LAG3","PDCD1","CTLA4","TIGIT","BTLA","LAYN","EOMES",
                         "TOX","NFATC1","SLA","ETV1","CD200","ARNT","ENTPD1","TOX2",
                         "PAG1","KIR2DL4","CD200R","CD276","CD96","CD47","VSIR","VSIG4",
                         "LGALS9","CD274","PDCD1LG2","NECTIN2","IDO1","HOPX","ARNT")
    
    ## myeloid related
    Pro_inflammatory_gene <- c("IL1B","CCL2","CCL3","CCL4","TNF","CXCL2","CXCL9",
                               "CXCL10","SOCS3","CCL5","CCL7","CCL8","CCL17","CCL22")
    Anti_inflammatory_gene <- c("CD163","MRC1","MSR1","CCL13","CCL18","IL10","FOLR2",
                                "STAB1","CD163L1","SELENOP","F13A1","MERTK","AXL","CCL22",
                                "IL1RN","GPNMB","CHI3L1","IL4","IL11","IL13","TGFB1",
                                "TNFRSF1A","TNFRSF1B","IL1R2","IL18BP")
    #血管生成
    Pro_angiogenic_gene <- c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2",
                             "FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2",
                             "MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6",
                             "TYMP","VAV2","VCAN","VEGFA")
    Antigen_presenting_gene <- c("CD74","TAP1","TAP2","TAPBP","CD36","HLA-DRA",
                                 "HLA-DRB1","HLA-DPB1","HLA-DOA","HLA-DOB","HLA-DRB5",
                                 "HLA-DRB6","HLA-A","HLA-B","HLA-C")
    M1_marker <- c("IL23","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A",
                   "IL1B","IL6","CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7")
    M2_marker <- c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24",
                   "LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC",
                   "CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B",
                   "FASL","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4")
    
    ## Function
    JAK_STAT_gene <- c("JAK2","JAK3","LYN","PTPN1","STAT3","STAT5A","STAT5B","STAT6")
    TGFB_KEGG_gene <- c("CDKN2B","FST","LEFTY1","ACVR1C","COMP","CREBBP","GDF7","DCN",
                        "E2F4","E2F5","EP300","AMH","AMHR2","ID1","ID2","ID3","ID4",
                        "IFNG","BMP8A","INHBA","INHBB","INHBC","RHOA","GDF6","LTBP1",
                        "SMAD1","SMAD2","SMAD3","SMAD4","SMAD5","SMAD6","SMAD7","SMAD9",
                        "MYC","NODAL","PITX2","PPP2CA","PPP2CB","PPP2R1A","PPP2R1B",
                        "MAPK1","MAPK3","SMURF1","RBL1","RBL2","ROCK1","RPS6KB1","RPS6KB2",
                        "SMURF2","BMP2","SKP1","BMP4","BMP5","BMP6","BMP7","BMP8B",
                        "BMPR1A","BMPR1B","BMPR2","SP1","TFDP1","TGFB1","TGFB2","TGFB3",
                        "LEFTY2","TGFBR1","TGFBR2","THBS1","THBS2","THBS3","THBS4","TNF",
                        "SKP1P2","GDF5","INHBE","CUL1","CHRD","ACVR1","ACVR2A","NOG",
                        "ACVR2B","ZFYVE9","ACVRL1","ROCK2","ZFYVE16","RBX1")
    
    #TNF, 从msigdb的hallmark里获取
    # gse_H <- msigdbr(species = "Homo sapiens",category = "H")
    # save(gse_H,file = file.path(local_path,'msigDB_Hallmark.rda'))
    # load(file.path(local_path,'msigDB_Hallmark.rda'))
    # use_tmp <- gse_H %>% dplyr::select(gs_name,gene_symbol) %>% data.table()
    # TNF_gene <- use_tmp[grep('TNF',gs_name),gene_symbol]
    # 
    # 
    # apoptosis_gene <- use_tmp[grep('apoptosis',gs_name,ignore.case = T),gene_symbol]
    
    
  }
  
  ##Cancer Cell TableS2
  {
    CD8_exhaustion <- c("LAG3","HAVCR2","PDCD1","PTMS","FAM3C","IFNG","AKAP5","CD7","PHLDA1",
                        "ENTPD1","SNAP47","TNS3","CXCL13","RDH10","DGKH","KIR2DL4","LYST",
                        "MIR155HG","RAB27A","CSF1","CTLA4","TNFRSF9","CD27","CCL3","ITGAE",
                        "PAG1","TNFRSF1B","GALNT1","GBP2","MYO7A")
    Progenitor_score <- c("XCL1","MS4A4C","ID3","CXCL10","SLAMF6","TCF7","IL7R","TNFSF8",
                          "CTLA2A","LTB","CD9","SOCS3","JUN","GPR183","TRAF1","ITGB1","GM8730","EMB","GM10073")
  }
  
  ##老师结直肠癌工作功能marker
  {
    ##老师结直肠癌工作
    # activation.genes <- c("GZMA","GZMB","GZMH","GNLY","NKG7","IFNG","PRF1")
    # exhaustion.genes <- c("CTLA4","PDCD1","TIGIT","HAVCR2","TNFRSF9","ENTPD1",
    #                       "CD27","CXCL13","LAYN","LAG3")
    # 
    # Exh.markers <- c('LAG3','HAVCR2','PDCD1','PTMS','FAM3C','IFNG','AKAP5','CD7','PHLDA1','ENTPD1',
    #                  'SNAP47','TNS3','CXCL13','RDH10','DGKH','KIR2DL4','LYST','MIR155HG','RAB27A',
    #                  'CSF1','CTLA4','TNFRSF9','CD27','CCL3','ITGAE','PAG1','TNFRSF1B','GALNT1','GBP2','MYO7A')
    
    NK.markers <- c('PIK3R1','PIK3CA','LAT','LAT','LAT','MAPK3','MAPK3','KLRD1',
                    'PAK1','SYK','SYK','SYK','PIK3R1','HLA-A','IL18','VAV1','VAV1',
                    'IL18','HLA-A','ITGB1','KLRC1','KLRC2','KLRC3','KLRD1','PAK1',
                    'MAPK3','MAP2K1','PTPN6','SYK','B2M','PTK2B','VAV1','PIK3CA',
                    'RAC1','KLRC1','KLRC3','KLRD1','KLRC4','LAT','RAC1','ITGB1','PTPN6',
                    'PTPN6','ITGB1','PTK2B','PTK2B','PTK2B','PIK3R1','PIK3R1','PIK3R1','KLRC1','KLRC1')
    Exh.markers <- c('LAG3','HAVCR2','PDCD1','PTMS','FAM3C','IFNG','AKAP5','CD7','PHLDA1','ENTPD1',
                     'SNAP47','TNS3','CXCL13','RDH10','DGKH','KIR2DL4','LYST','MIR155HG','RAB27A',
                     'CSF1','CTLA4','TNFRSF9','CD27','CCL3','ITGAE','PAG1','TNFRSF1B','GALNT1','GBP2','MYO7A')
    Prog.markers <- c('XCL1','MS4A4C','ID3','CXCL10','SLAMF6','TCF7','IL7R','TNFSF8','CTLA2A',
                      'LTB','CD9','SOCS3','JUN','GPR183','TRAF1','ITGB1','GM8730','EMB','GM10073')
  }
  
  ##初砚硕的marker
  {
    cys_markser <- fread(file.path(afterN_path,'初砚硕_cd8_signature_score.txt'))
    Cytotoxicity_cys <- cys_markser$Cytotoxicity[cys_markser$Cytotoxicity!='']
    Activation_cys <- cys_markser$`Activation:Effector function`[cys_markser$`Activation:Effector function`!='']
    Exhaustion_cys <- cys_markser$Exhaustion[cys_markser$Exhaustion!='']
    Pro_apoptosis_cys <- cys_markser$`Pro-apoptosis`[cys_markser$`Pro-apoptosis`!='']
    
    cys_markser_cd4 <- fread(file.path(afterN_path,'初砚硕_cd4_signature_score.txt'),header = T)
    Treg_gene_cys <- cys_markser_cd4$`Treg signature`[cys_markser_cd4$`Treg signature`!='']
  }
  
  ##张泽民 pancancer 髓系
  {
    
    zhang_table <- read.xlsx(file.path(cloud_path,'文献/张泽民单细胞/zhang2021_pancancer_髓样细胞','1-s2.0-S0092867421000106-mmc5.xlsx'),startRow = 2)
    sig_M1_zhang <- zhang_table$M1 %>% na.omit() %>% as.character()
    sig_M1_zhang <- gsub("\\s+$", "", sig_M1_zhang)
    sig_M1_zhang <- gsub("IL23", "IL23A", sig_M1_zhang) #IL23A
    
    sig_M2_zhang <- zhang_table$M2 %>% na.omit() %>% as.character()
    sig_M2_zhang <- gsub("\\s+$", "", sig_M2_zhang)
    sig_M2_zhang <- gsub("FASL", "FASLG", sig_M2_zhang) #IL23A
    
    sig_Angiogenesis_zhang <- zhang_table$Angiogenesis %>% na.omit() %>% as.character()
    sig_Angiogenesis_zhang <- gsub("\\s+$", "", sig_Angiogenesis_zhang)
    sig_Phagocytosis_zhang <- zhang_table$Phagocytosis %>% na.omit() %>% as.character()
    sig_Phagocytosis_zhang <- gsub("\\s+$", "", sig_Phagocytosis_zhang)
  }
  
  ##science TCR分析 marker表格
  {
    science_table <- fread('global_data/1-s2.0-S153561082300082X-mmc3.xlsx')
    science_table <- data.table(read.xlsx('global_data/1-s2.0-S153561082300082X-mmc3.xlsx','Transcriptional signatures',startRow = 2))
    science_table <- science_table[-1,]
    
    sig_science_NeoTCR_cd8 <- science_table[,`NeoTCR-CD8.(243.genes)`]
    sig_science_NeoTCR_cd4 <- science_table[!is.na(`NeoTCR-CD4.(40.genes)`),`NeoTCR-CD4.(40.genes)`]
    sig_science_CD8_exhaustion <- science_table[!is.na(`CD8+.exhaustion`),`CD8+.exhaustion`]
    sig_science_Virus_Specific <- science_table[!is.na(`Virus.Specific`),`Virus.Specific`]
    sig_science_Tumor_Specific <- science_table[!is.na(`Tumor.Specific`),`Tumor.Specific`]
    sig_science_CD8_Tumor_reactivity <- science_table[!is.na(`CD8.Tumor-reactivity`),`CD8.Tumor-reactivity`]
    sig_science_Influenza_TIL <- science_table[!is.na(`Influenza.TIL`),`Influenza.TIL`]
    sig_science_Progenitor_score <- science_table[!is.na(`Progenitor.score`),`Progenitor.score`]
    sig_science_MANA_TIL <- science_table[!is.na(`MANA-TIL`),`MANA-TIL`]
    
  }
  
  ##张泽民 pancancer NK文章里的marker
  {
    sig_zhangNK_cytotoxicity <- c('GZMA', 'GZMB', 'GZMH', 'GZMM', 'GZMK', 'GNLY', 'PRF1', 'CTSW')
    sig_zhangNK_inflammatory <- c('CCL2','CCL3','CCL4','CCL5','CXCL10','CXCL9','IL1B','IL6','IL7','IL15','IL18')
    sig_zhangNK_stress <- c("BAG3","CALU","DNAJB1","DUSP1","EGR1","FOS","FOSB","HIF1A","HSP90AA1","HSP90AB1","HSP90B1",
                            "HSPA1A","HSPA1B","HSPA6","HSPB1","HSPH1","IER2","JUN,","JUNB","NFKBIA","NFKBIZ","RGS2",
                            "SLC2A3","SOCS3","UBC","ZFAND2A","ZFP36","ZFP36L1")
    
  }

}

## 1.4 全局数据（KEGG GO MisgDB）
{
  ## KEGG 通路基因集合
  {
    load(file.path(global_data_path,'KEGGSets.rds'))
    load(file.path(global_data_path,'msigDB_Hallmark.rda'))
    # KEGGSets
    # msigdb_set 
    # names(KEGGSets)
    # names(msigdb_set)
    
    ##免疫相关基因
    #hallmark
    # Hallmark_immune <- unlist(msigdb_set[c(  "Hallmark_interferon_gamma_response",
    #                                          "Hallmark_interferon_alpha_response",
    #                                          "Hallmark_tnfa_signaling_via_nfkb",
    #                                          "Hallmark_il2_stat5_signaling")]) %>%  unique()
    # KEGG_immune <- unlist(KEGGSets[c(  "T cell receptor signaling pathway - Homo sapiens (human)",
    #                                    "Th1 and Th2 cell differentiation - Homo sapiens (human)",  
    #                                    "Th17 cell differentiation - Homo sapiens (human)",
    #                                    "TNF signaling pathway - Homo sapiens (human)",
    #                                    "Toll-like receptor signaling pathway - Homo sapiens (human)")] ) %>%  unique()
    # pathway_immune <- unique(c(Hallmark_immune,KEGG_immune))
    
    # intersect(sig_TCRsignal,KEGGSets[["T cell receptor signaling pathway - Homo sapiens (human)"]])
    sig_KEGG_TCR <- KEGGSets[["T cell receptor signaling pathway - Homo sapiens (human)"]]
    # names(KEGGSets)
    # grep('mTOR',names(KEGGSets),value = T)
    sig_KEGG_JAK <- KEGGSets[["JAK-STAT signaling pathway - Homo sapiens (human)"]]
    sig_KEGG_PI3K <- KEGGSets[["PI3K-Akt signaling pathway - Homo sapiens (human)"]]
    sig_KEGG_mTOR <- KEGGSets[["mTOR signaling pathway - Homo sapiens (human)"]]
    sig_KEGG_TNF <- KEGGSets[["TNF signaling pathway - Homo sapiens (human)"]]
    sig_KEGG_NFKB <- KEGGSets[["NF-kappa B signaling pathway - Homo sapiens (human)"]]
    sig_KEGG_Toll <- KEGGSets[["Toll-like receptor signaling pathway - Homo sapiens (human)"]]
    sig_KEGG_BCR <- KEGGSets[['B cell receptor signaling pathway - Homo sapiens (human)']]
    sig_KEGG_Antigen_processing <- KEGGSets[["Antigen processing and presentation - Homo sapiens (human)"]]
    
  }
  
  ## MsigDB canonical_pathway
  {
    #集合数量与网站上还是有一些差别
    # [1] "CP:BIOCARTA : 292"
    # [1] "CP:KEGG : 186"
    # [1] "CP:PID : 196"
    # [1] "CP:REACTOME : 1615"
    # [1] "CP:WIKIPATHWAYS : 664"
    
    # library(msigdbr)
    # msigdb_cp_path <- data.table()
    # for (i in c('CP:BIOCARTA','CP:KEGG','CP:PID','CP:REACTOME','CP:WIKIPATHWAYS')) {
    #   tmp <- msigdbr(species = "Homo sapiens", category = 'C2',subcategory = i) 
    #   print(paste0(i,' : ',tmp$gs_name %>% unique() %>% length()))
    #   tmp %>% select(gs_name,gene_symbol)
    #   msigdb_cp_path <- rbind(msigdb_cp_path,tmp %>% select(gs_name,gene_symbol))
    # }
    # print(paste0('总计 : ',msigdb_cp_path$gs_name %>% unique() %>% length()))
    # msigdb_cp_path <- unique(msigdb_cp_path)
    # save(msigdb_cp_path,file = file.path(global_data_path,'MsigDB_canonical_pathway.rda'))
    load(file.path(global_data_path,'MsigDB_canonical_pathway.rda'))
    #msigdb_cp_path
  }
  
  ## MsigDB Go
  {
    # msigdb_go_bp <- msigdbr(species = "Homo sapiens", category = 'C5',subcategory = 'GO:BP')
    # msigdb_go_bp <- msigdb_go_bp %>% select(gs_name,gene_symbol) %>% data.table()
    # save(msigdb_go_bp,file = file.path(global_data_path,'msigdb_go_bp.rda'))
    
    load(file.path(global_data_path,'msigdb_go_bp.rda'))
    #msigdb_go_bp 通路基因数目跟网站不对应
    # table(msigdb_go_bp[gs_name=='GOBP_B_CELL_ACTIVATION',])
  }
  
  {
    load('global_data/msigDB_Hallmark.rda')
    # names(msigdb_set)
    # grep('interferon',names(msigdb_set),value = T)
    
    sig_msigdb_IFNA <- msigdb_set[['Hallmark_interferon_alpha_response']]
    sig_msigdb_IFNG <- msigdb_set[['Hallmark_interferon_gamma_response']]
    
  }
  ## immport数据库
  {
    immport_table <- fread('global_data/GeneList_immPort.txt')
    # immport_table$Category %>% unique()
    
    sig_TCRsignal <- immport_table[Category=='TCRsignalingPathway',Symbol]
    sig_Cytotoxicity <- immport_table[Category=='NaturalKiller_Cell_Cytotoxicity',Symbol]
    
    
  }
}



