gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline_SeuratV5"
outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "deconvolution_from_WGS_data"
source(file.path(path.to.main.src, "processes_src", "import_libraries.R"))
source(file.path(path.to.main.src, "processes_src", "helper_functions.R"))

path.to.main.input <- "/media/hieunguyen/GSHD_HN01/raw_data/Tabula_Sapiens/rds/all_cells.rds"
path.to.main.output <- file.path(outdir, PROJECT)
path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.01.output, "wilcox_test"), showWarnings = FALSE, recursive = TRUE)

if ("presto" %in% installed.packages()){
  devtools::install_github("immunogenomics/presto")  
}
library("presto")

if (file.exists(file.path(path.to.01.output, "mat.rds")) == FALSE){
  s.obj <- readRDS(path.to.main.input)
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(use.celltype = sprintf("%s_%s", cell_type, tissue_in_publication))
  
  sample.list <- list()
  for (celltype in unique(meta.data$use.celltype)){
    sample.list[[celltype]] <- subset(meta.data, meta.data$use.celltype == celltype)$barcode
  }
  
  mat <- GetAssayData(object = s.obj, assay = "RNA", slot = "data")
  
  saveRDS(mat, file.path(path.to.01.output, "mat.rds"))  
  saveRDS(meta.data, file.path(path.to.01.output, "metadata.rds"))
  saveRDS(sample.list, file.path(path.to.01.output, "sample_list.rds"))
} else {
  mat <- readRDS(file.path(path.to.01.output, "mat.rds"))
  meta.data <- readRDS(file.path(path.to.01.output, "metadata.rds"))
  sample.list <- readRDS(file.path(path.to.01.output, "sample_list.rds"))
}

test.res <- list()
for (select.celltype in unique(meta.data$use.celltype)){
  if (file.exists(file.path(path.to.01.output, "wilcox_test", sprintf("%s.test.rds", select.celltype))) == FALSE){
    print(sprintf("Working on %s", select.celltype))
    labels <- unlist(lapply(colnames(mat), function(x){
      if (x %in% sample.list[[select.celltype]]){
        return(1)
      } else {
        return(0)
      }
    }))
    test.res <- wilcoxauc(mat, labels)
    saveRDS(test.res, file.path(path.to.01.output, "wilcox_test", sprintf("%s.test.rds", select.celltype)))    
  } else {
    print(sprintf("Reading in data from cell type %s", select.celltype))
    test.res[[select.celltype]] <- readRDS(file.path(path.to.01.output, "wilcox_test", sprintf("%s.test.rds", select.celltype)))
  }
}

top.genes <- data.frame()
for (celltype in names(test.res)){
  tmp <- test.res[[celltype]]
  tmp.sig <- subset(tmp, tmp$padj <= 0.05 & abs(tmp$logFC) >= 1) %>%
    rowwise() %>%
    mutate(abs.logFC = abs(logFC)) %>%
    mutate(group = ifelse(logFC > 0, "up", "down"))
  tmp.sig$celltype <- celltype
  top100.up <- subset(tmp.sig, tmp.sig$group == "up") %>% arrange(desc(abs.logFC)) %>% head(100)
  top100.down <- subset(tmp.sig, tmp.sig$group == "down") %>% arrange(desc(abs.logFC)) %>% head(100)
  top.genes <- rbind(top.genes, rbind(top100.up, top100.down))
}

unique.genedf <- data.frame(ENSEMBL = unique(top.genes$feature))
read.csv(file.path(path.to.01.output, "convertdf.csv")) %>%
  subset(select = -c(X))
unique.genedf <- merge(unique.genedf, convertdf, by.x = "ENSEMBL", by.y = "ENSEMBL")