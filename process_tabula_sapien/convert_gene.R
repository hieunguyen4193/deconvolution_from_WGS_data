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

library(org.Hs.eg.db)
if (packageVersion("clusterProfiler") != "4.13.3"){
  remove.packages("DOSE")
  remove.packages("GOSemSim")
  remove.packages("yulab.utils")
  remove.packages("clusterProfiler")
  remotes::install_github("GuangchuangYu/GOSemSim", upgrade = "never")
  remotes::install_github("GuangchuangYu/clusterProfiler", upgrade = "never")
}

library("clusterProfiler")
if (file.exists(file.path(path.to.01.output, "convertdf.csv")) == FALSE){
  s.obj <- readRDS(path.to.main.input)
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(use.celltype = sprintf("%s_%s", cell_type, tissue_in_publication))
  convertdf <- bitr(row.names(s.obj), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  library("biomaRt")
  mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",
                              mart    = useMart("ENSEMBL_MART_ENSEMBL",
                                                host    = "https://grch37.ensembl.org"))
  resultTable <- biomaRt::getBM(attributes = c("start_position",
                                               "end_position",
                                               "chromosome_name",
                                               "hgnc_symbol"),       
                                values     = convertdf$SYMBOL,         
                                mart       = mart)     
  colnames(resultTable) <- c("gene_start", "gene_end", "chrom", "symbol")
  resultTable <- resultTable[, c("symbol", "chrom", "gene_start", "gene_end")]
  resultTable <- subset(resultTable, resultTable$chrom %in% seq(1,22)) %>%
    subset(symbol != "")
  
  convertdf <- merge(resultTable, convertdf[, c("ENSEMBL", "SYMBOL")], by.x = "symbol", by.y = "SYMBOL")
  convertdf <- convertdf %>% arrange(desc(chrom))
  convertdf <- convertdf[!duplicated(convertdf$symbol),]
  write.csv(convertdf, file.path(path.to.01.output, "convertdf.csv"))
} else {
  convertdf <- read.csv(file.path(path.to.01.output, "convertdf.csv")) %>%
    subset(select = -c(X))
}


