cat("\nAffymetrix pipeline selected\n")


### Import Raw Data

# for (file in opt$files){
#   if (grepl("\\.gz$", file)) {
#     R.utils::gunzip(filename = file, remove = TRUE)
#     opt$files<-list.files(file.path(tempin,"00-RawData"))}
# }
# rm(file)
cat("Extracted runsheet files: ",opt$files)
str(opt$files)

workdir <- opt$out



raw <- oligo::read.celfiles(opt$files)
cat("\nAffymetrix platform subtype: ",class(raw),"\n")

# try({raw <- oligo::read.celfiles(opt$files, pkgname=database)}) #try forcing probe database from options on raw import
# if(!exists(raw)){
#   raw <- oligo::read.celfiles(opt$files)
# }

### Copy Raw Files to output directory

dir.create(file.path(workdir,"Processed_Data"), showWarnings = FALSE)
dir.create(file.path(workdir,"Processed_Data",opt$glds), showWarnings = FALSE)
dir.create(file.path(workdir,"Processed_Data",opt$glds,"00-RawData"), showWarnings = FALSE)
file.copy(from = opt$files, to = file.path(workdir,"Processed_Data",opt$glds,"00-RawData"), overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)

### Create Checksum file
checksums <- tools::md5sum(opt$files)
names(checksums) <- basename(opt$files)
write.table(checksums, file.path(workdir,"Processed_Data",opt$glds,"00-RawData","md5sum.txt"),quote = FALSE)



### Generate Raw Data QA HTML Report
if(opt$reports == TRUE){
  rmarkdown::render(file.path(codebase_dir,"qa_summary_raw.rmd"),"html_document", output_file="raw_qa",output_dir=file.path(workdir,"Processed_Data",opt$glds,"00-RawData"))
}


### Background Correction and Normalization
if (class(raw)=="ExonFeatureSet" || class(raw)=="GeneFeatureSet"){
  data <- oligo::rma(raw, target = "core", background=TRUE, normalize=TRUE)
  data.bgonly <- oligo::rma(raw, target = "core", background=TRUE, normalize=FALSE)
  cat("RMA background correction and quantile normalization performed with gene level summarization.\n")
}

if (class(raw)=="ExpressionFeatureSet"){
  data <- oligo::rma(raw, normalize = TRUE, background = TRUE)
  data.bgonly <- oligo::rma(raw, normalize = FALSE, background = TRUE)
  cat("RMA background correction and quantile normalization performed.\n")
}

### Generate Normalized Data QA HTML Report

if(opt$reports == TRUE){
  rmarkdown::render(file.path(codebase_dir,"qa_summary_normalized.rmd"),"html_document", output_file="normalized_qa",output_dir=file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))
}


###  Write out the expression values
dir.create(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"), showWarnings = FALSE)
setwd(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))
expression <- data.frame(Biobase::exprs(data))
write.table(expression,"normalized.txt",quote=FALSE, append=FALSE, sep = "\t", col.names=NA)

### Import Probe Data
if (length(opt$probe >= 1)){
  options(connectionObserver = NULL)
  database <- sub('\\.annotation.tar.gz$', '', basename(opt$probe)) 
  cat("\nLoading local probe annotation database: ",database,"\n")
  if(!require(database, character.only=TRUE)) {
    BiocManager::install(database, ask = FALSE)
  }
  install.packages(opt$probe,repos = NULL, verbose = FALSE, quiet = TRUE)
  library(database, character.only=TRUE)
  
  
}else {
  package <- raw@annotation
  package <- gsub("pd.","",package)
  package <- gsub(".v1","transcriptcluster",package)
  package <- gsub("[.]","",package)
  package <- paste0(package,".db")
  database <- package
  cat("\nSearch for package: ",database)
  if(!require(database, character.only=TRUE)) {
    BiocManager::install(database, ask = FALSE)
  }
  library(database, character.only=TRUE)
}
keytype<-"PROBEID"
keys<-rownames(expression)


### Map assay database annotations
annotation <- data.frame(REFSEQ=mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "REFSEQ",multiVals = "first"))

try(annotation$ENSEMBL<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "ENSEMBL",multiVals = "first"))
try(annotation$SYMBOL<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "SYMBOL",multiVals = "first"))
try(annotation$GENENAME<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "GENENAME",multiVals = "first"))
try(annotation$ENTREZID<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "ENTREZID",multiVals = "first"))
try(annotation$TAIR<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "TAIR",multiVals = "first"))
try(annotation$GOSLIM_IDS<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "GO",multiVals = "first"))


### Map STRING annotations
try({
  string_db <- STRINGdb::STRINGdb$new( version="11", species=organism_table$taxon[organism_table$species == opt$species],score_threshold=0)
  string_map<-string_db$map(annotation,"ENTREZID",removeUnmappedRows = TRUE, takeFirst = TRUE)
  string_cols <-string_map[,c("ENTREZID","STRING_id")]
  string_cols <- string_cols[!duplicated(string_cols$ENTREZID),]
  annotation <- dplyr::left_join(annotation,string_cols,by="ENTREZID")

  rm(string_map,string_db)
})
rm(keytype,keys)

### Generate normalized annotated expression text file
cat("\nGenerating normalized-annotated.txt file\n")
setwd(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))
expression <- cbind(annotation,expression)
write.table(expression,"normalized-annotated.txt",quote=FALSE, append = FALSE, row.names = FALSE, sep = "\t")
write.table(annotation,"probe_annotations.txt",quote=FALSE, append = FALSE, row.names = FALSE, sep = "\t")
### Rename assay samples with ISAtab sample names
index<-sapply(targets$t1$SampleName,function(x){grep(x,sampleNames(data))})
sampleNames(data)<-targets$t1$SampleName[order(index)]
### Sort assay data to be in same order as ISAtab metadata
targets$t1 <- targets$t1[order(index),]
targets$t2 <- targets$t2[order(index),]
targets$t3 <- targets$t3[order(index),]
rm(index)

### Annotate the expression set object and save as a file
cat("\nGenerating normalized-annotated.rda file\n")
setwd(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))
fData(data)<-annotation
save(data,file = "normalized-annotated.rda")
fData(data.bgonly)<-annotation
#save(data.bgonly,file = "uncorrected-annotated.rda")

### Gene level estimation by maximum interquartile range
cat("\nPerforming gene level estimation by max interquartile range")
data.filt <- data
#data.filt$genes <- annotation

# annotation_stats$gene_level_features <-dim(data.filt)[1]
# annotation_stats$numDupsRemoved <- 0
# annotation_stats$numRemoved.ACCNUM <- 0
# try({data.filt <- genefilter::nsFilter(data, require.entrez=TRUE,
#                                        remove.dupEntrez=TRUE, var.func=IQR,
#                                        var.cutoff=0.5, var.filter=TRUE,
#                                        filterByQuantile=TRUE, feature.exclude="^AFFX")
# cat("Probe Aggregation Summary \n")
# str(data.filt$filter.log)
# filter.log <- data.filt$filter.log
# 
# data.filt <- data.filt$eset
# annotation_stats$gene_level_features <-dim(data.filt@featureData@data)[1]
# annotation_stats$numDupsRemoved <- filter.log$numDupsRemoved
# annotation_stats$numRemoved.ACCNUM <- filter.log$numRemoved.ACCNUM
#})

### Basic linear model fit
cat("\nConstructing linear model\n")
library(limma)


group__ <- factor(targets$t3$Group, levels = unique(targets$t3$Group))
design <- model.matrix(~ 0 + group__)
colnames(design)<-gsub("group__","",colnames(design)) #remove design name formatting
fit <- lmFit(data.filt, design)

if (is.fullrank(design) == FALSE){
  cat("The following groups are non estimable:",nonEstimable(design))
}

fit.groups <- colnames(fit$design)[which(fit$assign == 1)]
fit.index <-  which(levels(group__) %in% fit.groups)
fit.group.names <- unique(targets$t2$Group)

# ### Remove low expression probes
# cat("\nRemoving low expression probes for DGE\n")
# CutOff <- quantile(as.matrix(data),probs=.33)
# 
# hist_res <- graphics::hist(as.matrix(data.filt), 100, col = "cornsilk", freq = FALSE,
#                            main = "Probe Filtering Intensity Cutoff",
#                            border = "antiquewhite4",
#                            xlab = "Median intensities")
# 
# abline(v = CutOff, col = "coral4", lwd = 2)
# keep <- fit$Amean > CutOff
# fit <- fit[keep,] # filter out probes below cutoff expression level
# annotation_stats$expressed_genes <- dim(fit$genes)[1]
# path <- file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData","QC_Repports")
# setwd(path)
# 
# rm(hist_res,keep,CutOff,path)

### Create Contrast Model
cat("\nCreating contrast model\n")
combos<-combn(fit.groups,2) # generate matrix of pairwise group combinations for comparison
combos.names<-combn(fit.group.names,2)
contrasts<-c(paste(combos[1,],combos[2,],sep = "-"),paste(combos[2,],combos[1,],sep = "-")) # format combinations for limma:makeContrasts
contrast.names <-c(paste(combos.names[1,],combos.names[2,],sep = "v"),paste(combos.names[2,],combos.names[1,],sep = "v")) # format combinations for output table file names


cont.matrix <- makeContrasts(contrasts = contrasts,levels=design)

contrast.fit <- contrasts.fit(fit, cont.matrix)
contrast.fit <- eBayes(contrast.fit)
results<-decideTests(contrast.fit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0.5) # FDR .05
# try({
#   colnames(results@.Data) <- contrast.names
#   summary <- as.data.frame(summary(results))
#   summary <- summary[,c(2,1,3)]
#   colnames(summary)<-c("CONTRAST","REGULATION","GENE COUNT SIG")
#   DT::datatable(summary, caption = "Summary of Differentially Regulated Genes (P<=05)")
# })
rm(combos,combos.names,cont.matrix)

### Construct DGE Output Tables
cat("Building DGE tables\n")
dir.create(file.path(workdir,"Processed_Data",opt$glds,"02-Limma_DGE"), showWarnings = FALSE)
setwd(file.path(workdir,"Processed_Data",opt$glds,"02-Limma_DGE"))

output_table <- fit$genes
reduced_output_table <- fit$genes
cat("\nDim of fit$genes: ",dim(output_table),"\n")
#try(expr <- as.data.frame(data.filt$E[keep,]))
#try(expr <- as.data.frame(data.filt@assayData$exprs[rownames(data.filt@assayData$exprs) %in% rownames(fit$genes),]))
expr <- as.data.frame(data.filt@assayData$exprs)
cat("\nDim of expr: ",dim(expr),"\n")
cat("\nDim of data.filt.exprs: ",dim(data.filt@assayData$exprs),"\n")

output_table <- cbind(output_table,expr)
reduced_output_table <- cbind(reduced_output_table,expr)

# add all sample mean column
output_table$All.mean <- fit$Amean
reduced_output_table$All.mean <- fit$Amean
# add all sample stdev column
output_table$All.stdev <- contrast.fit$s2.post
reduced_output_table$All.stdev <- contrast.fit$s2.post
# add F statistic p-value (similar to ANOVA p-value) column
output_table$F.p.value <- contrast.fit$F.p.value
reduced_output_table$F.p.value <- contrast.fit$F.p.value
uu<- unique(targets$t2$Group)
# Add group mean values
group_means<-fit$coefficients
colnames(group_means)<-paste0("Group.Mean_",uu)
output_table<-cbind(output_table,group_means)
reduced_output_table<-cbind(reduced_output_table,group_means)
rm(group_means)
# add group stdev columns
group_stdev<-fit$stdev.unscaled * fit$coefficients
colnames(group_stdev)<-paste0("Group.Stdev_",uu)
output_table<-cbind(output_table,group_stdev)
reduced_output_table<-cbind(reduced_output_table,group_stdev)
rm(group_stdev)
# iterate through contrasts
for (i in 1:length(contrasts)){
  top <- topTable(contrast.fit, coef = i, number = Inf, genelist = contrast.fit$genes$ID, adjust.method = "BH", sort.by = "none")
  table <- top[,c(1,4,5)] # Pull columns for Log2fc, P.value, Adj.p.value
  colnames(table)<- c("Log2fc","P.value","Adj.p.value")
  table.reduced <- table
  table$Updown <- sign(top$logFC)
  table$Sig.1 <- top$adj.P.Val<=0.1
  table$Sig.05 <- top$adj.P.Val<=0.05
  table$Log2_P.value <- log2(top$P.Value) # For volcano plot
  table$Log2_Adj.p.value <- log2(top$adj.P.Val) # For volcano plot
  colnames(table.reduced)<-paste(colnames(table.reduced),contrast.names[i],sep = "_")
  colnames(table)<-paste(colnames(table),contrast.names[i],sep = "_")
  output_table<-cbind(output_table,table)
  reduced_output_table<-cbind(reduced_output_table,table.reduced)
}
rm(i,top,table,table.reduced)

### Export DGE Output Data Tables
setwd(file.path(workdir,"Processed_Data",opt$glds,"02-Limma_DGE"))
write.csv(reduced_output_table,"differential_expression.csv", row.names = FALSE)
write.csv(output_table,"visualization_output_table.csv", row.names = FALSE)
contrast.output <- contrast.fit$contrasts
row.names(contrast.output)<-uu
contrast.order <- order(match(contrasts,colnames(contrast.fit$contrasts)))

colnames(contrast.output)<-contrast.names
write.csv(contrast.output,"contrasts.csv")

rm (uu,group__,fit.index,fit.groups,fit.group.names,contrasts,contrast.names)

### Export Metadata files
dir.create(file.path(workdir,"Processed_Data",opt$glds,"Metadata"), showWarnings = FALSE)
path<-file.path(workdir,"Processed_Data",opt$glds,"Metadata")
setwd(path)
file.copy(from = opt$isa, to = file.path(path,basename(opt$isa)),overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
try(file.copy(from = Sys.glob(opt$probe), to = file.path(path,basename(opt$probe)),overwrite = FALSE, recursive = FALSE, copy.mode = FALSE))
file.copy(from = opt$runsheet, to = file.path(path,basename(opt$runsheet)), overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
rm(path)

cat("All data files have been written to:  ",file.path(workdir,"Processed_Data",opt$glds))
