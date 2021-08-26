### Agilent One Channel Processing

### Import Raw Data
cat("\nStarting Agilent 1-channel Processing Pipeline\n")

filetypes <- unique(toupper(tools::file_ext(list.files(file.path(tempin,"00-RawData")))))
if (length(filetypes)>1){
  cat("\nFiletypes\n",filetypes)
  stop("Execution halted due to multiple filetypes present for raw data")
}

cat("\n Filetype: ",filetypes, "\n")

workdir <- opt$out
library(limma)

if (filetypes=="GPR"){
  raw <- read.maimages(basename(opt$files), source="genepix", path=file.path(tempin,"00-RawData"),wt.fun=wtflags(weight=0,cutoff=-50), green.only=TRUE, names = targets$t1$SampleName)
}else if (filetypes=="TXT" || filetypes=="RAW.TXT"){
  raw <- read.maimages(basename(opt$files), source="agilent", path=file.path(tempin,"00-RawData"), green.only=TRUE, names = targets$t1$SampleName)
}else{
  stop("Unsupported file extension.")
}

### If data imports properly, start building output files

dir.create(file.path(workdir,"Processed_Data"), showWarnings = FALSE)
dir.create(file.path(workdir,"Processed_Data",opt$glds), showWarnings = FALSE)
dir.create(file.path(workdir,"Processed_Data",opt$glds,"00-RawData"), showWarnings = FALSE)
path <- file.path(workdir,"Processed_Data",opt$glds,"00-RawData")
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

data <- backgroundCorrect(raw, method="normexp", offset=50)
data.bgonly <- backgroundCorrect(raw, method="normexp", offset=50)
data <- normalizeBetweenArrays(data, method="quantile")
cat("\nNormexp background correction and Quantile normalization performed.\n")

### Normalized QA Report
if(opt$reports == TRUE){
  rmarkdown::render(file.path(codebase_dir,"qa_summary_normalized.rmd"),"html_document", output_file="normalized_qa",output_dir=file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))
}

### Begin tracking feature annotation stats
# annotation_stats <- list()
# annotation_stats$total_features <- dim(data$genes)[1]

###  Write out the expression values
dir.create(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"), showWarnings = FALSE)
setwd(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))
data <- data[data$genes$ControlType == 0,]

expression<-cbind(data$genes$ProbeUID,data$E)
cat("\n Length expression features: ",dim(expression),"\n")
write.table(expression,"normalized.txt",quote=FALSE, append=FALSE, sep = "\t")

### Import Probe Annotation

database<-ann.dbi

cat("\nKey headers: ",colnames(data$genes),"\n")
keys<-data$genes$SystematicName
keys <- sub("\\.\\d+", "", keys) #remove any version suffixes from IDs
keys <- toupper(keys)

keytype <- NULL
testout <- NULL
testkeys <- keys[100:200]
cat("\nTesting keys against available keytypes.\n")
for(annkey in columns(eval(parse(text = ann.dbi)))){
  try(testout<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = testkeys,keytype = annkey, column = "ENTREZID",multiVals = "first"))
  if(length(testout)>10){
    keytype<-annkey
    cat("\nKeytype identified: ",keytype,"\n")
    break;
  }
}

### Map assay database annotations
annotation <- data.frame(ID=keys)
try(annotation$REFSEQ<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "REFSEQ",multiVals = "first"))
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
cat("\n Length annotation features: ",dim(annotation),"\n")

### Generate normalized annotated expression text file
setwd(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))
expression <- cbind(annotation,expression)
write.table(expression,"normalized-annotated.txt",quote=FALSE, append = FALSE, row.names = FALSE, sep = "\t")
write.table(annotation,"probe_annotations.txt",quote=FALSE, append = FALSE, row.names = FALSE, sep = "\t")
### Annotate the expression set object and save as a file
setwd(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))
data$annotation<-annotation
save(data,file = "normalized-annotated.rda")

data.filt <- data

### Basic linear model fit

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

### Create Contrast Model
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
# rm(combos,combos.names,cont.matrix)

### Construct DGE Output Tables

dir.create(file.path(workdir,"Processed_Data",opt$glds,"02-Limma_DGE"), showWarnings = FALSE)
setwd(file.path(workdir,"Processed_Data",opt$glds,"02-Limma_DGE"))

output_table <- data.filt$annotation
reduced_output_table <- data.filt$annotation
cat("\n 176: annotation: ",dim(annotation), colnames(annotation),"\n")
cat("\n 177: fit$genes:",dim(fit$genes),colnames(fit$genes),"\n")
cat("\n 178: data.filt$annotation:",dim(data.filt$annotation),colnames(data.filt$annotation),"\n")
try(expr <- as.data.frame(data.filt$E))
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

# Add group mean columns
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
try(file.copy(from = opt$probe, to = file.path(path,basename(opt$probe)),overwrite = FALSE, recursive = FALSE, copy.mode = FALSE))
file.copy(from = opt$runsheet, to = file.path(path,basename(opt$runsheet)), overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
rm(path)
cat("All data files have been written to:  ",file.path(workdir,"Processed_Data",opt$glds))
