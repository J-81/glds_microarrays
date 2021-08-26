### Agilent Two Channel Microarray Processing

### Import Raw Data

cat("\nStarting Agilent 2-channel Processing Pipeline\n")

filetypes <- unique(toupper(tools::file_ext(list.files(file.path(tempin,"00-RawData")))))
if (length(filetypes)>1){
  cat("\nFiletypes\n",filetypes)
  stop("Execution halted due to multiple filetypes present for raw data")
}

cat("\n Filetype: ",filetypes, "\n")

workdir <- opt$out
library(limma)

if (filetypes=="GPR"){
  raw <- read.maimages(basename(opt$files), source="genepix", path=file.path(tempin,"00-RawData"),wt.fun=wtflags(weight=0,cutoff=-50))
}else if (filetypes=="TXT" || filetypes=="RAW.TXT"){
  raw <- read.maimages(basename(opt$files), source="agilent", path=file.path(tempin,"00-RawData"))
}else{
  stop("Unsupported file extension.")
}

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
  try(rmarkdown::render(file.path(codebase_dir,"qa_summary_raw.rmd"),"html_document", output_file="raw_qa",output_dir=file.path(workdir,"Processed_Data",opt$glds,"00-RawData")))
}

### Import Probe Annotation

# if (grepl("\\.soft$", opt$probe) || grepl("\\.gpl.txt$", opt$probe) || grepl("GPL",opt$probe)){
#     probedata<- GEOquery:::parseGPL(opt$probe)
#     probedata<- probedata@dataTable@table
#     probedata <- probedata[with(probedata, order(ROW,COL)),] # sort by row & col to match RG order
#     cat("\nLocal GPL probe annotation imported\n")
#     RG$genes <- dplyr::left_join(RG$genes,probedata,by = c("Name" = "NAME"))
# } else if (grepl("\\.adf.txt$", opt$probe)){
#   probedata<- ArrayExpress:::readFeatures(opt$probefile$name, dirname(opt$probefile$datapath),ArrayExpress:::skipADFheader(opt$probefile$name, dirname(opt$probefile$datapath),proc=F))
#   cat("\nLocal ADF probe annotation imported\n")
# }else{
#   cat("\nAnnotation will be attempted with raw file ID keys\n")
# }


### Map Annotations from Primary Keytype
try(annotation<-data.frame(ID=raw$genes$ID)) #assumes Array ID is primary annotation
try(annotation<-data.frame(ID=raw$genes$GeneName))
try(annotation<-data.frame(ID=raw$genes$Name)) #assumes Array ID is primary annotation
annotation$ID <- sub("\\.\\d+", "", annotation$ID) #remove any version suffixes from IDs

keytype <- NULL
testout <- NULL
testkeys <- annotation$ID[100:200]
cat("\nTesting keys against available keytypes.\n")
for(annkey in columns(eval(parse(text = ann.dbi)))){
  try(testout<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = testkeys,keytype = annkey, column = "ENTREZID",multiVals = "first"))
  if(length(testout)>10){
    keytype<-annkey
    cat("\nKeytype identified: ",keytype,"\n")
    break;
  }
}
try(annotation$REFSEQ<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "REFSEQ",multiVals = "first")))
try(annotation$ENSEMBL<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "ENSEMBL",multiVals = "first")))
try(annotation$SYMBOL<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "SYMBOL",multiVals = "first"))) #GENENAME
try(annotation$GENENAME<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "GENENAME",multiVals = "first"))) #DESCRIPTION
try(annotation$TAIR<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "TAIR",multiVals = "first")))
try(annotation$ORF<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "ORF",multiVals = "first")))
try(annotation$ENTREZID<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "ENTREZID",multiVals = "first")))
try(annotation$GOSLIM_IDS<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "GO",multiVals = "first")))

annotation <- data.frame(annotation,stringsAsFactors = FALSE)
is.na(annotation) <- annotation == "NULL"


### Map STRING DB annotations
try({
  string_db <- STRINGdb::STRINGdb$new( version="11", species=organism_table$taxon[organism_table$species == opt$species],score_threshold=0)
  string_map<-string_db$map(annotation,"REFSEQ",removeUnmappedRows = FALSE, takeFirst = TRUE)
  annotation$STRING_ID <- string_map$STRING_id
  rm(string_map,string_db)
})

annotation.subset <- annotation
raw.annotated <- raw


### Background Correction and Normalization

raw.annotated <- backgroundCorrect(raw.annotated, method="normexp", offset=50,normexp.method="saddle")
cat("\nBackground correction by NormExp\n")
MA <- normalizeWithinArrays(raw.annotated, method="loess", weights = raw.annotated$weights)
cat("\nWithin Array Normalization by Loess\n")
MA <- normalizeBetweenArrays(MA, method="Aquantile")
cat("\nBetween array normalization by Aquantile\n")
MA.summarized <- MA
MA.summarized$genes<-annotation.subset


cat("\n Initial feature count: ",dim(MA.summarized$genes),"\n")
### Filter out Control Probes
cat("\n MA columns",colnames(MA.summarized$genes))
try({
    filter <- which(!(MA.summarized$genes$ID %in% c("DarkCorner","GE_BrightCorner","NegativeControl")))
    MA.summarized<- MA.summarized[filter,]
    annotation.subset <- annotation.subset[filter,]})
cat("\n Features after control probes removed: ",dim(MA.summarized$genes),"\n")


### Filter out Non-annotated Probes

try({
    MA.summarized <- MA.summarized[!is.na(MA.summarized$genes$ENTREZID),]
    annotation.subset <- MA.summarized$genes
    #annotation_stats$annotated_genes <- dim(MA.summarized)[1]
  })

cat("\n Features after unannotated probes removed: ",dim(MA.summarized$genes),"\n")


### Generate Gene Summary Values

# MA.summarized <- avereps(MA.summarized, ID = MA.summarized$genes$ID)
# annotation.subset <- MA.summarized$genes
#annotation_stats$unique_probes <- dim(MA.summarized)[1]

# MA.summarized$genes <- MA.summarized$genes[,-c("ID")]
# annotation.subset <- annotation.subset[,-c("ID")]

###  Generate Expression Value  Files
dir.create(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"), showWarnings = FALSE)
setwd(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))

### Normalized QA Report
if(opt$reports == TRUE){
  rmarkdown::render(file.path(codebase_dir,"qa_summary_normalized.rmd"),"html_document", output_file="normalized_qa",output_dir=file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))
}

#annotation.subset <- dplyr::select(MA.summarized$genes, contains(c("ID")))
#colnames(annotation.subset)[1]<-opt$primary
#annotation.subset <-apply(annotation.subset, 1, unlist)
expression <- as.data.frame(MA.summarized$M)
write.table(expression,"normalized.txt",quote=FALSE, append=FALSE, sep = "\t", col.names=NA)

### Write out annotated expression values

expression <- cbind(annotation.subset,expression)
write.table(expression,"normalized-annotated.txt",quote=FALSE, append=FALSE, sep = "\t", row.names = FALSE)
write.table(annotation,"probe_annotations.txt",quote=FALSE, append = FALSE, row.names = FALSE, sep = "\t")
### Write out annotated MA file
save(MA.summarized,file = "normalized-annotated.rda")

### Build Differential Gene Expression Model

library(statmod)

if (targets$design == "Replicate Array"){
    cat("\nReplicate Array Design Detected\n")
    uu<-uniqueTargets(targets$t2)
    uuu<-uniqueTargets(targets$t3)
    fff <- factor(targets$t3$Cy3, levels=uuu)
    fit <- lmFit(MA.summarized, ref=uuu[1])
    contrasts<-c(paste(uuu[1],uuu[2],sep = "-"),paste(uuu[2],uuu[1],sep = "-"))
    fit <- eBayes(fit)
    contrast.names<-c(paste(uu[1],uu[2],sep = "v"),paste(uu[2],uu[1],sep = "v"))
    results<-decideTests(fit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0.5)
    colnames(results)<-contrast.names[2]
    cat("\nContrast Names: ",contrasts,"\n")
}

if (targets$design == "Common Reference"){

    uu<-unique(targets2[,2])
    uuu1<-unique(targets3[,2])
    uuu2<-unique(targets3[,3]) #confirm reference column
    design <- modelMatrix(targets3, ref = uuu2)
    fit <- lmFit(MA.summarized, design)
    
    f <- factor(targets3[,2], levels=unique(uuu1)) #confirm non reference column
    
    # Create Contrast Model
    combos<-combn(levels(f),2) # generate matrix of pairwise group combinations for comparison
    ff <- factor(targets2[,2], levels=unique(uu))
    combos.names<-combn(levels(ff),2)
    contrasts<-c(paste(combos[1,],combos[2,],sep = "-"),paste(combos[2,],combos[1,],sep = "-")) # format combinations for limma:makeContrasts
    contrast.names <-c(paste(combos.names[1,],combos.names[2,],sep = "v"),paste(combos.names[2,],combos.names[1,],sep = "v")) # format combinations for output table file names
    
    
    cont.matrix <- makeContrasts(contrasts = contrasts,levels=design)
    contrast.fit <- contrasts.fit(fit, cont.matrix)
    contrast.fit <- eBayes(contrast.fit)
    results<-decideTests(contrast.fit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0.5) # FDR .05
}
  
if (targets$design == "Separate Channels"){
    targets$t2 <- targetsA2C(targets$t2, channel.codes = c(1,2), channel.columns = list(Target=c("Cy3","Cy5")),grep = FALSE)
    targets$t3 <- targetsA2C(targets$t3, channel.codes = c(1,2), channel.columns = list(Target=c("Cy3","Cy5")),grep = FALSE)
    
    uuu <- unique(targets$t3$Target)
    uu <- unique(targets$t2$Target)
    fff <- factor(targets$t3$Target, levels=uuu)
    design <- model.matrix(~0+fff)
    colnames(design) <- uuu
    corfit <- intraspotCorrelation(MA.summarized, design)
    fit <- lmscFit(MA.summarized, design, correlation=corfit$consensus)
    # if (opt$filter_nonexpressed){
    #   
    #   CutOff <- quantile(fit$Amean,probs=.33)
    #   
    #   hist_res <- graphics::hist(as.matrix(fit$coefficients), 100, col = "cornsilk", freq = FALSE, 
    #                              main = "Probe Filtering Intensity Cutoff",
    #                              border = "antiquewhite4",
    #                              xlab = "Median intensities")
    #   
    #   abline(v = CutOff, col = "coral4", lwd = 2)
    #   
    #   keep <- fit$Amean > CutOff
    #   fit <- fit[keep,] # filter out probes below cutoff expression level
    #   annotation.subset <- annotation.subset[keep,]
    #   MA.summarized <- MA.summarized[keep,]
    # }
    
    # Create Contrast Model
    combos<-combn(levels(fff),2) # generate matrix of pairwise group combinations for comparison
    uu <- unique(targets$t2$Target)
    ff <- factor(targets$t2$Target, levels=uu)
    combos.names<-combn(levels(ff),2)
    contrasts<-c(paste(combos[1,],combos[2,],sep = "-"),paste(combos[2,],combos[1,],sep = "-")) # format combinations for limma:makeContrasts
    contrast.names <-c(paste(combos.names[1,],combos.names[2,],sep = "v"),paste(combos.names[2,],combos.names[1,],sep = "v")) # format combinations for output table file names
    cont.matrix <- makeContrasts(contrasts = contrasts,levels=design)
    contrast.fit <- contrasts.fit(fit, cont.matrix)
    contrast.fit <- eBayes(contrast.fit, trend = TRUE)
    results<-decideTests(contrast.fit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0.5)
    
}

### Generate DGE Output Files

dir.create(file.path(workdir,"Processed_Data",opt$glds,"02-Limma_DGE"), showWarnings = FALSE)
setwd(file.path(workdir,"Processed_Data",opt$glds,"02-Limma_DGE"))

########## Construct Output Data Tables
# output_table <- annotation.subset[,-c(1)]
# reduced_output_table <- annotation.subset[,-c(1)]
output_table <- annotation.subset
reduced_output_table <- annotation.subset
output_table <- cbind(output_table,as.data.frame(MA.summarized$A))
reduced_output_table <- cbind(reduced_output_table,as.data.frame(MA.summarized$A))

if (targets$design == "Replicate Array"){
  cat("\nReplicate array DGE block 1\n")
  output_table$All.mean <- fit$Amean
  reduced_output_table$All.mean <- fit$Amean
  
  output_table$All.stdev <- fit$s2.post
  reduced_output_table$All.stdev <- fit$s2.post
  
  output_table$F.p.value <- fit$F.p.value
  reduced_output_table$F.p.value <- fit$F.p.value
  
}else {
  output_table$All.mean <- contrast.fit$Amean
  reduced_output_table$All.mean <- contrast.fit$Amean
  
  output_table$All.stdev <- contrast.fit$s2.post
  reduced_output_table$All.stdev <- contrast.fit$s2.post
  
  output_table$F.p.value <- contrast.fit$F.p.value
  reduced_output_table$F.p.value <- contrast.fit$F.p.value
  
}


if (targets$design == "Separate Channels"){
  
  ########## Add Group Mean Values
  group_means<-as.data.frame(fit$coefficients)
  colnames(group_means)<-paste0("Group.Mean_",uu)
  output_table<-cbind(output_table,group_means)
  reduced_output_table<-cbind(reduced_output_table,group_means)
  rm(group_means)
  
  # add group stdev columns
  group_stdev<-as.data.frame(fit$stdev.unscaled * fit$coefficients)
  colnames(group_stdev)<-paste0("Group.Stdev_",uu)
  output_table<-cbind(output_table,group_stdev)
  reduced_output_table<-cbind(reduced_output_table,group_stdev)
  rm(group_stdev)  
  
}

if (targets$design == "Replicate Array"){
  cat("\nReplicate array DGE block 2\n")
  uu<-uu[1]
  ########## Add Group Mean Values
  group_means<-as.data.frame(fit$coefficients)
  colnames(group_means)<-paste0("Group.Mean_",uu)
  output_table<-cbind(output_table,group_means)
  reduced_output_table<-cbind(reduced_output_table,group_means)
  rm(group_means)
  
  # add group stdev columns
  group_stdev<-as.data.frame(fit$stdev.unscaled * fit$coefficients)
  colnames(group_stdev)<-paste0("Group.Stdev_",uu)
  output_table<-cbind(output_table,group_stdev)
  reduced_output_table<-cbind(reduced_output_table,group_stdev)
  rm(group_stdev)  
  
}


if (targets$design == "Replicate Array"){
  # Contrast 1
  cat("\nReplicate array DGE block 3\n")
  
  top <- topTable(fit, coef = 1, number = Inf, genelist = fit$genes$ID, adjust.method = "BH", sort.by = "none")
  table <- top[,c(2,5,6)] # Pull columns for Log2fc, P.value, Adj.p.value
  colnames(table)<- c("Log2fc","P.value","Adj.p.value")
  table.reduced <- table
  table$Updown <- sign(top$logFC)
  table$Sig.1 <- top$adj.P.Val<=0.1
  table$Sig.05 <- top$adj.P.Val<=0.05
  table$Log2_P.value <- log2(top$P.Value) # For volcano plot
  table$Log2_Adj.p.value <- log2(top$adj.P.Val) # For volcano plot
  colnames(table.reduced)<-paste(colnames(table.reduced),contrast.names[1],sep = "_")
  colnames(table)<-paste(colnames(table),contrast.names[1],sep = "_")
  output_table<-cbind(output_table,table)
  reduced_output_table<-cbind(reduced_output_table,table.reduced)
  
  # Repeat for contrast 2
  table <- top[,c(2,5,6)] # Pull columns for Log2fc, P.value, Adj.p.value
  colnames(table)<- c("Log2fc","P.value","Adj.p.value")
  table$Log2fc <- -1*table$Log2fc
  table.reduced <- table
  table$Updown <- sign(table$Log2fc)
  table$Sig.1 <- top$adj.P.Val<=0.1
  table$Sig.05 <- top$adj.P.Val<=0.05
  table$Log2_P.value <- log2(top$P.Value) # For volcano plot
  table$Log2_Adj.p.value <- log2(top$adj.P.Val) # For volcano plot
  
  colnames(table.reduced)<-paste(colnames(table.reduced),contrast.names[2],sep = "_")
  colnames(table)<-paste(colnames(table),contrast.names[2],sep = "_")
  output_table<-cbind(output_table,table)
  reduced_output_table<-cbind(reduced_output_table,table.reduced)
}

if (targets$design == "Separate Channels" || targets$design == "Common Reference"){
  
  # iterate through contrasts
  for (i in 1:length(contrasts)){
    top <- topTable(contrast.fit, coef = i, number = Inf, genelist = contrast.fit$genes$ID, adjust.method = "BH", sort.by = "none")
    table <- top[,c(2,5,6)] # Pull columns for Log2fc, P.value, Adj.p.value
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
}

### Export DGE Output Data Tables


write.csv(reduced_output_table,"differential_expression.csv", row.names = FALSE)
write.csv(output_table,"visualization_output_table.csv", row.names = FALSE)


if (targets$design == "Replicate Array"){
  cat("\nReplicate array DGE block 4\n")
  
  #contrast.names<-c(paste(uu[1],uu[2],sep = "v"),paste(uu[2],uu[1],sep = "v"))
  write.csv(contrast.names,"contrasts.csv")
} else {
  contrast.output <- contrast.fit$contrasts
  row.names(contrast.output)<-uu
  colnames(contrast.output)<-contrast.names
  write.csv(contrast.output,"contrasts.csv")
}

### Export Metadata files
dir.create(file.path(workdir,"Processed_Data",opt$glds,"Metadata"), showWarnings = FALSE)
path<-file.path(workdir,"Processed_Data",opt$glds,"Metadata")
setwd(path)
file.copy(from = opt$isa, to = file.path(path,basename(opt$isa)),overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
try(file.copy(from = Sys.glob(opt$probe), to = file.path(path,basename(opt$probe)),overwrite = FALSE, recursive = FALSE, copy.mode = FALSE))
file.copy(from = opt$runsheet, to = file.path(path,basename(opt$runsheet)), overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
rm(path)

cat("All data files have been written to:  ",file.path(workdir,"Processed_Data",opt$glds))
