process LOAD_RUNSHEET {
  conda = "${projectDir}/envs/minimal.yml"

  input:
    path(runsheet)
    val(gldsAccession)
    val(meta)
  
  output:
    path("runsheet.RData")
  
  script:
    """
    #! /usr/bin/env Rscript
    #2

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
    
    source(file.path(codebase_dir, "microarray_functions.R"))

    targets <- stagingTargets("${ runsheet }")  
    
    targets\$glds <- "${ gldsAccession }"
    targets\$platform <- "${ meta.platform }"
    targets\$organism <- "${ meta.organism_sci }"
    print(targets)
    
    save(targets, file="runsheet.RData") 
    """
}

process QA_NORMALIZED {
  conda = "${projectDir}/envs/minimal.yml"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/01-NormalizedData",
    mode: params.publish_dir_mode

  input:
    path(runsheet_RData)
    path(normalized_files_RData)
  
  output:
    path("normalized_qa.html")
    path("visualization_PCA_table.csv")
  
  script:
    """
    #! /usr/bin/env Rscript
    #2

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
 
    load("${ runsheet_RData }") # named 'targets'
    load("${ normalized_files_RData }") # named 'normalized'

    print(class(normalized))
    
    
    if ((targets\$platform == "Affymetrix") | (targets\$platform == "Nimblegen 1-channel")){
      plots <- "OLIGO"
    }else if((targets\$platform == "Illumina Expression") | (targets\$platform == "Agilent 1-channel") | (targets\$platform == "Agilent 2-channel") | (targets\$platform == "Nimblegen 2-channel")){
      plots <- "LIMMA"
    }else {
      plots <- NULL
    }
 
    rmarkdown::render(file.path(codebase_dir,"qa_summary_normalized.rmd"),
                      "html_document", 
                      output_file="normalized_qa",
                      output_dir=getwd()) # generates PCA_raw
    write.csv(PCA_raw\$x, file = "visualization_PCA_table.csv")
    """
}


process QA_RAW {
  conda = "${projectDir}/envs/minimal.yml"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/00-RawData",
    mode: params.publish_dir_mode

  input:
    path(runsheet_RData)
    path(raw_files_RData)
  
  output:
    path("raw_qa.html")
    path("raw_visualization_PCA_table.csv")
  
  script:
    """
    #! /usr/bin/env Rscript
    # 3i 

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
 
    load("${ runsheet_RData }") # named 'targets'
    load("${ raw_files_RData }") # named 'raw'
   
    
    if ((targets\$platform == "Affymetrix") | (targets\$platform == "Nimblegen 1-channel")){
      plots <- "OLIGO"
    }else if((targets\$platform == "Illumina Expression") | (targets\$platform == "Agilent 1-channel") | (targets\$platform == "Agilent 2-channel") | (targets\$platform == "Nimblegen 2-channel")){
      plots <- "LIMMA"
    }else {
      plots <- NULL
    }
 
    rmarkdown::render(file.path(codebase_dir,"qa_summary_raw.rmd"),
                      "html_document", 
                      output_file="raw_qa",
                      output_dir=getwd()) # generates PCA_raw

    write.csv(PCA_raw\$x, file = "raw_visualization_PCA_table.csv")
    """
}

process NORMALIZE {
  conda = "${projectDir}/envs/minimal.yml"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/01-NormalizedData",
    mode: params.publish_dir_mode,
    pattern: "normalized.txt"

  input:
    path(raw_files_RData)
    path(runsheet_RData)
  
  output:
    path("normalized.RData"), emit: rdata
    path("normalized.txt")
  
  script:
    """
    #! /usr/bin/env Rscript

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
 
    load("${ raw_files_RData }") # named 'raw'
    
    ### Background Correction and Normalization
    if (class(raw)=="ExonFeatureSet" || class(raw)=="GeneFeatureSet"){
      normalized <- oligo::rma(raw, target = "core", background=TRUE, normalize=TRUE)
      normalized.bgonly <- oligo::rma(raw, target = "core", background=TRUE, normalize=FALSE)
      cat("RMA background correction and quantile normalization performed with gene level summarization.\n")
    }
    
    if (class(raw)=="ExpressionFeatureSet"){
      normalized <- oligo::rma(raw, normalize = TRUE, background = TRUE)
      normalized.bgonly <- oligo::rma(raw, normalize = FALSE, background = TRUE)
      cat("RMA background correction and quantile normalization performed.\n")
    }


    ### extract expression table and save to csv
    expression_df <- data.frame(Biobase::exprs(normalized))
    save(normalized, expression_df, file="normalized.RData")
    write.table(expression_df,"normalized.txt",quote=FALSE, append=FALSE, sep = "\t", col.names=NA)
    """
}

process READ_RAW {
  conda = "${projectDir}/envs/minimal.yml"

  input:
    path(raw_files)
    path(runsheet_RData)
  
  output:
    path("raw.Rdata")
  
  script:
    raw_file_extension = raw_files[0].extension
    if ( raw_file_extension == "CEL")
    """
    #! /usr/bin/env Rscript
    library(oligo)
    library(stringr)
    
    load("${ runsheet_RData }") # named 'targets'

    # remove potential .gz from FileNames
    # decompression happens during staging
    target_filenames <- str_replace_all(targets\$t3\$FileName, ".gz\$", "")
    
    raw <- oligo::read.celfiles(target_filenames, sampleNames=targets\$t3\$SampleName)

    save(raw, file="raw.Rdata")
    print("Saved to 'raw.Rdata'")
    """
    else
    """
    echo "No read in for this file extension is known"
    """
}

process DETERMINE_PLATFORM_DESIGN {
  conda = "${projectDir}/envs/minimal.yml"

  input:
    path(runsheet)
  
  output:
    path("platform_design.txt")
  
  script:
    """
    #! /usr/bin/env Rscript
    library(oligo)
    
    celFiles <- Sys.glob(file.path("raw_files","*.CEL"))
    print(celFiles)
    raw <- oligo::read.celfiles(celFiles)

    save(raw, file="raw.Rdata")
    print("Saved to 'raw.Rdata'")
    """
}

process IMPORT_PROBE {
  conda = "${projectDir}/envs/minimal.yml"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/01-NormalizedData",
    mode: params.publish_dir_mode

  input:
    path(raw_files_RData)
    path(runsheet_RData)
    path(normalized_RData)
  
  output:
    path("normalized-annotated.txt")
    path("probe_annotations.txt")
    path("normalized-annotated.rda"), emit: annotated_rdata
  
  script:
    """
    #! /usr/bin/env Rscript

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
    organism_table <- read.csv(file = file.path(codebase_dir,"organisms.csv"), header = TRUE, stringsAsFactors = FALSE)

    load("${ raw_files_RData }") # named 'raw'
    load("${ runsheet_RData }") # named 'targets'
    load("${ normalized_RData }") # named 'expression'
    
    ### Import Probe Data
    if (length(targets\$probe >= 1)){
      options(connectionObserver = NULL)
      database <- sub('\\\\.annotation.tar.gz\$', '', basename(targets\$probe)) 
      cat("\nLoading local probe annotation database: ",database,"\n")
      if(!require(database, character.only=TRUE)) {
        BiocManager::install(database, ask = FALSE)
      }
      install.packages(targets\$probe,repos = NULL, verbose = FALSE, quiet = TRUE)
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
    keys<-rownames(expression_df)
   
    print(keytype)
    print(keys)
 
    ### Map assay database annotations
    annotation <- data.frame(REFSEQ=mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "REFSEQ",multiVals = "first"))
    
    try(annotation\$ENSEMBL<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "ENSEMBL",multiVals = "first"))
    try(annotation\$SYMBOL<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "SYMBOL",multiVals = "first"))
    try(annotation\$GENENAME<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "GENENAME",multiVals = "first"))
    try(annotation\$ENTREZID<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "ENTREZID",multiVals = "first"))
    try(annotation\$TAIR<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "TAIR",multiVals = "first"))
    try(annotation\$GOSLIM_IDS<-mapIds(eval(parse(text = database),env=.GlobalEnv),keys = keys,keytype = keytype, column = "GO",multiVals = "first"))
    
    
    ### Map STRING annotations
    try({
      string_db <- STRINGdb::STRINGdb\$new( version="11", species=organism_table\$taxon[organism_table\$species == targets\$species],score_threshold=0)
      string_map<-string_db\$map(annotation,"ENTREZID",removeUnmappedRows = TRUE, takeFirst = TRUE)
      string_cols <-string_map[,c("ENTREZID","STRING_id")]
      string_cols <- string_cols[!duplicated(string_cols\$ENTREZID),]
      annotation <- dplyr::left_join(annotation,string_cols,by="ENTREZID")
    
      rm(string_map,string_db)
    })
    rm(keytype,keys)
    
    ### Generate normalized annotated expression text file
    cat("\nGenerating normalized-annotated.txt file\n")
    expression_df <- cbind(annotation,expression_df)
    write.table(expression_df,"normalized-annotated.txt",quote=FALSE, append = FALSE, row.names = FALSE, sep = "\t")
    write.table(annotation,"probe_annotations.txt",quote=FALSE, append = FALSE, row.names = FALSE, sep = "\t")

    ### Annotate the expression set object and save as a file
    cat("\nGenerating normalized-annotated.rda file\n")

    fData(normalized)<-annotation
    save(normalized, file = "normalized-annotated.rda")
    """
}

process DGE {
  conda = "${projectDir}/envs/minimal.yml"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/02-DGE",
    mode: params.publish_dir_mode

  input:
    path(normalized_annotated_RData)
    path(runsheet_RData)
  
  output:
   path("contrasts.csv")
   path("visualization_output_table.csv")
   path("differential_expression.csv")  

  script:
    """
    #! /usr/bin/env Rscript
    #2

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
    
    source(file.path(codebase_dir, "microarray_functions.R"))

    load("${ normalized_annotated_RData }") # named 'normalized'
    load("${ runsheet_RData }") # named 'targets'

    
    ### Gene level estimation by maximum interquartile range
    cat("\nPerforming gene level estimation by max interquartile range")
    normalized.filt <- normalized
    ### Basic linear model fit
    cat("\nConstructing linear model\n")
    library(limma)
    
    
    group__ <- factor(targets\$t3\$Group, levels = unique(targets\$t3\$Group))
    design <- model.matrix(~ 0 + group__)
    colnames(design)<-gsub("group__","",colnames(design)) #remove design name formatting
    fit <- lmFit(normalized.filt, design)
    
    if (is.fullrank(design) == FALSE){
      cat("The following groups are non estimable:",nonEstimable(design))
    }
    
    fit.groups <- colnames(fit\$design)[which(fit\$assign == 1)]
    fit.index <-  which(levels(group__) %in% fit.groups)
    fit.group.names <- unique(targets\$t2\$Group)



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
        
    ### Construct DGE Output Tables
    cat("Building DGE tables\n")
    output_table <- fit\$genes
    reduced_output_table <- fit\$genes
    cat("\nDim of fit\$genes: ",dim(output_table),"\n")
    expr <- as.data.frame(normalized.filt@assayData\$exprs)
    cat("\nDim of expr: ",dim(expr),"\n")
    cat("\nDim of normalized.filt.exprs: ",dim(normalized.filt@assayData\$exprs),"\n")
    
    output_table <- cbind(output_table,expr)
    reduced_output_table <- cbind(reduced_output_table,expr)
    
    # add all sample mean column
    output_table\$All.mean <- fit\$Amean
    reduced_output_table\$All.mean <- fit\$Amean
    # add all sample stdev column
    output_table\$All.stdev <- contrast.fit\$s2.post
    reduced_output_table\$All.stdev <- contrast.fit\$s2.post
    # add F statistic p-value (similar to ANOVA p-value) column
    output_table\$F.p.value <- contrast.fit\$F.p.value
    reduced_output_table\$F.p.value <- contrast.fit\$F.p.value
    uu<- unique(targets\$t2\$Group)
    # Add group mean values
    group_means<-fit\$coefficients
    colnames(group_means)<-paste0("Group.Mean_",uu)
    output_table<-cbind(output_table,group_means)
    reduced_output_table<-cbind(reduced_output_table,group_means)
    rm(group_means)
    # add group stdev columns
    group_stdev<-fit\$stdev.unscaled * fit\$coefficients
    colnames(group_stdev)<-paste0("Group.Stdev_",uu)
    output_table<-cbind(output_table,group_stdev)
    reduced_output_table<-cbind(reduced_output_table,group_stdev)
    rm(group_stdev)
    # iterate through contrasts
    for (i in 1:length(contrasts)){
      top <- topTable(contrast.fit, coef = i, number = Inf, genelist = contrast.fit\$genes\$ID, adjust.method = "BH", sort.by = "none")
      table <- top[,c(1,4,5)] # Pull columns for Log2fc, P.value, Adj.p.value
      colnames(table)<- c("Log2fc","P.value","Adj.p.value")
      table.reduced <- table
      table\$Updown <- sign(top\$logFC)
      table\$Sig.1 <- top\$adj.P.Val<=0.1
      table\$Sig.05 <- top\$adj.P.Val<=0.05
      table\$Log2_P.value <- log2(top\$P.Value) # For volcano plot
      table\$Log2_Adj.p.value <- log2(top\$adj.P.Val) # For volcano plot
      colnames(table.reduced)<-paste(colnames(table.reduced),contrast.names[i],sep = "_")
      colnames(table)<-paste(colnames(table),contrast.names[i],sep = "_")
      output_table<-cbind(output_table,table)
      reduced_output_table<-cbind(reduced_output_table,table.reduced)
    }
    
    ### Export DGE Output Data Tables
    write.csv(reduced_output_table,"differential_expression.csv", row.names = FALSE)
    write.csv(output_table,"visualization_output_table.csv", row.names = FALSE)
    contrast.output <- contrast.fit\$contrasts
    row.names(contrast.output)<-uu
    contrast.order <- order(match(contrasts,colnames(contrast.fit\$contrasts)))
    
    colnames(contrast.output)<-contrast.names
    write.csv(contrast.output,"contrasts.csv")
    
    
    """
}
