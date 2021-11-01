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
  
  script:
    """
    #! /usr/bin/env Rscript

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
                      output_dir=getwd())
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
  
  script:
    """
    #! /usr/bin/env Rscript

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
 
    load("${ runsheet_RData }") # named 'targets'
    load("${ raw_files_RData }") # named 'raw'
   
    # 4 
    
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
                      output_dir=getwd())
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

    save(normalized, file="normalized.RData")

    ### extract expression table and save to csv
    expression <- data.frame(Biobase::exprs(normalized))
    write.table(expression,"normalized.txt",quote=FALSE, append=FALSE, sep = "\t", col.names=NA)
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

