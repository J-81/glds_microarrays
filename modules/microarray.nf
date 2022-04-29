process LOAD_RUNSHEET {

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
  publishDir "${ params.outputDir }/${ params.gldsAccession }/01-NormalizedData",
    mode: params.publish_dir_mode

  input:
    path("normalized.RData")
    path("targets.RData")
  
  output:
    path("normalized_qa.html")
    path("visualization_PCA_table.csv"), optional: true // for two channel processing, sample wise is not generated
  
  script:
    """
    #! /usr/bin/env Rscript
    #2

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
    # replace with .rmd files in channels
    file.copy(file.path(codebase_dir,"qa_summary_normalized.rmd"), ".", overwrite=TRUE)
 
    
    rmarkdown::render("qa_summary_normalized.rmd",
                      "html_document", 
                      output_file="normalized_qa",
                      output_dir=getwd()) # generates PCA_raw
    """
}


process QA_RAW {
  publishDir "${ params.outputDir }/${ params.gldsAccession }/00-RawData",
    mode: params.publish_dir_mode

  input:
    path("raw_data.RData")
    path("targets.RData")
  
  output:
    path("raw_qa.html")
    path("raw_visualization_PCA_table.csv")
  
  script:
    """
    #! /usr/bin/env Rscript
    # 4i 

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
    # replace with .rmd files in channels
    file.copy(file.path(codebase_dir,"qa_summary_raw.rmd"), ".", overwrite=TRUE)
 
    rmarkdown::render("qa_summary_raw.rmd",
                      "html_document", 
                      output_file="raw_qa",
                      output_dir=getwd()) # generates PCA_raw

    """
}

process NORMALIZE {
  publishDir "${ params.outputDir }/${ params.gldsAccession }/01-NormalizedData",
    mode: params.publish_dir_mode,
    pattern: "normalize*.{txt,html}"

  input:
    path("raw_data.RData")
    path("targets.RData")
  
  output:
    path("normalized.RData"), emit: rdata
    path("normalized.txt")
    path("normalize.html")
  
  script:
    """
    #! /usr/bin/env Rscript
    #3

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
    # replace with .rmd files in channels
    file.copy(file.path(codebase_dir,"normalize.rmd"), ".", overwrite=TRUE)
 
    rmarkdown::render("normalize.rmd",
                      "html_document", 
                      output_file="normalize",
                      output_dir=getwd())
 
    """
}

process READ_RAW {
  publishDir "${ params.outputDir }/${ params.gldsAccession }/00-RawData",
    mode: params.publish_dir_mode,
    pattern: "load_raw.html"

  input:
    path(raw_files)
    path(runsheet_RData)
  
  output:
    path("raw_data.Rdata"), emit: rdata
    path("load_raw.html")
  
  script:
    """
    #! /usr/bin/env Rscript
    # 3i 

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
 
    load("${ runsheet_RData }") # named 'targets'
  
    # replace with .rmd files in channels
    file.copy(file.path(codebase_dir,"load_raw.rmd"), ".", overwrite=TRUE)
 
    rmarkdown::render("load_raw.rmd",
                      "html_document", 
                      output_file="load_raw",
                      output_dir=getwd())
    """
}


process IMPORT_PROBE {
  publishDir "${ params.outputDir }/${ params.gldsAccession }/01-NormalizedData",
    mode: params.publish_dir_mode

  input:
    path(raw_data_RData)
    path(runsheet_RData)
    path(normalized_RData)
  
  output:
    path("normalized-annotated.txt")
    path("probe_annotations.txt")
    path("normalized-annotated.rda"), emit: annotated_rdata
    path("annotations.html")
  
  script:
    """
    #! /usr/bin/env Rscript
    #12a

    codebase_dir <- file.path("${ workflow.projectDir }","bin")

    # replace with .rmd files in channels
    file.copy(file.path(codebase_dir,"annotations.rmd"), ".", overwrite=TRUE)
 
    rmarkdown::render("annotations.rmd",
                      "html_document", 
                      output_file="annotations",
                      output_dir=getwd())

    """
}

process DGE {
  publishDir "${ params.outputDir }/${ params.gldsAccession }/02-DGE",
    mode: params.publish_dir_mode

  input:
    path(normalized_annotated_RData)
    path(runsheet_RData)
  
  output:
   path("contrasts.csv")
   path("visualization_output_table.csv")
   path("differential_expression.csv")  
   path("dge.html")  

  script:
    """
    #! /usr/bin/env Rscript
    #6b

    codebase_dir <- file.path("${ workflow.projectDir }","bin")
    
    # replace with .rmd files in channels
    file.copy(file.path(codebase_dir,"dge.rmd"), ".", overwrite=TRUE)
 
    rmarkdown::render("dge.rmd",
                      "html_document", 
                      output_file="dge",
                      output_dir=getwd())

    """
}
