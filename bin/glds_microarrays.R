#!/usr/bin/env Rscript


### Get path of launched main script
args <- commandArgs(trailingOnly = FALSE)
codebase_dir <- dirname(sub("--file=","",args[grep("--file=", args)]))

### Enable better error traceback messaging
### Disable for production as this will cause R to not return an error code back to the interpreter
#options(error = quote({
#  dump.frames(to.file=T, dumpto='last.dump')
#  load('last.dump.rda')
#  print(last.dump)
#  q()
#}))



library("optparse")

option_list = list(
  make_option(c("-f", "--files"), type="character", default=NULL, 
              help="Dataset file paths", metavar="character"),
  make_option(c("-u", "--runsheet"), type="character", default=NULL, 
              help="Runsheet file path", metavar="character"),
  make_option(c("-r", "--reports"), action="store_true", default=FALSE, 
              help="Generate QA HTML Reports"),
  make_option(c("-i", "--isa"), type="character", default=NULL, 
              help="Study ISAtab.zip file path", metavar="character"),
  make_option(c("-t", "--staging"), type="character", default=NULL, 
              help="Use API staging CSV file"),
  make_option(c("-p", "--probe"), type="character", default=NULL, 
              help="Probe annotation file path", metavar="character"),
  make_option(c("-s", "--species"), type="character", default=NULL, 
              help="Specify organism species (e.g. Arabidopsis thaliana)", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL, 
              help="Organism annotation file path", metavar="character"),
  make_option(c("-g", "--glds"), type="character", default=NULL, 
              help="GeneLab study identifier (e.g. GLDS-121)", metavar="character"),
  make_option(c("-m", "--platform"), type="character", default=NULL, 
              help="Microarray Platform", metavar="character"), 
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Curated input directory path", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=getwd(), 
              help="output parent directory path [default= %default]", metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #####Testing
# opt<-list()
# opt$files <- "demo_datasets/GLDS-22/00-RawData/*.txt"
# opt$isa <- "demo_datasets/GLDS-22/Metadata/*ISA.zip"
# opt$glds <- "GLDS-22"
# opt$species <- "Arabidopsis thaliana"
# opt$platform <- "Agilent 1-channel"
# opt$out <- "~/Documents/glds_projects/glds_microarrays"
# setwd("~/Documents/glds_projects/glds_microarrays")
# #####

if (is.null(opt$glds)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (GLDS)", call.=FALSE)
}

source(file.path(codebase_dir, "microarray_functions.R"))

### Get Study Files
if (length(opt$files) > 0){
  cat("\nImporting files for:",opt$glds, "\n")
  opt$files <- Sys.glob(file.path(opt$files))
}

if (length(opt$runsheet >= 1)){
  tempstage <- tempdir()
  unlink(list.files(tempstage, full.names = TRUE))
  cat("Staging table path: ",opt$runsheet,"\n")
  table <- read.csv(opt$runsheet,header = TRUE, stringsAsFactors = FALSE)
  cat("Staging headers: ",colnames(table),"\n")
  opt$species <- table$Characteristics.Organism.[1]
  cat("\nParsed organism: ",opt$species,"\n")
  opt$platform <- table$Study.Assay.Technology.Platform[1]
  labels <- table$Label
  if(length(unique(labels))==1 && opt$platform == 'Agilent'){
    opt$platform <- "Agilent 1-channel"
  }else if(length(unique(labels))==2 && opt$platform == 'Agilent'){
    opt$platform <- "Agilent 2-channel"
  }
  cat("\nParsed platform: ",opt$platform,"\n")
  cat("\nStaging file list: ",table$array_data_file_path,"\n")

  opt$files <- table$array_data_file_path
  
 # for (file in 1:length(opt$files)){ # test whether file path exists from runsheet, if not, tries to find it using local directory structure
 #   if (file.exists(opt$files[file])){
 #     cat("The file exists: ",opt$files[file],"\n")
 #   } else{
 #     altfile <- paste0(dirname(dirname(opt$runsheet)),"/00-RawData/",basename(opt$files[file]))
 #     if (file.exists(altfile)) {
 #       
 #       cat("The file exists: ",altfile,"\n")
 #       opt$files[file] <- altfile
 #     } else {
 #       
 #       cat("The files do not exist: ",altfile,"  ",opt$files[file],"\n")
 #     }
 #   }
 # }
 # 
 # cat("\nExtracted file list: ",opt$files,"\n")
}

if (length(opt$staging >= 1)){
  tempstage <- tempdir()
  unlink(list.files(tempstage, full.names = TRUE))
  cat("Staging table path: ",opt$staging,"\n")
  table <- read.csv(opt$staging,header = TRUE, stringsAsFactors = FALSE)
  cat("Staging headers: ",colnames(table),"\n")
  opt$species <- table$Characteristics.Organism.[1]
  cat("\nParsed organism: ",opt$species,"\n")
  opt$platform <- table$Study.Assay.Technology.Platform[1]
  labels <- table$Label
  if(length(unique(labels))==1 && opt$platform == 'Agilent'){
    opt$platform <- "Agilent 1-channel"
  }else if(length(unique(labels))==2 && opt$platform == 'Agilent'){
    opt$platform <- "Agilent 2-channel"
  }
  cat("\nParsed platform: ",opt$platform,"\n")
  cat("\nStaging file list: ",table$array_data_file_path,"\n")
  
  for (file in 1:length(table)){
    utils::download.file(table$array_data_file_path[file], destfile = file.path(tempstage,table$array_data_file[file]),quiet = FALSE)
  }

  opt$files <- list.files(tempstage,full.names = TRUE)
  cat("\nExtracted file list: ",opt$files,"\n")
  
}
### Get organism annotation package
options(connectionObserver = NULL)
organism_table <- read.csv(file = file.path(codebase_dir,"organisms.csv"), header = TRUE, stringsAsFactors = FALSE)
ann.dbi <- organism_table$annotations[organism_table$species == opt$species] # Organism specific gene annotation database
ann.dbi=as.character(ann.dbi)
if(!require(ann.dbi, character.only=TRUE)) {
  BiocManager::install(ann.dbi, ask = FALSE)
  library(ann.dbi, character.only=TRUE)
}
#require(ann.dbi, character.only = TRUE, quietly = TRUE)
cat("\nOrganism annotation set loaded: ",ann.dbi,"\n")




### Get Files from Curated Directory
if (!is.null(opt$dir)){
  cat("Building from local curated directory\n")
  opt$isa <- list.files(file.path(getwd(),opt$dir,"Metadata"),pattern = "*ISA.zip", full.names = TRUE)
  opt$probe <- list.files(file.path(getwd(),opt$dir,"Metadata"),pattern = "*annotation*", full.names = TRUE)
  opt$files <- list.files(file.path(getwd(),opt$dir,"00-RawData"),full.names = TRUE)
  str(opt)
}

### Copy files to temp directory
tempin <- tempdir()
#unlink(list.files(tempin, full.names = TRUE))
dir.create(file.path(tempin,"00-RawData"), showWarnings = FALSE)
cat("files: ",opt$files)
file.copy(from = opt$files, to = file.path(tempin,"00-RawData"), overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
opt$files <- list.files(file.path(tempin,"00-RawData"),full.names = TRUE)
cat("files: ",opt$files)
dir.create(file.path(tempin,"Metadata"),showWarnings = FALSE)
file.copy(from = opt$isa, to = file.path(tempin,"Metadata"), overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
opt$isa <- list.files(file.path(tempin,"Metadata"),pattern = "*ISA.zip", full.names = TRUE)
file.copy(from = opt$probe, to = file.path(tempin,"Metadata"), overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
opt$probe <- list.files(file.path(tempin,"Metadata"),pattern = "*annotation*|GPL*", full.names = TRUE)
cat("Temp Probe file: ",opt$probe, "\n")
### Determine Platform

### Selects appropriate package for plots

if ((opt$platform == "Affymetrix") | (opt$platform == "Nimblegen 1-channel")){
  plots <- "OLIGO"
}else if((opt$platform == "Illumina Expression") | (opt$platform == "Agilent 1-channel") | (opt$platform == "Agilent 2-channel") | (opt$platform == "Nimblegen 2-channel")){
  plots <- "LIMMA"
}else {
  plots <- NULL
}
cat("\nPlots generate by package:",plots,"\n")

### Parse ISA


# all_targets<-buildTargets(opt)
# if (length(all_targets) == 1){
#   targets <- all_targets[[1]]
# }else{
#   targets <- all_targets[[which_assay]]
# }
targets <- stagingTargets(opt)
cat("\nTargets\n")
str(targets)



### Call appropriate platform processing script

if(opt$platform=="Affymetrix"){
  source(file.path(codebase_dir,"platforms","affymetrix.R"))
}else if(opt$platform=="Agilent 1-channel"){
  source(file.path(codebase_dir,"platforms","agilent_one_channel.R"))
}else if(opt$platform=="Agilent 2-channel"){
  source(file.path(codebase_dir,"platforms","agilent_two_channel.R"))
}else{
  print_help(opt_parser)
  stop("Platform not currently supported.n", call.=FALSE)
}
