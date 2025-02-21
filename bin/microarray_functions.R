
headers<-list(
  ae1=c("metaColumn","metaRow","row","column"),
  genepix=c("Block","Column","Row","X","Y"),
  arrayvision=c("Primary","Secondary"),
  agilent=c("Row","Col","PositionX","PositionY"),
  scanalyze=c("GRID","COL","ROW","LEFT","TOP","RIGHT","BOT"),
  scanarray=c('Array Column','Array Row','Spot Column','Spot Row','X','Y'),
  quantarray=c('Array Column', 'Array Row', 'Column', 'Row'),
  spotfinder=c("MC","MR","SC","SR","C","R"),
  mev=c("MC","MR","C","R","UID"),
  codelink=c("Logical_row","Logical_col","Center_X","Center_Y"),
  bluefuse=c("COL","ROW","SUBGRIDCOL","SUBGRIDROW"),
  UCSFSpot=c("Arr-colx","Arr-rowy","Spot-colx","Spot-rowy"),
  NimbleScanFeature=c("X","Y","PROBE_ID","X_PIXEL","Y_PIXEL"),
  NimblegenNASA=c("X_BC","Y_BC","Feature_ID","ProbID_BC"),
  imagene=c('Meta Column', 'Meta Row', 'Column', 'Row', 'Field', 'Gene ID'),
  ImaGene3=c("Meta_col","Meta_row","Sub_col","Sub_row","Name","Selected"),
  ImaGene7=c("Block","Column","Row","Ch1 XCoord","Ch1 YCoord", "Ch2 XCoord", "Ch2 YCoord"),
  ImaGeneFields=c("Field","Column","Row","XCoord","YCoord"),
  CSIRO_Spot=c("grid_c","grid_r","spot_c","spot_r","indexs"))

### Use API staging script
library(stringr)
stagingTargets <- function(runsheet_path) {

    na_strings <- c('NA','na','null','NULL','Null','')
    target<-list()
    table <- read.csv(runsheet_path,header = TRUE, stringsAsFactors = FALSE)
    study_samples <- table$sample_name
    target$channels = length(unique(table$Label))
    labelnames <- unique(table$Label)
    files <- table$array_data_file
    target$raw_files <- str_remove(files,".gz$")
    filetypes <- unique(toupper(tools::file_ext(target$raw_files)))
    if (length(filetypes)>1){
      cat("\nFiletypes\n",filetypes)
      stop("Execution halted due to multiple filetypes present for raw data")
    }
    target$filetype <- filetypes
    target$seperate_channel_files <- eval((length(unique(files)) == 2*length(unique(table$Hybridization.Assay.Name))))
    factorcols <- sapply(colnames(table),function(x) {stringr::str_detect(x,"Factor")})
   
    factors <- as.data.frame(table[,factorcols==TRUE])
    cat("Parsed factor headers: ",colnames(factors),"\n")
    colnames(factors)<-paste("factor",1:dim(factors)[2], sep = "_")
    is.na(factors) <- Reduce("|", lapply(na_strings, "==", factors)) # avoids whitespace issues for group strings
    group <- apply(factors,1,function(x) {paste(na.omit(x),collapse=" & ")}) # groups are concatenations of non empty factor values per sample
    table$Hybridization.Assay.Name <- sub("\\ \\d+", "", table$Hybridization.Assay.Name)
    table$Source.Name <- sub("\\ \\d+", "", table$Source.Name)
    print(colnames(table))
    
    if ((target$channels == 1) && (target$seperate_channel_files == FALSE)){
      cat("One color microarray design...","\n")
      target$t1 <- data.frame(SampleName=table$sample_name, Group=group, ArrayName=table$Hybridization.Assay.Name, FileName=files)
      target$t1$SourceName <- table$Source.Name
      target$t1 <- cbind(target$t1,factors)
      target$paired_samples <- (length(unique(target$t1$SourceName)) < length(unique(target$t1$SampleName)))
    }else if ((target$channels == 2) && (target$seperate_channel_files == FALSE)){
      cat("Two color microarray design...","\n")
      target$t1 <- data.frame(Label=table$Label, Group=group, ArrayName=table$Hybridization.Assay.Name, FileName=files)
      
      target$t1 <- tidyr::pivot_wider(target$t1, names_from = Label, values_from = Group)
    }else if ((target$channels == 1) && (target$seperate_channel_files == TRUE)){
      cat("One color microarray design with seperate channel files...","\n")
      target$t1 <- data.frame(SampleName=table$sample_name, Group=group, ArrayName=table$Hybridization.Assay.Name, FileName=files)      
      target$t1$SourceName <- table$Source.Name
      target$t1 <- target$t1[order(target$t1$SampleName,target$t1$FileName),]
      target$t1$Label <- rep(c("FileName_1","FileName_2"),dim(table)[1]/2)
      target$t1 <- tidyr::pivot_wider(target$t1, names_from = Label, values_from = FileName)
      target$t1 <- cbind(target$t1,factors)
      target$paired_samples <- (length(unique(target$t1$SourceName)) < length(unique(target$t1$SampleName)))
    }else if ((target$channels == 2) && (target$seperate_channel_files == TRUE)){
      cat("Two color microarray design with seperate channel files...","\n")
      target$t1 <- data.frame(Label=table$Label, Group=group, ArrayName=table$Hybridization.Assay.Name, FileName=files)
      target$t1 <- target$t1[order(target$t1$FileName,target$t1$Label),]
      target$t1$FileLabel <- paste0("FileName",target$t1$Label)
      t1a <- tidyr::pivot_wider(target$t1[,c("Label","Group","ArrayName")], names_from = Label, values_from = Group)
      t1b <- tidyr::pivot_wider(target$t1[,c("FileLabel","FileName","ArrayName")], names_from = FileLabel, values_from = FileName)
      target$t1 <- dplyr::left_join(t1a,t1b,by = "ArrayName")
    }else {
      cat("Error detecting labeling pattern", "\n")
    }
    
    target$technical_replicates <- (length(unique(target$t1$ArrayName)) < dim(target$t1)[1])
    
    if (target$channels == 2){
      group_pairs <- target$t1[,labelnames]
      Ref <- as.data.frame(t(apply(group_pairs,1,function(x){unique(group) %in% x})))
      Ref <- apply(Ref,2,all)
      Ref_group <- unique(group)[Ref] # groups common to all arrays
      sample_pairs <- table[,which(colnames(table) %in% c("Label","sample_name", "Hybridization.Assay.Name"))]
      sample_pairs <- tidyr::pivot_wider(sample_pairs, names_from = Label, values_from = `sample_name`, id_cols = `Hybridization.Assay.Name`)
      sample_pairs <- sample_pairs[,-c(1)]
      common_samples <- as.data.frame(t(apply(sample_pairs,1,function(x){unique(table$sample_name) %in% x})))
      common_samples <- unique(table$sample_name)[which(apply(common_samples,2,all))]
      
      contrasts <- unique(rbind(group_pairs,group_pairs[,c(2,1)],deparse.level = 0)) # available contrast group pairs
      contrast_space<-t(combn(unique(group),2)) # set of all contrast group pairs
      
      if(length(Ref_group) == 2){
        cat("Replicate array design ...", "\n")
        target$design <- "Replicate Array"
      }else if (length(common_samples) == 1){
        cat("Commnon reference design ...", "\n")
        target$design <- "Common Reference"
        replace_group <- target$t1[,labelnames]
        replace_group[sample_pairs == common_samples] <- "Ref"
        target$t1[,labelnames] <- replace_group
      }else if ((length(Ref_group) <= 1) && (all(apply(contrast_space,1,function(x){x %in% contrasts})))){
        cat("Direct two color design ...", "\n")
        target$design <- "Direct Two-Color"
      }else if ((length(Ref_group) <= 1) && (!all(apply(contrast_space,1,function(x){x %in% contrasts})))){
        cat("Seperate channel unconnected design ...", "\n")
        target$design <- "Separate Channels"
      }
    }
    
    label_cols <- which(colnames(target$t1) %in% labelnames)
    target$t2 <- target$t1
    try({target$t2$Group <- sapply(target$t1$Group,function(x){paste0("(",x,")")})},silent = TRUE)
    try({target$t2[,label_cols[1]] <- sapply(target$t1[,label_cols[1]],function(x){paste0("(",x,")")})},silent = TRUE)
    try({target$t2[,label_cols[2]] <- sapply(target$t1[,label_cols[2]],function(x){paste0("(",x,")")})},silent = TRUE)
    
    target$t3 <- target$t1
    try({target$t3$Group <- sapply(target$t1$Group,function(x){make.names(x,unique = FALSE, allow_ = TRUE)})},silent = TRUE) #Make syntactically valid group names
    try({target$t3[,label_cols[1]] <- sapply(target$t1[,label_cols[1]],function(x){make.names(x,unique = FALSE, allow_ = TRUE)})},silent = TRUE)
    try({target$t3[,label_cols[2]] <- sapply(target$t1[,label_cols[2]],function(x){make.names(x,unique = FALSE, allow_ = TRUE)})},silent = TRUE)
    
  
  return(target)

}

