nextflow.enable.dsl=2
// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("┇  Microarray Processing Pipeline: $workflow.manifest.version  ┇")
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("usage: nextflow run J-81/glds_microarrays -r $workflow.revision --gldsAccession GLDS-000 [--isaArchive] [--skipVV] [--outputDir] [--stageLocal]") 
  println()
  println("required arguments:")
  println("  --gldsAccession GLDS-000")
  println("                        the GLDS accession id to process through the RNASeq Concensus Pipeline.")
  println("optional arguments:")
  println("  --help                show this help message and exit")
  println("  --isaArchive          supplies an isa archive file instead of retrieving the isa archive from the GeneLab repository")
  println("  --runsheet            supplies the path to a runsheet instead of generating the runsheet based on the GeneLab repository")
  println("  --skipVV              skip automated V&V processes. Default: false")
  println("  --outputDir           directory to save staged raw files and processed files. Default: <launch directory>")
  println("  --stageLocal          download the raw reads files for the supplied GLDS accession id.  Set to false to disable raw read download and processing.  Default: true")
  exit 0
  }

println "PARAMS: $params"

/**************************************************
* CHECK REQUIRED PARAMS AND LOAD  *****************
**************************************************/
// Get all params sourced data into channels
// Set up channel containing glds accession number
if ( params.gldsAccession ) {ch_glds_accession = Channel.from( params.gldsAccession )} else { exit 1, "Missing Required Parameter: gldsAccession. Example for setting on CLI: --gldsAccession GLDS-194"}

if ( !params.outputDir ) {  params.outputDir = "$workflow.launchDir" }

/**************************************************
* DEBUG WARNING  **********************************
**************************************************/
if ( null ) {
  println("${c_back_bright_red}WARNING WARNING: DEBUG OPTIONS ENABLED!")
  println("WARNING WARNING: DEBUG OPTIONS ENABLED!${c_reset}")
} else {
}

/**************************************************
* WORKFLOW SPECIFIC PRINTOUTS  ********************
**************************************************/

include { RUN } from "./modules/genelab.nf"
include { QA_RAW;
          QA_NORMALIZED;
          READ_RAW;
          NORMALIZE;
          LOAD_RUNSHEET } from "./modules/microarray.nf"


include { staging as STAGING } from './stage_analysis.nf'

workflow {
  main:
    STAGING( params.gldsAccession )
    
    STAGING.out.raw_files | map {it[1]} | unique { it.name } | collect | set{ ch_raw_files } 

    STAGING.out.raw_files | map {it[0]} | toSortedList | map {it[0]} |  set { ch_meta }

    ch_raw_files | READ_RAW | NORMALIZE

    LOAD_RUNSHEET( STAGING.out.runsheet, params.gldsAccession, ch_meta )

    QA_RAW( READ_RAW.out, LOAD_RUNSHEET.out )
    QA_NORMALIZED( NORMALIZE.out.rdata, LOAD_RUNSHEET.out )

    // RUN( STAGING.out.runsheet, params.gldsAccession, ch_raw_files )
}
