/*
* Workflow that accepts a GLDS accession and generates the following:
* 1. Download ISA.zip and generates Microarray Samplesheet
* 2. Downloads Raw Files
* 3. For legacy datasets, denest data-set wide zipped raw files
*/

// This ensures DSL2 syntax and process imports
nextflow.enable.dsl=2


// Import process(es) and function(s) from separate module file
include { DOWNLOAD_ISA;
          RUNSHEET_FROM_GLDS;
          STAGE_RAW_FILES;
          get_runsheet_paths } from "./modules/genelab.nf"

/**************************************************
* ACTUAL WORKFLOW  ********************************
**************************************************/
workflow staging{
  take:
    ch_glds_accession
  main:
    // ISA FROM GLDS REPO
    if( !params.runsheet && !params.isaArchive ) {
      DOWNLOAD_ISA( params.gldsAccession )
      RUNSHEET_FROM_GLDS( tuple(params.gldsAccession, "microarray", params.dp_tools_config_version),
                          DOWNLOAD_ISA.out.isazip ) 
                          | set { ch_runsheet }
    // LOCAL ISA SUPPLIED
    } else if( !params.runsheet ) {
    	RUNSHEET_FROM_GLDS( tuple(params.gldsAccession, "microarray", params.dp_tools_config_version), \
                          params.isaArchive ) 
                          | set { ch_runsheet }
    } else {
    // LOCAL RUNSHEET SUPPLIED
    ch_runsheet = Channel.fromPath( params.runsheet )
    }

    ch_runsheet | splitCsv(header: true)
                | map{ row -> get_runsheet_paths(row) }
                | set{ ch_samples }

    ch_samples | first | view
    STAGE_RAW_FILES( ch_samples )

    emit:
      raw_files = params.stageLocal ? STAGE_RAW_FILES.out : null
      runsheet = ch_runsheet
}
