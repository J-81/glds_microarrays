/*
* Workflow that accepts a GLDS accession and generates the following:
* 1. Download ISA.zip and generates Microarray Samplesheet
* 2. Downloads Raw Files
* 3. For legacy datasets, denest data-set wide zipped raw files
*/

// This ensures DSL2 syntax and process imports
nextflow.enable.dsl=2


// Import process from separate module file
include { RUNSHEET_FROM_GLDS } from "./modules/genelab.nf" addParams(runsheet: "microarray")
include { STAGE_RAW_FILES;
          get_runsheet_paths } from'./modules/genelab.nf'

/**************************************************
* ACTUAL WORKFLOW  ********************************
**************************************************/
workflow staging{
  take:
    ch_glds_accession
  main:

    isa_ch = Channel.fromPath( params.isaArchive )
    RUNSHEET_FROM_GLDS( params.gldsAccession, isa_ch )

    RUNSHEET_FROM_GLDS.out.runsheet | splitCsv(header: true)
                                    | map{ row -> get_runsheet_paths(row) }
                                    | set{ ch_samples }
    ch_samples | first | view
    STAGE_RAW_FILES( ch_samples )

    emit:
      raw_files = params.stageLocal ? STAGE_RAW_FILES.out : null
      runsheet = RUNSHEET_FROM_GLDS.out.runsheet
}
