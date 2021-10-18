nextflow.enable.dsl=2


process RUN {
  conda = "${projectDir}/envs/minimal.yml"
  publishDir = "${ params.outputDir }/${ params.gldsAccession }"

  input:
    path(runsheet)
    val(gldsAccession)

  output:
    path("Processed_Data/*")

  script:
    """
	glds_microarrays.R --glds ${gldsAccession} --reports --runsheet ${runsheet}
    """

}


/* Processes dealing with retrieving data from GeneLab
*/
process MICROARRAY_RUNSHEET_FROM_GLDS {
  // Downloads isazip and creates run sheets using GeneLab API
  tag "${ glds_accession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/Metadata",
    mode: params.publish_dir_mode

  input:
    val(glds_accession)
    val(isa_archive)

  output:
    path("AST_autogen_*_${ glds_accession }_*.csv"), emit: runsheet

  script:
    if ( !isa_archive )
    """
    retrieve_isa_from_genelab.py --accession ${ glds_accession }\
                                 --alternate-url \
                                 --to-Microarray-runsheet
    """
    else
    """
    retrieve_isa_from_genelab.py --accession ${ glds_accession }\
                                 --local-isa-zip ${ isa_archive }\
                                 --to-Microarray-runsheet 
    """
}


workflow {
  main:
    MICROARRAY_RUNSHEET_FROM_GLDS( params.gldsAccession, params.isa_archive )
    RUN( MICROARRAY_RUNSHEET_FROM_GLDS.out.runsheet, params.gldsAccession )
}
