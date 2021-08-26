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

workflow {
  main:
    println(params)
    RUN(params.runsheet, params.gldsAccession)
}
