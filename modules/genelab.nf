/* Processes dealing with retrieving data from GeneLab
*/
process RUNSHEET_FROM_GLDS {
  // Downloads isazip and creates run sheets using GeneLab API
  tag "${ glds_accession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/Metadata",
    mode: params.publish_dir_mode

  input:
    val(glds_accession)
    path(isa_archive)

  output:
    path("AST_autogen_*_${ glds_accession }_*.csv"), emit: runsheet

  script:
    runsheet_arg = ( params.runsheet ==  "microarray" ? "--to-Microarray-runsheet" : "" )
    if ( isa_archive.name == "NO_FILE" )
    """
    retrieve_isa_from_genelab.py --accession ${ glds_accession }\
                                 --alternate-url \
                                 ${ runsheet_arg }
    """
    else
    """
    retrieve_isa_from_genelab.py --accession ${ glds_accession }\
                                 --local-isa-zip ${ isa_archive }\
                                 ${ runsheet_arg }
    """
}


process RUN {
  conda = "${projectDir}/envs/minimal.yml"
  publishDir = "${ params.outputDir }/${ params.gldsAccession }"

  input:
    path(runsheet)
    val(gldsAccession)
    path("raw_files/*")

  output:
    path("Processed_Data/*")

  script:
    """
	glds_microarrays.R --glds ${gldsAccession} --reports --runsheet ${runsheet} --files raw_files
    """

}

process STAGE_RAW_FILES {
  // Stages the raw files into appropriate publish directory
  tag "${ meta.id }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/${ meta.raw_file_root_dir }",
    mode: params.publish_dir_mode

  input:
    val(meta)

  output:
    tuple val(meta), path(resolved_sample_filename)

  script:
      dl_filename = file(meta.raw_file_path).name.split('\\?version')[0]
      resolved_sample_filename = meta.raw_file_comment != "N/A" ?  meta.raw_file_comment : meta.raw_file.replaceAll(".gz","")
      if ( meta.raw_file_path.startsWith("https") )
      """
      echo ${meta}

      wget -O ${dl_filename}  ${meta.raw_file_path}

      # unzip (only applies to legacy dataset wide zipped files)
      unzip ${dl_filename}  || true

      # decompress if gzipped
      gunzip *.gz || true
      """
      else
      """
      cp ${meta.raw_file_path} .  

      # decompress if gzipped
      gunzip *.gz || true
      """
}


def get_runsheet_paths(LinkedHashMap row) {
    def ORGANISMS = ["mus_musculus":"MOUSE",
                     "danio_rerio":"ZEBRAFISH",
                     "rattus_norvegicus":"RAT",
                     "homo_sapiens":"HUMAN",
                     "drosophila_melanogaster":"FLY",
                     "caenorhabditis_elegans":"WORM",
                     "arabidopsis_thaliana":"ARABIDOPSIS"]

    def meta = [:]
    meta.id                         = row.sample_name
    meta.organism_sci               = row.organism.replaceAll(" ","_").toLowerCase()
    meta.organism_non_sci           = ORGANISMS[meta.organism_sci]
    meta.raw_file_comment           = row["Comment[Array Data File Name]"]
    meta.raw_file                   = row["array_data_file"]
    meta.raw_file_path              = row.array_data_file_path
    meta.platform                   = row["Study Assay Technology Platform"]
    meta.channels                   = row.Channels

    meta.raw_file_root_dir          = row.raw_data
    
    return meta
}
