/* Processes dealing with retrieving data from GeneLab
*/
process DOWNLOAD_ISA {
  // Downloads ISA archive from genelab API
  publishDir "${ params.outputDir }/${ params.gldsAccession }/Metadata",
    mode: params.publish_dir_mode

  input:
    val(glds_accession)

  output:
    path("*.zip"), emit: isazip

  script:
    """
    dpt-get-isa-archive --accession ${ glds_accession }\
      --alternate-url
    """
}

process RUNSHEET_FROM_GLDS {
  // Creates runsheet from ISA archive
  publishDir "${ params.outputDir }/${ params.gldsAccession }/Metadata",
    mode: params.publish_dir_mode

  input:
    tuple val(glds_accession), val(config_type), val(config_version)
    path(isazip)

  output:
    // e.g. 'GLDS-205_microarray_v0_runsheet.csv'
    path("${ glds_accession }_${ config_type }_v${ config_version }_runsheet.csv"), emit: runsheet

  script:
    """
    dpt-isa-to-runsheet --accession ${ glds_accession } \
      --isa-archive ${ isazip } --config-type ${ config_type } \
      --config-version ${ config_version }
    """
}


process RUN {
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
  publishDir "${ params.outputDir }/${ params.gldsAccession }/00-RawData",
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
                     "bacillus_subtilis":"BACSU",
                     "arabidopsis_thaliana":"ARABIDOPSIS"]

    def meta = [:]
    meta.id                         = row["Sample Name"]
    meta.organism_sci               = row.organism.replaceAll(" ","_").toLowerCase()
    meta.organism_non_sci           = ORGANISMS[meta.organism_sci]
    meta.raw_file_comment           = row["Comment[Array Data File Name]"]
    meta.raw_file                   = row["Array Data File Name"]
    meta.raw_file_path              = row["Array Data File Path"]
    meta.platform                   = row["Study Assay Technology Platform"]
    
    return meta
}
