
Channel
          .fromPath(params.inputTable)
          .splitCsv(sep: "\t")
          .map { row -> [ row[0], row[1]]}
          .set { block_input}


mapFile = file("data/mapLDref.tsv.gz")


process makeLDfiles {

  container 'yocra3/hail_s3a:1.0'

  memory '32 GB'
  cpus 12

  input:
  tuple val(chr), val(block) from block_input
  file(mapFile)

  output:
  file("*.csv.bgz") into blockLD

  script:
  """
  makeLDfiles.py $chr $block
  """
}

process makeRDSfiles {

  container 'yocra3/rsession_chd_marato:release-1.2.5'

  publishDir "results/rds_files/", mode: 'copy'

  memory '30 GB'
  cpus 1

  input:
  file(mapFile) from blockLD

  output:
  file("*.rds") into blockRDS

  script:
  """
  generateLD_RDSfiles.R $mapFile
  """
}
