options(timeout=300)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

download.file(
  url = 'https://github.com/hammerlab/multi-omic-urothelial-anti-pdl1/raw/master/data_clinical.csv', 
  destfile = file.path(work_dir, 'data_clinical.csv')
)
download.file(
  url = 'https://github.com/hammerlab/multi-omic-urothelial-anti-pdl1/raw/master/data/2850417_Neoantigen_RNA_bams.csv',
  destfile = file.path(work_dir, '2850417_Neoantigen_RNA_bams.csv')
)
download.file(
  url = 'https://zenodo.org/record/546110/files/kallisto.zip?download=1',
  destfile = file.path(work_dir, 'kallisto.zip')
)