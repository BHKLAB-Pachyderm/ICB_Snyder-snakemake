library(data.table)
library(readxl) 
library(stringr)
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

# download.file(
#   url = 'https://github.com/hammerlab/multi-omic-urothelial-anti-pdl1/raw/master/data_effects.csv', 
#   destfile = file.path(work_dir, 'data_effects.csv')
# )
# download.file(
#   url = 'https://github.com/hammerlab/multi-omic-urothelial-anti-pdl1/raw/master/data_variants.csv', 
#   destfile = file.path(work_dir, 'data_variants.csv')
# )


# CLIN.txt
clin <- read.csv(file.path(work_dir, 'data_clinical.csv'), header=FALSE, stringsAsFactors=FALSE, fileEncoding="latin1")
clin[1, ] <- str_replace_all(clin[1, ], '\\W', '.')
clin <- clin[, clin[1, ] != 'NA']
colnames(clin) <- clin[1, ]
clin <- clin[-1, ]
clin <- clin[, c('patient_id', colnames(clin)[colnames(clin) != 'patient_id'])]

# EXPR.txt
counts <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(counts) <- c('patient_id', 'target_id', 'est_counts')
filemap <- read.csv(file.path(work_dir, '2850417_Neoantigen_RNA_bams.csv'), header=FALSE, sep='\t')

for(patient in filemap$V1){
  patient_id <- str_extract(patient, "(?<=-)(\\d+)(?=-)")
  print(patient_id)
  dir <- filemap$V2[filemap$V1 == patient]
  dir <- unlist(str_split(dir, "/"))
  dir <- str_replace(dir[6], '.bam', '')
  count <- read.csv(file.path(work_dir, 'kallisto', str_c(dir, '-kallisto'), 'abundance.tsv'), sep='\t')
  count$patient_id <- patient_id
  counts <- rbind(counts, count[, c('patient_id', 'target_id', 'est_counts')])
}

expr <- data.frame(matrix(ncol = length(unique(counts$patient_id)), nrow = length(unique(counts$target_id))))
colnames(expr) <- unique(counts$patient_id)
rownames(expr) <- unique(counts$target_id)
for(patient in colnames(expr)){
  print(patient)
  df <- counts[counts$patient_id == patient, c('target_id', 'est_counts')]
  expr[patient] <- unlist(lapply(rownames(expr), function(target_id){
    return(df[df$target_id == target_id, 'est_counts'])
  }))
}
saveRDS(expr, file=file.path(work_dir, 'expr.rds'))
load(file.path(work_dir, "Gencode.v40.annotation.RData"))

# SNV.txt
# TO DO Curate SNV data