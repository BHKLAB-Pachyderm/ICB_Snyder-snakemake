library(data.table)
library(readxl) 
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

# CLIN.txt
clin <- read.csv(file.path(work_dir, 'data_clinical.csv'), header=FALSE, stringsAsFactors=FALSE, fileEncoding="latin1")
clin[1, ] <- str_replace_all(clin[1, ], '\\W', '.')
clin <- clin[, clin[1, ] != 'NA']
colnames(clin) <- clin[1, ]
clin <- clin[-1, ]
clin <- clin[, c('patient_id', colnames(clin)[colnames(clin) != 'patient_id'])]
clin$patient_id <- paste0('P', clin$patient_id)
clin <- clin[, colnames(clin)[!str_detect(colnames(clin), 'Unnamed')]]

write.table( clin , file=file.path(work_dir, "CLIN.txt") , quote=FALSE , sep="\t" , col.names=TRUE , row.names=TRUE )

# EXPR.txt
unzip(file.path(work_dir, 'kallisto.zip'), exdir=file.path(work_dir))

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

# expr <- readRDS(file.path(work_dir, 'expr.rds'))

colnames(expr) <- unlist(lapply(colnames(expr), function(patient){
  if(nchar(patient) < 4){
    return(paste0(strrep('0', 4 - nchar(patient)), patient))
  }
  return(patient)
}))

colnames(expr) <- paste0('P', colnames(expr))

write.table( expr , file=file.path(work_dir, "EXPR.txt") , quote=FALSE , sep="\t" , col.names=TRUE , row.names=TRUE )

file.remove(file.path(work_dir, 'data_clinical.csv'))
file.remove(file.path(work_dir, '2850417_Neoantigen_RNA_bams.csv'))
file.remove(file.path(work_dir, 'kallisto.zip'))
unlink(file.path(work_dir, "kallisto"), recursive = TRUE)

# SNV.txt
# TO DO Curate SNV data