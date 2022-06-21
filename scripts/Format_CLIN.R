args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")

clin = read.csv( file.path(input_dir, "CLIN.txt") , stringsAsFactors=FALSE , sep="\t")

clin = cbind( clin[ , c( "patient_id" , "Sex" ,  "Age" , "Best.Response.RECIST.1.1" , "is_deceased" , "os" , "pfs" , "progressed" ) ] , "Ureteral" , "PD-1/PD-L1" , NA , NA, NA , NA , NA, NA )
colnames(clin) = c( "patient" , "sex" , "age"  ,"recist"  , "os" , "t.os" , "t.pfs" , "pfs" , "primary" , "drug_type" , "response.other.info" , "response" , "histo" , "stage" , "dna" , "rna")

clin$patient = paste( "P" , clin$patient , sep="" )
clin$rna = "tpm"
clin$dna = "wes"

clin$t.pfs = as.numeric( as.character( clin$t.pfs ) ) / 30.5
clin$t.os = as.numeric( as.character( clin$t.os ) ) / 30.5

clin$pfs = ifelse( as.character( clin$pfs ) %in% "True" , 1 , 0 )
clin$os = ifelse( as.character( clin$os ) %in% "True" , 1 , 0 )

clin$response = Get_Response( data=clin )		

clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

