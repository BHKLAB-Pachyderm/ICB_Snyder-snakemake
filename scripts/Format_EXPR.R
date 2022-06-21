library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

exp = read.csv( file.path(input_dir, "EXPR.txt") , stringsAsFactors=FALSE , sep="\t" )
exp$patient_id = paste( "P" , exp$patient_id , sep="" )

gene = sort( unique(  exp$gene_name ) )
patient = sort( unique(  exp$patient_id ) )

expr =  matrix( nrow=length(gene) , ncol= length(patient) )
colnames(expr) = patient
rownames(expr) = gene

for( i in 1:length(patient) ){
	p = exp[ exp$patient_id %in% patient[i] , ]
	expr[ exp$gene_name , patient[i] ] = p$est_counts
}


##################################
## Get Gene Length

genes <- rownames(expr)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="hgnc_symbol",values=genes, mart=human)

size=gene_coords$end_position - gene_coords$start_position
names(size) = gene_coords[,"hgnc_symbol"]

gene = intersect( names(size) , rownames(expr) )

size = size[ gene ]
expr = expr[ gene , ]

##################################
##################################
## Compute TPM data

GetTPM <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

tpm = log2( GetTPM(expr,size) + 1 )


##################################
##################################

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
tpm = tpm[ , colnames(tpm) %in% case[ case$expr %in% 1 , ]$patient ]

write.table( tpm , file=file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
