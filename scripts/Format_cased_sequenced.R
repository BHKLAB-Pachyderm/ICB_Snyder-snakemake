args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

clin = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" )

case = as.data.frame( cbind( 
  clin$patient_id , 
  rep( 0 , length(clin$patient_id ) ) , 
  rep( 0 , length(clin$patient_id ) ) , 
  rep( 0 , length(clin$patient_id ) ) 
) )
colnames(case) = c( "patient" , "snv" , "cna" , "expr" )
rownames(case) = clin$patient_id

expr <- read.csv( file.path(input_dir, "EXPR.txt") , stringsAsFactors=FALSE , sep="\t" )
case$expr <- lapply(case$patient, function(p){
  if(p %in% colnames(expr)){
    return(1)
  }else{
    return(0)
  }
})

case$snv = as.numeric( as.character( case$snv ) )
case$cna = as.numeric( as.character( case$cna ) )
case$expr = as.numeric( as.character( case$expr ) )

write.table( case , file=file.path(output_dir, "cased_sequenced.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
