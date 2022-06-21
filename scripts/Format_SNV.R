library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

snv = read.csv( file.path(input_dir, "SNV.txt") , stringsAsFactors=FALSE , sep="\t" )

data = snv[ , c( "patient" , "Chr" , "Pos" , "Ref" , "Alt" ,  "Gene" , "Effect" ) ] 
colnames(data) = c( "Sample" , "Chr" , "Pos" , "Ref" , "Alt" ,  "Gene" , "Effect" )


data$Ref = ifelse( data$Ref %in% "-" , "" , data$Ref )
data$Alt = ifelse( data$Alt %in% "-" , "" , data$Alt )

data = cbind( data ,
				apply( data[ , c( "Ref", "Alt" ) ] , 1 , function(x){ ifelse( nchar(x[1]) != nchar(x[2]) , "INDEL", "SNV") } )
			)

colnames(data) = c( "Sample" , "Chr" , "Pos" , "Ref" , "Alt" ,  "Gene" , "Effect" , "MutType"  )

data$Chr = paste( "chr" , data$Chr , sep="" ) 

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
data = data[ data$Sample %in% case[ case$snv %in% 1 , ]$patient , c( "Sample" , "Gene" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "MutType" ) ]

data$Effect = ifelse( data$Effect %in% "DE_NOVO_START_IN_FRAME" , "De_novo_Start_InFrame" , data$Effect )                
data$Effect = ifelse( data$Effect %in% "DE_NOVO_START_OUT_FRAME" , "De_novo_Start_OutOfFrame" , data$Effect )                

write.table( data , file=file.path(output_dir, "SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
