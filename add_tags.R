## Edit & Uncomment - set working directory for your machine
#setwd('zf_cds')

## Uncomment next 2 lines if you need to install biomaRt package
# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")

## Uncomment if you need to install data.table or stringr packages
# install.packages("data.table")
# install.packages("stringr")

library('data.table')
library('biomaRt')
library('stringr')

#read in tab-separated file (Comment lines removed: grep -v "^###$" filename.gff3 > filename.tsv)
## Uncomment and edit filename
#infile <- fread('filename.tsv', stringsAsFactors=F)

#adds column for nonversioned accession
infile[,'ID' := str_match(V9, "ID=([:alnum:]+)")[,2]]
#This adds to every row (unnecessary in our case)
#infile[, 'ID' := str_match(V9, "(?:ID|Parent)=([:alnum:]+)")[,2]]

#Gets list of accessions for biomaRt to find
list <- infile$ID[!is.na(infile$ID)]

#Uses the Danio mart from ENSEMBL
mart = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
danioMart<- useDataset("drerio_gene_ensembl",mart=useMart("ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org"))

#Gets attributes using transcript id from biomaRt
attr<- getBM(filters="ensembl_transcript_id", attributes=c("ensembl_transcript_id", "external_gene_name", "description"), values=list, mart=danioMart)


#Function parses Description string extra info into valid GFF tags
processDescStr <- function(str) {
  rexp <- "(.*?) \\[(.*?)\\]"
  parts <- str_match(str, rexp)[-1]
  parts[1]<-paste('Description', parts[1], sep="=")
  parts[1]<-gsub(",", "", parts[1])
  parts[2]<-gsub(':', '=', parts[2])
  return(paste(parts, sep="", collapse=";"))
}

#Creates new columns for descriptions and gene symbols
attr$description <- sapply(attr$description, processDescStr)
attr$external_gene_name <- paste("Symbol", attr$external_gene_name, sep="=")

#Merges new information with GFF table
#Important: IDs in infile are only on the line of the parent feature for each entry.
dt <- merge(infile,attr, by.x="ID", by.y="ensembl_transcript_id", all.x=T, sort=F)

###
#http://stackoverflow.com/questions/13673894/suppress-nas-in-paste
paste3 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}

dt$V9 <- paste3(dt$V9, dt$external_gene_name, dt$description, sep=";")
dt$ID <- NULL
dt$external_gene_name <- NULL
dt$description <- NULL

## Uncomment and edit filename
#write.table(dt, file="filename.gff3", quote=F, sep="\t", row.names=F, col.names=F)

#Just added dividers the old fashoned way in bash:
#sed -i '/\tgene\t/i\###' filename.gff3

