#!/usr/bin/env Rscript

#library("annotate")

get_array2gid <- function (array_name, tax_id){
    cat("processing:", array_name, "\n");
    biocLite(paste0(array_name, ".db"))
    library(paste0(array_name, ".db"), character.only=TRUE)
    if (array_name == "hgu95av2")
        array2gid_map <- hgu95av2ENTREZID
    else if (array_name == "hgu95a")
        array2gid_map <- hgu95aENTREZID
    else if (array_name == "hgu133a")
        array2gid_map <- hgu133aENTREZID
    else if (array_name == "hgu133a2")
        array2gid_map <- hgu133a2ENTREZID
    else if (array_name == "mgu74a")
        array2gid_map <- mgu74aENTREZID        
    else if (array_name == "mgu74av2")
        array2gid_map <- mgu74av2ENTREZID                
    else if (array_name == "mgu74b")
        array2gid_map <- mgu74bENTREZID        
    else if (array_name == "mgu74bv2")
        array2gid_map <- mgu74bv2ENTREZID                        
    else if (array_name == "mgu74c")
        array2gid_map <- mgu74cENTREZID        
    else if (array_name == "mgu74cv2")
        array2gid_map <- mgu74cv2ENTREZID                        
    else if (array_name == "rgu34a")
        array2gid_map <- rgu34aENTREZID                        
    else if (array_name == "rgu34b")
        array2gid_map <- rgu34bENTREZID                        
    else if (array_name == "rgu34c")
        array2gid_map <- rgu34cENTREZID                        
    else if (array_name == "zebrafish")
        array2gid_map <- zebrafishENTREZID
    else if (array_name == "chicken")
        array2gid_map <- chickenENTREZID
    else if (array_name == "drosgenome1")
        array2gid_map <- drosgenome1ENTREZID
    else if (array_name == "drosophila2")
        array2gid_map <- drosophila2ENTREZID
    else if (array_name == "celegans")
        array2gid_map <- celegansENTREZID
    else
        return (data.frame());
         
    all_probes <- ls(array2gid_map)
    gid <- unlist(mget(all_probes, array2gid_map))
    data.frame(names(gid), unname(gid), rep(array_name, times=length(gid)),rep(tax_id, times=length(gid)))
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
    stop("Please provide a file name");
}
source("http://bioconductor.org/biocLite.R");

df = rbind(get_array2gid ("hgu95av2", "9606"), get_array2gid ("hgu95a", "9606"), get_array2gid ("hgu133a", "9606"),get_array2gid ("hgu133a2", "9606"));
df = rbind(df, get_array2gid ("mgu74a", "10090"), get_array2gid ("mgu74av2", "10090"), get_array2gid ("mgu74b", "10090"), get_array2gid ("mgu74bv2", "10090"), get_array2gid ("mgu74c", "10090"), get_array2gid ("mgu74cv2", "10090"));
df = rbind(df, get_array2gid ("rgu34a", "10116"), get_array2gid ("rgu34b", "10116"), get_array2gid ("rgu34c", "10116"))
df = rbind(df, get_array2gid ("zebrafish", "7955"))
df = rbind(df, get_array2gid ("chicken", "9031"))
df = rbind(df, get_array2gid ("drosgenome1", "7227"), get_array2gid ("drosophila2", "7227"))
df = rbind(df, get_array2gid ("celegans", "6239"))

colnames(df)=c("array_id", "gene_id", "array_name", "tax_id");
df = df[complete.cases(df),] #removing rows with NA values
df = df[!duplicated(df[,c("array_id","gene_id","tax_id")]),] #removing duplicates.

write.table(df, file=args[1], quote=FALSE, sep=",", col.names=TRUE, row.names=FALSE);
