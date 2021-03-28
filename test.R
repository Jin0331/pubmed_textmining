#!/usr/bin/env Rscript

# command line input
arg <- commandArgs(trailingOnly=TRUE)

if(length(arg) == 0){
    stop("At least one argument must be supplied(TCGA barcode).", call.=FALSE) 
}

print(arg[1])