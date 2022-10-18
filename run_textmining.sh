#!/bin/bash

while getopts f:t: flag
do
    case "${flag}" in
        f) foldername=${OPTARG};;
        t) abb=${OPTARG};;
    esac
done
 
while IFS= read -r line || [ -n "$line" ]; do

    docker stop textmining & docker rm textmining

    echo $line
    docker run -t --name textmining sempre813/textmining:latest Rscript R/pubmed_apriori_docker.R ${line}
    docker cp textmining:/home/rstudio/RAW_DATA/. ${foldername}

done < $abb
