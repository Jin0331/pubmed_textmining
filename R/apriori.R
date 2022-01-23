#!/usr/bin/env Rscript

# command line input
arg <- commandArgs(trailingOnly=TRUE)

if(length(arg) == 0){

    stop("At least one argument must be supplied(TCGA barcode).", call.=FALSE) 

} else {
    suppressPackageStartupMessages({
        library(arules)
        library(RMariaDB)
        library(tidyverse)
        library(data.table)})


    # db table load
    con_textmining <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.86", port = 3306, user = "root", password = "sempre813!", dbname = "Textmining")
    con_textmining_result <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.86", port = 3306, user = "root", password = "sempre813!", dbname = "Textmining_Result")

    # [cancerType_gene_disease_pair]
    cancer_type <- arg[1]
    item_table <- tbl(con_textmining, paste0(cancer_type, "_gene_disease_pair")) %>% collect()

    # only human gene
    not_human_pmid <- item_table %>% filter(is.na(HGNC), type == "Gene") %>% pull(1) %>% unique()
        item_table_conversion_raw <- item_table %>% filter(pmid %in% not_human_pmid == F) %>%
        select(-MEDIC) %>% filter(type == "Gene") %>% rename(mapping = HGNC) %>%  ### gene 부분 
        bind_rows(., item_table %>% select(-HGNC) %>% filter(type == "Disease") %>% rename(mapping = MEDIC)) %>% 
        arrange(pmid) %>%  as.data.frame()

    # Transaction 생성
    pmid_list <- split(item_table_conversion_raw$mapping , item_table_conversion_raw$pmid)

    suppressWarnings({
        pmid_trans <- as(pmid_list, "transactions")  
    })

    # parameter create
    parameter <- new("APparameter", support = 10e-9, confidence = 0.1, target = "rules", minlen = 2, maxlen = 2)

    # gene / main disease symbol filter
    gene_symbol <- item_table_conversion_raw %>% filter(type == "Gene") %>% .$mapping %>% unique()
    main_terms <- tbl(con_textmining, "Term_dict_MeSH") %>% collect() %>% 
        filter(Cancer_Type == cancer_type) %>% select(-Cancer_Type) %>% 
        as.character() %>% .[!is.na(.)] # main term index 1, other 2:n

    # run apriori, lhs == gene & rhs == main_terms
    suppressWarnings({
    apriori_result <- apriori(data = pmid_trans, parameter = parameter) %>% 
        arules::subset(x = ., subset = lhs %in% gene_symbol) %>% 
        arules::subset(subset = (rhs %in% main_terms[1]))
    })


    # fisher's exact test

    apriori_result_DF <- apriori_result %>% DATAFRAME() %>% as_tibble() %>% 
        dplyr::select(-coverage) %>% 
        mutate(LHS = str_remove_all(string = LHS, pattern = "(\\{)|(\\})"),
               RHS = str_remove_all(string = RHS, pattern = "(\\{)|(\\})")) %>% 
        arrange(desc(count)) %>% 
        dplyr::rename(Name = LHS, Cancer_type = RHS, SUPPORT = support, CONFIDENCE = confidence, LIFT = lift, COUNT = count)


    # apriori result db import
    apriori_result_DF %>% 
        copy_to(dest = con_textmining_result, df = ., name = paste0(cancer_type, "_Result"), overwrite = T, temporary = F, indexes = list("Name"))

    print(paste0(cancer_type, " is done!@!@!"))

    dbDisconnect(con_textmining)
    dbDisconnect(con_textmining_result)
}