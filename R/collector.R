#!/usr/bin/env Rscript

# command line input
arg <- commandArgs(trailingOnly=TRUE)

if(length(arg) == 0){

    stop("At least one argument must be supplied(TCGA barcode).", call.=FALSE) 

} else {

    # library load
    suppressPackageStartupMessages({
        library(httr)
        library(xml2)
        library(rentrez)
        library(tidyverse)
        library(RMariaDB)
        library(org.Hs.eg.db)
        library(AnnotationDbi)})

    con_textmining <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.86", port = 3306, user = "root", password = "-",
                                 dbname = "Textmining")
                                 
    # search_terms : terms of specific cancer type 
    cancer_type_search <- function(search_terms){
    #https://pubmed.ncbi.nlm.nih.gov/help/#searching-for-a-phrase
    if(length(search_terms) == 0)
        return(stop("no input terms!!"))
    else {
        cancer_terms <- c()
        for(term in search_terms){
        temp <- entrez_search(db="pubmed", term = paste0('"', term, '"', " AND 1900:2020[PDAT]"), retmax=999999999)
        print(temp$QueryTranslation)
        cancer_terms <- c(cancer_terms, temp$ids)
        }
    }
        if(is.list(cancer_terms)){
            cancer_terms %>% unlist() %>% unique() %>% return()
        } else{
            cancer_terms %>% unique() %>% return()
        }
    }

    term <- tbl(con_textmining, "Term_dict") %>% collect() %>% 
        filter(Cancer_Type == arg[1]) %>% dplyr::select(-Cancer_Type) %>%
        as.character() %>% .[!is.na(.)]
    
    print(term)

    pmid_search <- cancer_type_search(search_terms = term)

    # Pubtator --> title, abstract, gene/disease keyword extractions
    url <- "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/"
    n <- 0 # splice count
    index <- 1
    title_abstract <- list()
    gene_disease <- list()
    total <- pmid_search %>% length()

    print(paste0("pmid : ", total))

    # Pubtator API run
    while(n <= total){ #length(pmid)
        # NA remove
        pmid <- pmid_search[(n+1):(n+100)] %>% .[!is.na(.)]  %>% paste(collapse = ",")
        
        # GET 방식 limit 100 request
        re <- FALSE
        tryCatch(
            expr = {
            all_node <- GET(paste0(url, "biocxml?pmids=", pmid, "&concepts=gene,disease")) %>%
                httr::content(encoding = "UTF-8") %>% xml_find_all(".//document")
            },
            error = function(e) { print(e);re <<- TRUE}
        )

        # time out 
        if(re) {
            print("time out")
            next             
        } else {
            print(all_node)
        }
        
        pubtator_node <- mclapply(X = all_node, FUN = function(pubtator_node){
            # pmid extraction
            pmid <- pubtator_node %>% xml_find_all(., ".//id") %>% xml_text()
            
            # search node
            xml_start <-  xml_find_all(pubtator_node, ".//annotation") # title, abstract의 모든 annotation
            xml_title <- xml_find_all(pubtator_node, ".//passage")[1] # title area
            xml_abtstart <- xml_find_all(pubtator_node, ".//passage")[2] # abstract area
            
            ## title
            title_node <- xml_children(xml_title)
            col_name <- title_node %>% xml_name() 
            col_name_compare <- title_node %>% xml_attrs() %>% as.character()
            
            # xml_anme, xml_attrs 비교 
            for(index in seq(col_name_compare)){
                if(col_name_compare[index] == "character(0)"){
                    col_name_compare[index] <- col_name[index]
                }
            }
            temp_vec <- title_node %>% xml_text()
            names(temp_vec) <- col_name_compare %>% map_chr(.x = ., tolower)
            title_ab_DF <- bind_cols(tibble(pmid = pmid), as.data.frame(t(temp_vec))) %>% as_tibble()
            
            ## abstract
            abstract <- xml_abtstart %>% 
                xml_find_first(., ".//text") %>% 
                xml_text() %>% 
                tibble(abstract = .)
            
            title_ab_DF <-  title_ab_DF %>% bind_cols(., abstract) 
            
            # %>% 
            #   select(pmid, year, authors, title = text, abstract) # try exception
            
            # annotation gene/disease to DF
            gene_disease_DF <- map(.x = xml_start, function(value){
            temp <- xml_children(value) 
            vec_name <- temp %>% xml_attrs() %>% as.character()
            vec <- temp %>% xml_text()
            names(vec) <- vec_name %>% map_chr(.x = ., tolower)
            
            tibble(pmid = pmid, type = vec["type"], NCBI_homologene = vec["ncbi homologene"], 
                    identifier = vec["identifier"], Symbol = vec["character(0)"]) %>% return()
            }) %>% bind_rows()
            
            return(list(title = title_ab_DF, info = gene_disease_DF))
        }, mc.cores = 10) 
        
        # 회차별 모든 info 정보 combine, PK = pmid
        temp_gene_disease <- map(pubtator_node, function(value){
            value[["info"]]
        }) %>% bind_rows()
        # gene_disease <- gene_disease %>% bind_rows(., temp_gene_disease)
        gene_disease[[index]] <- temp_gene_disease
        
        # 회차별 모든 title 정보 combine, PK = pmid
        temp_title_abstract <- map(pubtator_node, function(value){
            value[["title"]]
         }) %>% bind_rows() %>% 
            dplyr::select(pmid, year, journal, authors, title = text, abstract)
        # title_abstract <- title_abstract %>% bind_rows(., temp_title_abstract)
        title_abstract[[index]] <- temp_title_abstract
        
        n <- n + 100 # get size
        index <- index + 1 # list 저장 index
        
        print(paste0(n, " is done!@!"))
        rm(all_node, pubtator_node, temp_gene_disease, temp_title_abstract)
        Sys.sleep(1.3) # van 방지 
    }

        gene_disease <- gene_disease %>% bind_rows()
        title_abstract <- title_abstract %>% bind_rows()

    # gene-disease pair
    # ; 분리 후 행 추가
    gene_disease_gene <- gene_disease %>% filter(type == "Gene") %>% 
        separate_rows(identifier, convert = T) %>% 
        filter(identifier != "None")  %>% 
        mutate(identifier = as.character(identifier)) # None remove

    # Entrez to HGNC(official), dict
    cols <- c("ENTREZID", "SYMBOL", "ENSEMBL", "GENENAME")
    ent2hgnc <- gene_disease_gene$identifier %>% unique() %>%  as.character() %>%   ## mapping dict
        AnnotationDbi::select(org.Hs.eg.db, keys = ., columns = cols, keytype="ENTREZID") %>% 
        dplyr::select(identifier = ENTREZID, HGNC = SYMBOL) %>% distinct() %>% as_tibble()
    gene_disease_gene <- gene_disease_gene %>% left_join(x = ., y = ent2hgnc, by = "identifier")

    # disease filter
    ## mesh id mapping
    disease_dict <- tbl(con_textmining, "Disease_dict_medic") %>% collect() %>% 
        dplyr::select(identifier = DiseaseID, MEDIC)
    gene_disease_disease <- gene_disease %>% 
        filter(type == "Disease", identifier != "None") %>% 
        left_join(x = ., y = disease_dict, by = "identifier") %>% 
        mutate(MEDIC = ifelse(is.na(MEDIC), identifier, MEDIC))
    
    # gene-disease bind_rows
    gene_disease_filter <- bind_rows(gene_disease_disease, gene_disease_gene) %>% arrange(pmid)

    # title_abstract import
    title_abstract %>% dplyr::select(pmid, year, journal, authors,title, abstract) %>% 
        copy_to(dest = con_textmining, df = ., name = paste0(arg[1], "_pubmed"), overwrite = T, temporary = F, indexes = list("pmid"))

    # gene-disease import
    gene_disease_filter %>% 
        copy_to(dest = con_textmining, df = ., name = paste0(arg[1], "_gene_disease_pair"), overwrite = T, temporary = F, indexes = list("pmid"))

    print(paste0(arg[1], " is succeed!!!"))
    dbDisconnect(con_textmining)
}
