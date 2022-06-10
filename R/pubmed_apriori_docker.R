#!/usr/bin/env Rscript

# command line input
arg <- commandArgs(trailingOnly=TRUE)

if(length(arg) == 0){
  
  stop("At least one argument must be supplied(TCGA barcode).", call.=FALSE) 
  
} else {
  
  # function
  cancer_type_search <- function(search_terms){
    #https://pubmed.ncbi.nlm.nih.gov/help/#searching-for-a-phrase
    if(length(search_terms) == 0)
      return(stop("no input terms!!"))
    else {
      cancer_terms <- c()
      for(term in search_terms){
        temp <- entrez_search(db="pubmed", term = paste0('"', term, '"[title/abstract]', 'AND 1900:2022[PDAT]'), retmax=999999999)
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
  
  # library load
  suppressPackageStartupMessages({
    library(httr)
    library(xml2)
    library(rentrez)
    library(tidyverse)
    library(RMariaDB)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(parallel)})
  
  con_textmining <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.91", port = 3306, user = "root", password = "sempre813!",
                                   dbname = "Textmining")  
  # date
  run_date <- Sys.time() %>% str_split(pattern = " ") %>% unlist() %>% .[1]
  
  # RAW data directory
  dir.create("RAW_DATA", showWarnings = FALSE)
  
  # search_terms : terms of specific cancer type 
  # pubmed collect
  if(file.exists(paste0("RAW_DATA/", arg[1], "_", run_date,"_gene_disease_pair.txt")) == FALSE){
    print("Pubmed and Pubtator Search!!")
    term <- read_delim(file = "dict/term_dict.txt", delim = "\t", show_col_types = FALSE) %>% 
      dplyr::filter(Cancer_Type == arg[1]) %>% 
      dplyr::select(-Cancer_Type) %>%
      as.character() %>% .[!is.na(.)] %>% 
      tolower()
    
    
    print(term)
    pmid_search <- cancer_type_search(search_terms = term)
    
    # Pubtator --> title, abstract, gene/disease keyword extractions
    url <- "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml"
    n <- 0 # splice count
    index <- 1
    title_abstract <- list()
    gene_disease <- list()
    total <- pmid_search %>% length()
    
    print(paste0("pmid : ", total))
    
    # Pubtator API run
    while(n <= total){ #length(pmid)
      # NA remove
      pmid <- pmid_search[(n+1):(n+1000)] %>% .[!is.na(.)]
      
      # GET 방식 limit 100 request
      re <- FALSE
      tryCatch(
        expr = {
          all_node <- POST(url, body = list(pmids = pmid, concepts = c("gene","disease")), encode = "json") %>% 
            httr::content(encoding = "UTF-8") %>% xml_find_all(".//document")
        },
        error = function(e) { print(e);re <<- TRUE}
      )
      
      # time out 
      if(re) {
        print("time out")
        next             
      } else {
        # print(all_node)
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
      
      n <- n + 1000 # get size
      index <- index + 1 # list 저장 index
      
      print(paste0(n, " is done!@!"))
      rm(all_node, pubtator_node, temp_gene_disease, temp_title_abstract)
      Sys.sleep(1.3) # van 방지 
    }
    
    gene_disease <- gene_disease %>% bind_rows()
    title_abstract <- title_abstract %>% bind_rows()
    
    # gene-disease pair
    # ; 분리 후 행 추가
    gene_disease_gene <- gene_disease %>% 
      filter(type == "Gene") %>% 
      separate_rows(identifier, convert = T) %>% 
      filter(identifier != "None")  %>% 
      mutate(identifier = as.character(identifier)) # None remove
    
    # Entrez to HGNC(official), dict
    cols <- c("ENTREZID", "SYMBOL", "ENSEMBL", "GENENAME")
    ent2hgnc <- gene_disease_gene$identifier %>% unique() %>%  as.character() %>%   ## mapping dict
      AnnotationDbi::select(org.Hs.eg.db, keys = ., columns = cols, keytype="ENTREZID") %>%  ## only homosapience if mouse gene, NA
      dplyr::select(identifier = ENTREZID, HGNC = SYMBOL) %>% distinct() %>% as_tibble()
    gene_disease_gene <- gene_disease_gene %>% left_join(x = ., y = ent2hgnc, by = "identifier")
    
    # disease filter
    ## mesh id mapping
    disease_dict <- read_delim(file = "dict/ctd_dict_medic.txt", delim = "\t", show_col_types = FALSE) %>% 
      dplyr::select(identifier = DiseaseID, MEDIC)
    gene_disease_disease <- gene_disease %>% 
      filter(type == "Disease", identifier != "None") %>% 
      left_join(x = ., y = disease_dict, by = "identifier") %>% 
      mutate(MEDIC = ifelse(is.na(MEDIC), identifier, MEDIC))
    
    # gene-disease bind_rows
    gene_disease_filter <- bind_rows(gene_disease_disease, gene_disease_gene) %>% arrange(pmid)
    
    # title_abstract import
    title_abstract %>% dplyr::select(pmid, year, journal, authors,title, abstract) %>% 
      write_delim(file = paste0("RAW_DATA/", arg[1], "_", run_date, "_pubmed.txt"), delim = "\t")
    
    # gene-disease import
    gene_disease_filter %>% 
      write_delim(file = paste0("RAW_DATA/", arg[1], "_", run_date,"_gene_disease_pair.txt"), delim = "\t")
    
    print(paste0(arg[1], "'s pubmed collect is succeed!!!"))
  }
  
  ## RUN Apriori ====
  suppressPackageStartupMessages({
    library(arules)
  })
  
  
  # db table load
  con_textmining <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.91", port = 3306, user = "root", password = "sempre813!", dbname = "Textmining")

  # [cancerType_gene_disease_pair]
  cancer_type <- arg[1]
  item_table <- read_delim(file = paste0("RAW_DATA/", cancer_type, "_", run_date,"_gene_disease_pair.txt"), delim = "\t",
                           show_col_types = FALSE)
  
  # only human gene
  not_human_pmid <- item_table %>% 
    filter(is.na(HGNC), type == "Gene") %>% 
    pull(1) %>% 
    unique()
  
  item_table_conversion_raw <- item_table %>% 
    filter(pmid %in% not_human_pmid == F) %>%
    dplyr::select(-MEDIC) %>% 
    filter(type == "Gene") %>% 
    dplyr::rename(mapping = HGNC) %>%  ### gene 부분 
    bind_rows(., item_table %>% dplyr::select(-HGNC) %>% filter(type == "Disease") %>% dplyr::rename(mapping = MEDIC)) %>% 
    arrange(pmid) %>%  
    as.data.frame()
  
  # Transaction 생성
  pmid_list <- split(item_table_conversion_raw$mapping , item_table_conversion_raw$pmid)
  
  suppressWarnings({
    pmid_trans <- as(pmid_list, "transactions")  
  })
  
  # parameter create
  parameter <- new("APparameter", support = 10e-9, confidence = 0.1, target = "rules", minlen = 2, maxlen = 2)
  
  # gene / main disease symbol filter
  gene_symbol <- item_table_conversion_raw %>% 
    dplyr::filter(type == "Gene") %>% .$mapping %>%
    unique()
  main_terms <- read_delim(file = "dict/term_dict_mesh.txt", delim = "\t", show_col_types = FALSE) %>% 
    dplyr::filter(Cancer_Type == cancer_type) %>% 
    dplyr::select(-Cancer_Type) %>% 
    as.character() %>% .[!is.na(.)] %>% trimws(which = "both")# main term index 1, other 2:n
  
  # run apriori, lhs == gene & rhs == main_terms

  suppressWarnings({
    apriori_result <- apriori(data = pmid_trans, parameter = parameter) 
  })
  
  # fisher's exact test
  
  apriori_result_DF <- apriori_result %>% 
    DATAFRAME() %>% 
    as_tibble() %>% 
    dplyr::select(-coverage) %>% 
    mutate(LHS = str_remove_all(string = LHS, pattern = "(\\{)|(\\})"),
                                 RHS = str_remove_all(string = RHS, pattern = "(\\{)|(\\})")) %>% 
    filter(RHS %in% main_terms) %>% 
    filter(LHS %in% gene_symbol) %>% 
    arrange(desc(count)) %>% 
    dplyr::rename(gene = LHS, type = RHS, SUPPORT = support, CONFIDENCE = confidence, LIFT = lift, COUNT = count)
  
  
  
  # apriori result db import
  apriori_result_DF %>% 
    copy_to(dest = con_textmining, df = ., name = cancer_type, overwrite = T, temporary = F, indexes = list("gene"))
  
  print(paste0(cancer_type, " is done!@!@!"))
  
  dbDisconnect(con_textmining)
}
