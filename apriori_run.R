library(arules)
library(RMariaDB)
library(tidyverse)
library(data.table)


# db table load
con_textmining <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.86", port = 3306, user = "root", password = "sempre813!", dbname = "Textmining")
# [cancerType_gene_disease_pair]
cancer_type <- "COAD"
item_table <- tbl(con_textmining, paste0(cancer_type, "_gene_disease_pair")) %>% collect()

# only human gene
not_human_pmid <- item_table %>% filter(is.na(HGNC), type == "Gene") %>% pull(1) %>% unique()
item_table_conversion_raw <- item_table %>% filter(pmid %in% not_human_pmid == F) %>%
  select(-MEDIC) %>% filter(type == "Gene") %>% rename(mapping = HGNC) %>%  ### gene 부분 
  bind_rows(., 
            item_table %>% select(-HGNC) %>% filter(type == "Disease") %>% rename(mapping = MEDIC) ### diasese 부분
            ) %>%  ###### 이부분 
  arrange(pmid) %>%  as.data.frame()

# item_table_select <- item_table_conversion_raw %>% 
#   mutate(mapping = str_remove(string = mapping, pattern = "(MESH:)|(OMIM:)")) %>%  ##### 특수문자 error 였던듯 함! #### symbol로 할경우 제외
#   select(pmid, item = mapping)

# basket 생성
pmid_list <- split(item_table_conversion_raw$mapping , item_table_conversion_raw$pmid)
pmid_trans <- as(pmid_list, "transactions")

# parameter create
parameter <- new("APparameter", support = 0.001, confidence = 0.8, target = "rules", minlen = 2, maxlen = 5) 
apriori_result <- apriori(data = pmid_trans, parameter = parameter)
