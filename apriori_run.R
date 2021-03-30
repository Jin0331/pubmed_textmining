suppressPackageStartupMessages({
  library(arules)
  library(RMariaDB)
  library(tidyverse)
  library(data.table)
})


# db table load
con_textmining <- DBI::dbConnect(drv = MariaDB(), host = "192.168.0.86", port = 3306, user = "root", password = "sempre813!", dbname = "Textmining")

# [cancerType_gene_disease_pair]
cancer_type <- "THCA"
item_table <- tbl(con_textmining, paste0(cancer_type, "_gene_disease_pair")) %>% collect()

# only human gene
not_human_pmid <- item_table %>% filter(is.na(HGNC), type == "Gene") %>% pull(1) %>% unique()
item_table_conversion_raw <- item_table %>% filter(pmid %in% not_human_pmid == F) %>%
  select(-MEDIC) %>% filter(type == "Gene") %>% rename(mapping = HGNC) %>%  ### gene 부분 
  bind_rows(., 
            item_table %>% select(-HGNC) %>% filter(type == "Disease") %>% rename(mapping = MEDIC) ### diasese 부분
            ) %>%  mutate(mapping = ifelse(type == "Gene", paste0("G_", mapping), paste0("D_", mapping))) %>% ###### 이부분 
  arrange(pmid) %>%  as.data.frame()

# item_table_select <- item_table_conversion_raw %>% 
#   mutate(mapping = str_remove(string = mapping, pattern = "(MESH:)|(OMIM:)")) %>%  ##### 특수문자 error 였던듯 함! #### symbol로 할경우 제외
#   select(pmid, item = mapping)

# basket 생성
pmid_list <- split(item_table_conversion_raw$mapping , item_table_conversion_raw$pmid)
pmid_trans <- as(pmid_list, "transactions")

# parameter create
# run apriori
parameter <- new("APparameter", support = 10e-6, confidence = 0.2, target = "rules", minlen = 2, maxlen = 2)  # 0.001
apriori_result <- apriori(data = pmid_trans, parameter = parameter)

# fishersExactTest
quality(apriori_result) <- bind_cols(quality(apriori_result), 
                                     p_value = interestMeasure(apriori_result, measure = "fishersExactTest",
                                                              complements = T, transactions = pmid_trans, reuse = F))
# gene symbol
gene_symbol <- item_table_conversion_raw %>% filter(type == "Gene") %>% .$mapping %>% unique()

# disease symbol
disease_dict <- tbl(con_textmining, "Disease_dict_medic") %>% collect() 

term <- tbl(con_textmining, "Term_dict") %>% collect() %>% 
  filter(Cancer_Type == cancer_type) %>% dplyr::select(-Cancer_Type) %>%
  as.character() %>% .[!is.na(.)]


apriori_result_filter <- apriori_result %>% DATAFRAME() %>% as_tibble() %>% 
  select(-coverage) %>% 
  filter(str_detect(string = LHS, pattern = "\\{G")) 


apriori_result_filter$RHS %>% table() %>% sort(decreasing = T) %>% View()
  # filter(str_detect(string = RHS, pattern = "\\{D_Colorectal Neoplasms\\}")) %>% 
  # arrange(p_value)

fdr_value <- apriori_result_filter$p_value %>% p.adjust(p = ., method = "fdr", n = length(.)) %>% tibble(FDR = .)
apriori_result_filter %>% bind_cols(., fdr_value) %>% arrange(FDR) %>% #View()
  write_delim(file = "COAD_textming(apriori)_dummy.txt", delim = "\t")
