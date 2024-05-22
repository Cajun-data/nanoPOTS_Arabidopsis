## James Fulcher & Kyle Nadeau
## Pacific Northwest National Laboratory
## Single Cell Proteomics

## Imports
library(tidyverse)
library(gprofiler2)

#############################################################
############ Gene Ontology Summarizing function #############
#############################################################


clustGO <- function(x,organism = "mmusculus", Bfactor = 0.5,
                    S_type = c("Survey", "Summary"),
                    bg = bg){

  a <- c(1:length(x))
  
  temp <- data.frame(Algorithm_Test = numeric(0),
                    Clusters = numeric(0),
                    Annotated_Genes = numeric(0),
                    GO_Clusters = numeric(0),
                    GO_Score = numeric(0))  
  
for(i in a) {
  
clust_order <- data.frame((x[[i]]))
  
  clust_order <- clust_order %>%
    distinct(Cluster, Gene)
  
  list_data <- split(subset(clust_order, select = -Cluster), clust_order$Cluster)

  gostres <- gost(query = list_data,
                  organism = organism,
                  ordered_query = F, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = F, 
                  measure_underrepresentation = FALSE, evcodes = TRUE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "annotated", custom_bg = bg, 
                  numeric_ns = "", sources = c("GO"), as_short_link = FALSE)
  if(is.null(gostres)) {
    temp2 <- data.frame(Algorithm_Test = names(x[i]),
                       Clusters = NA,
                       Annotated_Genes = NA,
                       GO_Clusters = NA,
                       GO_Score = NA)
    
    temp[nrow(temp)+1,] <- temp2
  }
  if(!is.null(gostres)) {
  sig_genes <- as.data.frame(gostres$result) %>%
    mutate(intersection = gsub(",NA","",intersection)) %>%
    separate_rows(intersection,sep = ",") %>%
    distinct(intersection,.keep_all = T) %>%
    filter(query_size  < (nrow(clust_order)/3)) %>%
    tally()
  
  GO_Cluster_terms <- as.data.frame(gostres$result) %>%
    mutate(intersection = gsub(",NA","",intersection)) %>%
    mutate(FBscore = (1+(Bfactor^2))*(precision*recall)/((precision*(Bfactor^2)) + recall)) %>%
    group_by(query) %>%
    slice_max(order_by = FBscore, n =1, with_ties = F) %>%
    group_by(term_name) %>%
    summarise(FBscore = mean(FBscore))
  
  CF1 <- nrow(GO_Cluster_terms)/length(unique(clust_order$Cluster))
  CF2 <- sig_genes$n/nrow(clust_order)
  temp2 <- data.frame(Algorithm_Test = names(x[i]),
                    Clusters = length(unique(clust_order$Cluster)),
                    Annotated_Genes = sig_genes$n,
                    GO_Clusters = nrow(GO_Cluster_terms),
                    GO_Score = sum(GO_Cluster_terms$FBscore)*CF1*CF2)
  
  temp[nrow(temp)+1,] <- temp2
  }
}
  if(S_type == "Summary") {
    FscoreMax <- data.frame()
  FscoreMax <- temp %>%
    filter(GO_Score == max(GO_Score))
  
  x <- x[FscoreMax$Algorithm_Test]
  
  clust_order <- data.frame((x[[1]]))
  
  clust_order <- clust_order %>%
    distinct(Cluster, Gene,.keep_all = T)
  
  list_data <- split(subset(clust_order, select = -Cluster), clust_order$Cluster)
  
  gostres <- gost(query = list_data,
                  organism = organism,
                  ordered_query = F, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = F, 
                  measure_underrepresentation = FALSE, evcodes = TRUE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "annotated", custom_bg = bg, 
                  numeric_ns = "", sources = c("GO"), as_short_link = FALSE)
  
  All_GO_terms <- as.data.frame(gostres$result)%>%
    mutate(intersection = gsub(",NA","",intersection)) %>%
    mutate(FBscore = (1+(Bfactor^2))*(precision*recall)/((precision*(Bfactor^2)) + recall)) %>%
    rename(GO_Cluster = query) %>%
    dplyr::select(-effective_domain_size, -source_order, -parents,
                  -evidence_codes, -significant)
    
  Rep_GO_terms <- as.data.frame(gostres$result) %>%
      mutate(intersection = gsub(",NA","",intersection)) %>%
      mutate(FBscore = (1+(Bfactor^2))*(precision*recall)/((precision*(Bfactor^2)) + recall)) %>%
      group_by(query) %>%
      slice_min(order_by = p_value, n =1, with_ties = T) %>%
      dplyr::select(-effective_domain_size, -source_order, -parents,
                    -evidence_codes, -significant) %>%
      mutate(term_name = str_c(term_name, collapse = ","),
             term_id = str_c(term_id, collapse = ","),
             source = str_c(source, collapse = ",")) %>%
      rename(GO_Cluster = query) %>%
      distinct() %>%
    ungroup()
    
    
    GO_Cluster_terms <- as.data.frame(gostres$result) %>%
      mutate(intersection = gsub(",NA","",intersection)) %>%
      mutate(FBscore = (1+(Bfactor^2))*(precision*recall)/((precision*(Bfactor^2)) + recall)) %>%
      group_by(query) %>%
      slice_min(order_by = p_value, n =1, with_ties = T) %>%
      mutate(term_name = str_c(term_name, collapse = ","),
             term_id = str_c(term_id, collapse = ","),
             source = str_c(source, collapse = ","),
             term_size = str_c(term_size, collapse = ",")) %>%
      separate_rows(intersection,sep = ",") %>%
      distinct(intersection, .keep_all = T) %>%
      ungroup() %>%
      dplyr::select(-effective_domain_size, -source_order, -parents,
                    -evidence_codes, -significant, -query,
                    -term_size)
    
    GO_Cluster_terms2 <- as.data.frame(gostres$result) %>%
      mutate(intersection = gsub(",NA","",intersection)) %>%
      separate_rows(intersection,sep = ",") %>%
      filter(!intersection %in% GO_Cluster_terms$intersection) %>%
      group_by(term_name, query) %>%
      add_count(name = "intersection_size2") %>%
      ungroup %>%
      mutate(precision = intersection_size2/query_size) %>%
      mutate(FBscore = (1+(Bfactor^2))*(precision*recall)/((precision*(Bfactor^2)) + recall)) %>%
      group_by(query) %>%
      slice_min(order_by = p_value, n =1, with_ties = T) %>%
      ungroup() %>%
      dplyr::select(-effective_domain_size, -source_order, -parents,
                    -evidence_codes, -significant, -intersection_size2, -query,
                    -term_size)
    combined_df <- full_join(GO_Cluster_terms, GO_Cluster_terms2)
    
    complete_results <- full_join(combined_df, clust_order,
                                  by = c("intersection" = "Gene"),
                                 relationship = "many-to-many") %>%
      distinct(intersection, .keep_all = T) %>%
      rename(Gene = intersection,
             GO_term_pvalue = p_value) %>%
      dplyr::select(intersect(c('Cluster',
                                       'uncertainty',
                                   #    'Entry.Name',
                                    #   "Protein.ID",
                                       "Gene",
                                       "term_name",
                                       "term_id",
                                       "GO_term_pvalue",
                                       "source",
                                       "FBscore"
                                       ), names(.)))

    
    listN <- c("Cluster_GO_Comparison",
               paste0(FscoreMax$Algorithm_Test, "_ALL_GO_terms", sep =""),
               paste0(FscoreMax$Algorithm_Test, "_REP_GO_terms", sep =""),
               paste0(FscoreMax$Algorithm_Test, "_Clusters_Meta"))
    
      result <- list(temp,
                     All_GO_terms,
                     Rep_GO_terms,
                     complete_results
                     )
      names(result) <- listN
      return(result)
  }
  return(temp)
}

