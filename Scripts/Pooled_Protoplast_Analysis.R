library(tidyverse)
library(ggpubr)
library(ggfortify)
library(corrplot)
library(proDA)
library(data.table)
library(ggbreak)


mod_pep <- read_tsv("combined_modified_peptide.tsv") %>%
  setnames(., gsub("nanoPOTs_Protoplasts_", "", names(.))) %>%
  setnames(., gsub(" Intensity", "", names(.))) 

meta <- mod_pep[,1:15]

meta2 <- meta %>%
  select(`Modified Sequence`, Protein)


allcell <- mod_pep %>%
  filter(grepl("ARATH", Protein)) %>%
  select(`Modified Sequence` | contains("cell"))


allcell_long <- allcell %>%
  pivot_longer(!`Modified Sequence`, names_to = "SampleID", values_to = "Intensity")


Cell10 <- mod_pep %>%
  filter(grepl("ARATH", Protein)) %>%
  select(`Modified Sequence` | contains("10cell"))

Cell10 <- Cell10 %>% remove_rownames %>% column_to_rownames(var="Modified Sequence")
Cell10[Cell10 == 0] <- NA
Cell10[,1:4] <- log2(Cell10[,1:4])
Cell10 <- as.data.frame(median_normalization(as.matrix(Cell10)))
Cell10_long <- Cell10 %>%
  rownames_to_column(., "Modified Sequence") %>%
  pivot_longer(!`Modified Sequence`, names_to = "SampleID", values_to = "Intensity")

Cell3 <- mod_pep %>%
  filter(grepl("ARATH", Protein)) %>%
  select(`Modified Sequence` | contains("3cell")) 

Cell3 <- Cell3 %>% remove_rownames %>% column_to_rownames(var="Modified Sequence")
Cell3[Cell3 == 0] <- NA
Cell3[,1:4] <- log2(Cell3[,1:4])
Cell3 <- as.data.frame(median_normalization(as.matrix(Cell3)))
Cell3_long <- Cell3 %>%
  rownames_to_column(., "Modified Sequence") %>%
  pivot_longer(!`Modified Sequence`, names_to = "SampleID", values_to = "Intensity")

Cell1 <- mod_pep %>%
  filter(grepl("ARATH", Protein)) %>%
  select(`Modified Sequence` | contains("1cell")) 


Cell1 <- Cell1 %>% remove_rownames %>% column_to_rownames(var="Modified Sequence")
Cell1[Cell1 == 0] <- NA
Cell1[,1:8] <- log2(Cell1[,1:8])
Cell1 <- as.data.frame(median_normalization(as.matrix(Cell1)))
Cell1_long <- Cell1 %>%
  rownames_to_column(., "Modified Sequence") %>%
  pivot_longer(!`Modified Sequence`, names_to = "SampleID", values_to = "Intensity")


x_long <- rbind(Cell10_long,Cell3_long,Cell1_long) %>%
  mutate(Cells = case_when(grepl("10cell", SampleID) ~ "10",
                           grepl("3cell", SampleID) ~ "3",
                           grepl("1cell", SampleID) ~ "1")) %>%
  mutate(Intensity = 2^Intensity) %>%
  filter(!is.na(Intensity))

xmed <- x_long %>%
  filter(!is.na(Intensity)) %>%
  group_by(Cells, `Modified Sequence`) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 2) %>%
  mutate(CV = sd(Intensity)/mean(Intensity)) %>%
  distinct(Cells, `Modified Sequence`,CV) %>%
  group_by(Cells) %>%
  mutate(med = median(CV)) %>%
  distinct(Cells, med) %>%
  ungroup()


# Figure 1B
# plot <- x_long %>%
#   filter(!is.na(Intensity)) %>%
#   group_by(Cells, `Modified Sequence`) %>%
#   add_count(name = "Obs") %>%
#   filter(Obs >= 2) %>%
#   mutate(CV = sd(Intensity)/mean(Intensity)) %>%
#   distinct(Cells, `Modified Sequence`,CV) %>%
#   ggplot()+
#   aes(x =fct_relevel(Cells, "10", "3","1"), y = CV, fill = fct_relevel(Cells, "10", "3","1"))+
#   geom_violin(alpha = 0.5, show.legend = FALSE, size = 0)+
#   scale_y_continuous(limits = c(0,1.6))+
#   theme_bw(base_size = 20) +
#   theme(legend.position = "none",
#         panel.background = element_rect(fill= 'white'),
#         axis.text.y=element_text(color = 'black', size = 20),
#         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line()) +
#   #xlab("Condition")+
#   ylab("Coefficient of Variation (CV)") +
#   xlab("Number of Protoplasts")+
#   scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
#   geom_text(data = xmed, aes(x = fct_relevel(Cells, "10", "3","1"), y = med, label = paste(round(med, digits = 2))), 
#             size = 9, vjust = -2.2,hjust = -0.25, show.legend = FALSE) +
#   scale_fill_manual(values = c("#88419d", "#8c96c6", "#b3cde3"))
# 
# data_summary <- function(x) {
#   m <- median(x)
#   ymin <- m-sd(x)
#   ymax <- m+sd(x)
#   return(c(y=m,ymin=ymin,ymax=ymax))
# }
# plot +
#   stat_summary(fun.data=data_summary)+
#   ggsave("CVs.png", width = 8, height = 4.7)


x_long %>%
  group_by(Cells) %>%
  summarize(log2(median(Intensity))) %>%
  view()


########## Rollup

protein <- x_long %>%
  inner_join(., meta2) %>%
  distinct(Protein, `Modified Sequence`, SampleID, Intensity) %>%
  filter(!is.na(Intensity)) %>%
  mutate(Intensity = log2(Intensity)) %>%
  dplyr::rename(protein_list = Protein,
                sample_list = SampleID,
                id = `Modified Sequence`,
                quant = Intensity) %>%
  as.list() %>%
  iq::fast_MaxLFQ(.) %>%
  .[[1]]




  
  
protein_long <- protein %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein") %>%
  pivot_longer(-Protein, names_to = "SampleID",
               values_to = "Intensity") %>%
    mutate(Cells = case_when(grepl("10cell", SampleID) ~ "10",
                             grepl("3cell", SampleID) ~ "3",
                             grepl("1cell", SampleID) ~ "1")) %>%
    filter(!is.na(Intensity)) %>%
  mutate(Intensity = 2^Intensity)

protein_long %>%
  group_by(Cells) %>%
  summarize(log2(median(Intensity))) %>%
  view()


### Figure 1C
# protein_long %>%
#   ggplot()+
#   aes(x = fct_relevel(Cells, "10", "3", "1"),
#       y = log2(Intensity), fill = fct_relevel(Cells, "10", "3", "1"))+
#   geom_boxplot(alpha = 1, show.legend = FALSE, size = 1)+
#   #scale_y_continuous(limits = c(0,20))+
#   theme_bw(base_size = 20) +
#   theme(panel.background = element_rect(fill= 'white'),
#         axis.text.y=element_text(color = 'black', size = 20),
#         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line()) +
#   #xlab("Condition")+
#   labs(y = expression("log"["2"]~"(Intensity)"))+
#   xlab("Number of Protoplasts")+
#   scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
#   scale_fill_manual(values = c( "#88419d", "#8c96c6", "#b3cde3",
#                                 "#edf8fb"))+
#   geom_hline(yintercept = 15.76, linetype = "dashed", color = "black")+
#   geom_hline(yintercept = 15.76+log2(3), linetype = "dashed", color = "red")+
#   geom_hline(yintercept = 15.76+log2(10), linetype = 'dashed', color = "red")+
#   ggsave(file = "log2_LFQensities_FragPipe.png", width = 8, height = 4.7)
#   
  x_peptide <- x_long %>%
    distinct(SampleID, `Modified Sequence`, Cells) %>%
    group_by(SampleID) %>%
    add_count(name = "n") %>%
    group_by(Cells) %>%
    summarise_at("n", funs(mean, sd)) %>%
    mutate(Measure = "Peptides")


  x_pepgene <- protein_long %>%
    filter(!is.na(Intensity)) %>%
    distinct(SampleID,Protein, Cells) %>%
    group_by(SampleID) %>%
    add_count(name = "n") %>%
    distinct(SampleID,Cells,n) %>%
    group_by(Cells) %>%
    summarise_at("n", funs(mean, sd)) %>%
    mutate(Measure =  "Proteins") %>%
    full_join(.,x_peptide)


### FIgure 1A
  # x_pepgene %>%
  #   ggplot()+
  #   aes(x = fct_relevel(Cells, "10", "3","1"), y = mean, fill = fct_relevel(Cells, "10", "3","1"),
  #       alpha = Measure)+
  #   geom_bar(stat= "identity", position = position_dodge(width = 0.95), show.legend = FALSE, color = "black")+
  #   theme_bw(base_size = 20) +
  #   theme(panel.background = element_rect(fill= 'white'),
  #         legend.position = "none",
  #         axis.text.y=element_text(color = 'black', size = 20),
  #         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
  #         panel.grid.minor = element_blank(),
  #         panel.grid.major = element_blank(),
  #        # panel.border = element_blank(),
  #         axis.line = element_line()) +
  #   ylab("Total Identifications (n)")+
  #   scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  #   scale_fill_manual(values = c("#88419d", "#8c96c6", "#b3cde3", "#edf8fb"))+
  #   scale_alpha_manual(values=c(0.5,0.5)) +
  #   xlab("Number of Protoplasts")+
  #   geom_errorbar(position=position_dodge(width =0.95), aes(x=fct_relevel(Cells, "10", "3","1"),
  #                      ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black",  size=1.2) +
  #   scale_y_break(c(3300,16000), space = 0.1, scales = "free")+
  #   ggsave("proteins_peptide.png", width = 8, height = 4.7)

  

  
  
  corx <- protein_long %>%
    filter(!is.na(Intensity)) %>%
    mutate(Intensity = log2(Intensity)) %>%
    select(-Cells) %>%
    spread(SampleID, Intensity) %>%
    ungroup() 
  
  corx2 <- corx %>% remove_rownames %>% column_to_rownames(var="Protein")
  
  corTest <- cor(corx2, method = "pearson", use = "pairwise.complete.obs") 
  
  
  ### Supplementary Figure S3
  # pdf(file = "Pearson_clustering_protein.pdf")
  # 
  # corrplot(corTest, method = 'shade', order = 'hclust',
  #          hclust.method = 'mcquitty',
  #          is.corr = FALSE, tl.srt = 45,
  #          col.lim = c(0,1),
  #          diag = FALSE,tl.col = "black")
  # 
  # dev.off()

  
  
