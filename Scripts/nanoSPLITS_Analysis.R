library(tidyverse)
library(ggpubr)
library(ggfortify)
library(corrplot)
library(proDA)
library(data.table)
library(ggbreak)
library(biomaRt)
library(VennDiagram)


mod_pep <- read_tsv("combined_modified_peptide.tsv") %>%
  setnames(., gsub("nanoPOTs_Protoplasts_", "", names(.))) %>%
  setnames(., gsub(" Intensity", "", names(.))) %>%
  setnames(., gsub("_Protoplasts_spintest_", "", names(.)))

meta <- mod_pep[,1:15]

meta2 <- meta %>%
  dplyr:select(`Modified Sequence`, Protein)


allcell <- mod_pep %>%
  filter(grepl("ARATH", Protein)) %>%
  dplyr::select(`Modified Sequence` | contains("cell"))


allcell_long <- allcell %>%
  pivot_longer(!`Modified Sequence`, names_to = "SampleID", values_to = "Intensity")


Cell10 <- mod_pep %>%
  filter(grepl("ARATH", Protein)) %>%
  dplyr::select(`Modified Sequence` | contains("10cell"))

Cell10 <- Cell10 %>% remove_rownames %>% column_to_rownames(var="Modified Sequence")
Cell10[Cell10 == 0] <- NA
Cell10[,1:7] <- log2(Cell10[,1:7])
Cell10 <- as.data.frame(median_normalization(as.matrix(Cell10)))
Cell10_long <- Cell10 %>%
  rownames_to_column(., "Modified Sequence") %>%
  pivot_longer(!`Modified Sequence`, names_to = "SampleID", values_to = "Intensity")

Cell3 <- mod_pep %>%
  filter(grepl("ARATH", Protein)) %>%
  dplyr::select(`Modified Sequence` | contains("3cell")) 

Cell3 <- Cell3 %>% remove_rownames %>% column_to_rownames(var="Modified Sequence")
Cell3[Cell3 == 0] <- NA
Cell3[,1:7] <- log2(Cell3[,1:7])
Cell3 <- as.data.frame(median_normalization(as.matrix(Cell3)))
Cell3_long <- Cell3 %>%
  rownames_to_column(., "Modified Sequence") %>%
  pivot_longer(!`Modified Sequence`, names_to = "SampleID", values_to = "Intensity")

Cell1 <- mod_pep %>%
  filter(grepl("ARATH", Protein)) %>%
  dplyr::select(`Modified Sequence` | contains("1cell")) 


Cell1 <- Cell1 %>% remove_rownames %>% column_to_rownames(var="Modified Sequence")
Cell1[Cell1 == 0] <- NA
Cell1[,1:12] <- log2(Cell1[,1:12])
Cell1 <- as.data.frame(median_normalization(as.matrix(Cell1)))
Cell1_long <- Cell1 %>%
  rownames_to_column(., "Modified Sequence") %>%
  pivot_longer(!`Modified Sequence`, names_to = "SampleID", values_to = "Intensity")


x_long <- rbind(Cell10_long,Cell3_long,Cell1_long) %>%
  mutate(Cells = case_when(grepl("10cell", SampleID) ~ "10",
                           grepl("3cell", SampleID) ~ "3",
                           grepl("1cell", SampleID) ~ "1")) %>%
  mutate(Intensity = 2^Intensity) %>%
  filter(!is.na(Intensity)) %>%
  mutate(isnanoSPLITS = case_when(grepl("nanoSPLITs", SampleID) ~ "Yes",
                                  T ~ "No")) %>%
  filter(isnanoSPLITS == "Yes") %>%
  mutate(SampleID = gsub("nanoSPLITs", "", SampleID)) %>%
  dplyr::select(-isnanoSPLITS)

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


## Supplementary Figure S5B
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
    filter(!is.na(Intensity)) 
  
  
  
  
  
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


### Supplementary Figure S5A
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
  #         #panel.border = element_blank(),
  #         axis.line = element_line()) +
  #   ylab("Total (n)")+
  #   scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  #   scale_fill_manual(values = c("#88419d", "#8c96c6", "#b3cde3", "#edf8fb"))+
  #   scale_alpha_manual(values=c(0.5,0.5)) +
  #   xlab("Number of Protoplasts")+
  #   geom_errorbar(position=position_dodge(width =0.95), aes(x=fct_relevel(Cells, "10", "3","1"),
  #                      ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black",  size=1.2) +
  #   scale_y_break(c(3300,14000), space = 0.1, scales = "free")+
  #   ggsave("proteins_peptide.png", width = 8, height = 4.7)

  
  
  
  ################## RNAseq analysis
  
  
  RNAx_TPM <- read.delim("tpmCounts.txt", check.names = F)  
  
  gene_list2  <- read.delim("idmapping_2024_04_18.txt", check.names = F) 
  
  listDatasets(useMart(biomart="plants_mart",host="plants.ensembl.org"))
  
  mart <- useMart(biomart="plants_mart",host="plants.ensembl.org", dataset = "athaliana_eg_gene")
  
  gene_list <- getBM(attributes = c(#"tair_locus", 
                                    "tair_locus_model",
                                    #"uniprot_gn_trans_name", 
                                   # "uniprot_gn_id",
                                   "uniprotswissprot",
                                    "uniprotsptrembl"),
                    #filters = "ensembl_peptide_id",
                   #  values = c(prot_quant$Protein, unique(meta$Mapped.Proteins)),
                     mart = mart) 
  
  gene_list <- gene_list %>%
    rename(Gene = tair_locus_model) %>%
    mutate(Gene = gsub("\\..*","", Gene )) %>%
    mutate(`Protein ID` = case_when(uniprotswissprot != "" ~ uniprotswissprot,
                                    uniprotswissprot == "" ~ uniprotsptrembl)) %>%
    filter(Gene != "") %>%
    filter(`Protein ID` != "") %>%
    distinct(Gene, `Protein ID`)%>%
    rbind(., gene_list2)
 
  
  RNAx2 <- RNAx_TPM %>%
    pivot_longer(!Gene, names_to = "SampleID", values_to = "TPM") %>%
    mutate(Cells = case_when(grepl("10cell", SampleID) ~ "10",
                             grepl("3cell", SampleID) ~ "3",
                             grepl("1cell", SampleID) ~ "1"))   %>%
    filter(TPM >= 10) %>%
    mutate(TPM = log10(TPM)) %>%
    filter(!is.na(TPM)) %>%
    filter(SampleID %in% unique(x_long$SampleID))
  
  
  ## Figure 2A
  # 
  # RNAx2 %>%
  #   distinct(SampleID,Gene, Cells) %>%
  #   group_by(SampleID) %>%
  #   add_count(name = "n") %>%
  #   ungroup() %>%
  #   distinct(SampleID,Cells,n) %>%
  #   group_by(Cells) %>%
  #   summarise_at("n", funs(mean, sd)) %>%
  #   ggplot()+
  #   aes(x = fct_relevel(Cells, "10", "3","1"), y = mean, fill = fct_relevel(Cells, "10", "3","1"))+
  #   geom_bar(stat= "identity",show.legend = FALSE, alpha = 0.5, color = "black")+
  #   theme_bw(base_size = 20) +
  #   theme(panel.background = element_rect(fill= 'white'),
  #         axis.text.y=element_text(color = 'black', size = 20),
  #         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
  #         panel.grid.minor = element_blank(),
  #         panel.grid.major = element_blank(),
  #         panel.border = element_blank(),
  #         axis.line = element_line(),
  #         text=element_text(family="Helvetica")) +
  #   ylab("Genes Detected (n)")+
  #   scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  #   scale_fill_manual(values = c("#88419d", "#8c96c6", "#b3cde3"))+
  #   geom_errorbar( aes(x=fct_relevel(Cells, "10", "3","1"),
  #                      ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=1, size=1.2) +
  #   xlab("Number of Protoplasts")+
  #   ggsave("genes.png", width = 8, height = 4.7)
  
  
 ################
  
RNA_sc <- RNAx2 %>%
    filter(Cells == 1) %>%
    group_by(`Gene`) %>%
    summarize(TPM = mean(TPM)) %>%
    left_join(., gene_list) %>%
    filter(!is.na(`Protein ID`))

  
  
meta3 <- meta %>%
  distinct(Protein, `Protein ID`)

protein_sc <- protein_long %>%
  filter(Cells == 1) %>%
  left_join(.,meta3) %>%
  group_by(`Protein ID`) %>%
  summarize(Int = mean(Intensity)) 


combined <- full_join(RNA_sc, protein_sc) %>%
  group_by(TPM) %>%
  add_count(name= "n") %>%
  mutate(Var = case_when(is.na(Gene) ~ "Keep",
                        !is.na(Int)~ "Keep",
                         n > 1 & is.na(Int) ~"Keep",
                         n > 1 & !is.na(Int) ~ "Discard",
                         n == 1 ~ "Keep")) %>%
  filter(Var == "Keep") %>%
  group_by(TPM) %>%
  mutate(n2 = case_when(n >1 & n < 250 ~ rank(Int,
                                              ties.method = "random"),
                       T ~ 1)) %>%
  filter(n2 == 1)

mRNA_all <- combined %>%
  filter(!is.na(TPM)) %>%
  pull(`Protein ID`)

prot_all <- combined %>%
  filter(!is.na(Int)) %>%
  pull(`Protein ID`)

## Figure 2B
# venn.diagram(x = list(prot_all, mRNA_all),
#              category.names = c("Protein" , "mRNA"),
#              height = 4.7 , 
#              width = 6 , 
#              units = "in",
#              resolution = 300,
#              disable.logging = T,
#              sigdigs = 0,
#              col=c("darkred", 'lightcoral'),
#              fill = c(alpha("darkred",0.3), alpha('lightcoral',0.4)),
#              cex = 1,
#              ext.text = F,
#              cat.cex = 0,
#              ext.line.lwd = 0.000001,
#              filename = "venn.png",
#              output = TRUE)




## Figure 2C
# png("Gene_protein_overlap.png", height = 4.7, width = 8,
#     units = "in", res = 300)
# combined %>%
#   mutate(Type = case_when(!is.na(TPM) & !is.na(Int) ~ "Detected as protein and transcript",
#                            is.na(TPM) ~ "Detected as protein only",
#                            is.na(Int) ~ "Detected as transcript only")) %>%
#   ggplot()+
#   aes(x = TPM, fill = fct_relevel(Type, rev))+
#   geom_histogram(alpha = 0.7,
#                  position = "identity",
#                  bins = 100)+
#   xlab("log10(TPM)")+
#   ylab("Genes (n)")+
#   theme_minimal(base_size = 20)+
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line())+
#   # scale_fill_brewer(name = "", palette = 4,type = "qual",
#   #                  direction = -1)+
#   theme(legend.position = "bottom",
#         legend.direction = "vertical",
#         legend.box = "horizontal")+
#   labs(fill = "" )
# dev.off()



RNAx_Counts <- read.delim("rawCounts.txt", check.names = F)  

# RNAx_Counts %>%
#   pivot_longer(-Gene,
#                names_to = "SampleID",
#                values_to = "Counts") %>%
#   mutate(Type = case_when(grepl("ATM", Gene) ~ "Mitochondrial",
#                           TRUE ~ "Other")) %>%
#   group_by(SampleID, Type) %>%
#   summarize(Sum1 = sum(Counts)) %>%
#   ungroup() %>%
#   pivot_wider(
#               names_from = "Type",
#               values_from = "Sum1") %>%
#   mutate(PercentMito = 100*Mitochondrial/(Mitochondrial+Other)) %>%
#   mutate(Cells = case_when(grepl("10cell", SampleID) ~ "10",
#                            grepl("3cell", SampleID) ~ "3",
#                            grepl("1cell", SampleID) ~ "1"))   %>%
#     ggplot()+
#     aes(x = fct_relevel(Cells, "10", "3","1"), y = PercentMito, color = fct_relevel(Cells, "10", "3","1"),
#         size = 4)+
#     geom_jitter(width = 0.2, height = 0, show.legend = FALSE)+
#     theme_bw(base_size = 20) +
#     theme(panel.background = element_rect(fill= 'white'),
#           legend.position = "none",
#           axis.text.y=element_text(color = 'black', size = 20),
#           axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.border = element_blank(),
#           axis.line = element_line()) +
#     ylab("Percent Mitochondrial \n Reads")+
#     scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
#     scale_color_manual(values = c("#88419d", "#8c96c6", "#b3cde3", "#edf8fb"))+
#   scale_y_continuous(limits = c(0,5))+
#   xlab("Number of Protoplasts")+
#   ggsave("PercentMito.png", width = 8, height = 4.7)
