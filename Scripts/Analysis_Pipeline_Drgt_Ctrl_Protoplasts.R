#### Import packages
library(tidyverse)
library(circlize)
library(proDA)
library(DreamAI)
library(mclust)
library(impute)
library(ggpubr)
library(sva)
library(dendextend)
library(rstatix)
library(gt)
library(qualpalr)
library(DEP)
library(stats)
library(gprofiler2)
library(ComplexHeatmap)
library(PCAtools)
library(EnhancedVolcano)
library(scales)
library(factoextra)
library(RColorBrewer)
library(dialects)
library(biomaRt)
library(ggpointdensity)
library(viridis)
library(VennDiagram)

############# Functions

digestFasta2 <- function (fasta.df = NULL) 
{
  oldw <- getOption("warn")
  options(warn = -1)
  if (missing(fasta.df)) 
    stop("ERROR: Need to specify fasta data frame")
  fasta.test <- c("accession", "name", "sequence")
  fasta.colnames <- colnames(fasta.df)
  if (!(all(fasta.colnames == fasta.test))) 
    stop("ERROR: SRL format not recognised.")
  fasta.dt <- data.table::data.table(fasta.df)
  x <- gsub("((?<=[K])(?=[^P]))|((?<=[R])(?=[^P]))", ",", 
            fasta.dt$sequence, perl = TRUE)
  fasta.dt$sequence <- x
  digest.dt <- splitstackshape::cSplit(fasta.dt, "sequence", 
                                       sep = ",", direction = "long", drop = FALSE, fixed = FALSE)
  digest.dt$sequence <- as.character(digest.dt$sequence)
  digest.2.dt <- digest.dt[!(nchar(digest.dt$sequence) < 5 | 
                               nchar(digest.dt$sequence) > 51), ]
  # digest.3.dt <- unique(digest.2.dt, by = "sequence")
  digest.df <- data.frame(digest.2.dt)
  options(warn = oldw)
  return(digest.df)
}


rotate <- function(x) t(apply(x, 2, rev))

reorder_mat <- function(mat, order){
  n <- length(order)
  if(!inherits(mat, "matrix")){
    stop("'mat' must be a matrix")
  } else if (!inherits(order, "character")){
    stop("'order' must be a character vector")
  } else if(!(isSymmetric(mat))){
    stop("The matrix 'mat' must be symmetric")
  } else if (n != length(colnames(mat))){
    stop("'order' must have as many elements as there are rows and
         columns in 'mat'")
  } else if(length(which(colnames(mat) %in% order)) != n){
    stop("The column names of the matrix you want to reorder must
         be present in the vector 'order'")
  } else if (length(which(row.names(mat) %in% order)) != n){
    print("The row names of the matrix you want to reorder must
          be present in the vector 'order'")
  } else {
    
    mat2 <- mat[order, order]
    
    return(mat2)
  }
}


#################### Import Data

mod_pep <- read_tsv("combined_modified_peptide.tsv") %>%
  setNames(., gsub("nanoPOTs_Protoplasts_", "", names(.))) %>%
  setNames(., gsub(" Intensity", "", names(.))) 

meta <- mod_pep[,1:15]

meta2 <- meta %>%
  dplyr::select(`Modified Sequence`, Protein)

######## Clean Up MetaData and import TAIR IDs

meta_prot <- meta %>%
  distinct( Protein, .keep_all = T) %>%
  mutate(name = `Protein ID`) %>%
  rename(ID = `Protein ID`) %>%
  filter(grepl("ARATH", Protein)) %>%
  dplyr::select(Protein, ID, `Entry Name`, name,Gene, `Protein Description`)

missing_IDs <- read_delim("missingIDs.txt")
listDatasets(useMart(biomart="plants_mart",host="plants.ensembl.org"))

mart <- useMart(biomart="plants_mart",host="plants.ensembl.org", dataset = "athaliana_eg_gene")

gene_list <- getBM(attributes = c(
  "tair_locus_model",
  "uniprotswissprot",
  "uniprotsptrembl"),
  mart = mart) 


gene_list <- gene_list %>%
  rename(TAIR = tair_locus_model) %>%
  mutate(ID= case_when(uniprotswissprot != "" ~ uniprotswissprot,
                                  uniprotswissprot == "" ~ uniprotsptrembl)) %>%
  filter(TAIR != "") %>%
  filter(ID != "") %>%
  distinct(TAIR, ID) %>%
  full_join(., missing_IDs)

IDcheck1 <- meta_prot %>%
  left_join(., gene_list) %>%
  mutate(name = case_when(is.na(TAIR) ~ ID,
                          TRUE ~ TAIR)) %>%
  filter(!is.na(`Protein Description`)) %>%
  group_by(name) %>%
  add_count(name  = "n") %>%
  filter(n > 1) %>%
  dplyr::select(-n) %>%
  mutate(name = paste(TAIR, row_number(), sep =(".")))
  
IDcheck2 <- meta_prot %>%
  left_join(., gene_list) %>%
  mutate(name = case_when(is.na(TAIR) ~ ID,
                          TRUE ~ TAIR)) %>%
  filter(!is.na(`Protein Description`)) %>%
  group_by(name) %>%
  add_count(name  = "n") %>%
  filter(n == 1) %>%
  group_by(ID) %>%
  add_count(name  = "n") %>%
  filter( n > 1) %>%
  distinct(ID, .keep_all = T) %>%
  dplyr::select(-n)

meta_prot <- meta_prot %>%
  left_join(., gene_list) %>%
  mutate(name = case_when(is.na(TAIR) ~ ID,
                          TRUE ~ TAIR)) %>%
  filter(!is.na(`Protein Description`)) %>%
  group_by(name) %>%
  add_count(name  = "n") %>%
  filter(n == 1) %>%
  group_by(ID) %>%
  add_count(name  = "n2") %>%
  filter( n2 == 1) %>%
  dplyr::select(-n, -n2) %>%
  full_join(., IDcheck1) %>%
  full_join(., IDcheck2) %>%
  distinct(ID, .keep_all = T)

bg <- meta_prot %>%
  pull(name)

################################ Median normalization
allcell <- mod_pep %>%
  filter(grepl("ARATH", Protein)) %>%
  dplyr::select(`Modified Sequence` | contains("nanoPOTs"))


allcell <- allcell %>% remove_rownames %>% column_to_rownames(var="Modified Sequence")
allcell[allcell == 0] <- NA
allcell[,1:117] <- log2(allcell[,1:117])
allcell <- as.data.frame(median_normalization(as.matrix(allcell)))

x_long <- allcell %>%
  rownames_to_column(var="Modified Sequence") %>%
  pivot_longer(!`Modified Sequence`, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Intensity = 2^Intensity) %>%
  filter(!is.na(Intensity)) 

########### Roll peptide values up to protein level
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
  rownames_to_column(var="Protein") %>%
  pivot_longer(!`Protein`, names_to = "SampleID", values_to = "Intensity") %>%
  full_join(., meta_prot) %>%
  filter(!is.na(Intensity)) %>%
  dplyr::select(name, SampleID, Intensity)


### Figure 3A
# protein_long %>%
#   distinct(name, SampleID) %>%
#   mutate(Group1 = case_when(grepl("Ctrl_Day0", SampleID) ~ "Ctrl-50%",
#                             grepl("Drgt_Day0", SampleID)~ "WD-50%",
#                             grepl("Ctrl_Day7", SampleID) ~ "Ctrl-30%",
#                             grepl("Drgt_Day7", SampleID) ~ "WD-30%"
#   )) %>%
#   group_by(SampleID) %>%
#   add_count(name = "n") %>%
#   distinct(SampleID, n, Group1) %>%
#   ungroup()  %>% 
#   mutate(Group1 = factor(Group1, c("Ctrl-50%","WD-50%",
#                                    "Ctrl-30%","WD-30%"))) %>%
#   ggplot()+
#   aes(y = n, x = Group1, color = Group1)+
#   geom_boxplot(outliers = F )+
#   geom_jitter(position = position_jitter(w = 0.1, h = 0), size = 3)+
#   theme_bw(base_size = 20) +
#   theme(legend.position = "none",
#         panel.background = element_rect(fill= 'white'),
#         axis.text.y=element_text(color = 'black', size = 22),
#         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 22),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         text=element_text(family="Helvetica")) +
#   ylab("Proteins Identified (n)")+
#   xlab("")+
#   scale_y_continuous(limits = c(2000,3300))+
#   geom_hline(yintercept = 2873, linetype = "dashed")
# ggsave("QC_proteins.png", width = 8, height = 4.7)



Sample_Metadata <-protein_long %>%
  distinct(name, SampleID) %>%
  mutate(Group1 = case_when(grepl("Ctrl_Day0", SampleID) ~ "Ctrl-50%",
                            grepl("Drgt_Day0", SampleID)~ "WD-50%",
                            grepl("Ctrl_Day7", SampleID) ~ "Ctrl-30%",
                            grepl("Drgt_Day7", SampleID) ~ "WD-30%"
  )) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  distinct(SampleID, n, Group1) %>%
  mutate(Group1 = factor(Group1, c("Ctrl-50%","WD-50%",
                                   "Ctrl-30%","WD-30%"))) %>%
  ungroup()

######### Filter for proteins with >=90% observations
protein_wide <- protein_long %>%
  filter(SampleID %in% Sample_Metadata$SampleID) %>%
  group_by(name) %>% 
  add_count(name = "n") %>%
  filter(n >= 117*0.9) %>% 
  dplyr::select(-n) %>%
  spread(SampleID, Intensity) %>% 
  column_to_rownames(var= "name")

protoP_impute <-  DreamAI(protein_wide, k = 10, maxiter_MF = 10, ntree = 100, 
                         maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2), 
                         gamma_ADMIN = NA, gamma = 50, CV = FALSE, 
                         fillmethod = "row_mean", maxiter_RegImpute = 10, 
                         conv_nrmse = 1e-06, iter_SpectroFM = 40, method = "KNN", 
                         out = c("KNN")) 


protoP_impute <- as.data.frame(protoP_impute$Ensemble)

pca_j<- PCAtools::pca(protoP_impute, scale = TRUE, center = T)


Group <- Sample_Metadata$Group1


############ Supplementary Figure S6
# PCAtools::biplot(pca_j, x = "PC1", y =  "PC2", lab = Group,
#                  showLoadings = F,
#                  ntopLoadings = 2,
#                  labSize = 0, drawConnectors = FALSE,
#                  legendPosition = "right",
#                  encircle = TRUE
# )# ggsave("PCA_precluster.png", height = 6, width = 8)




############# Identify  and remove fragmented protoplasts
PCA_transform <- as.data.frame(pca_j$rotated[,1:2])
fviz_nbclust(PCA_transform, kmeans, method = 'wss')
set.seed(25)
PCA_Km<- kmeans(PCA_transform, centers = 2, nstart = 50)

######## Supplementary Figure S7A
# png("fragmented_protoplasts_cluster.png", height = 4.7, width = 8, 
#     units = "in", res = 300)
# fviz_cluster(PCA_Km, data = PCA_transform,
#              labelsize = NA,
# )+
#   theme_bw(base_size = 20) +
#   theme(legend.position = "bottom",
#         panel.background = element_rect(fill= 'white'),
#         axis.text.y=element_text(color = 'black', size = 20),
#         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         text=element_text(family="Helvetica")) 
# dev.off()



meta2 <- as.data.frame(PCA_Km$cluster) %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(Cluster1 = `PCA_Km$cluster`) %>%
  dplyr::select(-`PCA_Km$cluster`)

meta3 <- meta2 %>% 
           filter(Cluster1 != 2)
meta4 <- meta2 %>% 
  filter(Cluster1 == 2)



protein_long2 <- protein_long %>%
  filter(!SampleID %in% meta4$SampleID) 


############### Supplementary Figure S7B
# stat.test2 <-  protein_long %>%
#   distinct(name, SampleID) %>%
#   mutate(Group2 = case_when(SampleID %in% meta4$SampleID ~ "Fragmented",
#                             T ~ "Intact")) %>%
#   group_by(SampleID) %>%
#   add_count(name = "n")  %>%
#   distinct(SampleID, n, Group2) %>%
#   ungroup() %>%
#   t_test(n ~ Group2) %>%
#   add_xy_position(x = "Group2", fun ="median_iqr") %>%
#   add_significance()
# 
# stat.test2 <- stat.test2 %>%
#   mutate(y.position = y.position + 150)
# 
# protein_long %>%
#   distinct(name, SampleID) %>%
#   mutate(Group1 = case_when(grepl("Ctrl_Day0", SampleID) ~ "Ctrl-50%",
#                             grepl("Drgt_Day0", SampleID)~ "WD-50%",
#                             grepl("Ctrl_Day7", SampleID) ~ "Ctrl-30%",
#                             grepl("Drgt_Day7", SampleID) ~ "WD-30%"
#   )) %>%
#   group_by(SampleID) %>%
#   add_count(name = "n") %>%
#   distinct(SampleID, n, Group1) %>%
#   ungroup()  %>%
#   mutate(Group1 = factor(Group1, c("Ctrl-50%","WD-50%",
#                                    "Ctrl-30%","WD-30%"))) %>%
#   mutate(Group2 = case_when(SampleID %in% meta4$SampleID ~ "Cluster 2",
#                             T ~ "Cluster 1")) %>%
#   ggplot()+
#   aes(y = n, x = Group2, color = Group2)+
#   geom_boxplot(outliers = F )+
#   geom_jitter(position = position_jitter(w = 0.1, h = 0), size = 3)+
#   theme_bw(base_size = 20) +
#   theme(legend.position = "none",
#         panel.background = element_rect(fill= 'white'),
#         axis.text.y=element_text(color = 'black', size = 20),
#         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         text=element_text(family="Helvetica")) +
#   ylab("Proteins Identified (n)")+
#   xlab("")+
#   stat_pvalue_manual(stat.test2)+
#   scale_y_continuous(limits = c(2000,3400))
# ggsave("QC_proteins.png", width = 8, height = 4.7)

############

fasta <- importFasta(fasta.file = "./2024-03-11-decoys-contam-uniprotkb_proteome_UP000006548_2024_03_12.fasta")
fasta2 <- digestFasta2(fasta) %>%
  filter(!grepl("rev_", accession))  %>%
  group_by(accession) %>%
  tally() %>%
  rename(Protein = accession)

iBAQ <- x_long %>%
  full_join(., meta) %>%
  filter(!is.na(Intensity)) %>%
  filter(SampleID %in% meta3$SampleID) %>%
  left_join(., fasta2) %>%
  filter(`Prev AA` %in% c("K", "R", "-", "M")) %>%
  dplyr::select(SampleID, Protein, `Protein ID`, Gene, Intensity, `Modified Sequence`,n) %>%
  group_by(SampleID, Protein) %>%
  mutate(iBAQ = sum(Intensity)/n) %>%
  distinct(SampleID, iBAQ, Protein, `Protein ID`, Gene) %>%
  filter(grepl("Ctrl", SampleID)) %>%
  group_by(SampleID) %>%
  mutate(Total = sum(iBAQ)) %>%
  ungroup() %>%
  mutate(riBAQ = iBAQ/Total)  %>%
  group_by(Protein) %>%
  mutate(Copy_Number = mean(riBAQ)*25800000000) %>%
  distinct(Protein, Copy_Number) %>%
  ungroup() %>%
  arrange(-Copy_Number) %>%
  mutate(Rank1 = row_number(-Copy_Number)) %>%
  left_join(., meta_prot) %>%
  distinct(Rank1 ,.keep_all = T)



markers <- read.csv(file = "042424_Protein_Complexes.csv")


gene_list2 <- gene_list %>%
  mutate(TAIR = gsub("\\..*","", TAIR )) %>%
  right_join(.,markers)


marker_prot <- meta_prot %>%
  mutate(TAIR = gsub("\\..*","", TAIR )) %>%
  inner_join(., gene_list2)  %>%
  dplyr::select(name, ID, Type, Gene) %>%
  full_join(., iBAQ) %>%
  distinct(Rank1, .keep_all = T) %>%
  filter(!is.na(Rank1))


marker_prot2 <- marker_prot %>%
  filter(!is.na(Type))

## Figure 3B
# marker_prot %>%
#   ggplot()+
#   aes(x = Rank1, y = Copy_Number,
#       color = Type)+
#   geom_point(size =ifelse(!is.na(marker_prot$Type),0,1))+
#   geom_point(size =ifelse(!is.na(marker_prot$Type),3,NA))+
#   geom_label_repel(data = marker_prot2, aes(label = Gene),
#                    size = 4,
#                    nudge_x = 40,
#                    label.padding = 0.1,
#                    force = 5,
#                    force_pull = 0.1)+
#   theme_bw(base_size = 20) +
#   theme(panel.background = element_rect(fill= 'white'),
#         legend.position = "none",
#         axis.text.y=element_text(color = 'black', size = 20),
#         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line())+
# 
#   xlab("Rank")+
#   ylab("log(Copy Number)")+
#   yscale("log10", .format = TRUE)+
# ggsave("LHC_RuBisco.png", height = 4.7, width = 8)

#####################################



comparison <- read.csv("Hildebrandt.csv")


marker_prot <- iBAQ %>%
  mutate(TAIR = gsub("\\..*","", TAIR )) %>%
  inner_join(., comparison)



### Supplementary Figure S9B
# marker_prot %>%
#   mutate(Mol_Fraction2 = Mol_Fraction2*25800000000) %>%
#   ggscatter(x = "Copy_Number",
#             y = "Mol_Fraction2",
#             add = "reg.line",
#             color = NA,
#             conf.int = T,
#             add.params = list(color = "black")
#   )+
#   stat_cor(method="pearson")+
#   yscale("log10", .format = TRUE)+
#   xscale("log10", .format = TRUE)+
#   geom_pointdensity(size = 2,
#                     show.legend= F)+
#   scale_color_viridis()+
#   xlab("log(Copy Number)\n Present study")+
#   ylab("log(Copy Number)\n Heinemann et al., 2021")+
# ggsave("hildebrandt.png", height = 4.7, width = 8)





##########################


iBAQ2 <- x_long %>%
  full_join(., meta) %>%
  filter(!is.na(Intensity)) %>%
  filter(SampleID %in% meta3$SampleID) %>%
  left_join(., fasta2) %>%
  filter(`Prev AA` %in% c("K", "R", "-", "M")) %>%
  dplyr::select(SampleID, Protein, `Protein ID`, Gene, Intensity, `Modified Sequence`,n) %>%
  group_by(SampleID, Protein) %>%
  mutate(iBAQ = sum(Intensity)/n) %>%
  distinct(SampleID, iBAQ, Protein, `Protein ID`, Gene) %>%
  filter(grepl("Ctrl", SampleID)) %>%
  group_by(SampleID) %>%
  mutate(Total = sum(iBAQ)) %>%
  ungroup() %>%
  mutate(riBAQ = iBAQ/Total)  %>%
  group_by(Protein) %>%
  mutate(Copy_Number = riBAQ*25800000000) %>%
  distinct(SampleID, Protein, Copy_Number) %>%
  left_join(., meta_prot) %>%
  filter(grepl("LHCA", Gene)) %>%
  mutate(Type = case_when(Gene == "LHCA5" | Gene == "LHCA6" ~ "LHCA5/6",
                          TRUE ~ "LHCA1-4"))






# ## Supplementary Figure S8
# iBAQ2 %>%
#   ggplot()+
#   aes(x = Type, y = Copy_Number)+
#   geom_boxplot(outliers = F )+
#   geom_jitter(width = 0.1, height = 0)+
#   yscale("log10", .format = TRUE)+
#   theme_bw(base_size = 20)+
#   theme(#legend.position = "none",
#     panel.grid.major = element_blank())+
#   xlab(NULL)+
#   ylab("log(Copy Number)")+
#   ggsave("LHC_comparison.png", height = 4.7, width = 8)



presentStudy_all <- iBAQ %>%
  mutate(name = gsub("\\..*","", name )) %>%
  pull(`name`)

Hild_all <- comparison %>%
  pull(`TAIR`)



### Supplementary Figure S9A
# venn.diagram(x = list(presentStudy_all, Hild_all),
#              category.names = c("Protein" , "mRNA"),
#              height = 4.7 ,
#              width = 6 ,
#              units = "in",
#              resolution = 300,
#              disable.logging = T,
#              sigdigs = 0,
#              col=c("darkblue", 'lightblue'),
#              fill = c(alpha("darkblue",0.3), alpha('lightblue',0.4)),
#              cex = 1,
#              ext.text = F,
#              cat.cex = 0,
#              ext.line.lwd = 0.000001,
#              filename = "venn_compare.png",
#              output = TRUE)

# Supp_Table <- iBAQ %>%
#   mutate(TAIR = gsub("\\..*","", TAIR )) %>%
#   mutate(Hein = case_when(TAIR %in% Hild_all ~ "Present",
#                           TRUE ~ "Not Present"))
# 
# 
# ### Supplementary Table
# write.csv(Supp_Table, file = "Supp_Table.csv")




###### Figure 3A
# protein_long2 %>%
#   distinct(name, SampleID) %>%
#   mutate(Group1 = case_when(grepl("Ctrl_Day0", SampleID) ~ "Ctrl-50%",
#                             grepl("Drgt_Day0", SampleID)~ "WD-50%",
#                             grepl("Ctrl_Day7", SampleID) ~ "Ctrl-30%",
#                             grepl("Drgt_Day7", SampleID) ~ "WD-30%"
#   )) %>%
#   group_by(SampleID) %>%
#   add_count(name = "n") %>%
#   distinct(SampleID, n, Group1) %>%
#   ungroup()  %>% 
#   mutate(Group1 = factor(Group1, c("Ctrl-50%","WD-50%",
#                                    "Ctrl-30%","WD-30%"))) %>%
#   ggplot()+
#   aes(y = n, x = Group1, color = Group1)+
#   geom_boxplot(outliers = F )+
#   geom_jitter(position = position_jitter(w = 0.1, h = 0), size = 3)+
#   theme_bw(base_size = 20) +
#   theme(legend.position = "none",
#         panel.background = element_rect(fill= 'white'),
#         axis.text.y=element_text(color = 'black', size = 22),
#         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 22),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         text=element_text(family="Helvetica")) +
#   ylab("Proteins Identified (n)")+
#   xlab("")+
#   scale_y_continuous(limits = c(2000,3300))+
#   #  geom_vline(xintercept = 2.5)+
#   geom_hline(yintercept = 2873, linetype = "dashed")+
#   ggsave("QC_proteins_postcluster.png", width = 8, height = 4.7)




meta_final <- inner_join(Sample_Metadata, meta3) %>%
  arrange(SampleID)



protein_wide <- protein_long2 %>%
  spread(SampleID, Intensity) %>% 
  column_to_rownames(var= "name")




###############  Perform differential abundance analysis


protoP2 <- protein_wide %>%
  rownames_to_column(var = "name")



protoP2[,2:88] <- 2^(protoP2[,2:88])



protoP3 <- full_join(protoP2, meta_prot)

LFQ_columns <- grep("_",colnames(protoP3)) # get LFQ column numbers
experimental_design <- meta_final %>%
  mutate(label = SampleID) %>%
  mutate(condition = Group1) %>%
  group_by(condition) %>%
  mutate(replicate = dense_rank(desc(label))) %>%
  ungroup() %>%
  dplyr::select(label, condition, replicate) %>%
  arrange(desc(label))

data_se <- make_se(protoP3, LFQ_columns, experimental_design)

data_diff <- DEP::test_diff(data_se, type = "manual",
                            test = c("WD.50._vs_Ctrl.50.",
                                     "WD.30._vs_Ctrl.30."))

data_diff@elementMetadata@listData <- data_diff@elementMetadata@listData %>%
  as.data.frame() %>%
  mutate_at(vars(matches("p.val")), ~ p.adjust(.x, method = "BH")) %>%
  dplyr::select(-contains("p.adj"),
                -contains("CI.L"),
                -contains("CI.R")) %>%
  rename_with(.fn = ~ str_replace(., "p.val", "p.adj"),
              .cols = contains("p.val")) %>%
  as.list() 

dep <- add_rejections(data_diff, alpha = 0.05, lfc = 0.5)


data_results <- as.data.frame(dep@elementMetadata@listData) 

### Supplementary Table 1
 # write.csv(data_results, file= "Supp_Table1.csv")



################## Perform gene ontological enrichment based on DA data
Drgt_Down <- data_results %>%
  filter(WD.30._vs_Ctrl.30._diff <= -0.5) %>%
  filter(WD.30._vs_Ctrl.30._p.adj <= 0.05) %>%
  arrange(WD.30._vs_Ctrl.30._diff) %>%
  pull(name)

Drgt_Up <- data_results %>%
  filter(WD.30._vs_Ctrl.30._diff >= 0.5) %>%
  filter(WD.30._vs_Ctrl.30._p.adj <= 0.05) %>%
  arrange(WD.30._vs_Ctrl.30._diff) %>%
  pull(TAIR)


###################### First protein descreased in abundance

gostres <- gost(query = list(#"Drgt_Up" = Drgt_Up
  "Drgt_Down" = Drgt_Down
),
organism = "athaliana", ordered_query = F, 
multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
measure_underrepresentation = FALSE, evcodes = TRUE, 
user_threshold = 0.05, correction_method = "fdr", 
domain_scope = "annotated", custom_bg = bg, 
numeric_ns = "", as_short_link = FALSE)

p<- gostplot(gostres, capped = F, interactive = F)

downfilter <- c("GO:0003735", 
                "GO:0005840", "KEGG:03010",
                "GO:0006412")

##### Supplementary Table S3
#TableS3 <- gostres$result %>%
#  dplyr::select(-parents)
#write.csv(TableS3, file =  "Table_S3.csv")

terms <- gostres$result %>%
  filter(term_id %in% downfilter )


###### Figure 4D
#publish_gostplot(p,highlight_terms = terms$term_id,
#                 height = 6, width = 7.2, filename = "gostplot_down.png" )


term_down <- gostres$result %>%
  group_by(query) %>%
  separate_rows(intersection) %>%
  dplyr::distinct(term_name, intersection) %>%
  distinct(intersection, .keep_all = T) %>%
  inner_join(.,meta_prot, by = c("intersection" = "TAIR"))

############### Proteins increased in abundance
gostres <- gost(query = list("Drgt_Up" = Drgt_Up
                           #  "Drgt_Down" = Drgt_Down
                             ),
                organism = "athaliana", ordered_query = F, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = bg, 
                numeric_ns = "", as_short_link = FALSE)

p<- gostplot(gostres, capped = F, interactive = F)

upfilter <- c("GO:0009414", 
              "GO:0000325", "GO:1901700", "GO:0015250",
              "GO:0000325",
              "GO:0009737")

#### Supplementary Table 2
#TableS2 <- gostres$result %>%
 # dplyr::select(-parents)
#write.csv(TableS2, file =  "Table_S2.csv")


terms <- gostres$result %>%
  filter(term_id %in% upfilter )


######## Figure 4C
#publish_gostplot(p,highlight_terms = terms$term_id,
#                 height = 8.2, width = 7.4, filename = "gostplot_up.png" )

############### PCA, post-removal of fragmented protoplasts

protoP_wide <- protein_long2 %>%
  group_by(name) %>% 
  add_count(name = "n") %>%
  filter(n >= 87*0.9) %>% 
  dplyr::select(-n) %>%
  spread(SampleID, Intensity) %>% 
  column_to_rownames(var= "name")

protoP_impute <-  DreamAI(protoP_wide, k = 10, maxiter_MF = 10, ntree = 100, 
                        maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2), 
                        gamma_ADMIN = NA, gamma = 50, CV = FALSE, 
                        fillmethod = "row_mean", maxiter_RegImpute = 10, 
                        conv_nrmse = 1e-06, iter_SpectroFM = 40, method = "KNN", 
                        out = c("KNN")) 

protoP_impute <- as.data.frame(protoP_impute$Ensemble)

pca_j<- PCAtools::pca(protoP_impute, scale = TRUE, center = T)

Group <- meta_final$Group1

############ Figure 4A
# PCAtools::biplot(pca_j, x = "PC1", y =  "PC2", lab = Group,
#                  showLoadings = F,
#                  ntopLoadings = 5,
#                  labSize = 0, drawConnectors = FALSE,
#                  legendPosition = "right",
#                  encircle = TRUE
# )+
#   theme_bw(base_size = 12)+
#   ggsave("PCA_postcluster.png", height = 6, width = 7)
# 

################ Volcano plot, Figure 4B and Supplementary Figure S10
term_up <- gostres$result %>%
  filter(term_id == "GO:0009414") %>%
  separate_rows(intersection) %>%
  dplyr::select(query, intersection)  %>%
  inner_join(.,meta_prot, by = c("intersection" = "TAIR")) %>%
  filter(query == "Drgt_Up")

# 
# cs_x <- data_results %>%
#   dplyr::select(Protein, name ,Gene, contains("WD.30")) %>%
#   mutate(adj.pvalue = WD.30._vs_Ctrl.30._p.adj) %>%
#   mutate(log2FC =  WD.30._vs_Ctrl.30._diff) %>%
#   dplyr::select(Protein,Gene, adj.pvalue, log2FC) %>%
#   filter(!is.na(adj.pvalue))  %>%
#   mutate(Color = case_when(Protein %in% term_up$Protein ~ "#984ea3",
#                            TRUE ~ "NaN")) %>%
#   mutate(Color = case_when(Color == "NaN" & adj.pvalue <= 0.05 & !between(log2FC, -0.5, 0.5) ~ "red2",
#                            Color == "NaN" & adj.pvalue > 0.05 & !between(log2FC, -0.5, 0.5) ~ "forestgreen",
#                            Color == "NaN" & adj.pvalue > 0.05 ~ "grey30",
#                            Color == "NaN" & between(log2FC, -0.5, 0.5)~ "grey30",
#                            Color == "NaN" & adj.pvalue <= 0.05 ~ "royalblue1",
#                            TRUE ~ Color)) %>%
#   arrange(desc(Color)) 
# 
# 
# rownames(cs_x) <- cs_x$Protein 
# 
# keyvals <- cs_x$Color
# 
# keyvals[is.na(keyvals)] <- 'NA'
# names(keyvals)[keyvals == 'forestgreen'] <- 'Log2FC'
# names(keyvals)[keyvals == 'red2'] <- 'p-value and log2FC'
# names(keyvals)[keyvals == '#984ea3'] <- 'Response to\nWater Deprivation'
# names(keyvals)[keyvals == 'royalblue1'] <- 'p-value'
# names(keyvals)[keyvals == 'grey30'] <- 'NS'
# names(keyvals)[keyvals == 'NA'] <- ' '
# 
# 
# 
# png('Arab_volc_plot2.png', width = 7.2, height = 6,
#     units = "in",
#     res =300)
# EnhancedVolcano(cs_x,
#                 lab = cs_x$Gene,
#                 selectLab = ifelse(cs_x$Color == "#984ea3", cs_x$Gene,NA),
#                 labSize = 3.5,
#                 drawConnectors = T,
#                 title = NULL,
#                 subtitle = NULL,
#                 caption = NULL,
#                 boxedLabels = T,
#                 arrowheads = F,
#                 x = 'log2FC',
#                 y = 'adj.pvalue',
#                 pCutoff = 0.05,
#                 FCcutoff = 0.5,
#                 axisLabSize = 14,
#                 xlim = c(-2.5,3.7),
#                 ylim = c(0,17),
#                 colCustom = keyvals,
#                 colAlpha = c(ifelse(cs_x$Color != '#984ea3', 0.3, 1)),
#                 pointSize = c(ifelse(cs_x$Color != "grey30", 3, 3)),
#                 legendPosition = 'bottom',
#                 legendLabSize = 12,
#                 legendIconSize = 4.0)+
#   theme(panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.line = element_line())
# dev.off()




######## supplementary Figure S10
# cs_x <- data_results %>%
#   dplyr::select(name,Protein, contains("WD.50")) %>%
#   mutate(adj.pvalue = WD.50._vs_Ctrl.50._p.adj) %>%
#   mutate(log2FC =  WD.50._vs_Ctrl.50._diff) %>%
#   dplyr::select(Protein, adj.pvalue, log2FC) %>%
#   filter(!is.na(adj.pvalue))  %>%
#   mutate(Color = case_when(TRUE ~ "NaN")) %>%
#   mutate(Color = case_when(Color == "NaN" & adj.pvalue <= 0.05 & !between(log2FC, -0.5, 0.5) ~ "red2",
#                            Color == "NaN" & adj.pvalue > 0.05 & !between(log2FC, -0.5, 0.5) ~ "forestgreen",
#                            Color == "NaN" & adj.pvalue > 0.05 ~ "grey30",
#                            Color == "NaN" & between(log2FC, -0.5, 0.5)~ "grey30",
#                            Color == "NaN" & adj.pvalue <= 0.05 ~ "royalblue1",
#                            TRUE ~ Color)) %>%
#   arrange(desc(Color)) 
# 
# 
# rownames(cs_x) <- cs_x$Protein 
# 
# keyvals <- cs_x$Color
# 
# keyvals[is.na(keyvals)] <- 'NA'
# names(keyvals)[keyvals == 'forestgreen'] <- 'Log2FC'
# names(keyvals)[keyvals == 'red2'] <- 'p-value and log2FC'
# names(keyvals)[keyvals == 'royalblue1'] <- 'p-value'
# names(keyvals)[keyvals == 'grey30'] <- 'NS'
# names(keyvals)[keyvals == 'NA'] <- ' '
# 
# 
# 
# png('Arab_volc_plot_moderateDr.png', width = 7.4, height = 8.2,
#     units = "in",
#     res =300)
# EnhancedVolcano(cs_x,
#                 lab = NA,
#                 labSize = 3.5,
#                 drawConnectors = T,
#                 title = NULL,
#                 subtitle = NULL,
#                 caption = NULL,
#                 boxedLabels = T,
#                 arrowheads = F,
#                 x = 'log2FC',
#                 y = 'adj.pvalue',
#                 pCutoff = 0.05,
#                 FCcutoff = 0.5,
#                 axisLabSize = 14,
#                 xlim = c(-2.5,3.7),
#                 ylim = c(0,5),
#                 colCustom = keyvals,
#                 colAlpha = c(ifelse(cs_x$Color != '#984ea3', 0.3, 1)),
#                 pointSize = c(ifelse(cs_x$Color != "grey30", 3, 3)),
#                 legendPosition = 'bottom',
#                 legendLabSize = 12,
#                 legendIconSize = 4.0)+
#   theme(panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.line = element_line())
# dev.off()








####################################### Functional networks


final_table <- protein_wide %>%
  as.data.frame() %>%
  mutate(Gene1 = rowSums(is.na(.))) %>%
  filter(Gene1 <= max(Gene1)*0.1) %>% 
  dplyr::select(-Gene1) %>%
  as.matrix() %>%
  DreamAI(., k = 10, maxiter_MF = 10, ntree = 100, 
          maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2), 
          gamma_ADMIN = NA, gamma = 50, CV = FALSE, 
          fillmethod = "row_mean", maxiter_RegImpute = 10, 
          conv_nrmse = 1e-06, iter_SpectroFM = 40, method = "KNN", 
          out = c("KNN"))

final_table <- final_table$KNN


VB <- final_table


VB_cor <- fast_cor(t(VB), method = "pearson")

#### 19 is optimal, VII

clustVB_summary <- Clust_compare(VB,
                                 clusterNumbers= c(18:20),
                                 nameAlgorithm = c("mclust"
                                 ),
                                 models = "VII",
                                 shrinkage_values = c(0.01,0))


clustVB_summary <- clustGO(clustVB_summary, organism = "athaliana", Bfactor = 0.5,
                            S_type = "Summary",
                            bg = bg)
############

clust_order <- clustVB_summary[[4]] %>%
  mutate(term_name = case_when(Cluster == 15 ~ "ATPase activity",
                   TRUE ~ term_name))

contam_filter <- clust_order %>%
  group_by(Cluster) %>%
  add_count(name = "Total_Genes") %>%
  mutate(Unnannotated = sum(is.na(term_id))) %>%
  mutate(Proportion = Unnannotated/Total_Genes) %>%
  filter(Proportion > 0.9 | is.na(Cluster)|
           Total_Genes == 1 ) %>%
  distinct(Cluster,Gene) %>%
  distinct(Cluster)

contam_filter <- contam_filter %>%
  left_join(.,clust_order)


VB <- VB %>%
  as.data.frame  %>%
  filter(!rownames(.) %in% contam_filter$Gene) %>%
  as.matrix()

clust_order<- clust_order %>%
  filter(!Gene %in% contam_filter$Gene)


clust_names <- clust_order %>%
  dplyr::select(Cluster,Gene, uncertainty) %>%
  mutate(Cluster = as.character(Cluster)) %>%
  left_join(., clustVB_summary[[3]], by = c("Cluster" = "GO_Cluster")) %>%
  distinct(Cluster, term_name) %>%
  group_by(term_name) %>%
  mutate(term_name = paste0(str_extract(term_name,"^([^,])+"), " "),
         term_name = make.unique(term_name, sep = "_")) %>%
  arrange(as.numeric(Cluster)) %>%
  ungroup() 

clust_order <- clust_order %>%
  dplyr::select(Cluster,Gene, uncertainty, term_name) %>%
  mutate(Cluster = as.character(Cluster)) %>%
  arrange(Gene)


VB <- VB[ order(row.names(VB)), ]


All_cors <- VB_cor[[1]]


######## Cluster arrangment manually entered here but determined by
###### hclust with Euclidean distance and method = "complete"
###### on cluster correlation averages

clust_order_f <- All_cors %>%
  inner_join(., clust_order, by = c("Gene1" = "Gene")) %>%
  rename(Cluster1 = Cluster) %>%
  dplyr::select(Gene1, Gene2, cor, p, FDR, Cluster1) %>%
  inner_join(., clust_order, by = c("Gene2" = "Gene")) %>%
  rename(Cluster2 = Cluster) %>%
  filter(Cluster1 != Cluster2) %>%
  group_by(Cluster1, Cluster2) %>%
  mutate(Med1 = median(cor)) %>%
  distinct(Cluster1, Cluster2, Med1) %>%
  ungroup()%>%
  group_by(Cluster2) %>%
  slice_max(., order_by = Med1, n = 1, with_ties = F) %>%
  arrange(desc(Med1)) %>%
  dplyr::select(Cluster2) %>%
  rename(Cluster = Cluster2) %>%
  mutate(Cluster = factor(Cluster,
                          levels = c("16", "11","7",
                                     "12", "2", "14",
                                     "18", "5", "17",
                                     "8","19",  "9","13",
                                     "4", "6", "1", "3",
                                     "15"))) %>%
  ungroup() %>%
  arrange(Cluster)

clust_order2 <- clust_order %>%
  ungroup() %>%
  arrange(desc(factor(Cluster,
                      levels = rev(clust_order_f$Cluster))),
          uncertainty)

VB_cor_mat <- VB_cor[[2]]

include_list <- clust_order2$Gene

VB_cor_mat <- VB_cor_mat[include_list,include_list]


VB_cor_mat <- rotate(reorder_mat(VB_cor_mat,
                                 order = clust_order2$Gene))


clust_order2 <- clust_order2 %>%
  mutate(DifferentialA = case_when(Gene %in% Drgt_Up ~ "Up",
                                   Gene %in% Drgt_Down ~ "Down",
                                   TRUE ~ "NoDiff"))

index1 <- clust_order2 %>%
  dplyr::select(DifferentialA)

index2 <- clust_order2 %>%
  dplyr::select(Cluster)

chr_level = as.factor(unique(index1$DifferentialA))
chr_level2 = as.factor(unique(index2$Cluster))


set.seed(290)
colrs = qualpal(length(unique(clust_order2$DifferentialA)),
                colorspace = "pretty_dark")
colrs <- colrs$hex

set.seed(290)
colrs2 = qualpal(length(unique(clust_order2$Cluster)),
                 colorspace = "pretty")
colrs2 <- colrs2$hex

names(colrs) <- unique(clust_order2$DifferentialA)
colrs[names(colrs) == "NoDiff"] <- NA
colrs[names(colrs) == "Down"] <- "#1E1E52"
colrs[names(colrs)== "Up"] <- "#36872a"
names(colrs2) <- unique(clust_order2$Cluster)


column_ha = HeatmapAnnotation(Module = rev(index2$Cluster),
                              AnnotationG = rev(index1$DifferentialA),
                              col = list(AnnotationG = colrs,
                                         Module = colrs2
                              ),
                              show_legend = F,
                              show_annotation_name = F,
                              border =F,
                              simple_anno_size = unit(4, "cm"))

row_ha = rowAnnotation(AnnotationG = index1$DifferentialA,
                       Module = index2$Cluster,
                       col = list(AnnotationG = colrs,
                                  Module = colrs2),
                       show_annotation_name = F,
                       border = F,
                       simple_anno_size = unit(4, "cm"),
                       simple_anno_size_adjust	= F,
                       show_legend = F,
                       #  gp = gpar(col  = list(AnnotationG = colrs,
                       #                            Module = colrs2)),
                       annotation_legend_param = list( AnnotationG = list(
                         ncol = 2,
                         title = "Differentially Abundant",
                         title_position = "topcenter",
                         at = chr_level
                       ),
                       Module = list(
                         ncol = 2,
                         title = "Cluster",
                         title_position = "topcenter",
                         at = chr_level2
                       )))


###### Figure 5A
# png("Arab_heatmap_updated.png", width = 2500, height = 2500)
# Heatmap(VB_cor_mat, cluster_rows = F,
#         cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         use_raster = F,
#         show_heatmap_legend = F,
#         # heatmap_legend_param = list(title = "Pearson R"),
#         top_annotation =column_ha,
#         right_annotation =row_ha
# )
# dev.off()



leaf_table <- clustVB_summary2[[3]]

cluster_size <- clust_order2 %>%
  group_by(Cluster) %>%
  add_count(name = "n1") %>%
  filter(!is.na(term_name)) %>%
  add_count(name = "n2") %>%
  distinct(Cluster, n1, n2)


# Figure 5B
# leaf_table %>%
#   mutate(term_name = case_when(GO_Cluster == 15 ~ "ATPase activity",
#                                TRUE ~ term_name)) %>%
#   rename(Cluster = GO_Cluster) %>%
#   left_join(., cluster_size) %>%
#   # mutate(`Percent Annotated` = intersection_size/n*100) %>%
#   #mutate(precision = round(precision, digits = 2)) %>%
#   mutate(p_value = round(-log10(p_value),digits = 2)) %>%
#   dplyr::select(Cluster, n1, intersection_size, term_name,
#                 p_value) %>%
#   rename(`Cluster Size` = n1) %>%
#   rename(`Gene Ontological Term` = term_name) %>%
#   rename(`Proteins Annotated` = intersection_size) %>%
#   rename(`-log10(p-value)` = p_value) %>%
#   arrange(desc(factor(Cluster,
#                       levels = rev(clust_order_f$Cluster)))) %>%
#   gt() %>%
#   data_color(columns =  `-log10(p-value)`,
#              palette = "viridis") %>%
#   data_color(columns = `Cluster`,
#              palette = colrs2,
#              ordered = T) %>%
#   tab_style(
#     style = cell_borders(sides = "all", 
#                          color = "#000000",
#                          style = "solid",
#                          weight = px(1)),
#     locations = cells_body()) %>%
#   cols_align(
#     align = "center",
#     columns = everything()
#   ) %>%
#   gtsave("VB_Cell_Cluster_table2.png",
#          expand = 30)


clust_order3 <- colrs2 %>%
  as.data.frame() %>%
  rownames_to_column(var = "Cluster") %>%
  full_join(., clust_order2)

## Supplementary Table S5
#write.csv(clust_order3, file = "TableS5.csv")


########################


GO_Figure <- clustVB_summary[["VII_clusters_19_ALL_GO_terms"]] %>% 
  filter(GO_Cluster %in% c("11","7","12")) %>%
  filter(term_id %in% c("GO:0022626", "GO:0009532",
                        "GO:0009547", "GO:0042644",
                        "GO:0010467","GO:0031981", "GO:0006413")) %>%
  #   distinct(p_value, .keep_all = TRUE) %>%
  # group_by(GO_Cluster) %>% 
  #  slice_min(p_value, n= 3, with_ties = F) %>%
  mutate(GO_Cluster = factor(GO_Cluster,
                             levels = c("12",
                                        "7",
                                        "11"))) %>%
  ungroup() %>%
  arrange(GO_Cluster)



######## Figure 5C
# png("GO_Figure.png", height = 2, width = 5,
#     units = "in",
#     res = 300)
# GO_Figure %>%
#   ggplot()+
#   aes(y = GO_Cluster, x = term_name,
#       color = -log10(`p_value`))+
#   geom_point(size = 8)+
#   scale_color_viridis_c()+
#   ylab("")+
#   xlab(NULL)+
#   theme_bw(base_size = 8)+
#   theme(legend.position = "bottom")
# # ggplot2::facet_grid(source~.,
# #                      space = "free_y",
# #                     scales = "free_y", switch = "y")
# dev.off()











########################

### Supplementary Table S4
# TableS4 <- clustVB_summary[[2]] %>%
#   rename(Cluster = GO_Cluster) %>%
#   arrange(desc(factor(Cluster,
#                       levels = rev(clust_order_f$Cluster))),
#           p_value)
# 
# TableS4 <- colrs2 %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "Cluster") %>%
#   full_join(., TableS4)
# 
# write.csv(TableS4, file = "TableS4.csv")

######## Control 30% correlations

Ctrl_samples <- meta_final %>%
  filter(grepl("Ctrl-30%", Group1)) %>%
  pull(SampleID)

Ctrl_prots <- final_table[include_list,Ctrl_samples]

Cor1 <- fast_cor(t(Ctrl_prots), method = "pearson")

Ctrl_cor1 <- Cor1[[1]] %>%
  dplyr::select(-df)

###### WD 30% correlations
Drgt_samples <- meta_final %>%
  filter(grepl("WD-30%", Group1)) %>%
  pull(SampleID)

Drgt_prots <- final_table[include_list,Drgt_samples]

Cor2 <- fast_cor(t(Drgt_prots), method = "pearson")

Drgt_cor <- Cor2[[1]] %>%
  rename(Drgt_p = p) %>%
  rename(Drgt_FDR = FDR) %>%
  rename(Drgt_cor = cor) %>%
  dplyr::select(-df)


Dist_mat <- Cor2[[2]] - Cor1[[2]]



Dist_compare  <- rotate(reorder_mat(Dist_mat,
                                 order = clust_order2$Gene))


col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))



## Supplementary Figure S12 
# png("DistanceMatrix_heatmap.png", width = 3200, height = 3200)
# Heatmap(Dist_compare, cluster_rows = F,
#         cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         use_raster = F,
#         show_heatmap_legend = F,
#         top_annotation =column_ha,
#         right_annotation =row_ha,
#         col = col_fun
# )
# dev.off()


Drgt_dist <- Dist_compare %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene1") %>%
  pivot_longer(-Gene1,
               names_to = "Gene2",
               values_to = "Pearson_Distance")








Ctrl_mat  <- rotate(reorder_mat(Cor1[[2]],
                                order = clust_order2$Gene))


## Supplementary Figure S12 
# png("Control_Matrix_heatmap.png", width = 3200, height = 3200)
# Heatmap(Ctrl_mat, cluster_rows = F,
#         cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         use_raster = F,
#         show_heatmap_legend = F,
#         heatmap_legend_param = list(title = "Pearson R"),
#         top_annotation =column_ha,
#         right_annotation =row_ha,
# )
# dev.off()



Drgt_mat  <- rotate(reorder_mat(Cor2[[2]],
                                order = clust_order2$Gene))


## Supplementary Figure S12 
# png("Drought_Matrix_heatmap.png", width = 3200, height = 3200)
# Heatmap(Drgt_mat, cluster_rows = F,
#         cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         use_raster = F,
#         show_heatmap_legend = F,
#         heatmap_legend_param = list(title = "Pearson R"),
#         top_annotation =column_ha,
#         right_annotation =row_ha,
# )
# dev.off()




###################### Correlation Average matrices

  Full_compare <- full_join(Ctrl_cor1, Drgt_cor) %>%
    filter(cor != 1) %>%
    dplyr::select(-p, -Drgt_p) %>%
    mutate(Dist_Cor = Drgt_cor-cor) %>%
    left_join(., clust_order2, by = c("Gene1" = "Gene")) %>%
    rename(Cluster1 = Cluster) %>%
    left_join(., clust_order2, by = c("Gene2" = "Gene")) %>%
    rename(Cluster2 = Cluster) %>%
    dplyr::select(Gene1, Gene2, cor, FDR, Drgt_cor,
                  Drgt_FDR, Dist_Cor, Cluster1, Cluster2) %>%
  mutate(cs2=pmap_chr(list(`Cluster1`,`Cluster2`), ~paste(sort(c(...)), collapse = "_"))) %>%
  distinct(Cluster1, Cluster2, cor, Drgt_cor, .keep_all = T) 
  
  Full_mat <- Full_compare %>%
    group_by(Cluster1, Cluster2) %>%
    summarize(Ctrl_Avg = mean(cor),
              Drgt_Avg = mean(Drgt_cor),
              Dist_Avg = mean(Dist_Cor))
  
 Ctrl_Avg<-  Full_mat %>%
    dplyr::select(Cluster1, Cluster2, Ctrl_Avg) %>%
    pivot_wider(names_from  = "Cluster1",
                values_from= "Ctrl_Avg") %>%
    column_to_rownames(var = "Cluster2") %>%
    as.matrix()
 
 
Ctrl_Avg <-  rotate(reorder_mat(Ctrl_Avg,
                    order = as.character(clust_order_f$Cluster)))




column_ha3 = HeatmapAnnotation(Module = rev(clust_order_f$Cluster),
                               foo = anno_empty(border = F),
                               col = list(Module = colrs2),
                               show_legend = F,
                               show_annotation_name = F,
                               border =T,
                               gap = unit(c(5,10), "mm"),
                               simple_anno_size = unit(2, "cm"))

row_ha3 = rowAnnotation(foo = anno_empty(border = F),
                        Module = clust_order_f$Cluster,
                        col = list(Module = colrs2),
                        show_annotation_name = F,
                        border = T,
                        simple_anno_size = unit(2, "cm"),
                        gap = unit(c(5,10), "mm"),
                        simple_anno_size_adjust	= F,
                        show_legend = F
)


col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))


###### Figure 6A
# png("Ctrl_Avg.png", width = 1200, height = 1200)
# Heatmap(Ctrl_Avg, cluster_rows = F,
#         cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         use_raster = F,
#         show_heatmap_legend = F,
#         col = col_fun,
#        top_annotation =column_ha3,
#         right_annotation =row_ha3,
# )
# dev.off()

 

Drgt_Avg<-  Full_mat %>%
  dplyr::select(Cluster1, Cluster2, Drgt_Avg) %>%
  pivot_wider(names_from  = "Cluster1",
              values_from= "Drgt_Avg") %>%
  column_to_rownames(var = "Cluster2") %>%
  as.matrix()

Drgt_Avg <-  rotate(reorder_mat(Drgt_Avg,
                                order = as.character(clust_order_f$Cluster)))



###### Figure 6A
# png("Drgt_Avg.png", width = 1200, height = 1200)
# Heatmap(Drgt_Avg, cluster_rows = F,
#         cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         use_raster = F,
#         show_heatmap_legend = F,
#         col = col_fun,
#         top_annotation =column_ha3,
#         right_annotation =row_ha3,
# )
# dev.off()


Dist_Avg<-  Full_mat %>%
  dplyr::select(Cluster1, Cluster2, Dist_Avg) %>%
  pivot_wider(names_from  = "Cluster1",
              values_from= "Dist_Avg") %>%
  column_to_rownames(var = "Cluster2") %>%
  as.matrix()


Dist_Avg <-  rotate(reorder_mat(Dist_Avg,
                                order = as.character(clust_order_f$Cluster)))

###### Figure 6A
# png("Dist_Avg2.png", width = 1200, height = 1200)
# Heatmap(Dist_Avg, cluster_rows = F,
#         cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         use_raster = F,
#         show_heatmap_legend = F,
#         col = col_fun,
#         top_annotation =column_ha3,
#         right_annotation =row_ha3,
# )
# dev.off()

####### Distribution comparisons
Sig_subset <- Full_compare %>%
  pivot_longer(c(-Gene1, -Gene2, -Cluster1, -Cluster2, -cs2, -FDR,
                 -Drgt_FDR),
               names_to = "Type",
               values_to = "cor") %>%
  filter(cs2 == "6_6") %>%
  filter(Type != "Dist_Cor")

Sig_subset %>%
  group_by(Type) %>%
  summarize(Avg = mean(cor))


  
########## Figures 6B and 6C  
# Sig_subset %>%
#   ggplot()+
#   aes(x = cor, fill = Type) %>%
#   geom_histogram(position = "identity", alpha = 0.5, bins = 100)+
#   theme_bw(base_size = 12)+
#   theme(legend.position = "none",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.y=element_text(color = 'black', size = 12),
#     axis.text.x=element_text(color='black',size = 12))+
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) 
 # geom_vline(xintercept = 0.12, linetype = "dashed", color = "black")+
 # geom_vline(xintercept = 0.55, linetype = "dashed", color = "black")
 # ggsave("DeltaPearson_66.png", height = 5, width = 5)
  

################

log2FC <- inner_join(clust_order2, data_results, by = c("Gene" = "name"))

cs_x <- log2FC  %>%
  filter(Cluster %in% c("6", "11", "12")) %>%
  dplyr::select(Gene,Protein, contains("WD.30"), Cluster) %>%
  mutate(adj.pvalue = WD.30._vs_Ctrl.30._p.adj) %>%
  mutate(log2FC =  WD.30._vs_Ctrl.30._diff) %>%
  dplyr::select(Protein, adj.pvalue, log2FC, Cluster) %>%
  filter(!is.na(adj.pvalue))  %>%
  mutate(Color = case_when(Cluster == "6" ~ "grey30",
                           Cluster == "11"~ "#B66CCC",
                           Cluster == "12" ~ "#75B1C6",
                           TRUE ~ "NaN")) %>%
  mutate(Shape = case_when(Cluster == "6" ~ 15,
                           Cluster == "11"~ 16,
                           Cluster == "12" ~ 17)) %>%
  arrange(desc(Color)) 


rownames(cs_x) <- cs_x$Protein 

keyvals <- cs_x$Color

keyvals[is.na(keyvals)] <- 'NA'
names(keyvals)[keyvals == 'grey30'] <- 'Cluster 6'
names(keyvals)[keyvals == '#B66CCC'] <- 'Cluster 11'
names(keyvals)[keyvals == '#75B1C6'] <- 'Cluster 12'
names(keyvals)[keyvals == 'NA'] <- ' '

 ######## Figure 6D
# png('Cluster_Volc_plot.png', width = 7, height = 7,
#     units = "in",
#     res =300)
# EnhancedVolcano(cs_x,
#                 lab = NA,
#                 labSize = 3.5,
#                 drawConnectors = T,
#                 title = NULL,
#                 subtitle = NULL,
#                 caption = NULL,
#                 boxedLabels = T,
#                 arrowheads = F,
#                 x = 'log2FC',
#                 y = 'adj.pvalue',
#                 pCutoff = 0.05,
#                 FCcutoff = 0.5,
#                 axisLabSize = 14,
#                 xlim = c(-2,0),
#                 ylim = c(0,12),
#                 colCustom = keyvals,
#                 pointSize = c(ifelse(cs_x$Color != "grey30", 3, 3)),
#                 legendPosition = 'bottom',
#                 legendLabSize = 12,
#                 legendIconSize = 4.0)+
#   theme(panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.line = element_line())
# dev.off()



############################### Cluster 11 example

include_list <- clust_order2 %>%
  filter(Cluster == "11") %>%
  pull(Gene)

Mat_subset <- VB_cor_mat[include_list,include_list]


Mat_subset <- rotate(reorder_mat(Mat_subset,
                                 order = include_list))


########## Supplementary Figure S11
# png("Cluster11.png", width = 600, height = 600)
# Heatmap(Mat_subset, cluster_rows = F,
#         cluster_columns = F,
#         show_row_names = T,
#         show_column_names = T,
#         use_raster = F,
#         show_heatmap_legend = F,
#         col = col_fun,
# )
# dev.off()
