library(tidyverse)
library(SingleCellExperiment)
library(scp)
library(ggplot2)


mod_pep <- read_tsv("combined_modified_peptide.tsv") %>%
  setNames(., gsub("nanoPOTs_Protoplasts_", "", names(.))) %>%
  setNames(., gsub(" Intensity", "", names(.))) %>%
  setNames(., gsub("_Protoplasts_spintest_", "", names(.)))

meta <- mod_pep[,1:15]

meta2 <- meta %>%
  dplyr::select(`Modified Sequence`, Protein)


peptide <- mod_pep %>%
  filter(grepl("ARATH", Protein)) %>%
  dplyr::select(`Modified Sequence` | contains("_")) %>%
  column_to_rownames(var="Modified Sequence") %>%
  dplyr::select(-contains("10cell") | contains("3cell"))

x_long <- mod_pep %>%
  filter(grepl("ARATH", Protein)) %>%
  dplyr::select(`Modified Sequence` | contains("_")) %>%
  column_to_rownames(var="Modified Sequence") %>%
  dplyr::select(-contains("10cell") | contains("3cell")) %>%
  rownames_to_column(var="Modified Sequence") %>%
  pivot_longer(!`Modified Sequence`, names_to = "SampleID", values_to = "Intensity") %>%
  filter(Intensity  > 0) 

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


protein <- !is.na(protein)
peptide <- peptide != 0

sampleGroups <- data.frame(SampleID = colnames(protein)) %>%
  mutate(Group = "Protoplasts")

sce <- SingleCellExperiment(assays = List(protein),)
sce2 <- SingleCellExperiment(assays = List(peptide),)
sim <- QFeatures(experiments = List(Assay1 = sce,
                                    Assay2 = sce2))
sim$SampleType <- sampleGroups$Group

protein_Missing <- reportMissingValues(sim, "Assay1", by = sim$SampleType)%>%
  mutate(Type = "Protein-Level")


csc1 <- cumulativeSensitivityCurve(
  sim, "Assay1", by = sim$SampleType,
) %>%
  mutate(Type = "Protein-Level")


ji1 <- jaccardIndex(sim, "Assay1") %>%
  mutate(Type = "Protein-Level")


peptide_Missing <- reportMissingValues(sim, "Assay2", by = sim$SampleType)%>%
  mutate(Type = "Peptide-Level")



csc2 <- cumulativeSensitivityCurve(
  sim, "Assay2", by = sim$SampleType,
) %>%
  mutate(Type = "Peptide-Level")

csc <- full_join(csc1, csc2)


csc %>%
  ggplot()+
  aes(x = SampleSize, y = Sensitivity, colour = Type) +
  geom_point() +
  geom_hline(yintercept = peptide_Missing$LocalSensitivityMean,
             linetype = "dashed",
             color = "#F8766D")+
  geom_hline(yintercept = peptide_Missing$TotalSensitivity,
             color = "#F8766D")+
  geom_hline(yintercept = protein_Missing$LocalSensitivityMean,
             linetype = "dashed",
             color = "#00BFC4")+
  geom_hline(yintercept = protein_Missing$TotalSensitivity,
             color = "#00BFC4")+
  scale_y_continuous(limits = c(0,30000))+
  theme_bw(base_size = 20) +
  theme(#legend.position = "none",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = 20),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica"))+
  ylab("Sensitivity (n)")+
  xlab("Number of cells")+
  ggsave("sensitivity_index.png", width = 8, height = 4.7)

ji2 <- jaccardIndex(sim, "Assay2") %>%
  mutate(Type = "Peptide-Level")

ji <- full_join(ji1, ji2)

ji %>%
  ggplot()+
  aes(y = jaccard, x = Type, fill = Type) +
  geom_violin()+
  scale_y_continuous(limits = c(0,1))+
    theme_bw(base_size = 20) +
    theme(legend.position = "none",
          panel.background = element_rect(fill= 'white'),
          axis.text.y=element_text(color = 'black', size = 20),
          axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          text=element_text(family="Helvetica"))+
  ylab("Jaccard Index")+
  xlab("")+
  ggsave("jaccard_index.png", width = 8, height = 4.7)



data_completedness <- full_join(peptide_Missing, protein_Missing) %>%
  mutate(Completeness = Completeness * 100)


 data_completedness %>%
  ggplot()+
  aes(x = Completeness, y = LocalSensitivityMean, color = Type,
      group = Type)+
  geom_point(size =10,
             shape = 5)+
  geom_point()+
  theme_bw()+
  scale_x_continuous(limits = c(0,90))+
  scale_y_continuous(limits = c(0,20000))+
  geom_pointrange(aes(ymin=LocalSensitivityMean-LocalSensitivitySd,
                    ymax=LocalSensitivityMean+LocalSensitivitySd),
                position=position_dodge(0.05))+
   theme_bw(base_size = 20) +
   theme(#legend.position = "none",
         panel.background = element_rect(fill= 'white'),
         axis.text.y=element_text(color = 'black', size = 20),
         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(),
         text=element_text(family="Helvetica"))+
   ylab("Local Sensitivty")+
   xlab("Data completeness (%)")+
   ggsave("data_completeness.png", width = 8, height = 4.7)

