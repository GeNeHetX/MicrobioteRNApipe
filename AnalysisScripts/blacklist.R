library(ggplot2)
library(tidyverse)

library(dplyr)

required_packages <- c("tidyverse", "moments", "pROC")

missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if (length(missing_packages) > 0) {
  message("Installation des packages manquants : ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

invisible(lapply(required_packages, library, character.only = TRUE))

library(tidyr)
library(purrr)
library(moments)

library(pROC)


tab_species_TRIM <- read.csv("path_to_extraction_table.csv", header = TRUE, strip.white = TRUE)

### chemin pour base de donne BacDive
setwd("path_to/data_db_bact")
fichiers <- list.files(pattern = "*.csv")

all_data <- fichiers %>%
  purrr::map_df(~read.csv(.x))

df_classif_clean <- all_data %>%
  fill(Species, .direction = "down") %>%
  group_by(Species) %>%
  summarise(
    Categories_Globales = paste(unique(na.omit(c(`Category.1`, `Category.2`, `Category.3`))), collapse = "; "),
    Source_Isolation    = dplyr::first(na.omit(`Isolation.source`)), 
    .groups = "drop"
  )


setwd("path_to_work_direct")

df_final <- tab_species_TRIM %>%
    left_join(df_classif_clean, by = "Species")

df_final_TRIM <- df_final %>%
  mutate(Source_Microbiome = case_when(
    grepl("Homo sapiens", Species) ~ "Hôte",
    #grepl("Gastrointestinal|Oral|Gut|Fecal|feces|swab|intestine|stomach|mouth|colon|duodenum|jejunum|saliva|rectum|oropharynx|tongue|caecum|rumen|abdomen", Categories_Globales, ignore.case = TRUE) ~ "Digestif/Oral",
    #grepl("Gastrointestinal|Oral|Gut|Fecal|feces|swab|intestine|stomach|mouth|colon|duodenum|jejunum|saliva|rectum|oropharynx|tongue|caecum|rumen|abdomen", Source_Isolation, ignore.case = TRUE) ~ "Digestif/Oral",

    grepl("duodenum|jejunum|stomach|gastric|small intestine|ileum", Categories_Globales, ignore.case = TRUE) ~ "Microbiote Digestif/Oral",
    grepl("duodenum|jejunum|stomach|gastric|small intestine|ileum", Source_Isolation, ignore.case = TRUE) ~ "Microbiote Digestif/Oral",

    grepl("colon|colorectal|rectum|feces|fecal|stool|caecum|appendix|gut|Gastrointestinal|intestine|abdomen|rumen", Categories_Globales, ignore.case = TRUE) ~ "Microbiote Digestif/Oral",
    grepl("colon|colorectal|rectum|feces|fecal|stool|caecum|appendix|gut|Gastrointestinal|intestine|abdomen|rumen", Source_Isolation, ignore.case = TRUE) ~ "Microbiote Digestif/Oral",

    grepl("Oral|mouth|saliva|oropharynx|tongue|dental|teeth|gingival|periodontal", Categories_Globales, ignore.case = TRUE) ~ "Microbiote Digestif/Oral",
    grepl("Oral|mouth|saliva|oropharynx|tongue|dental|teeth|gingival|periodontal", Source_Isolation, ignore.case = TRUE) ~ "Microbiote Digestif/Oral",

    grepl("skin|hand|elbow|foot|neck|epidermis|integument|nare|nose|nasal|acne", Categories_Globales, ignore.case = TRUE) ~ "Autre Microbiote Humain/Cutané",
    grepl("skin|hand|elbow|foot|neck|epidermis|integument|nare|nose|nasal|acne", Source_Isolation, ignore.case = TRUE) ~ "Autre Microbiote Humain/Cutané",

    grepl("Human", Categories_Globales, ignore.case = TRUE) ~ "Autre Microbiote Humain/Cutané",
    grepl("Human", Source_Isolation, ignore.case = TRUE) ~ "Autre Microbiote Humain/Cutané",

    grepl("Environmental|Condition|Soil|Water|Rhizosphere|forest|marine|plant|lake|river|mud|sediment|sludge|pond|stream|aquifer|cave|air|dust|compost|field|beach|wood|root|leaf|flower|spore|tree|grass", Categories_Globales, ignore.case = TRUE) ~ "Contaminant Env",
    grepl("Environmental|Condition|Soil|Water|Rhizosphere|forest|marine|plant|lake|river|mud|sediment|sludge|pond|stream|aquifer|cave|air|dust|compost|field|beach|wood|root|leaf|flower|spore|tree|grass", Source_Isolation, ignore.case = TRUE) ~ "Contaminant Env",

    grepl("hospital|clinic|medical|patient|intensive care|surgery|prosthesis|abscess|wound|sputum|lesion|ulcer", Categories_Globales) ~ "Autre",
    grepl("hospital|clinic|medical|patient|intensive care|surgery|prosthesis|abscess|wound|sputum|lesion|ulcer", Source_Isolation) ~ "Autre",
    grepl("Infection|Pathogen", Categories_Globales, ignore.case = TRUE) ~ "Autre",
    grepl("irus|hage", Species, ignore.case = TRUE) ~ "Autre",
    grepl("andidatus", Species) ~ "Autre",
    is.na(Categories_Globales) ~ "Autre",
    TRUE ~ "Autre"
  )) %>%
  select(-Categories_Globales, -Source_Isolation)

df_final_TRIM_pauline <- df_filtre %>% filter(ORIGINE == "stomach_cancer") %>% select(-Categories_Globales, -Source_Isolation)
write.csv(df_final_TRIM_pauline, "Resultat_Pipe_estomac.csv", row.names = FALSE)


df_filtre <- df_final_TRIM %>%
  filter(Categorie =="Bacterie")


df_filtre <- df_filtre %>%
  #filter(Reads>5)%>%
  group_by(Sample, Species) %>%
  mutate(Reads_Somme = sum(Reads, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(Total_Reads_Sample = sum(Reads, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    CPM = if_else(Total_Reads_Sample > 0, (Reads_Somme / Total_Reads_Sample) * 1e6, 0),
    Score = as.numeric(as.character(Score))
  )




df_filtre_CPM_TRIM <- df_filtre %>%
  filter(ORIGINE %in% c("Colon", "Duodenum", "Jejunum", "CHC", "PERMAJI", "BESANCON", "IPMN", "RETINE", "VESSIE", "stomach_cancer", "slice_pancreas_tumor_TT", "pancreas_cancerPanNET", "ovary"),Score==0.2)%>%
  mutate(Uniq_Minimizers_Clean = if_else(is.na(Uniq_Minimizers), 0, as.numeric(Uniq_Minimizers))) %>%
  mutate(
    Minimizers_Per_Million = if_else(Total_Reads_Sample > 0, (Uniq_Minimizers_Clean * 1e6) / Total_Reads_Sample, 0)
  )%>%
  group_by(Sample, Species) %>%
  mutate(Reads_Somme = sum(Reads, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(CPM = if_else(Total_Reads_Sample > 0, (Reads_Somme / Total_Reads_Sample) * 1e6, 0))





source("path_to/fonctions_matrice.R")
ma_super_matrice_especes_complets <- calculer_grosses_stats_especes_moyenne_mediane(df_filtre_CPM_TRIM , scores_a_tester =  c(0.2))


ma_super_matrice_especes_complets_test <- ma_super_matrice_especes_complets %>%
  mutate(Groupe_Test = case_when(
    Species %in% c(especes_strictes_colon) ~ "Colon",
    Species %in% especes_contaminants_env_200 ~ "Contaminant",
    TRUE ~ "Autre" 
  ))


cols_controles_sans_IG <- c(
  "CHC_S02_Q05_R", 
  "stomach_cancer_S02_Q05_R", 
  "slice_pancreas_tumor_TT_S02_Q05_R", 
  "pancreas_cancerPanNET_S02_Q05_R",
  "ovary_S02_Q05_R")

matrice_enrichie <- ma_super_matrice_especes_complets_test %>%
  rowwise() %>%                                                                                                                                                       
  mutate(
    Presence_Controles_sans_IG = mean(c_across(all_of(cols_controles_sans_IG)) > 0, na.rm = TRUE)
  ) %>%
  ungroup()

roc_presence <- roc(
  response = matrice_enrichie$Groupe_Test, 
  predictor = matrice_enrichie$Presence_Controles_sans_IG, 
  levels = c("Colon", "Contaminant"),
  direction = ">",
  quiet = TRUE
)

auc(roc_presence)
seuil_presence_dure <- coords(roc_presence, "best", best.method = "youden")$threshold

roc_minimizers <- roc(
  response = matrice_enrichie$Groupe_Test, 
  predictor = matrice_enrichie$Colon_S02_Moy_M, 
  levels = c("Contaminant", "Colon"),
  direction = "<",
  quiet = TRUE
)
seuil_minimizers_opti <- coords(roc_minimizers, "best", best.method = "youden")$threshold


seuil_presence_dure <- seuil_presence_dure[is.finite(seuil_presence_dure)][1]
seuil_minimizers_opti <- seuil_minimizers_opti[is.finite(seuil_minimizers_opti)][1]

blackliste_especes <- matrice_enrichie %>%
  filter(
    Presence_Controles_sans_IG > seuil_presence_dure,
    Colon_S02_Moy_M < seuil_minimizers_opti
    #Species %in% c(especes_contaminants_env_200, especes_strictes_colon)
  ) %>%
  select(Species, Groupe_Test, Source_Microbiome, Presence_Controles_sans_IG, Colon_S02_Moy_M)%>% arrange(desc(Colon_S02_Moy_M))

print(seuil_presence_dure)
print(seuil_minimizers_opti)
print(nrow(ma_super_matrice_especes_complets)-nrow(blackliste_especes))
print(table(blackliste_especes$Source_Microbiome))
print(table(ma_super_matrice_especes_complets$Source_Microbiome))



write.csv(blackliste_especes, "Blacklist.csv", row.names = FALSE)






