calculer_grosses_stats_especes_moyenne_mediane <- function(data, scores_a_tester = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)) {
  library(dplyr)
  library(tidyr)
  library(moments) 
  library(purrr)
  
  df_propre <- data %>%
    mutate(confidence_score_num = round(as.numeric(as.character(Score)), 2)) %>%
    filter(!is.na(confidence_score_num))

  seuils_numeriques <- round(as.numeric(scores_a_tester), 2)
  labels_seuils <- paste0("S", gsub("\\.", "", as.character(seuils_numeriques)))

  df_long_stats <- map2_df(seuils_numeriques, labels_seuils, function(s, nom_label) {
    df_propre %>%
      filter(confidence_score_num == s) %>%
      
    group_by(Species, ORIGINE, Source_Microbiome) %>%
    summarise(
      ### READS BRUTS
      Moy_R = mean(Reads, na.rm = TRUE),
      Med_R = median(Reads, na.rm = TRUE),
      Q75_R = quantile(Reads, probs = 0.75, na.rm = TRUE),
      Q95_R = quantile(Reads, probs = 0.95, na.rm = TRUE),
      Q05_R = quantile(Reads, probs = 0.05, na.rm = TRUE),
      Q10_R = quantile(Reads, probs = 0.1, na.rm = TRUE),
      Freq_R = mean(Reads > 0, na.rm = TRUE) * 100, 
      
      ### MINIMIZERS BRUTS
      Moy_M = mean(Uniq_Minimizers, na.rm = TRUE),
      Med_M = median(Uniq_Minimizers, na.rm = TRUE),
      Q75_M = quantile(Uniq_Minimizers, probs = 0.75, na.rm = TRUE),
      Q95_M = quantile(Uniq_Minimizers, probs = 0.95, na.rm = TRUE),
      Q05_M = quantile(Uniq_Minimizers, probs = 0.05, na.rm = TRUE),
      Q10_M = quantile(Uniq_Minimizers, probs = 0.1, na.rm = TRUE),
      Freq_M = mean(Uniq_Minimizers > 0, na.rm = TRUE) * 100, 
      
      ### CPM (Reads Normalisés charge bactérienne trimming)
      Moy_R_CPM = mean(CPM, na.rm = TRUE),
      Med_R_CPM = median(CPM, na.rm = TRUE),
      Q75_R_CPM = quantile(CPM, probs = 0.75, na.rm = TRUE),
      Q95_R_CPM = quantile(CPM, probs = 0.95, na.rm = TRUE),
      Q05_R_CPM = quantile(CPM, probs = 0.05, na.rm = TRUE),
      Q10_R_CPM = quantile(CPM, probs = 0.1, na.rm = TRUE),
      
      ### MPM (Minimizers Normalisés charge bactérienne trimming)
      Moy_M_MPM = mean(Minimizers_Per_Million, na.rm = TRUE),
      Med_M_MPM = median(Minimizers_Per_Million, na.rm = TRUE),
      Q75_M_MPM = quantile(Minimizers_Per_Million, probs = 0.75, na.rm = TRUE),
      Q95_M_MPM = quantile(Minimizers_Per_Million, probs = 0.95, na.rm = TRUE),
      Q05_M_MPM = quantile(Minimizers_Per_Million, probs = 0.05, na.rm = TRUE),
      Q10_M_MPM = quantile(Minimizers_Per_Million, probs = 0.1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~ if_else(is.na(.) | is.nan(.), 0, .))) %>%
    mutate(Seuil_Label = nom_label)
  })

  rapport_final_especes <- df_long_stats %>%
    pivot_wider(
      names_from = c(ORIGINE, Seuil_Label),
      values_from = c(
        Moy_R, Med_R, Q75_R, Q95_R, Freq_R,
        Moy_M, Med_M, Q75_M, Q95_M, Freq_M,
        Moy_R_CPM, Med_R_CPM, Q75_R_CPM, Q95_R_CPM,
        Moy_M_MPM, Med_M_MPM, Q75_M_MPM, Q95_M_MPM,
        Q05_R, Q10_R,Q05_M, Q10_M,Q05_R_CPM, Q10_R_CPM,Q05_M_MPM, Q10_M_MPM
      ),
      names_glue = "{ORIGINE}_{Seuil_Label}_{.value}",
      values_fill = 0
    ) %>%
    select(Species, Source_Microbiome, everything())
    
  return(rapport_final_especes)
}


especes_contaminants_env_200 <- c(
  "Burkholderia contaminans", "Caulobacter mirabilis", "Brevundimonas albigilva", 
  "Burkholderia lata", "Aromatoleum toluolicum", "Klebsiella variicola", 
  "Burkholderia pyrrocinia", "Agrobacterium pusense", "Ralstonia insidiosa", 
  "Xanthobacter autotrophicus", "Bacillus velezensis", "Streptomyces armeniacus", 
  "Agrobacterium tumefaciens", "Burkholderia stabilis", "Ilumatobacter coccineus", 
  "Tessaracoccus flavescens", "Paraburkholderia pallida", "Burkholderia arboris", 
  "Enterobacter bugandensis", "Burkholderia anthina", "Burkholderia ubonensis", 
  "Telmatobacter bradus", "Microlunatus phosphovorus", "Burkholderia glumae", 
  "Maricaulis maris", "Sphingomonas koreensis", "Urbifossiella limnaea", 
  "Paracoccus contaminans", "Vibrio harveyi", "Roseateles sp.", 
  "Burkholderia thailandensis", "Ralstonia solanacearum", "Micropruina glycogenica", 
  "Burkholderia territorii", "Dermacoccus nishinomiyaensis", "Staphylococcus capitis", 
  "Methyloversatilis sp.", "Hydrogenophilus thermoluteolus", "Kocuria rhizophila", 
  "Telmatocola sphagniphila", "Rothia aeria", "Rhodoluna lacicola", 
  "Actinomarinicola tropica", "Roseomonas gilardii", "Sulfuritortus calidifontis", 
  "Sulfurimicrobium lacus", "Tuwongella immobilis", "Asticcacaulis excentricus", 
  "Sphingomonas hankookensis", "Achromobacter marplatensis", "Sphingobium yanoikuyae", 
  "Polynucleobacter difficilis", "Aquabacterium olei", "Microvirga lotononidis", 
  "Aquihabitans daechungensis", "Verrucomicrobium spinosum", "Brevundimonas sp.", 
  "Gemmata obscuriglobus", "Burkholderia stagnalis", "Humibacillus xanthopallidus", 
  "Sandaracinus amylolyticus", "Aquirufa antheringensis", "Capillimicrobium parvum", 
  "Rhodoferax lithotrophicus", "Paracoccus kondratievae", "Xanthomonas campestris",

  "Thiobacillus denitrificans", "Limosilactobacillus pontis", "Thermoflavifilum sp.", 
  "Rhodococcus pyridinivorans", "Carnobacterium maltaromaticum", "Pseudobacter ginsenosidimutans", 
  "Affinirhizobium pseudoryzae", "Novosphingobium aromaticivorans", "Heliomicrobium modesticaldum", 
  "Streptomyces thermocarboxydus", "Cronobacter sakazakii", "Pseudonocardia saturnea", 
  "Alcanivorax sp.", "Bacteroides graminisolvens", "Methylobacterium aquaticum", 
  "Sinomonas atrocyanea", "Afipia carboxydohydrogena", "Brucella pituitosa", 
  "Pelolinea submarina", "Paracoccus versutus", "Pandoraea thiooxydans", 
  "Sphingomonas daechungensis", "Sphingomonas yabuuchiae", "Paracoccus aestuarii", 
  "Paroceanicella profunda", "Celeribacter indicus", "Flavobacterium commune", 
  "Microbacterium proteolyticum", "Tardiphaga alba", "Comamonas serinivorans", 
  "Modestobacter italicus", "Stappia sp.", "Arthrobacter agilis", 
  "Azonexus hydrophilus", "Shewanella xiamenensis", "Parvibaculum lavamentivorans", 
  "Kitasatospora purpeofusca", "Pigmentiphaga litoralis", "Silicimonas algicola", 
  "Hydrocarboniphaga effusa", "Nocardioides marinisabuli", "Janthinobacterium lividum", 
  "Marinospirillum sp.", "Methylobacterium radiotolerans", "Microbacterium azadirachtae", 
  "Flavisolibacter tropicus", "Leucobacter luti", "Nocardiopsis exhalans", 
  "Microbacterium lemovicicum", "Phycisphaera mikurensis", "Advenella alkanexedens", 
  "Devosia sp.", "Pseudonocardia dioxanivorans", "Saccharopolyspora erythraea", 
  "Acetobacter aceti", "Microbacterium terricola", "Pseudomonas azotoformans", 
  "Rhizobium bangladeshense", "Brevundimonas naejangsanensis", "Marisediminicola antarctica", 
  "Microbacterium testaceum", "Acidothermus cellulolyticus", "Kutzneria albida", 
  "Pikeienuella piscinae", "Rhodovastum atsumiense", "Allostella humosa", 
  "Listeria innocua",

  "Flavobacterium cerinum", "Janibacter cremeus", "Streptomyces malaysiensis", 
  "Flavobacterium magnum", "Oxalobacter vibrioformis", "Tessaracoccus flavus", 
  "Thiospirochaeta perfilievii", "Desulfosudis oleivorans", "Thermomonas brevis", 
  "Deinococcus actinosclerus", "Methylomusa anaerophila", "Pseudomonas lurida", 
  "Kocuria turfanensis", "Deinococcus metallilatus", "Agromyces flavus", 
  "Alteromonas macleodii", "Microbacterium invictum", "Pseudonocardia abyssalis", 
  "Thiobacillus sp.", "Agrobacterium rubi", "Labrys okinawensis", 
  "Thauera butanivorans", "Desulfobacca acetoxidans", "Paraburkholderia ginsengisoli", 
  "Blastochloris viridis", "Micromonospora zamorensis", "Sphingopyxis terrae", 
  "Bermanella marisrubri", "Kribbella flavida", "Larkinella insperata", 
  "Gordonia amarae", "Kocuria flava", "Pelosinus fermentans", 
  "Streptomyces virginiae", "Kocuria marina", "Orrella dioscoreae", 
  "Anaeromicropila herbilytica", "Dyella jiangningensis", "Salinibacter ruber", 
  "Alkalitalea saponilacus", "Corynebacterium bovis", "Pantoea anthophila", 
  "Xanthomonas axonopodis", "Advenella kashmirensis", "Bradyrhizobium elkanii", 
  "Isoptericola jiangsuensis", "Clostridium fermenticellae", "Flavobacterium sediminis", 
  "Sphingomicrobium flavum", "Xanthomonas arboricola", "Clostridium thermosuccinogenes", 
  "Corynebacterium humireducens", "Methylobacterium oryzae", "Heliorestis convoluta", 
  "Actinoplanes derwentensis", "Methylomicrobium lacus", "Oceanithermus profundus", 
  "Sphingopyxis lindanitolerans", "Actinocatenispora thailandica", "Alloyangia pacifica", 
  "Lysobacter soli", "Martelella mediterranea", "Streptomyces fradiae", 
  "Bradyrhizobium canariense", "Flavobacterium arcticum", "Sphingobium chlorophenolicum", 
  "Cupriavidus campinensis"
)

especes_strictes_colon <- c(
  "Faecalibacterium prausnitzii",   
  "Roseburia hominis",              
  "Roseburia faecis",               
  "Roseburia intestinalis",         
  "Anaerostipes hadrus",            
  "Anaerostipes caccae",            
  "Agathobacter rectalis",          
  
  "Phocaeicola vulgatus",           
  "Bacteroides thetaiotaomicron",   
  "Bacteroides uniformis",         
  "Bacteroides cellulosilyticus",   
  "Bacteroides ovatus",             
  "Bacteroides stercoris",          
  
  "Akkermansia muciniphila",        
  "Peptacetobacter hiranonis",      
  "Odoribacter splanchnicus",       
  "Phascolarctobacterium faecium",  
  "Phascolarctobacterium succinatutens", 
  
  "Subdoligranulum variabile",      
  "Dysosmobacter welbionis"
)

### couleurs pour les plots

eco_annot <- c(
  "Acidaminococcus"          = "Gut anaerobe",
  "Anaerococcus"             = "Gut anaerobe",
  "Anaerostipes"             = "Gut anaerobe",
  "Bacteroides"              = "Gut anaerobe",
  "Bifidobacterium"          = "Gut anaerobe",
  "Blautia"                  = "Gut anaerobe",
  "Clostridium"              = "Gut anaerobe",
  "Collinsella"              = "Gut anaerobe",
  "Coprococcus"              = "Gut anaerobe",
  "Dialister"                = "Gut anaerobe",
  "Dorea"                    = "Gut anaerobe",
  "Faecalibacterium"         = "Gut anaerobe",
  "Faecalimonas"             = "Gut anaerobe",
  "Finegoldia"               = "Gut anaerobe",
  "Lachnospira"              = "Gut anaerobe",
  "Lancefieldella"           = "Gut anaerobe",
  "Mediterraneibacter"       = "Gut anaerobe",
  "Parabacteroides"          = "Gut anaerobe",
  "Phascolarctobacterium"    = "Gut anaerobe",
  "Phocaeicola"              = "Gut anaerobe",
  "Prevotella"               = "Gut anaerobe",
  "Roseburia"                = "Gut anaerobe",
  "Ruminococcus"             = "Gut anaerobe",
  "Segatella"                = "Gut anaerobe",
  "Simiaoa"                  = "Gut anaerobe",
  "Succinivibrio"            = "Gut anaerobe",
  "Wujia"                    = "Gut anaerobe",

  "Aggregatibacter"          = "Oral",
  "Campylobacter"            = "Oral",
  "Capnocytophaga"           = "Oral",
  "Cetobacterium"            = "Oral",
  "Faucicola"                = "Oral",
  "Fusobacterium"            = "Oral",
  "Gemella"                  = "Oral",
  "Granulicatella"           = "Oral",
  "Haemophilus"              = "Oral",
  "Lautropia"                = "Oral",
  "Neisseria"                = "Oral",
  "Parvimonas"               = "Oral",
  "Porphyromonas"            = "Oral",
  "Rothia"                   = "Oral",
  "Schaalia"                 = "Oral",
  "Solobacterium"            = "Oral",
  "Streptococcus"            = "Oral",
  "Treponema"                = "Oral",
  "Veillonella"              = "Oral",

  "Buchnera"                 = "Enterobacteria",
  "Citrobacter"              = "Enterobacteria",
  "Enterobacter"             = "Enterobacteria",
  "Escherichia"              = "Enterobacteria",
  "Klebsiella"               = "Enterobacteria",
  "Kluyvera"                 = "Enterobacteria",
  "Morganella"               = "Enterobacteria",
  "Proteus"                  = "Enterobacteria",
  "Providencia"              = "Enterobacteria",
  "Salmonella"               = "Enterobacteria",
  "Serratia"                 = "Enterobacteria",
  "Shigella"                 = "Enterobacteria",
  "Yersinia"                 = "Enterobacteria",

  "Enterococcus"             = "Lactic",
  "Lactobacillus"            = "Lactic",
  "Lacticaseibacillus"       = "Lactic",
  "Lactiplantibacillus"      = "Lactic",
  "Lactococcus"              = "Lactic",
  "Leuconostoc"              = "Lactic",
  "Ligilactobacillus"        = "Lactic",
  "Limosilactobacillus"      = "Lactic",
  "Pediococcus"              = "Lactic",
  "Weissella"                = "Lactic",

  "Achromobacter"            = "Mucosal / Commensal",
  "Actinomyces"              = "Mucosal / Commensal",
  "Brevibacterium"           = "Mucosal / Commensal",
  "Corynebacterium"          = "Mucosal / Commensal",
  "Cutibacterium"            = "Mucosal / Commensal",
  "Dolosigranulum"           = "Mucosal / Commensal",
  "Gardnerella"              = "Mucosal / Commensal",
  "Kocuria"                  = "Mucosal / Commensal",
  "Listeria"                 = "Mucosal / Commensal",
  "Micrococcus"              = "Mucosal / Commensal",
  "Moraxella"                = "Mucosal / Commensal",
  "Nosocomiicoccus"          = "Mucosal / Commensal",
  "Staphylococcus"           = "Mucosal / Commensal",

  "Acinetobacter"="Contaminant", "Acidovorax"="Contaminant", "Actinoalloteichus"="Contaminant",
  "Asticcacaulis"="Contaminant", "Bacillus"="Contaminant", "Bordetella"="Contaminant",
  "Bradyrhizobium"="Contaminant", "Brevibacillus"="Contaminant", "Comamonas"="Contaminant",
  "Georhizobium"="Contaminant", "Kosakonia"="Contaminant", "Microbacterium"="Contaminant",
  "Nitratireductor"="Contaminant", "Paucibacter"="Contaminant", "Peribacillus"="Contaminant",
  "Pseudomonas"="Contaminant", "Rheinheimera"="Contaminant", "Shinella"="Contaminant",
  "Stutzerimonas"="Contaminant", "Vibrio"="Contaminant", "Williamsia"="Contaminant"
)