library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")

# Récupérer les arguments de la ligne de commande
args <- commandArgs(trailingOnly = TRUE)

# Initialiser la variable pour le répertoire de sortie
biom_data <- NULL

# Vérifier si l'option -tab est présente dans les arguments
if ("-tab" %in% args) {
  # Trouver l'index de l'option -tab
  indexTab <- which(args == "-tab")

  # Extraire la valeur associée à l'option -tab
  if (length(args) > (indexTab + 1)) {
    biom_data <- args[indexTab + 1]

  }
}

merged_metagenomes <- import_biom("all_sample_contig_table.biom")

merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

tab1 <- merged_metagenomes@tax_table@.Data
# Enregistrer en CSV
write.csv(tab1, file = "tax_table.csv")

tab2 <- merged_metagenomes@otu_table@.Data
write.csv(tab2, file = "otu_table.csv")