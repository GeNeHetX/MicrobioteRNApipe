library(phyloseq)
library(ggplot2)
library(vegan) 
library(openxlsx)
library(dplyr)
library(stringr)
library(optparse)

process_data <- function(otu_table_path, tax_table_path, sample_data_path) {
  seqtab <- as.matrix(t(read.csv(otu_table_path, sep = ",", header = TRUE, row.names = 1, check.names = FALSE)))
  clin <- read.xlsx(sample_data_path, rowNames = TRUE) 
  taxTab_species <- as.matrix(read.delim(tax_table_path, sep = ",", row.names = 1, header = TRUE)) 
  
  noms_long <- rownames(seqtab)
  noms_court <- rownames(clin) 
  match <- c()
  
  for(name in noms_court) {
    correspondance <- which(stringr::str_detect(noms_long, name))
    index <- ifelse(length(correspondance) == 0, NA, correspondance)
    match <- c(match, index)
  }
  
  names(match) <- noms_court
  
  rownames(seqtab) <- names(sort(match))
  
  ps_s <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE), sample_data(clin), tax_table(taxTab_species))
  
  ps_s_genus <- tax_glom(ps_s, "Genus", NArm = TRUE)
  
  alpha_div <- estimate_richness(ps_s_genus)
  write.csv(alpha_div, "alpha_div_results.csv")
  
  out.pcoa.log <- ordinate(ps_s_genus,  method = "MDS", distance = "bray")
  beta_div <- out.pcoa.log$values
  
  return(list(alpha_div = alpha_div, beta_div = beta_div, ps_s_genus = ps_s_genus))
}

generate_alpha_diversity_boxplot <- function(alpha_div) {
  nouveau_tableau <- data.frame(Sample_name = rownames(alpha_div), Observed = alpha_div[, 1])
  nouveau_tableau$Location <- str_split(nouveau_tableau$Sample_name, "_", simplify = TRUE)[, 1]
  nouveau_tableau$Sample_name <- gsub("\\.", "-", nouveau_tableau$Sample_name)
  
  clin2 <- clin
  clin2$Sample_name <- rownames(clin)
  clin2_subset <- clin2[clin2$Sample_name %in% nouveau_tableau$Sample_name, ]
  indices <- match(nouveau_tableau$Sample_name, clin2_subset$Sample_name)
  nouveau_tableau$Control <- clin2_subset$Control[indices]
  
  boxplot0 <- ggplot(nouveau_tableau, aes(x = Location, y = Observed, fill = Control)) +  
    geom_boxplot() + 
    scale_fill_brewer(palette = "Accent") + 
    labs(title = "Alpha Diversity Variation Based on Sample Origin", 
         x = "Locations", 
         y = "Observed Diversity") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(boxplot0)
}

generate_beta_diversity_plot <- function(ps_s_genus) {
  # Calculate beta diversity distances
  bray_dist <- distance(ps_s_genus, method = "bray")
  
  # Perform PCoA on beta diversity distances
  pcoa_result <- ordinate(ps_s_genus, method = "PCoA", distance = bray_dist)
  
  # Extract PCoA axes
  pcoa_axes <- pcoa_result$vectors
  
  # Create dataframe with PCoA axes and sample metadata
  pcoa_df <- as.data.frame(pcoa_axes)
  pcoa_df$Location <- ps_s_genus@sam_data$Location
  pcoa_df$Control <- ps_s_genus@sam_data$Control
  
  # Plot beta diversity using PCoA axes
  beta_plot <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = Location, shape = Control)) +
    geom_point() +
    labs(title = "Beta Diversity Plot") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(beta_plot)
}

# Argument parsing
option_list <- list(
    make_option(c("-o", "--otu_table"), type="character", default=NULL, 
                help="Path to the OTU table CSV file"),
    make_option(c("-t", "--tax_table"), type="character", default=NULL, 
                help="Path to the taxonomic table CSV file"),
    make_option(c("-s", "--sample_data"), type="character", default=NULL, 
                help="Path to the sample data Excel file")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Process data using provided file paths
if (!is.null(opt$otu_table) && !is.null(opt$tax_table) && !is.null(opt$sample_data)) {
    result <- process_data(opt$otu_table, opt$tax_table, opt$sample_data)
    
    # Generate alpha diversity boxplot
    alpha_plot <- generate_alpha_diversity_boxplot(result$alpha_div)
    ggsave("alpha_div_boxplot.png", plot = alpha_plot, width = 10, height = 8)
    
    # Generate beta diversity plot
    beta_plot <- generate_beta_diversity_plot(result$ps_s_genus)
    ggsave("beta_div_plot.png", plot = beta_plot, width = 10, height = 8)
} else {
    stop("Please provide paths to the OTU table, taxonomic table, and sample data files.")
}
