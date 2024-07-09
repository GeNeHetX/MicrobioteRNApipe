
process_data <- function(otu_table_path, tax_table_path, sample_data_path) {
  ## 
  seqtab <- as.matrix(t(read.csv(otu_table_path, sep = ",", header = TRUE, row.names = 1, check.names = FALSE)))
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
  
  ps_s_genus <- tax_glom(ps_s, "Genus", NArm = FALSE)
  
  # Calculate Shannon diversity index
  alpha_div <- estimate_richness(ps_s_genus, measures = "Observed")
  #, measures = "Shannon"
  write.csv(alpha_div, "alpha_div_results_shannon.csv")
  
  out.pcoa.log <- ordinate(ps_s_genus,  method = "MDS", distance = "bray")
  beta_div <- out.pcoa.log$values
  
  return(list(alpha_div = alpha_div, beta_div = beta_div, ps_s_genus = ps_s_genus))
}
generate_alpha_diversity_boxplot <- function(alpha_div) {
  nouveau_tableau <- data.frame(Sample_name = rownames(alpha_div), Shannon_index = alpha_div[, 1])
  nouveau_tableau$Location <- str_split(nouveau_tableau$Sample_name, "_", simplify = TRUE)[, 1]
  nouveau_tableau$Sample_name <- gsub("\\.", "-", nouveau_tableau$Sample_name)
  
  # Assign categories based on sample names
  nouveau_tableau$Category <- NA
  nouveau_tableau$Category[str_detect(nouveau_tableau$Sample_name, "^BPDAC")] <- "pancreas_microdissec"
  nouveau_tableau$Category[str_detect(nouveau_tableau$Sample_name, "^X")] <- "pancreas_microdissec"
  
  
  
  
  clin2 <- clin
  clin2$Sample_name <- rownames(clin)
  clin2_subset <- clin2[clin2$Sample_name %in% nouveau_tableau$Sample_name, ]
  indices <- match(nouveau_tableau$Sample_name, clin2_subset$Sample_name)
  nouveau_tableau$Control <- clin2_subset$Control[indices]
  nouveau_tableau <- nouveau_tableau[!(nouveau_tableau$Category %in% c("Pancreas_FFPE", "Pancreas_frozen", "liver", "bladder", "Biliary tract")), ]
  
  boxplot0 <- ggplot(nouveau_tableau, aes(x = Category, y = Shannon_index, fill = Control)) +  
    geom_boxplot(color = "black") +  # Make boxplot outlines black
    geom_jitter(position = position_jitter(width = 0.3), color = "black", alpha = 0.6) +  # Add jittered points with black color and lower alpha
    scale_fill_brewer(palette = "Accent") +
    labs(title = "Observed Diversity Variation Based on Sample Origin (Confidence score 0.1 + 16s + filter reads : 10)", 
         x = "Locations", 
         y = "Observed Diversity ") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis text size
      axis.text.y = element_text(size = 12),  # Adjust y-axis text size
      axis.title.x = element_text(size = 12),  # Adjust x-axis title size
      axis.title.y = element_text(size = 12),  # Adjust y-axis title size 
      legend.text = element_text(size = 12),
      # Adjust legend text size
    )
  return(boxplot0)
}

generate_beta_diversity_plot <- function(ps_s_genus) {
  # Calculate beta diversity distances
  bray_dist <- distance(ps_s_genus, method = "jaccard", binary = TRUE)
  #canberra
  
  # Perform PCoA on beta diversity distances
  pcoa_result <- ordinate(ps_s_genus, method = "PCoA", distance = bray_dist)
  
  # Extract PCoA axes
  pcoa_axes <- pcoa_result$vectors
  
  # Create dataframe with PCoA axes and sample metadata
  pcoa_df <- as.data.frame(pcoa_axes)
  pcoa_df$histological_aspect <- ps_s_genus@sam_data$histological_aspect
  pcoa_df$Location <- ps_s_genus@sam_data$Location
  pcoa_df$Control <- ps_s_genus@sam_data$Control
  
  # Filter out samples labeled as "pancreas_FFPE" and "pancreas_frozen"
  pcoa_df_filtered <- pcoa_df[!(pcoa_df$Location %in% c("pancreas_FFPE", "pancreas_frozen", "retine", "bladder", "Biliary tract")), ]
  
  # Exclude rows with NA values
  pcoa_df_filtered <- na.omit(pcoa_df_filtered)
  
  # Plot beta diversity using PCoA axes
  beta_plot <- ggplot(pcoa_df_filtered, aes(x = Axis.1, y = Axis.2, color = histological_aspect)) +
    geom_point() +
    labs(title = "Beta Diversity Plot based on histological aspect") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(beta_plot)
}



