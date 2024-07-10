## Authors

- [Ali YOUNCHA](https://github.com/MrAli1582)

  
# Introduction

This tutorial demonstrates how to directly modify the functions.R script to match
the output of the Alpha and Beta diversity analysis with your desired analysis.


# 1) Performing Alpha and Beta diversity on the species level :

If you desire performing your analysis on the `species` level, replace the following
line of the process_data function:

  `ps_s_genus <- tax_glom(ps_s, "Genus", NArm = FALSE)`

with:

  `ps_s_genus <- tax_glom(ps_s, "Species", NArm = FALSE)`

# 2) Changing the measure method used to perform and plot the alpha diversity :

For the alpha diversity function `process_data`, replace:

  `alpha_div <- estimate_richness(ps_s_genus, measures = "Shannon")`

with:

  `alpha_div <- estimate_richness(ps_s_genus, measures = "Observed") ## Or any other type of mesure you desire to use/test.`


# 3) Assigning categories for the sample names when plotting the boxplots of the alpha diversity:

If you want to categorize your samples when plotting the alpha diversity results using the `generate_alpha_diversity_boxplot`, replace :

  `nouveau_tableau$Category[str_detect(nouveau_tableau$Sample_name, "^BPDAC")] <- "pancreas_microdissec"`
  `nouveau_tableau$Category[str_detect(nouveau_tableau$Sample_name, "^X")] <- "pancreas_microdissec"`

with: (For example if you are working on Melanoma samples that start with TRZ in their name)

  `nouveau_tableau$Category[str_detect(nouveau_tableau$Sample_name, "^TRZ")] <- "Melanoma"`

# 4) Changing the distance calculation method when performing Beta diversity analysis:

The function `generate_beta_diversity_boxplot` uses by default `jaccard` method to calculate the beta diversity distances, if you want to change it, replace:

  `bray_dist <- distance(ps_s_genus, method = "jaccard", binary = TRUE)`

with:

  `bray_dist <- distance(ps_s_genus, method = "bray", binary = TRUE) ## or any other method that you desire to use/test.`

Also, by default PcoA (Principal Coordinate Analysis) is performed on the beta diversity distances. By replacing this line of the code:


  `pcoa_result <- ordinate(ps_s_genus, method = "PCoA", distance = bray_dist)`

with:

  `pcoa_result <- ordinate(ps_s_genus, method = "MDS", distance = bray_dist) ##or any other method that you desire to use/test. Here are the available methods for this function ("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")`


Moreover, changing the label that colors the dots on the Beta plot requires replacing the following line in the `generate_beta_diversity_boxplot` function:

  `pcoa_df$histological_aspect <- ps_s_genus@sam_data$histological_aspect`

with:

  `pcoa_df$sous_type_visuel <- ps_s_genus@sam_data$sous_type_visuel ## Or any other label that you want to use/test.`


Finally, if you want to test different combinations of beta axis. Change the following line of the code:

  `beta_plot <- ggplot(pcoa_df_filtered, aes(x = Axis.1, y = Axis.2, color = histological_aspect)) +
    geom_point() +
    labs(title = "Beta Diversity Plot based on histological aspect") +
    theme(plot.title = element_text(hjust = 0.5))`

with: (Based on the axis that you want to plot)

  `beta_plot <- ggplot(pcoa_df_filtered, aes(x = Axis.1, y = Axis.3, color = histological_aspect)) +
    geom_point() +
    labs(title = "Beta Diversity Plot based on histological aspect") +
    theme(plot.title = element_text(hjust = 0.5))`
