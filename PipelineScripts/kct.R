library(ggplot2)

# Check if arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) < 1) {
  stop("Please provide the folder path as an argument.")
}

# Get folder_path from command-line argument
folder_path <- commandArgs(trailingOnly = TRUE)[1]

# Get the list of files in the folder
files <- list.files(folder_path)

# Define the colors and symbols
point_color <- "darkblue"
point_symbol <- 19

# Create an empty data frame to store correlation test results
correlation_results <- data.frame(taxon = character(), spearman_correlation = numeric(), stringsAsFactors = FALSE)

# Create a vector to store unique taxons found
unique_taxons <- character()

# Loop through each file
for (file in files) {
  # Read the data from the current file
  data <- read.delim(file.path(folder_path, file), header = FALSE, stringsAsFactors = FALSE)
  
  # Extract rows with taxon "G"
  genus_data <- data[data[, 6] == "G",]
  
  # Extract unique taxons from the current file
  current_unique_taxons <- unique(genus_data[, 8])
  
  # Filter out taxons that were found before
  new_unique_taxons <- setdiff(current_unique_taxons, unique_taxons)
  
  # Update unique_taxons vector
  unique_taxons <- c(unique_taxons, new_unique_taxons)
  
  # Loop through each unique taxon in the current file
  for (taxon in new_unique_taxons) {
    # Remove "_" from taxon name
    taxon_cleaned <- gsub("_", "", taxon)
    
    # Initialize vectors to store data for the current taxon
    X <- numeric()
    Y <- numeric()
    
    # Loop through each subsequent file
    for (i in seq_along(files)) {
      # Read the data from the subsequent file
      subsequent_data <- read.delim(file.path(folder_path, files[i]), header = FALSE, stringsAsFactors = FALSE)
      
      # Extract rows with the same taxon as the current iteration
      taxon_data <- subsequent_data[subsequent_data[, 8] == taxon, ]
      
      # Check if taxon_data is empty
      if (nrow(taxon_data) > 0) {
        # Append Unique number of k-mers and Total number of k-mers to X and Y vectors
        X <- c(X, as.numeric(taxon_data[, 5]))  # Number of unique k-mers
        Y <- c(Y, as.numeric(taxon_data[, 4]))  # Total number of k-mers
      }
    }
    
    # Perform Spearman correlation test
    correlation <- try(cor.test(X, Y, method = "spearman"), silent = TRUE)
    print(correlation)
    
    # Check if Spearman correlation test was successful
    if (!inherits(correlation, "try-error")) {
      # Save the spearman correlation test result to the correlation_results data frame
      correlation_results <- rbind(correlation_results, data.frame(taxon = taxon_cleaned, spearman_correlation = correlation$estimate))
    }
  }
}

write.table(correlation_results, file = "correlation.csv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Correlation results saved as CSV\n")

scatterplot <- ggplot(correlation_results, aes(x = taxon, y = spearman_correlation)) +
  geom_point(color = "darkgreen", size = 1.5) +
  theme_minimal() +  # Use minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        panel.background = element_rect(fill = "white"),  # Set background color to white inside plot
        plot.background = element_rect(fill = "white")) +  # Set background color to white outside plot
  labs(x = "Taxon", y = "Spearman Correlation") +
  ggtitle("K-mer correlation tests") +
  scale_x_discrete(labels = NULL) +  # Hide taxon labels
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())  # Remove minor gridlines



ggsave("scatterplot.png", plot = scatterplot)
