process_data <- function(tax_table_path, otu_table_path) {
    data1 <- read.csv(tax_table_path)
    data2 <- read.csv(otu_table_path)
    
    # Identify rows to remove from the first file based on condition
    rows_to_remove <- which(data1$Kingdom %in% c("Eukaryota", "Viruses", "Archaea"))
    
    # Extract IDs to remove from the second file
    ids_to_remove <- data1$X[rows_to_remove]
    
    # Remove rows from the first file
    data1 <- data1[!data1$Kingdom %in% c("Eukaryota", "Viruses", "Archaea"), ]
    
    # Remove corresponding rows from the second file
    data2 <- data2[!data2$X %in% ids_to_remove, ]
    
    colnames(data2) <- gsub("\\.", "-", colnames(data2))
    
    # Write the filtered data back to original files, overwriting them
    write.csv(data1, tax_table_path, row.names = FALSE)
    write.csv(data2, otu_table_path, row.names = FALSE)
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    cat("Usage: Rscript script.R <tax_table_path> <otu_table_path>\n")
    quit(save = "no", status = 1)
}

tax_table_path <- args[1]
otu_table_path <- args[2]

# Call the function with provided arguments
process_data(tax_table_path, otu_table_path)
