# to check the data format and structure before reading in
# with help from ChatGPT
library(data.table)
library(Matrix)

con <- gzfile("data/GSE206785_scgex.txt.gz", "rt")  # Open compressed file
lines <- readLines(con, n = 10)  # Read first 10 lines
close(con)
cat(lines, sep = "\n")  # Print output
readLines("data/GSE206785_scgex.txt.gz", n = 5)
### this is csv, the data are seperated by comma

df <- fread("data/GSE206785_scgex.txt.gz", nrows = 5)  # Read only first 5 rows
str(df)  # Check structure
### data looks already

# determine the number of rows in the file
total_rows <- system("gunzip -c data/GSE206785_scgex.txt.gz | wc -l", intern = TRUE)
# total_rows <- as.numeric(total_rows) - 1  # Subtract 1 for the header
print(total_rows)

# Read in Chunks & Convert to Sparse Matrix
colnames <- fread("data/GSE206785_scgex.txt.gz", nrows = 1, header = F) %>% unlist()
chunk_size <- 10  # Adjust based on available memory
# num_chunks <- ceiling(total_rows / chunk_size)
num_chunks <- 2
sparse_list <- list()  # Store chunks

for (i in 0:(num_chunks - 1)) {
  skip_rows <- i * chunk_size + 1
  
  # Read chunk
  chunk <- fread("data/GSE206785_scgex.txt.gz", skip = skip_rows, nrows = chunk_size, header = F) %>% as.data.frame
  
  # Handle first column as row names
  setnames(chunk, colnames)  # Ensure column names match
  
  chunk <- column_to_rownames(chunk, var = colnames(chunk)[1])  # Set first column as row names
  
  # Convert to sparse matrix
  sparse_list[[i + 1]] <- as(Matrix(as.matrix(chunk), sparse = TRUE), "dgCMatrix")
  
  rm(chunk)  # Free memory
  gc()  # Run garbage collection
}

# Combine all chunks into one sparse matrix
# final_matrix <- do.call(rbind, sparse_list)
toy_mtx <- do.call(rbind, sparse_list)

# Save final result
saveRDS(final_matrix, "data/GSE206785_sparse_matrix.rds")
