# read in matrix
# 2026-03-21

library(tidyverse)
library(data.table)
library(Matrix)

# read in data
dt <- fread("data/GSE206785_scgex.txt.gz")
cell_ids <- dt$Cell
mtx <- as.matrix(dt[, -1, with = FALSE])
mtx <- Matrix(mtx, sparse = TRUE)
rownames(mtx) <- cell_ids
colnames(mtx) <- colnames(dt)[-1]

# data distribution
hist(mtx@x, breaks = 50)
### values are near 1 with some high values > might be log-transformed

# check the total sum of the counts
expm1(mtx) %>% rowSums() %>% tail()
### total sum of the counts are not same for all the cells, meaning the data is NOT normalized
### From here I assume the data is just log-transformed

# recover raw counts
raw_mtx <- expm1(mtx)
str(raw_mtx)
raw_mtx@x <- round(raw_mtx@x)
saveRDS(raw_mtx, "data/GSE206785_raw_matrix.rds")

## > For larger data, use Python and convert the matrix to sparse before reading in
