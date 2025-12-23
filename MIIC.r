library(readr)
library(miic)
library(dplyr)
library(JGL)
library(igraph)
library(ggraph)
library(ggplot2)
library(scales)


# Read the expression matrix
expr <- read_tsv("Data/expr_l1_clean.tsv")
meta <- read_tsv("Data/meta_l1_clean.tsv")

expr_df <- as.data.frame(expr)
expr_df[] <- lapply(expr_df, as.numeric)
rownames(expr_df) <- meta$sample

dim(expr_df)
head(expr_df[,1:5])

# Input of MIIC
miic_input <- cbind(
    expr_df,
    subgroup = meta$subgroup
)
head(miic_input)

# Test run MIIC on a subset of data
expr_test <- expr_df[, 1:50]
miic_test <- miic(
    input_data = expr_test,
    latent = "yes",
    n_shuffles = 10,
    conf_threshold = 0.001,
    n_threads = parallel::detectCores() - 1
)

if(require(igraph)) {
    plot(miic_test, method = "igraph")
}


# Run MIIC
set.seed(128)
miic_res <- miic(
    input_data = miic_input,
    latent = "yes",
    n_shuffles = 1,
    conf_threshold = 0.001,
    n_threads = max(1, parallel::detectCores() - 1)
)

if(require(igraph)) {
    plot(miic_res, method = "igraph")
}

# Save MIIC result
saveRDS(miic_res, "Results/miic_with_subgroup.rds")


