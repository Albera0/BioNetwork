library(readr)
library(dplyr)
library(JGL)
library(igraph)
library(ggraph)
library(ggplot2)
library(scales)


# Read the expression matrix
expr <- read_tsv("Data/expr_l1_clean.tsv")

expr_df <- as.data.frame(expr)
expr_df[] <- lapply(expr_df, as.numeric)
rownames(expr_df) <- expr$gen

dim(expr_df)
head(expr_df[,1:5])

mat <- as.matrix(expr_df)
dim(mat)


# Run JGL
jgl_fit <- JGL(
    Y = list(mat),
    penalty = "group",
    lambda1 = 0.1,
    lambda2 = 0
)

if(require(igraph)) {
    plot(jgl_fit, method = "igraph")
}


# # Save edge list
edge_df <- as_data_frame(g, what = "edges")
write_tsv(edge_df, "Results/jgl_edges_all_samples.tsv")


# Save the JGL model
saveRDS(jgl_fit, "Results/jgl_all_samples.rds")

