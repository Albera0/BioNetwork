library(tidyverse)
library(readxl)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(glmnet)
library(umap)

# Read the expression matrix
expr <- read_tsv("Data/GSE124814_HW_expr_matrix.tsv")

# turn the column into row and set gene names as row names
expr_mat <- expr %>% column_to_rownames(var = colnames(expr)[1])
expr_mat <- as.matrix(expr_mat)

dim(expr_mat)

# Read the sample information
meta <- read_xlsx("Data/GSE124814_sample_descriptions.xlsx")

# Read the group/subgroup information from metadata
meta$group <- meta[["source name"]]
meta$subgroup <- meta[["characteristics: subgroup relabeled"]]
meta$sample <- meta[["Sample name"]]

# Make sure the data is normalized
summary(as.vector(expr_mat))
expr_vec <- as.vector(expr_mat)

ggplot(data.frame(expr = expr_vec), aes(x = expr)) +
    geom_histogram(bins = 80, fill = "skyblue", color = "black") +
    theme_minimal() +
    labs(title = "Expression Value Distribution",
        x = "Expression value",
        y = "Count")

set.seed(128)

# Filter out lowly expressed genes
vars <- apply(expr_mat, 1, var)
topN <- 3000
expr_filt <- expr_mat[order(-vars)[1:topN], ]

dim(expr_mat)
dim(expr_filt)
length(meta$group)

# Lasso feature selection
X <- t(expr_filt)
y <- as.factor(meta$group)

cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1)
coef_sel <- coef(cvfit, s = "lambda.min")
selected_genes <- rownames(coef_sel)[which(coef_sel[, 1] != 0)][-1]
expr_l1 <- expr_filt[selected_genes, ]

dim(expr_l1)

# Calculate variance for each gene
vars_l1 <- apply(expr_l1, 1, var)

# Select top 50 most variable genes
top50_var <- names(sort(vars_l1, decreasing = TRUE))[1:50]

expr_l1_50 <- expr_l1[top50_var, ]

dim(expr_l1_50)

# PCA analysis
meta_l1 <- meta[match(colnames(expr_l1_50), meta$sample), ]

length(meta_l1$group)
ncol(expr_l1_50)
all(colnames(expr_l1_50) == meta_l1$sample)
meta_l1$group <- factor(meta_l1$group)
pca <- PCA(t(expr_l1_50), graph = FALSE)

dim(pca$ind$coord)
pc_scores <- pca$ind$coord[, 1:5]

fviz_pca_ind(pca, habillage = meta_l1$group, 
    addEllipses = TRUE, repel = TRUE, label = "none")

# UMAP analysis
umap_res <- umap(t(expr_l1_50))

umap_df <- data.frame(
    UMAP1 = umap_res$layout[,1],
    UMAP2 = umap_res$layout[,2],
    sample = meta_l1$sample,
    group = meta_l1$group
)

plot(umap_res$layout, col = as.factor(meta_l1$group), pch = 19,
    xlab = "UMAP1", ylab = "UMAP2")
legend("bottomleft", legend = levels(as.factor(meta_l1$group)),
    col = 1:length(levels(as.factor(meta_l1$group))), pch = 19)

# Clustering analysis
k <- length(unique(meta_l1$group))
km <- kmeans(pc_scores, centers = k)

# Visualize clustering results
plot(pc_scores, col = km$cluster, pch = 19,
    xlab = "PC1", ylab = "PC2")

# Remove outliers 
umap_df <- umap_df %>%
    mutate(
    z1 = scale(UMAP1),
    z2 = scale(UMAP2),
    is_outlier = (abs(z1) > 3) | (abs(z2) > 3)
    )

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = is_outlier)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = c("FALSE"="black", "TRUE"="red")) +
    theme_minimal() +
    labs(title = "Doublets / Outliers detection",
        color = "Outlier")

# Remove outliers from expression matrix and metadata
expr_l1_clean <- expr_l1_50[, !umap_df$is_outlier]
meta_l1_clean <- meta_l1[!umap_df$is_outlier, ]

dim(expr_l1_clean)
dim(meta_l1_clean)

# Prepare metadata for MIIC and only keep subgroup information
meta_miic <- meta_l1_clean %>%
    select(sample, subgroup)

meta_miic$subgroup <- as.factor(meta_miic$subgroup)

# MIIC Web input file
miic_input <- t(expr_l1_clean) %>%
    as.data.frame()

head(miic_input)

write_tsv(miic_input, "Data/miic_input_web.tsv")

# Save the processed expression matrix for JGL and MIIC
expr_l1_clean <- t(expr_l1_clean)

write_tsv(as.data.frame(expr_l1_clean), "Data/expr_l1_clean.tsv")
write_tsv(meta_miic, "Data/meta_l1_clean.tsv")

# Quick check of the saved metadata
meta_test <- read_tsv("Data/meta_l1_clean.tsv")

table(meta_test$subgroup)
sum(is.na(meta_test$subgroup))