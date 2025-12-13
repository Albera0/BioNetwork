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
miic_all <- miic(
    input_data = expr_df,
    latent = "yes",
    n_shuffles = 2,
    conf_threshold = 0.001,
    n_threads = parallel::detectCores() - 1
)

if(require(igraph)) {
    plot(miic_all, method = "igraph")
}

# Plot MIIC network
edges_miic <- miic_test$edges
edges_miic$weight <- as.numeric(as.character(edges_miic$weight))
edges_miic <- edges_miic[abs(edges_miic$weight) > 0.001, ]
g_miic <- graph_from_data_frame(edges_miic, directed = FALSE)

theta_diff_miic <- miic_test$theta_list$Medulloblastoma - miic_test$theta_list$Normal
diff_score <- apply(abs(theta_diff_miic), 1, sum)

common_genes <- intersect(V(g_miic)$name, names(diff_score))
V(g_miic)$diff_score <- NA
V(g_miic)$diff_score[match(common_genes, V(g_miic)$name)] <- diff_score[common_genes]

V(g_miic)$degree <- degree(g_miic)
V(g_miic)$size <- 3 + 5 * (V(g_miic)$degree / max(V(g_miic)$degree))

cancer_genes <- c("TP53","MYC","CDK6")

V(g_miic)$color <- NA
V(g_miic)$color[!is.na(V(g_miic)$diff_score)] <- col_numeric("RdBu", 
    domain = range(diff_score, na.rm = TRUE))(V(g_miic)$diff_score[!is.na(V(g_miic)$diff_score)])

V(g_miic)$color[V(g_miic)$name %in% cancer_genes] <- "red"

p_miic <- ggraph(g_sub, layout = "kk") +   # Kamada-Kawai 布局
    geom_edge_link(aes(width = weight_plot), alpha = 0.6, color = "grey50") +
    geom_node_point(aes(size = size, color = color)) +
    geom_node_text(aes(label = ifelse(name %in% cancer_genes, name, "")),
                    repel = TRUE, size = 3.5, fontface = "bold", segment.color = "grey50") +
    scale_size_continuous(range = c(2, 8)) +
    theme_void() +
    theme(legend.position = "right")

print(p_miic)


# Save MIIC network plot
ggsave("Results/MIIC_network_paper.png", plot = p_miic, width = 12, height = 12, dpi = 600)
ggsave("Results/MIIC_network_paper.pdf", plot = p_miic, width = 12, height = 12)






