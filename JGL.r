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

# JGL Network
meta$group <- factor(meta$group)
groups <- unique(meta$group)
print(groups)

data_list <- lapply(groups, function(g){
    mat <- expr_df[meta$group == g, ]
    as.matrix(mat)
})
names(data_list) <- groups
meta <- meta[match(rownames(expr_df), meta$sample), ]

lapply(data_list, dim)
all(rownames(expr_df) == meta$sample)


# Parameter set for JGL
lambda1 <- 0.1
lambda2 <- 0.05

# Run JGL
jgl_fit <- JGL(Y = data_list,
                penalty = "fused",
                lambda1 = lambda1,
                lambda2 = lambda2,
                return.whole.theta = TRUE)


if(require(igraph)) {
    plot(jgl_fit, method = "igraph")
}

# Extract the estimated precision matrices
theta_list <- jgl_fit$theta
length(theta_list)
lapply(theta_list, dim)

# Visualize JGL networks
edges_list <- lapply(1:length(theta_list), function(i){
    theta_mat <- theta_list[[i]]
    edges <- as.data.frame(as.table(theta_mat))
    edges <- edges %>% filter(abs(Freq) > 1e-3)
    colnames(edges) <- c("from", "to", "weight")
    edges$group <- names(theta_list)[i]
    edges
})

edges_all <- bind_rows(edges_list)
head(edges_all)

# Node
theta_mat <- theta_list[[1]]
g_jgl <- graph.adjacency(abs(theta_mat) > 1e-3, 
    mode = "undirected", diag = FALSE)

V(g_jgl)$degree <- degree(g_jgl)
V(g_jgl)$size <- 3 + 5 * (V(g_jgl)$degree / max(V(g_jgl)$degree))
cancer_genes <- character(0)
V(g_jgl)$color <- ifelse(V(g_jgl)$name %in% cancer_genes, "red", "skyblue")

edges <- get.edges(g_jgl, E(g_jgl))
E(g_jgl)$weight <- abs(theta_mat[edges])

# Visualize JGL network
p_jgl <- ggraph(g_jgl, layout = "fr") +
    geom_edge_link(aes(width = weight), alpha = 0.5, color = "grey50") +
    geom_node_point(aes(size = size, color = ifelse(name %in% cancer_genes, "red", "skyblue"))) +
    geom_node_text(aes(label = ifelse(name %in% cancer_genes, name, "")),
                    repel = TRUE, size = 3, fontface = "bold") +
    scale_edge_size_continuous(range = c(0.2, 2)) +
    scale_size_continuous(range = c(2,8)) +
    theme_void() +
    theme(legend.position = "none")

p_jgl

# Save JGL networks
ggsave("JGL_Tumor_network_paper.png", 
    plot = p_jgl, width = 12, height = 12, dpi = 600)

ggsave("JGL_Tumor_network_paper.pdf", plot = p_jgl, width = 12, height = 12)


# Difference network
theta_tumor  <- theta_list[[1]]
theta_normal <- theta_list[[2]]

theta_diff <- theta_tumor - theta_normal

thresh <- 1e-3
edge_idx <- which(abs(theta_diff) > thresh, arr.ind = TRUE)

edges_diff <- data.frame(
    from = rownames(theta_diff)[edge_idx[,1]],
    to   = colnames(theta_diff)[edge_idx[,2]],
    weight = theta_diff[edge_idx]
)

edges_diff <- edges_diff[edges_diff$from != edges_diff$to, ]
edges_diff$weight_plot <- abs(edges_diff$weight)
edges_diff$weight_plot[edges_diff$weight_plot == 0] <- 1e-5

g_diff <- graph_from_data_frame(edges_diff, directed = FALSE)
E(g_diff)$weight <- edges_diff$weight_plot
V(g_diff)$degree <- degree(g_diff)
V(g_diff)$size <- 2 + 6 * (V(g_diff)$degree / max(V(g_diff)$degree))

cancer_genes <- c("TP53", "MYC", "CDK6")
V(g_diff)$color <- ifelse(V(g_diff)$name %in% cancer_genes, "red", "skyblue")

# Visualize difference network
p_diff <- ggraph(g_diff, layout = "fr") +
    geom_edge_link(aes(width = abs(weight)), alpha = 0.6, color = "grey70") +
    geom_node_point(aes(size = size, color = color)) +
    geom_node_text(aes(label = ifelse(name %in% cancer_genes, name, "")),
                    repel = TRUE, size = 3.5, fontface = "bold", segment.color = "grey50") +
    scale_edge_size_continuous(range = c(0.2, 2)) +
    scale_size_continuous(range = c(2, 8)) +
    theme_void() +
    theme(legend.position = "none")

print(p_diff)


ggsave("JGL_diff_network_paper.png", 
    plot = p_diff, width = 12, height = 12, dpi = 600)

ggsave("JGL_diff_network_paper.pdf", plot = p_diff, width = 12, height = 12)