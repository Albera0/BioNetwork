library(readr)
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
rownames(expr_df) <- expr$gen

dim(expr_df)
head(expr_df[,1:5])

# Define cell types based on subgroup
meta$cell_type <- ifelse(
    is.na(meta$subgroup),
    "Normal",
    "Tumor"
)

dim(meta)
head(meta[,1:2])

mat_normal <- expr_df[meta$cell_type == "Normal", , drop = FALSE]
mat_tumor  <- expr_df[meta$cell_type == "Tumor", , drop = FALSE]

dim(mat_normal)
dim(mat_tumor)

head(mat_normal[,1:5])
head(mat_tumor[,1:5])

# Run JGL
jgl_fit <- JGL(
    Y = list(
        Normal = as.matrix(mat_normal),
        Tumor  = as.matrix(mat_tumor)
    ),
    penalty = "fused",
    lambda1 = 0.1,
    lambda2 = 0.1,
    return.whole.theta = TRUE
)

if(require(igraph)) {
    plot(jgl_fit, method = "igraph")
}


# Precision matrix
theta_normal <- jgl_fit$theta[[1]]
theta_tumor  <- jgl_fit$theta[[2]] 

# Create the edge list from the precision matrix
thr <- 1e-6
adj_normal <- abs(theta_normal) > thr
adj_tumor  <- abs(theta_tumor)  > thr
diag(adj_normal) <- 0
diag(adj_tumor)  <- 0

genes <- colnames(mat_normal)

# Create edge data frame
gene_idx <- setNames(seq_along(genes), genes)

edge_df <- expand.grid(
    from = genes,
    to   = genes,
    stringsAsFactors = FALSE
    ) %>%
    filter(from < to) %>%
    mutate(
        normal_idx = gene_idx[from],
        tumor_idx  = gene_idx[to],
        normal = adj_normal[cbind(normal_idx, tumor_idx)],
        tumor  = adj_tumor[cbind(normal_idx, tumor_idx)],
        type = case_when(
        normal & tumor  ~ "Shared",
        tumor  & !normal ~ "Tumor-specific",
        normal & !tumor  ~ "Normal-specific",
        TRUE             ~ NA_character_
        )
    ) %>%
    filter(!is.na(type)) %>%
    select(from, to, type)

# Create igraph object
g <- graph_from_data_frame(
    edge_df,
    directed = FALSE,
    vertices = data.frame(name = genes)
)

# Calculate node degree and set node size
V(g)$degree <- degree(g)
V(g)$size <- rescale(V(g)$degree, to = c(2.5, 7))

# Draw the network
set.seed(128)

p <- ggraph(g, layout = "fr") +
    geom_edge_link(
        aes(color = type),
        alpha = 0.6,
        width = 0.5
    ) +
    scale_edge_color_manual(
        values = c(
            "Shared"          = "grey70",
            "Tumor-specific"  = "#D73027",
            "Normal-specific" = "#4575B4"
        )
    ) +
    geom_node_point(
        aes(size = degree),
        color = "black"
    ) +
    geom_node_text(
        aes(label = name),
        size = 2.2,
        repel = TRUE,
        max.overlaps = 40
    ) +
    scale_size(range = c(2.5, 7)) +
    theme_void() +
    theme(
        legend.position = "right",
        legend.title = element_blank()
    )

# Save the plot
ggsave(
    "Results/JGL_Normal_Tumor_overlay_network.pdf",
    p,
    width = 9,
    height = 9,
    device = cairo_pdf
)

ggsave(
    filename = "Results/JGL_Normal_Tumor_overlay_network_flat.png",
    plot = p,
    width = 9,
    height = 9,
    dpi = 300
)

# Save edge list
edge_df <- as_data_frame(g, what = "edges")
write_tsv(edge_df, "Results/jgl_edges_all_samples.tsv")


# Save the JGL model
saveRDS(jgl_fit, "Results/jgl_all_samples.rds")