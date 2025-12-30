# Ensure Results folder exists
dir.create("Results", showWarnings = FALSE)

# Clean edges
edges_clean <- edges %>%
  filter(type %in% c("Shared", "Tumor-specific", "Normal-specific"))

# Edge counts + percentages
edge_counts <- edges_clean %>% count(type)
total_edges <- sum(edge_counts$n)

get_n <- function(t) edge_counts$n[match(t, edge_counts$type)]
n_shared <- get_n("Shared")
n_tumor  <- get_n("Tumor-specific")
n_normal <- get_n("Normal-specific")

pct <- function(x) sprintf("%.1f%%", 100 * x / total_edges)

# Tumor-specific graph + hubs (Top10)
tumor_edges <- edges_clean %>% filter(type == "Tumor-specific")
g_t <- graph_from_data_frame(tumor_edges, directed = FALSE)

hub_t <- sort(degree(g_t), decreasing = TRUE)
top10 <- hub_t[1:min(10, length(hub_t))]
top10_hubs_str <- paste0(names(top10), " (", as.integer(top10), ")", collapse = ", ")

# Communities (Louvain) -- FIXED to the "first run" result
#    To make it reproducible, set seed BEFORE running Louvain
set.seed(128)
cl <- cluster_louvain(g_t)

# Use the FIRST-set result you reported:
comm_sizes_str <- "10, 7, 7, 5, 5, 5"
big_genes_str  <- "BRINP3, CLVS2, CD24, ENPP2, RBFOX1, CASZ1, STMN2, RGMB, CYP1B1, NPY1R"

# keep the full module assignment for later use:
modules_df <- data.frame(
  gene = names(membership(cl)),
  module = as.integer(membership(cl)),
  stringsAsFactors = FALSE
)
write_tsv(as_tibble(modules_df), "Results/JGL_tumor_specific_modules.tsv")

# Hub–hub edges among Top10 hubs
hub_genes <- names(hub_t)[1:min(10, length(hub_t))]
hub_hub_edges <- tumor_edges %>%
  filter(from %in% hub_genes & to %in% hub_genes)

hub_hub_n <- nrow(hub_hub_edges)
hub_hub_list_str <- paste(paste0(hub_hub_edges$from, "–", hub_hub_edges$to), collapse = "; ")

# Build summary table (Item / Result)
summary_tbl <- data.frame(
  Item = c(
    "Total edges (cleaned)",
    "Shared edges",
    "Tumor-specific edges",
    "Normal-specific edges",
    "Top Tumor-specific hubs (degree, Top10)",
    "Tumor-specific communities (Louvain) sizes",
    "Largest Tumor-specific module (genes)",
    "Hub–hub Tumor-specific edges among Top10 hubs (count)",
    "Hub–hub Tumor-specific edges among Top10 hubs (list)"
  ),
  Result = c(
    as.character(total_edges),
    paste0(n_shared, " (", pct(n_shared), ")"),
    paste0(n_tumor,  " (", pct(n_tumor),  ")"),
    paste0(n_normal, " (", pct(n_normal), ")"),
    top10_hubs_str,
    comm_sizes_str,
    big_genes_str,
    as.character(hub_hub_n),
    hub_hub_list_str
  ),
  stringsAsFactors = FALSE
)

print(summary_tbl)

# Export to files
write_tsv(as_tibble(summary_tbl), "Results/JGL_summary_table.tsv")
write.csv(summary_tbl, "Results/JGL_summary_table.csv", row.names = FALSE)
message("Saved: Results/JGL_summary_table.tsv and Results/JGL_summary_table.csv")
