#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(pathview)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) {
    return(args[idx + 1])
  }
  default
}

de_table <- get_arg_value("--de-table", "results/real_processed/reprogramming/all_genes_failing_vs_nonfailing.tsv")
out_dir <- get_arg_value("--out-dir", "results/real_processed/pathview")
de_table <- normalizePath(de_table, mustWork = TRUE)
out_dir <- normalizePath(out_dir, mustWork = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
old_wd <- getwd()
setwd(out_dir)
on.exit(setwd(old_wd), add = TRUE)

message("Reading differential table: ", de_table)
de <- read.delim(de_table, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("gene_id", "mean_diff")
missing_cols <- setdiff(required_cols, colnames(de))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

de <- de[!is.na(de$gene_id) & !is.na(de$mean_diff), c("gene_id", "mean_diff", intersect(c("fdr", "pvalue"), colnames(de)))]
de <- de[order(-abs(de$mean_diff)), ]
de <- de[!duplicated(de$gene_id), ]

symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(de$gene_id),
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "SYMBOL"
)

symbol_map <- symbol_map[!is.na(symbol_map$ENTREZID) & !duplicated(symbol_map$SYMBOL), ]
write.table(symbol_map, file.path(out_dir, "symbol_to_entrez.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

merged <- merge(de, symbol_map, by.x = "gene_id", by.y = "SYMBOL", all.x = FALSE, all.y = FALSE)
merged <- merged[order(-abs(merged$mean_diff)), ]
merged <- merged[!duplicated(merged$ENTREZID), ]

gene_data <- merged$mean_diff
names(gene_data) <- merged$ENTREZID

write.table(
  data.frame(entrez_id = names(gene_data), mean_diff = unname(gene_data), gene_id = merged$gene_id),
  file.path(out_dir, "pathview_gene_input.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pathways <- data.frame(
  category = c("Inflammation", "Fibrosis", "ECM Remodeling", "ECM Remodeling"),
  pathway_id = c("hsa04060", "hsa04350", "hsa04512", "hsa04510"),
  pathway_name = c(
    "Cytokine-cytokine receptor interaction",
    "TGF-beta signaling pathway",
    "ECM-receptor interaction",
    "Focal adhesion"
  ),
  stringsAsFactors = FALSE
)

manifest_rows <- list()

for (i in seq_len(nrow(pathways))) {
  row <- pathways[i, ]
  suffix <- paste0(
    tolower(gsub("[^A-Za-z0-9]+", "_", row$category)),
    "_",
    row$pathway_id
  )
  message("Rendering ", row$pathway_id, " (", row$pathway_name, ")")

  pv <- pathview(
    gene.data = gene_data,
    pathway.id = row$pathway_id,
    species = "hsa",
    gene.idtype = "entrez",
    out.suffix = suffix,
    kegg.native = TRUE,
    multi.state = FALSE,
    same.layer = TRUE,
    low = list(gene = "#2166ac"),
    mid = list(gene = "#f7f7f7"),
    high = list(gene = "#b2182b")
  )

  png_candidate <- file.path(out_dir, paste0(row$pathway_id, ".", suffix, ".png"))
  xml_candidate <- file.path(out_dir, paste0(row$pathway_id, ".xml"))
  manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
    category = row$category,
    pathway_id = row$pathway_id,
    pathway_name = row$pathway_name,
    png_file = if (file.exists(png_candidate)) basename(png_candidate) else NA_character_,
    xml_file = if (file.exists(xml_candidate)) basename(xml_candidate) else NA_character_,
    mapped_genes = sum(merged$ENTREZID %in% names(gene_data)),
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest_rows)
write.table(manifest, file.path(out_dir, "pathview_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

message("Wrote pathview outputs to ", out_dir)