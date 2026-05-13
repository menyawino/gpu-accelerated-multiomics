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

de_table <- normalizePath(get_arg_value("--de-table"), mustWork = TRUE)
pathway_table <- normalizePath(get_arg_value("--pathway-table"), mustWork = TRUE)
out_dir <- normalizePath(get_arg_value("--out-dir", "results/real_processed/hf_vs_nonhf_pathways/topology"), mustWork = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
old_wd <- getwd()
setwd(out_dir)
on.exit(setwd(old_wd), add = TRUE)

de <- read.delim(de_table, check.names = FALSE, stringsAsFactors = FALSE)
pathways <- read.delim(pathway_table, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("gene_id", "mean_diff")
missing_cols <- setdiff(required_cols, colnames(de))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

de <- de[!is.na(de$gene_id) & !is.na(de$mean_diff), c("gene_id", "mean_diff")]
de <- de[order(-abs(de$mean_diff)), ]
de <- de[!duplicated(de$gene_id), ]

symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(de$gene_id),
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "SYMBOL"
)
symbol_map <- symbol_map[!is.na(symbol_map$ENTREZID) & !duplicated(symbol_map$SYMBOL), ]

merged <- merge(de, symbol_map, by.x = "gene_id", by.y = "SYMBOL", all.x = FALSE, all.y = FALSE)
merged <- merged[order(-abs(merged$mean_diff)), ]
merged <- merged[!duplicated(merged$ENTREZID), ]

gene_data <- merged$mean_diff
names(gene_data) <- merged$ENTREZID

manifest_rows <- list()

for (i in seq_len(nrow(pathways))) {
  row <- pathways[i, ]
  suffix <- paste0(
    tolower(gsub("[^A-Za-z0-9]+", "_", row$category)),
    "_",
    row$pathway_id
  )

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
    mapped_gene_count = sum(merged$ENTREZID %in% names(gene_data)),
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest_rows)
write.table(manifest, file.path(out_dir, "pathview_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

message("Wrote selected pathview outputs to ", out_dir)