# This parses the metadata files as provided by Wei Sun's team like
# v7_JAX_RNAseq1_Prod.xlsx
parse.jax.xlsx.metadata.file <- function(xlsx.file) {
  cell.line.info <- read.xlsx(xlsx.file, sheet="Clonal cell line", startRow = 4)
  cell.line.info <- cell.line.info[-1,]
  expression.alt.info <- read.xlsx(xlsx.file, sheet="Expression alteration", startRow=4)
  expression.alt.info <- expression.alt.info[-1,]
  diff.product.info <- read.xlsx(xlsx.file, sheet="Differentiated product", startRow=4)
  diff.product.info <- diff.product.info[-1,]
  lib.prep.info <- read.xlsx(xlsx.file, sheet="Library preparation", startRow=4)
  lib.prep.info <- lib.prep.info[-1,]
  seq.file.info <- read.xlsx(xlsx.file, sheet="Sequence file", startRow=4)
  seq.file.info <- seq.file.info[-1,]
  
  # In the new v7_JAX_RNAseq*_Prod.xlsx files, the Clonal cell line and Expression alteration
  # tabs are merged with the common column 'expression_alteration.label'
  # Note that WTs don't have alteration! Hence, we need all.x = TRUE
  old_nrows <- nrow(cell.line.info)
  common.col <- intersect(colnames(cell.line.info), colnames(expression.alt.info))
  stopifnot(length(common.col) == 1)
  stopifnot(all(expression.alt.info[, common.col] %in% cell.line.info[, common.col]))
  info <- merge(cell.line.info, expression.alt.info, all.x=TRUE)
  stopifnot(nrow(info) == old_nrows)
  
  # In the new v7_JAX_RNAseq*_Prod.xlsx files, the Clonal cell line and Differentiated product
  # tabs are merged with the common column 'clonal_cell_line.label'
  # Note that we might have the same cell line under different conditions (e.g., hypoxia vs normoxia)
  # As such we may have more entries in diff.product.info
  common.col <- intersect(colnames(info), colnames(diff.product.info))
  stopifnot(length(common.col) == 1)
  stopifnot(all(info[, common.col] %in% diff.product.info[, common.col]))
  stopifnot(all(diff.product.info[, common.col] %in% info[, common.col]))
  info <- merge(info, diff.product.info)
  stopifnot(nrow(info) == nrow(diff.product.info))
  
  # In the new v7_JAX_RNAseq*_Prod.xlsx files, the Differentiated product and Library preparation
  # tabs are merged with the common differentiated_product.label column.
  old_nrows <- nrow(info)
  common.col <- intersect(colnames(info), colnames(lib.prep.info))
  stopifnot(length(common.col) == 1)
  stopifnot(all(info[, common.col] %in% lib.prep.info[, common.col]))
  stopifnot(all(lib.prep.info[, common.col] %in% info[, common.col]))
  info <- merge(info, lib.prep.info)
  stopifnot(nrow(info) == nrow(lib.prep.info))
  
  # In the new v7_JAX_RNAseq*_Prod.xlsx files, the Library preparation and Sequence file
  # tabs are merged with the common library_preparation.label column.
  old_nrows <- nrow(info)
  common.col <- intersect(colnames(info), colnames(seq.file.info))
  stopifnot(length(common.col) == 1)
  
  # For our purposes, only keep 1 of the 2 FASTQ files
  seq.file.info <- subset(seq.file.info, sequence_file.read_index == "read1")
  
  missing <- info[, common.col][!(info[, common.col] %in% seq.file.info[, common.col])]
  if(length(missing) != 0) {
    print(missing)
    stop()
  }
  #stopifnot(all(info[, common.col] %in% seq.file.info[, common.col]))
  missing <- seq.file.info[, common.col][!(seq.file.info[, common.col] %in% info[, common.col])]
  if(length(missing) != 0) {
    print(missing)
    stop()
  }
  #stopifnot(all(seq.file.info[, common.col] %in% info[, common.col]))
  metadata <- merge(info, seq.file.info)
  stopifnot(nrow(metadata) == nrow(seq.file.info))
  
  metadata$sample <- gsub(metadata$sequence_file.label, pattern="_R1_001.fastq.gz", replacement="")
  metadata[,"model.system"] <- metadata$differentiated_product.model_system
  metadata[,"cell.type"] <- metadata$differentiated_product.model_system
  flag <- grepl(metadata[,"cell.type"], pattern="primitive syn", ignore.case=TRUE)
  metadata[flag, "cell.type"] <- "PrS"
  flag <- grepl(metadata[,"cell.type"], pattern="embryonic mesench", ignore.case=TRUE)
  metadata[flag, "cell.type"] <- "ExM"
  flag <- grepl(metadata[,"cell.type"], pattern="cortical brain", ignore.case=TRUE)
  metadata[flag, "cell.type"] <- "CBO"
  flag <- grepl(metadata[,"cell.type"], pattern="fibroblast", ignore.case=TRUE)
  metadata[flag, "cell.type"] <- "fibroblast"
  
  metadata[,"time.point"] <- metadata$differentiated_product.timepoint_value
  metadata[,"condition"] <- metadata$differentiated_product.wt_control_status
  #flag <- !is.na(metadata$expression_alteration.genes.altered_gene_symbol)
  # GRHL1 is sometimes mispelled as GHRL1
  flag <- !is.na(metadata$expression_alteration.genes.altered_gene_symbol) & ( metadata$expression_alteration.genes.altered_gene_symbol == "GHRL1" )
  metadata[flag,"expression_alteration.genes.altered_gene_symbol"] <- "GRHL1"
  flag <- metadata$condition == "KO"
  metadata[flag,"condition"] <- metadata[flag,"expression_alteration.genes.altered_gene_symbol"]
  metadata[,"strategy"] <- metadata$expression_alteration.genes.editing_strategy
  metadata[is.na(metadata$strategy),"strategy"] <- metadata[is.na(metadata$strategy),"differentiated_product.wt_control_status"]
  metadata[metadata$strategy == "Not applicable","strategy"] <- metadata[metadata$strategy == "Not applicable","differentiated_product.wt_control_status"]
  trans <- list("reversion of termination codon" = "REV PTC", "full coding length" = "Gene", "critical exon" = "CE", "termination codon" = "PTC", "WT" = "None")
  for(nm in names(trans)) {
    flag <- grepl(metadata$strategy, pattern=nm, ignore.case=TRUE)
    metadata[flag,"strategy"] <- trans[[nm]]
  }
  metadata[,"oxygen"] <- metadata$differentiated_product.treatment_condition
  trans <- list("hypoxia" = "Hypoxia", "normoxia" = "Normoxia", "Not applicable" = "Normoxia")
  for(nm in names(trans)) {
    flag <- grepl(metadata$oxygen, pattern=nm, ignore.case=TRUE)
    metadata[flag,"oxygen"] <- trans[[nm]]
  }
  
  metadata
}

# This accounts for the irregularities in the original metadata files like
# JAX_RNAseq_ExtraEmbryonic_metadata.xlsx
# metadata/JAX_RNAseq08.xlsx
process.metadata.file <- function(metadata.file) {
  sheet.names <- getSheetNames(metadata.file)
  seq.metadata <- read.xlsx(metadata.file, sheet = "Sequence file")
  seq.metadata <- seq.metadata[-c(1:4),]
  lib.metadata <- read.xlsx(metadata.file, sheet = "Library preparation")
  lib.metadata <- lib.metadata[-c(1:4),]
  #alt.metadata <- read.xlsx(metadata.file, sheet = "Expression alteration strategy")
  #alt.metadata <- alt.metadata[-c(1:4),]
  cell.line.sheet = sheet.names[sheet.names %in% c("Cell line ", "Clonal cell line")]
  cell.metadata <- read.xlsx(metadata.file, sheet = cell.line.sheet)
  cell.metadata <- cell.metadata[-c(1:4),]
  #cell.metadata <- merge(cell.metadata[, c("CELL.LINE.ID.(Required)", "GENE.EXPRESSION.ALTERATION.PROTOCOL.ID")],
  #                       alt.metadata[, c("GENE.EXPRESSION.ALTERATION.PROTOCOL.ID.(Required)", "ALTERED.GENE.SYMBOLS")],
  #                       by.x = c("GENE.EXPRESSION.ALTERATION.PROTOCOL.ID"),
  #                       by.y = c("GENE.EXPRESSION.ALTERATION.PROTOCOL.ID.(Required)"))
  
  diff.sheet.name <- sheet.names[sheet.names %in% c("Differentiated cell line", "Differentiated product")]
  diff.cell.line.metadata <- read.xlsx(metadata.file, sheet = diff.sheet.name)
  diff.cell.line.metadata <- diff.cell.line.metadata[-c(1:4),]
  
  left_on = colnames(seq.metadata)[grepl(x=colnames(seq.metadata), pattern="input.library.preparation", ignore.case=TRUE)]
  metadata <- merge(seq.metadata, lib.metadata, 
                    by.x = left_on,
                    by.y = c("LIBRARY.PREPARATION.ID.(Required)"))
  stopifnot(nrow(seq.metadata) == nrow(metadata))
  
  left_on = colnames(metadata)[grepl(x=colnames(metadata), pattern=".*differentiated.*id.*", ignore.case=TRUE)]
  right_on = colnames(diff.cell.line.metadata)[colnames(diff.cell.line.metadata) %in% c("DIFFERENTIATED.CELL.LINE.ID.(Required)", "DIFFERENTIATED.PRODUCT.ID.(Required)")]
  metadata <- merge(metadata, diff.cell.line.metadata, by.x = left_on, by.y = right_on)
  stopifnot(nrow(seq.metadata) == nrow(metadata))
  
  # The alteration protocol is not unique to a gene, so this will not work.
  # Instead, parse the KO out of the name
  parental.col <- colnames(cell.metadata)[colnames(cell.metadata) %in% c("PARENTAL.CELL.LINE.NAME.(Required)", "DERIVED.FROM.CELL.LINE.NAME.(Required)")]
  cell.line.col <- colnames(metadata)[colnames(metadata) %in% c("INPUT.CLONAL.CELL.LINE.ID.(Required)", "INPUT.CELL.LINE.ID.(Required)", "Input.CELL.LINE.ID.(Required)")]
  # "X1" for JAX_RNAseq09.xlsx
  clonal.cell.line.id.col <- colnames(cell.metadata)[colnames(cell.metadata) %in% c("X1", "CELL.LINE.ID.(Required)", "CLONAL.CELL.LINE.ID.(Required)")]
  metadata <- merge(metadata, cell.metadata[, c(clonal.cell.line.id.col, parental.col)],
                    by.x = cell.line.col,
                    by.y = clonal.cell.line.id.col)
  
  tbl <- metadata[metadata[,"READ.INDEX.(Required)"] == "read1",]
  
  diff.cell.line.description <- colnames(tbl)[colnames(tbl) %in% c("DIFFERENTIATED.CELL.LINE.DESCRIPTION", "DIFFERENTIATED.PRODUCT.DESCRIPTION")]
  tbl$description <- tbl[,diff.cell.line.description]
  tbl$condition <- tbl$description
  run.id.col <- "RUN.ID"
  if(!(run.id.col %in% colnames(tbl))) { run.id.col <- "LANE.INDEX" }
  tbl$batch <- tbl[, run.id.col]
  
  time.pt.col <- colnames(tbl)[colnames(tbl) %in% c("TIMEPOINT.VALUE", "TIMEPOINT.VALUE.(Required)")]
  model.system.col <- colnames(tbl)[colnames(tbl) %in% c("Model.System", "MODEL.SYSTEM")]

  tbl$cell.type <- "None"
  tbl$strategy <- "None"
  flag <- grepl(tbl$condition, pattern="primitive syncytium")
  tbl[flag, "cell.type"] <- "PrS"
  flag <- grepl(tbl$condition, pattern="extra-embryonic mesenchyme")
  tbl[flag, "cell.type"] <- "ExM"
  flag <- grepl(tbl$condition, pattern="premature termination codon")
  tbl[flag, "strategy"] <- "PTC"
  flag <- grepl(tbl$condition, pattern="full coding region")
  tbl[flag, "strategy"] <- "KO"
  flag <- grepl(tbl$condition, pattern="critical exon")
  tbl[flag, "strategy"] <- "CE"
  tbl$condition <- gsub(tbl$condition, pattern="KOLF2.2[J]? derived ", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern="KOLF2.2[J]? deleted ", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern="derived extra-embryonic mesenchyme[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern="derived primitive syncytium[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern=" by insertion of premature termination codon \\(PTC\\)[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern=" by deletion of full coding region \\(KO\\)[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern=" by deletion of critical exon \\(CE\\)[ ]*", replacement="")
  
  tbl$oxygen <- "Normoxia"
  flag <- grepl(tbl$condition, pattern="hypoxia")
  tbl[flag, "oxygen"] <- "Hypoxia"
  flag <- grepl(tbl$condition, pattern="normoxia")
  tbl[flag, "oxygen"] <- "Normoxia"
  tbl$condition <- gsub(tbl$condition, pattern="[ ]*in normoxia[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern="[ ]*in hypoxia[ ]*", replacement="")
  
  flag <- grepl(tbl$condition, pattern="primitive syncytium")
  tbl[flag, "condition"] <- "WT"
  flag <- grepl(tbl$condition, pattern="extra-embryonic mesenchyme")
  tbl[flag, "condition"] <- "WT"
  
  tbl$condition <- gsub(tbl$condition, pattern="[ ]*primitive syncytium[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern="[ ]*extra-embryonic mesenchyme[ ]*", replacement="")
  
  tbl <- tbl[,c("FILE.NAME.(Required)", "description", parental.col, "oxygen", "condition", "strategy", "cell.type", time.pt.col, model.system.col, "batch")]
  colnames(tbl) <- c("FILE.NAME.(Required)", "description", "genetic.background", "oxygen", "condition", "strategy", "cell.type", "time.point", "model.system", "batch")
  
  tbl
}

make.heatmap <- function(tbl, pops = NULL, zscore = FALSE, legend.name = "Deconvolved\nCell Type\nPercentage", fontsize = 12, plot.batch = TRUE, plot.cell.line = TRUE, ...) {
  suppressPackageStartupMessages(p_load(rcartocolor))
  suppressPackageStartupMessages(p_load(plyr))
  suppressPackageStartupMessages(p_load(Polychrome)) # for Glasbey / colors
  cols <- c("Prim", "EXM", "Early Progenitor", "TE", "ExM", "PrS", "CTB", "Pluripotency", "Peri-TB", "PrSyn", "ExMC", "Cycling-Cell")
  if(!is.null(pops)) {
    cols <- pops
  }
  cols <- cols[cols %in% colnames(tbl)]
  mat <- t(as.matrix(tbl[, cols]))
  if(zscore) mat <- t(scale(t(mat)))
  print(rowSums(mat))
  # Schematic representation of cell types derived from trophectoderm (TE): extra-embryonic mesenchyme (ExM), cytotrophoblast (CTB) and primitive syncytium (PrS). 
  
  genes <- sort(unique(tbl$condition))
  if("WT" %in% genes) genes <- c("WT", genes[genes != "WT"])
  strategies <- sort(unique(tbl$strategy))
  if("None" %in% strategies) strategies <- c("None", strategies[strategies != "None"])
  oxygens <- sort(unique(tbl$oxygen))
  lineages <- sort(unique(tbl$cell.type))
  batches <- c()
  cell.lines <- c()
  batch.column <- "sequence_file.run_id"
  cell.line.column <- "clonal_cell_line.parental_cell_line_name"
  if(plot.cell.line & (cell.line.column %in% colnames(tbl))) {
    cell.lines <- sort(unique(tbl[,cell.line.column]))
  }
  if(plot.batch & (batch.column %in% colnames(tbl))) {
    batches <- sort(unique(tbl[,batch.column]))
  }
  
  vals <- list()
  #vals[["Gene"]] = genes
  #vals[["Strategy"]] = strategies
  #vals[["Oxygen"]] = oxygens
  #vals[["Lineage"]] = lineages
  if (length(genes) > 1) vals[["Gene"]] = genes
  if (length(strategies) > 1) vals[["Strategy"]] = strategies
  if (length(oxygens) > 1) vals[["Oxygen"]] = oxygens
  if (length(lineages) > 1) vals[["Lineage"]] = lineages
  if (length(cell.lines) > 1) vals[["Cell Line"]] = cell.lines
  if (length(batches) > 1) vals[["Batch"]] = batches
  print(vals)
  cols <- llply(vals, 
                .fun = function(vec) {
                  if(length(vec) > 12) {
                    pal = glasbey.colors(length(vec))
                    names(pal) = vec
                  } else if(length(vec) > 2) {
                    pal = carto_pal(length(vec), "Safe")
                    names(pal) = vec
                  } else {
                    pal = c("red", "blue")
                    pal <- pal[1:length(vec)]
                    names(pal) <- vec
                  }
                  pal
                })
  #tiff(file, width = 3 * 480, height = 2 * 480)
  #sz <- 16
  #hm <- Heatmap(mt, name="Meta Program\nNES", row_km=row_kms[[nm]], column_km=column_kms[[nm]], use_raster=TRUE, raster_device="tiff", raster_by_magick = FALSE, raster_resize_mat = TRUE, column_names_rot = 45, column_names_gp = gpar(fontsize = sz), column_title_gp = gpar(fontsize = sz), row_title_gp = gpar(fontsize = sz), column_names_max_height = unit(8, "cm"), heatmap_legend_param = list(title_gp = gpar(fontsize = sz), labels_gp = gpar(fontsize = sz))) + rowAnnotation(sample = sample.labels, treated = treated.labels, col = list(sample = sample_cols, treated = treated_cols), annotation_name_gp= gpar(fontsize = sz), annotation_legend_param = list(sample = list(title_gp = gpar(fontsize = sz), labels_gp = gpar(fontsize = sz)), treated = list(title_gp = gpar(fontsize = sz), labels_gp = gpar(fontsize = sz))))
  
  if("Gene" %in% names(cols)) {
    if("WT" %in% names(cols[["Gene"]])) {
      cols[["Gene"]][["WT"]] <- "black"
    }
  }
  if("Strategy" %in% names(cols)) {
    if("None" %in% names(cols[["Strategy"]])) {
      cols[["Strategy"]][["None"]] <- "black"
    }
  }
  
  col_param = list(Gene = cols[["Gene"]], Strategy = cols[["Strategy"]])
  col_param = cols
  #  if (length(genes) > 2) vals[["Gene"]] = genes
  #  if (length(strategies) > 2) vals[["Strategy"]] = strategies
  #  if (length(oxygens) > 2) vals[["Oxygen"]] = oxygens
  #  if (length(lineages) > 2) vals[["Lineage"]] = lineages
  
  set.seed(1)
  if(length(cols) == 0) { col_param = NULL }
  #print(col_param)
  
  annotation_legend_param = list()
  for(nm in names(vals)) {
    print(vals[[nm]])
    annotation_legend_param[[nm]] = list(title_gp = gpar(fontsize = fontsize), at = vals[[nm]], labels_gp = gpar(fontsize = fontsize))
  }
  old.col.names <- c("condition", "strategy", "oxygen", "cell.type", batch.column, cell.line.column)
  new.col.names <- c("Gene", "Strategy", "Oxygen", "Lineage", "Batch", "Cell Line")
  flag <- new.col.names %in% names(vals)
  
  meta.df <- tbl[, old.col.names[flag], drop=FALSE]
  colnames(meta.df) <- new.col.names[flag]
  column_ha = HeatmapAnnotation(df = meta.df, col = col_param, 
                                annotation_name_gp = gpar(fontsize = fontsize),
                                annotation_legend_param = annotation_legend_param)
  
  Heatmap(mat, name = legend.name, bottom_annotation = column_ha, cluster_rows = FALSE, column_names_rot = 45, ...)
}

calculate.pairwise.degs <- function(obj, rds.file = NULL) {
  if(is.null(rds.file) || !file.exists(rds.file)) {
    
    Idents(obj)  <- "cell.type"
    mex <- calculate.mean.gene.expr.per.cell.type(obj)
    cell.types <- as.character(unique(obj[[]]$cell.type))
    names(cell.types) <- cell.types
    de_res <- 
      llply(cell.types,
            .parallel = FALSE,
            .fun = function(cell.type1) {
              print(cell.type1)
              res <- ldply(cell.types[cell.types != cell.type1],
                           .parallel = FALSE,
                           .fun = function(cell.type2) {
                             print(c(cell.type1, cell.type2))
                             mk <- FindMarkers(obj,
                                               slot = "data",
                                               assay = "SCT",
                                               ident.1 = cell.type1,
                                               ident.2 = cell.type2,
                                               only.pos = FALSE,
                                               min.pct = -Inf,
                                               logfc.threshold = -Inf)
                             mk <- as.data.frame(tibble::rownames_to_column(as.data.frame(mk), "GeneSymbol"))
                             mk               
                           })
              colnames(res)[1] <- "cell.type2"
              res
            })
    de_res_tbl <- ldply(de_res)
    # subset(de_res_tbl, GeneSymbol == "ABI2") ## why don't we have six entries
    colnames(de_res_tbl)[1] <- "cell.type1"
    mex <- mex[, c("cell.type", "Gene", "mean_expr", "relative_rank")]
    mex.tmp <- mex
    colnames(mex.tmp) <- c("cell.type", "Gene", "mean_expr.1", "relative_rank.1")
    de_res_tbl <- merge(de_res_tbl, mex.tmp, by.x = c("cell.type1", "GeneSymbol"), by.y = c("cell.type", "Gene"),
                        all.x = TRUE)
    colnames(mex.tmp) <- c("cell.type", "Gene", "mean_expr.2", "relative_rank.2")
    de_res_tbl <- merge(de_res_tbl, mex.tmp, by.x = c("cell.type2", "GeneSymbol"), by.y = c("cell.type", "Gene"),
                        all.x = TRUE)
    saveRDS(de_res_tbl, rds.file)
  }
  if(!is.null(rds.file)) {
    de_res_tbl <- readRDS(rds.file)
  }
  de_res_tbl
}

process.metadata.file.old <- function(metadata.file) {
  seq.metadata <- read.xlsx(metadata.file, sheet = "Sequence file")
  seq.metadata <- seq.metadata[-c(1:4),]
  lib.metadata <- read.xlsx(metadata.file, sheet = "Library preparation")
  lib.metadata <- lib.metadata[-c(1:4),]
  alt.metadata <- read.xlsx(metadata.file, sheet = "Expression alteration strategy")
  alt.metadata <- alt.metadata[-c(1:4),]
  cell.metadata <- read.xlsx(metadata.file, sheet = "Cell line ")
  cell.metadata <- cell.metadata[-c(1:4),]
  cell.metadata <- merge(cell.metadata[, c("CELL.LINE.ID.(Required)", "GENE.EXPRESSION.ALTERATION.PROTOCOL.ID")],
                         alt.metadata[, c("GENE.EXPRESSION.ALTERATION.PROTOCOL.ID.(Required)", "ALTERED.GENE.SYMBOLS")],
                         by.x = c("GENE.EXPRESSION.ALTERATION.PROTOCOL.ID"),
                         by.y = c("GENE.EXPRESSION.ALTERATION.PROTOCOL.ID.(Required)"))
  
  
  diff.cell.line.metadata <- read.xlsx(metadata.file, sheet = "Differentiated cell line")
  diff.cell.line.metadata <- diff.cell.line.metadata[-c(1:4),]
  
  metadata <- merge(seq.metadata, lib.metadata, 
                    by.x = c("INPUT.LIBRARY.PREPARATION.ID.(Required)"),
                    by.y = c("LIBRARY.PREPARATION.ID.(Required)"))
  stopifnot(nrow(seq.metadata) == nrow(metadata))
  
  metadata <- merge(metadata, diff.cell.line.metadata)
  stopifnot(nrow(seq.metadata) == nrow(metadata))
  
  tbl <- metadata[metadata[,"READ.INDEX.(Required)"] == "read1",]
  
  # The alteration protocol is not unique to a gene, so this will not work.
  # Instead, parse the KO out of the name
  #metadata <- merge(metadata, cell.metadata,
  #                  by.x = c("INPUT.CELL.LINE.ID.(Required)"),
  #                  by.y = c("CELL.LINE.ID.(Required)"))
  tbl$condition <- tbl$`DIFFERENTIATED.CELL.LINE.DESCRIPTION`
  tbl$batch <- tbl$`RUN.ID`
  
  tbl$cell.type <- "None"
  tbl$strategy <- "None"
  flag <- grepl(tbl$condition, pattern="primitive syncytium")
  tbl[flag, "cell.type"] <- "PrS"
  flag <- grepl(tbl$condition, pattern="extra-embryonic mesenchyme")
  tbl[flag, "cell.type"] <- "ExM"
  flag <- grepl(tbl$condition, pattern="premature termination codon")
  tbl[flag, "strategy"] <- "PTC"
  flag <- grepl(tbl$condition, pattern="full coding region")
  tbl[flag, "strategy"] <- "KO"
  flag <- grepl(tbl$condition, pattern="critical exon")
  tbl[flag, "strategy"] <- "CE"
  tbl$condition <- gsub(tbl$condition, pattern="KOLF2.2 derived ", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern="KOLF2.2 deleted ", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern="derived extra-embryonic mesenchyme[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern="derived primitive syncytium[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern=" by insertion of premature termination codon \\(PTC\\)[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern=" by deletion of full coding region \\(KO\\)[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern=" by deletion of critical exon \\(CE\\)[ ]*", replacement="")
  
  tbl$oxygen <- "Normoxia"
  flag <- grepl(tbl$condition, pattern="hypoxia")
  tbl[flag, "oxygen"] <- "Hypoxia"
  flag <- grepl(tbl$condition, pattern="normoxia")
  tbl[flag, "oxygen"] <- "Normoxia"
  tbl$condition <- gsub(tbl$condition, pattern="[ ]*in normoxia[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern="[ ]*in hypoxia[ ]*", replacement="")
  
  flag <- grepl(tbl$condition, pattern="primitive syncytium")
  tbl[flag, "condition"] <- "WT"
  flag <- grepl(tbl$condition, pattern="extra-embryonic mesenchyme")
  tbl[flag, "condition"] <- "WT"
  
  tbl$condition <- gsub(tbl$condition, pattern="[ ]*primitive syncytium[ ]*", replacement="")
  tbl$condition <- gsub(tbl$condition, pattern="[ ]*extra-embryonic mesenchyme[ ]*", replacement="")
  tbl
}
