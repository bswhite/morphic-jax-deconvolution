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

make.heatmap <- function(tbl, pops = NULL, zscore = FALSE, legend.name = "Deconvolved\nCell Type\nPercentage", fontsize = 12, ...) {
  p_load(rcartocolor)
  suppressPackageStartupMessages(p_load(plyr))
  cols <- c("Prim", "EXM", "Early Progenitor", "TE", "ExM", "PrS", "CTB", "Pluripotency", "Peri-TB", "PrSyn", "ExMC", "Cycling-Cell")
  if(!is.null(pops)) {
    cols <- pops
  }
  cols <- cols[cols %in% colnames(tbl)]
  mat <- t(as.matrix(tbl[, cols]))
  if(zscore) mat <- t(scale(t(mat)))
  print(rowSums(mat))
  # Schematic representation of cell types derived from trophectoderm (TE): extra-embryonic mesenchyme (ExM), cytotrophoblast (CTB) and primitive syncytium (PrS). 
  
  p_load(rcartocolor)
  genes <- sort(unique(tbl$condition))
  strategies <- sort(unique(tbl$strategy))
  if("None" %in% strategies) strategies <- c(strategies[strategies != "None"], "None")
  oxygens <- sort(unique(tbl$oxygen))
  lineages <- sort(unique(tbl$cell.type))
  
  vals <- list()
  #vals[["Gene"]] = genes
  #vals[["Strategy"]] = strategies
  #vals[["Oxygen"]] = oxygens
  #vals[["Lineage"]] = lineages
  if (length(genes) > 1) vals[["Gene"]] = genes
  if (length(strategies) > 1) vals[["Strategy"]] = strategies
  if (length(oxygens) > 1) vals[["Oxygen"]] = oxygens
  if (length(lineages) > 1) vals[["Lineage"]] = lineages
  cols <- llply(vals, 
                .fun = function(vec) {
                  if(length(vec) > 2) {
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
    annotation_legend_param[[nm]] = list(title_gp = gpar(fontsize = fontsize), labels_gp = gpar(fontsize = fontsize))
  }
  old.col.names <- c("condition", "strategy", "oxygen", "cell.type")
  new.col.names <- c("Gene", "Strategy", "Oxygen", "Lineage")
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
