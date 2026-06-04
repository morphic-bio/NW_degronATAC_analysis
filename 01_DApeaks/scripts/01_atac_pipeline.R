################################################################################################
# atac_pipeline.R
# Author: Amanda Everitt
#
# Usage:
  # Rscript 00_atac_pipeline.R \
  #   -meta ../00_inputs/YY1.metadata.tsv \
  #   -o ../01_outputs \
  #   --cores 4


################################################################################################

################################################################################################

library(optparse)

parser <- OptionParser()

parser <- add_option(parser, c("-m", "--meta"), type="character",  default=NULL, action="store",
                     help="Path to metadata csv",
                     dest="atac_meta")
parser <- add_option(parser, c("-o", "--outdir"), type="character",  default=NULL, action="store",
                     help="Path to output diretory",
                     dest="out_dir")
parser <- add_option(parser, c("--cores"), type="integer",  default=1, action="store",
                     help="number of cores for mclapply",
                     dest="ncores")
parser <- add_option(parser, c("--seed"), type="integer",  default=1, action="store",
                     help="set seed value for reproducibility",
                     dest="seed")
parser <- add_option(parser, c("--skipMotifEnrich"), action="store_true", default=FALSE,
                     help="Skip motif enrichment step",
                     dest="skip_motif")
parser <- add_option(parser, c("--skipGSEA"), action="store_true", default=FALSE,
                     help="Skip GSEA step",
                     dest="skip_gsea")
parser <- add_option(parser, c("--skipFixedPeaks"), action="store_true", default=FALSE,
                     help="Skip fixed-width peak centering step",
                     dest="skip_fixed_peaks")
parser <- add_option(parser, c("--saveOnFail"), action="store_true", default=FALSE,
                     help="Save workspace image on failure for debugging",
                     dest="save_on_fail")

opt <- parse_args(parser)
`%ni%` <- Negate(`%in%`)

################################################################################################
##amanda testing section.

# opt = list(
#   atac_meta = "../00_inputs/CCNC.metadata.tsv",
#   out_dir   = "../01_outputs",
#   ncores    = 1,
#   seed      = 42,
#   skip_fixed_peaks = TRUE,
#   skip_gsea = TRUE,
#   skip_motif = TRUE
#   )
# script_dir = ""
# source(file.path(script_dir, "01_atac_pipeline_functions.R"))

################################################################################################

# Locate script directory and source shared functions (assumed to be in same location)
args_cmd    <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_cmd[grep("--file=", args_cmd)])
script_dir  <- ifelse(length(script_path) && dir.exists(dirname(script_path)),
                      dirname(normalizePath(script_path)), "./")

source(file.path(script_dir, "01_atac_pipeline_functions.R"))

################################################################################################

#1. Load dependencies
pkgs <- c("optparse","dplyr", "tidyr", "ggplot2", "GenomicRanges", "SummarizedExperiment", "ComplexHeatmap",
          "rtracklayer", "stringr", "Rsubread", "DESeq2", "knitr", "ggrepel", "rmarkdown", "gridExtra",
          "GenomeInfoDb", "Rsamtools","TxDb.Hsapiens.UCSC.hg38.knownGene", "ChIPseeker", "org.Hs.eg.db",
          "reshape2", "BSgenome.Hsapiens.UCSC.hg38", "TFBSTools", "monaLisa","JASPAR2020", "parallel",
          "kableExtra", "clusterProfiler")

pipeline_complete <- FALSE
tryCatch({ #entire script wrapper for save_on_fail parameter.

timed_step("Loading dependencies", {
  check_dependencies(pkgs)
  invisible(lapply(pkgs, function(x)
    suppressPackageStartupMessages(library(x, character.only = TRUE))
  ))
})

set.seed(opt$seed)

##########################################################################

#2. Collect metadata
timed_step("Collecting metadata", {
  meta_df   <- read.table(opt$atac_meta, header = TRUE)
  check_meta(meta_df)
  gene_name <- unique(meta_df$gene)
  out_dir   <- create_unique_dir(file.path(opt$out_dir, gene_name))
  meta_df   <- addExtraInfo(meta_df, ncores = opt$ncores)
})

#3. Read in peaks, annotate, add summit info
timed_step("Reading in peak files", {
  tmp <- unique(meta_df[, c("timepoint", "peak_path")])
  tmp <- tmp[tmp$peak_path != "NULL", ]
  consensus_peaks <- list()
  
  #Read in condition peaks. 
  consensus_peaks[["condition_peaks"]] <- setNames(
    lapply(tmp$peak_path, function(p) {
      gr <- rtracklayer::import(p)
      keepStandardChromosomes(gr, pruning.mode = "coarse")
    }),
    tmp$timepoint
  )
  
  #Merge the union set. 
  consensus_peaks[["dataset_peaks"]] <- GenomicRanges::reduce(
    unlist(GRangesList(consensus_peaks[["condition_peaks"]])),
    ignore.strand = TRUE,
    min.gapwidth=10
  )
  
  #Annotate and retain info about all sample-specific summits locations. 
  consensus_peaks[["dataset_peaks"]] <- annotate_peaks(consensus_peaks[["dataset_peaks"]])
  consensus_peaks <- add_summitInfo(consensus_peaks)
})

#4. Get read counts per peak for each bam
timed_step("Generating read x peak matrix", {
  fc_object <- run_feature_counts(peaks = consensus_peaks[["dataset_peaks"]], 
                                  bams = meta_df$bam_path[meta_df$bam_path != "NULL"], 
                                  n_threads=opt$ncores)
})

#5. Filter to well covered peaks
timed_step("Filtering Peaks & Running DESEQ2 sweep", {
  deseq_sweep <- run_deseq_sweep(
    counts = fc_object$counts,
    metadata = meta_df,
    cpm_thresholds = c(0, 0.25, 0.5, 0.75, 1, 2, 3, 5, 10),
    replicate_num = 2, 
    default_thresh=1)

  consensus_peaks[["dataset_peaks_filt"]] <- subset(
    consensus_peaks[["dataset_peaks"]],
    peak_coord %in% rownames(deseq_sweep$count_matrix)
  )
  
  
  deseq_sweep$deseq_object$res_list_annot = list()
  for (i in names(deseq_sweep$deseq_object$res_list)){
    deseq_sweep$deseq_object$res_list_annot[[i]] = clean_de_results(deseq_df = deseq_sweep$deseq_object$res_list[[i]], 
                                                   peaks_annot = consensus_peaks[["dataset_peaks_filt"]])
  }
  deseq_object <- deseq_sweep$deseq_object
  
  rm(fc_object); gc()
})


#6. Add fixed-width peaks centered on summits
if (opt$skip_fixed_peaks) {
  message("##### Skipping fixed-width peaks (--skipFixedPeaks) #####")
} else {
  timed_step("Adding fixed-width peaks", {
    consensus_peaks[["dataset_peaks_filt"]] <- add_fixed_width(
      peaks                = consensus_peaks[["dataset_peaks_filt"]],
      bam_paths            = unique(meta_df$merged_bam_path),
      width                = 500,
      summit_spread_thresh = 10,
      ncores               = opt$ncores
    )
  })
}


#9. Motif enrichment
has_hits <- vapply(deseq_object$res_list, anyhits, logical(1))

if (opt$skip_motif) {
  message("##### Skipping motif enrichment (--skipMotifEnrich) #####")
  enrich_list <- list()
  
} else if (!any(has_hits)) {
  message("##### No datasets passed threshold — skipping motif enrichment #####")
  enrich_list <- list()
  
} else {
  enrich_list <- timed_step("Performing motif enrichment", {
    
    peaks_gr <- consensus_peaks[["dataset_peaks_filt"]]
    if ("peak_coord_500bp" %in% colnames(peaks_gr@elementMetadata)) {
      message("Fixed width peak column `peak_coord_500bp` found. Using summit-centered coordinates.")
      ranges(peaks_gr) <- ranges(GRanges(peaks_gr$peak_coord_500bp))
    } else {
      message("No fixed width peak column found. Standardizing to 500bp based on center.")
      peaks_gr <- resize(peaks_gr, width = 500, fix = "center")
    }

    input_info      <- collect_input_info(gr=peaks_gr, ncores=opt$ncores)
    conds_with_hits <- names(has_hits)[has_hits]

    setNames(lapply(conds_with_hits, function(cond) {
      df <- as.data.frame(deseq_object$res_list[[cond]])
      
      if (!identical(peaks_gr$peak_coord, rownames(df)))
        stop(sprintf("Row mismatch for %s: DESeq2 output not aligned with peaks", cond))
      
      perform_motif_enrichment(df, input_info, opt$seed)
    }), conds_with_hits)
    
    rm(input_info); rm(peaks_gr); gc()
  })
}


#10. GSEA
if (opt$skip_gsea) {
  message("##### Skipping GSEA (--skipGSEA) #####")
  gsea_list <- list()
  
} else {
  timed_step("Running GSEA", {
    gsea_list <- run_gsea(deseq_object$res_list_annot)
  })
}

##########################################################################

#8. Write outputs
timed_step("Writing outputs", {
  data.table::fwrite(as.data.frame(deseq_sweep$count_matrix),
                     paste0(out_dir, "/peak_counts.tsv.gz"),
                     sep = "\t", row.names = TRUE)
  
  data.table::fwrite(as.data.frame(consensus_peaks[["dataset_peaks_filt"]]),
                     paste0(out_dir, "/consensus_peaks.bed.gz"),
                     sep = "\t", row.names = FALSE)

  for (i in names(deseq_object$res_list_annot)){
    data.table::fwrite(deseq_object$res_list_annot[[i]], 
                       file = sprintf("%s/%s_vs_0h.tsv.gz", out_dir, i), 
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }

  for (cond in names(enrich_list)) {
    data.table::fwrite(enrich_list[[cond]]$enrich_df,
                       file = sprintf("%s/%s_vs_0h_DE_MotifEnrich.tsv.gz", out_dir, cond),
                       sep = "\t", quote = FALSE, row.names = TRUE)
  }
  
  for (i in names(gsea_list)) {
    for (j in names(gsea_list[[i]])) {
      data.table::fwrite(gsea_list[[i]][[j]],
                         file = sprintf("%s/%s_vs_0h_DE_GSEA_%s.tsv.gz", out_dir, j, i),
                         sep = "\t", quote = FALSE, row.names = TRUE)
    }
  }
  
  dracc_meta <- make_dracc_metadata(gene_name, deseq_object, has_hits, gsea_list)
  data.table::fwrite(dracc_meta, file = file.path(out_dir, "dracc_metadata.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  
  save(meta_df,
       deseq_sweep,
       consensus_peaks,
       enrich_list,
       gsea_list,
       file = paste0(out_dir, "/QC_data.RData"))
  
  output_gffs(gene_name, 
              consensus_peaks[["dataset_peaks_filt"]], 
              deseq_object$res_list, 
              out_dir, 
              padj_thresh = 0.05)
})
  
  
##########################################################################

rm(list = intersect(ls(), c("consensus_peaks", "deseq_sweep", "deseq_object",
                            "fc_object", "enrich_list", "gsea_list", "input_info")))
gc()

tryCatch({
  timed_step("Rendering QC report",
    rmarkdown::render(
      input       = file.path(script_dir, "QC_report.Rmd"),
      output_file = sprintf("%s/QC_report_%s.html", out_dir, gene_name),
      params      = list(
        gene       = gene_name,
        input_file = file.path(out_dir, "QC_data.RData"),
        out_dir    = out_dir)
    )
  )},
  error = function(e) warning(sprintf("QC report rendering failed: %s", e$message))
)

message(sprintf("##### Pipeline complete. Outputs: %s #####", out_dir))
pipeline_complete <- TRUE

}, error = function(e) {
  message(sprintf("ERROR: %s", e$message))
  if (isTRUE(opt$save_on_fail)) {
    save_dir  <- if (exists("out_dir", envir = .GlobalEnv) && dir.exists(out_dir)) out_dir else opt$out_dir
    save_path <- file.path(save_dir, "debug_dump.RData")
    message(sprintf("  --saveOnFail: writing debug image to %s", save_path))
    save.image(file = save_path)
  }
  if (!pipeline_complete) stop(e)
})

