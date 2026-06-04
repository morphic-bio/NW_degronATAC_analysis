#Author: Amanda Everitt
#in R

`%ni%` <- Negate(`%in%`)

################################################################################################

timed_step <- function(label, expr) {
  get_rss_kb <- function() as.numeric(gsub("[^0-9]", "", grep("VmRSS", readLines("/proc/self/status"), value = TRUE)))
  message(sprintf("\n##### %s #####", label))
  t0  <- proc.time()["elapsed"]
  rss_before <- get_rss_kb()
  result <- expr
  elapsed <- proc.time()["elapsed"] - t0
  rss_after <- get_rss_kb()
  delta_gb  <- (rss_after - rss_before) / 1024^2
  delta_str <- sprintf("%+.1f GB", delta_gb)
  message(sprintf("-> complete: %.1fs | RSS: %.1f GB (%s) #####\n", elapsed, rss_after / 1024^2, delta_str))
  invisible(result)
}


check_dependencies <- function(packages){
  missing <- packages[!packages %in% installed.packages()[, "Package"]]
  if (length(missing) > 0) {
    stop("Found missing R packages. Please install:\n", cat(missing))
  }
}


create_unique_dir <- function(dir_path) {
  if (!file.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Creating output directory")
    return(dir_path)
  }
  
  # Directory exists — find next available increment
  message("Directory exists. Creating an incremented directory")
  i = 1
  new_path <- paste0(dir_path, "_", i)
  while (file.exists(new_path)) {
    i <- i + 1
    new_path <- paste0(dir_path, "_", i)
  }
  dir.create(new_path, recursive = TRUE)
  return(new_path)
}


check_meta <- function(meta_df){
  required_cols <- c("gene", "sampleID", "timepoint", "bam_path", "merged_bam_path", "peak_path")
  missing_cols <- required_cols[required_cols %ni% colnames(meta_df)]
  if (length(missing_cols) > 0) {
    stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "),
         "\nRequired columns: ", paste(required_cols, collapse = ", "))
  }

  if (length(unique(meta_df$gene)) != 1) {
    stop("Metadata contains ", length(unique(meta_df$gene)), " unique genes: ",
         paste(unique(meta_df$gene), collapse = ", "),
         "\nThis pipeline expects metadata for a single gene only.")
  }

  zero_h_bams <- meta_df[meta_df$timepoint == "0h", ]$bam_path
  n_missing <- sum(is.na(zero_h_bams) | zero_h_bams == "")
  if (n_missing > 0) {
    stop("Found ", n_missing, " missing bam_path value(s) for timepoint '0h'.",
         "\nAll 0h samples require a valid bam_path for differential accessibility analysis.")
  }
  
}


################################################################################################

addExtraInfo <- function(meta_df, ncores = 1) {
  
  has_bam <- !is.na(meta_df$bam_path) & meta_df$bam_path != "NULL"
  
  stats_list <- mclapply(which(has_bam), function(i) {
    
    bam_file  <- meta_df$bam_path[i]
    peak_file <- meta_df$peak_path[i]
    
    bam_stats    <- idxstatsBam(bam_file)
    total_mapped <- sum(bam_stats$mapped)
    pct_mito     <- bam_stats$mapped[bam_stats$seqnames == "chrM"] / total_mapped * 100
    
    param <- ScanBamParam( what = "isize", flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isFirstMateRead = TRUE))
    frags          <- scanBam(bam_file, param = param)[[1]]$isize
    frags          <- frags[frags > 0]
    pct_nfr        <- mean(frags < 200) * 100
    pct_mono       <- mean(frags >= 200 & frags <= 400) * 100
    nfr_mono_ratio <- pct_nfr / pct_mono
    
    # frip_condition <- if (!is.na(peak_file) && peak_file != "NULL") {
    #   peaks <- rtracklayer::import(peak_file)
    #   sum(countBam(bam_file, param = ScanBamParam(which = reduce(peaks)))$records) / total_mapped
    # } else {
    #   NA_real_
    # }
    # 
    # frip_dataset <- sum(countBam(bam_file, param = ScanBamParam(which = consensus_gr))$records) / total_mapped
    
    data.frame(
      total_mapped   = total_mapped,
      pct_mito       = pct_mito,
      pct_nfr        = pct_nfr,
      pct_mono       = pct_mono,
      nfr_mono_ratio = nfr_mono_ratio
      #frip_condition = frip_condition,
      #frip_dataset   = frip_dataset
    )
  }, mc.cores = ncores)
  
  stats_list <- Filter(Negate(is.null), stats_list)
  stats_df <- do.call(rbind, stats_list)
  meta_df[which(has_bam), colnames(stats_df)] <- stats_df
  
  return(meta_df)
}



################################################################################################

add_summitInfo <- function(consensus_peaks){
  #Get all the summits
  get_summits <- lapply(consensus_peaks[["condition_peaks"]], function(gr) {
    GRanges(seqnames = seqnames(gr),
            ranges   = IRanges(start = start(gr) + gr$peak, width = 1),
            strand   = "*")
  })
  all_summits <- unlist(GRangesList(get_summits))

  #intersect the summits with the merged peaks.
  peaks <- GRanges(consensus_peaks[["dataset_peaks"]])
  hits  <- findOverlaps(peaks, all_summits, ignore.strand = TRUE)

  # factor over all peak indices so every peak gets an entry, then split in one pass
  f <- factor(queryHits(hits), levels = seq_along(peaks))
  summit_by_peak <- IntegerList(split(start(all_summits)[subjectHits(hits)], f))

  peaks$all_summits   <- summit_by_peak
  peaks$mean_summit   <- as.integer(mean(summit_by_peak))
  peaks$summit_spread <- as.integer(max(summit_by_peak) - min(summit_by_peak))

  consensus_peaks[["dataset_peaks"]] <- peaks
  return(consensus_peaks)
}

annotate_peaks <- function(gr){
  #for the purposes here, were going to hard code human specific resources.
  
  anno <- ChIPseeker::annotatePeak(
    gr,
    tssRegion = c(-3000, 3000),
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb = "org.Hs.eg.db"
  )
  
  
  #full_out_df <- dplyr::full_join(as.data.frame(anno), as.data.frame(gr), by=c("seqnames", "start","end", "width", "strand")) #if no formatting, much simpiler
  
  anno_df <- as.data.frame(anno) %>%
    mutate(gene_coord=paste0("chr", geneChr, ":", geneStart, "-", geneEnd)) %>%
    dplyr::select(-geneChr, -geneStart, -geneEnd)
  
  out_df <- dplyr::full_join(anno_df, as.data.frame(gr), by=c("seqnames", "start","end", "width", "strand"))
  
  full_out_df <- out_df %>%
    mutate(peak_coord=paste0(seqnames, ":", start, "-", end)) %>%
    dplyr::select("seqnames", "start", "end","width", "strand", "peak_coord", "gene_coord", "geneId", "annotation","ENSEMBL", "SYMBOL", "geneLength", "geneStrand","distanceToTSS") %>%
    mutate(annotation_class=sub("\\s*\\(.*", "", annotation))
  
  if (length(gr) != nrow(full_out_df)){
    message("When annotating peaks, there was a matrix mismatch. Please quality check this step.")
  }
  
  return(GRanges(full_out_df))
}


# add_fixed_width <- function(peaks, bam_paths, width = 500, summit_spread_thresh = 10, ncores = 1) {
# 
#   n_total <- length(peaks)
#   peaks$summit <- NA
# 
#   # For close summits, use the mean summit position
#   close_peaks <- which(peaks$summit_spread <= summit_spread_thresh)
#   peaks$summit[close_peaks] <- peaks$mean_summit[close_peaks]
#   message(sprintf("  %d/%d peaks: summit from mean (spread <= %d bp)", length(close_peaks), n_total, summit_spread_thresh))
# 
#   # For distant summits: single batched pileup across all peaks at once
#   needs_bam <- which(peaks$summit_spread > summit_spread_thresh)
#   if (length(needs_bam) > 0) {
#     message(sprintf("  %d/%d peaks: computing BAM-derived summits across %d BAMs...", length(needs_bam), n_total, length(bam_paths)))
# 
#     needs_bam_gr <- peaks[needs_bam]
#     peak_labels  <- paste0(seqnames(needs_bam_gr), ":", start(needs_bam_gr), "-", end(needs_bam_gr))
# 
#     sbp   <- Rsamtools::ScanBamParam(which = needs_bam_gr)
#     param <- Rsamtools::PileupParam(max_depth = 1e6, min_base_quality = 0,
#                                    min_mapq = 10, distinguish_strands = FALSE)
# 
#     pileup_all <- do.call(rbind, lapply(seq_along(bam_paths), function(i) {
#       message(sprintf("    [%d/%d] pileup: %s", i, length(bam_paths), basename(bam_paths[i])))
#       Rsamtools::pileup(bam_paths[i], scanBamParam = sbp, pileupParam = param)
#     }))
# 
#     # Sum counts across BAMs, then pick argmax position per peak region
#     bam_summits <- pileup_all %>%
#       group_by(which_label, pos) %>%
#       summarise(count = sum(count), .groups = "drop") %>%
#       group_by(which_label) %>%
#       summarise(summit = pos[which.max(count)], .groups = "drop")
# 
#     peaks$summit[needs_bam] <- bam_summits$summit[match(peak_labels, bam_summits$which_label)]
#     message("  BAM summits complete.")
#   }
# 
#   if (any(is.na(peaks$summit))) {
#     warning(sum(is.na(peaks$summit)), " peaks have no summit assigned — check summit_spread values.")
#   }
# 
#   # Fix the width around each summit
#   message(sprintf("  Resizing %d peaks to %d bp...", n_total, width))
#   half <- width %/% 2L
#   peak_start <- peaks$summit - half
#   peak_end   <- peaks$summit + half
# 
#   # Check chromosomal boundaries
#   if (any(peak_start < 1, na.rm = TRUE)) {
#     warning(sum(peak_start < 1, na.rm = TRUE), " peaks extend below chromosome start (pos < 1) and will be clipped.")
#     peak_start <- pmax(peak_start, 1L)
#   }
# 
#   chr_lens <- seqlengths(peaks)[as.character(seqnames(peaks))]
#   if (!all(is.na(chr_lens))) {
#     over_end <- !is.na(chr_lens) & peak_end > chr_lens
#     if (any(over_end)) {
#       warning(sum(over_end), " peaks extend beyond chromosome end and will be clipped.")
#       peak_end[over_end] <- chr_lens[over_end]
#     }
#   }
# 
#   peaks$peak_coord_500bp <- paste0(seqnames(peaks), ":", peak_start, "-", peak_end)
#   message(sprintf("  Done. %d peaks resized.", n_total))
#   return(peaks)
# }


pileup_bam_summits <- function(gr, bam_paths, ncores) {
  param <- Rsamtools::PileupParam(max_depth = 1e6, min_base_quality = 0, min_mapq = 10,
                                  distinguish_strands = FALSE, distinguish_nucleotides = FALSE)
  sbp   <- Rsamtools::ScanBamParam(which = gr)

  bam_results <- mclapply(seq_along(bam_paths), function(j) {
    message(sprintf("  [%d/%d] pileup: %s", j, length(bam_paths), basename(bam_paths[j])))
    Rsamtools::pileup(bam_paths[j], scanBamParam = sbp, pileupParam = param)
  }, mc.cores = min(length(bam_paths), ncores))

  pileup_combined <- do.call(rbind, bam_results)
  rm(bam_results); gc()

  summit_df <- pileup_combined %>%
    group_by(which_label, pos) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    group_by(which_label) %>%
    summarise(summit = pos[which.max(count)], .groups = "drop")

  setNames(summit_df$summit, summit_df$which_label)
}


resize_to_fixed_width <- function(peaks, width) {
  half       <- width %/% 2L
  peak_start <- peaks$summit - half
  peak_end   <- peaks$summit + half

  if (any(peak_start < 1, na.rm = TRUE)) {
    warning(sum(peak_start < 1, na.rm = TRUE), " peaks extend below chromosome start and will be clipped.")
    peak_start <- pmax(peak_start, 1L)
  }

  chr_lens <- seqlengths(peaks)[as.character(seqnames(peaks))]
  if (!all(is.na(chr_lens))) {
    over_end <- !is.na(chr_lens) & peak_end > chr_lens
    if (any(over_end)) {
      warning(sum(over_end), " peaks extend beyond chromosome end and will be clipped.")
      peak_end[over_end] <- chr_lens[over_end]
    }
  }

  peaks$peak_coord_500bp <- paste0(seqnames(peaks), ":", peak_start, "-", peak_end)
  message(sprintf("  Resized %d peaks to %d bp.", length(peaks), width))
  return(peaks)
}


add_fixed_width <- function(peaks, bam_paths, width = 500, summit_spread_thresh = 10, ncores = 1) {

  n_total      <- length(peaks)
  peaks$summit <- NA
  bam_paths <- bam_paths[!is.na(bam_paths) & bam_paths != "NULL" & nzchar(bam_paths)]

  close_peaks <- which(peaks$summit_spread <= summit_spread_thresh)
  peaks$summit[close_peaks] <- peaks$mean_summit[close_peaks]
  message(sprintf("  %d/%d peaks: summit from mean (spread <= %d bp)", length(close_peaks), n_total, summit_spread_thresh))

  needs_bam <- which(peaks$summit_spread > summit_spread_thresh)
  if (length(needs_bam) > 0) {
    needs_bam_gr <- peaks[needs_bam]
    peak_labels  <- paste0(seqnames(needs_bam_gr), ":", start(needs_bam_gr), "-", end(needs_bam_gr))
    message(sprintf("  %d/%d peaks: computing BAM-derived summits across %d BAMs...", length(needs_bam), n_total, length(bam_paths)))

    summits <- pileup_bam_summits(needs_bam_gr, bam_paths, ncores)
    peaks$summit[needs_bam] <- summits[peak_labels]
    message("  BAM summits complete.")
  }

  if (any(is.na(peaks$summit))) {
    warning(sum(is.na(peaks$summit)), " peaks have no summit assigned — check summit_spread values.")
  }

  peaks <- resize_to_fixed_width(peaks, width)
  return(peaks)
}


################################################################################################


run_feature_counts <- function(peaks, bams, n_threads=1){
  
  #reformat the GRanges object
  peaks_bed <- data.frame(
    Chr    = as.character(seqnames(peaks)),
    Start  = start(peaks),
    End    = end(peaks),
    Strand = as.character(strand(peaks))
  )
  peaks_bed$GeneID <- paste0(peaks_bed$Chr, ":", peaks_bed$Start, "-", peaks_bed$End)
  
  #Feature counts
  fc_result <- Rsubread::featureCounts(
    files = bams,
    annot.ext = peaks_bed,
    isPairedEnd = TRUE,
    countMultiMappingReads = FALSE,
    minMQS = 10,
    nthreads = n_threads
  )
  
  #counts <- fc_result$counts
  #stats <- fc_result$stat
  colnames(fc_result$counts) <- gsub(".bam", "", colnames(fc_result$counts))
  colnames(fc_result$stat) <- gsub(".bam", "", colnames(fc_result$stat))
  return(fc_result)
  #return(list("counts"=counts,"stats"=stats))
}



# filter_peaks <- function(counts, cpm_thresh=0.5, replicate_num=2){
#   
#   cpm <- apply(counts, 2, function(x) (x/sum(x))*1e6)
# 
#   
#   optimize_cpm_thresh <- function(cpm, thresholds = seq(0, 5, by = 0.1)) {
#     kept_counts <- sapply(thresholds, function(thresh) { sum(rowSums(cpm > thresh) >= replicate_num)})
#     
#     data.frame(cpm_thresh = thresholds, peaks_kept = kept_counts)
#   }
#   
#   results <- optimize_cpm_thresh(cpm)
#   
#   g1 <- ggplot(results, aes(x = cpm_thresh, y = peaks_kept)) +
#     geom_line(linewidth = 1) +
#     theme_classic() +
#     labs(
#       x = "CPM Threshold",
#       y = "Number of Peaks Retained",
#       title = "Peak Retention Across CPM Thresholds"
#     ) +
#     geom_hline(yintercept=150000) +
#     geom_vline(xintercept=cpm_thresh)
#   
# 
#   kept = rownames(cpm[rowSums(cpm > cpm_thresh) >= replicate_num,])
#   return(list("cpm_plot"=g1, "counts"=counts[rownames(counts) %in% kept, ]))
# }



################################################################################################


# identify_DE_peaks <- function(count_mat, sample_meta){
#   sample_meta$timepoint <- factor(sample_meta$timepoint)
#   
#   tp_allowed <- sample_meta %>%
#     dplyr::filter(timepoint != "0h") %>%
#     dplyr::count(timepoint) %>%
#     filter(n >= 2) %>% #require at least 2 replicates for any condition
#     pull(timepoint)
#   
#   dds <- DESeqDataSetFromMatrix(
#     countData = count_mat,
#     colData = sample_meta,
#     design = ~timepoint
#   )
#   
#   dds <- DESeq(dds, parallel = FALSE) #could change, but not a huge time sink.
#   ddsnorm <- counts(dds, normalized=TRUE)
#   
#   res_list = list()
#   for (cond in tp_allowed){
#     res <- results(dds, contrast = c("timepoint", cond, "0h"))
#     res_list[[cond]] = res
#   }
# 
#   return(list("res_list"=res_list, "ddsnorm"=ddsnorm))
# }


run_deseq_sweep <- function(counts, metadata, cpm_thresholds=c(), replicate_num = 2, default_thresh=1) {
  
  #a bit unneccessary, but this section requires 2 replicates for each condition are found.  
  viable_tp <- metadata %>% 
    filter(!is.null(sampleID) & sampleID != "NULL") %>% 
    dplyr::count(timepoint) %>%
    filter(n >= replicate_num) %>% 
    pull(timepoint)
  metadata <- metadata[metadata$timepoint %in% viable_tp, ]
  counts <- counts[, colnames(counts) %in% metadata$sampleID]
  conditions <- setdiff(viable_tp, "0h")

  # No non-0h timepoints with enough replicates — nothing to compare against.
  if ("0h" %ni% viable_tp) {
    message("No non-0h timepoints with sufficient replicates. Nothing to do — exiting.")
    quit(save = "no", status = 0)
  }

  cpm <- sweep(counts, 2, colSums(counts), "/") * 1e6
  metadata$timepoint <- relevel(factor(metadata$timepoint), ref = "0h")
  sweep_results <- vector("list", length(cpm_thresholds) * length(conditions))
  idx <- 1L
  
  for (thresh in unique(c(cpm_thresholds, default_thresh))) {
    
    #Filter peaks according to threshold. 
    keep <- rowSums(cpm > thresh) >= replicate_num
    counts_filt <- counts[keep, , drop = FALSE]
    n_peaks <- nrow(counts_filt)
    
    #Run DESEQ2
    dds <- DESeqDataSetFromMatrix(
      countData = counts_filt,
      colData = metadata,
      design = ~ timepoint
    )
    dds <- DESeq(dds, quiet = TRUE)
    
    deseq_list = list()
    for (cond in conditions) {
      res <- results(dds, contrast = c("timepoint", cond, "0h"))
      deseq_list[[cond]] = res
      
      sweep_results[[idx]] <- data.frame(
        cpm_threshold     = thresh,
        peaks_retained    = n_peaks,
        peaks_significant = sum(!is.na(res$padj) & res$padj < 0.05),
        timepoint         = cond
      )
      idx <- idx + 1L
    }
    
    #If this is the default threshold we want, return the objects. 
    if (thresh == default_thresh){
      #If these are the elements we want to return. 
      count_mat <- counts_filt
      deseq_object <- list("res_list"=deseq_list, "ddsnorm"=counts(dds, normalized=TRUE))
    }
  }
  
  sweep_df = do.call(rbind, sweep_results)
  
  return(list("sweep_results"=sweep_df, 
              "count_matrix"=count_mat,
              "deseq_object"=deseq_object))
}


clean_de_results <- function(deseq_df, peaks_annot) {
  cols_to_remove = c("seqnames", "start", "end", "width", "strand", "all_summits", "mean_summit", "summit_spread")
  
  df <- data.frame(deseq_df) %>%
      mutate(peak_coord = rownames(.)) %>%
      dplyr::left_join(as.data.frame(peaks_annot), by = "peak_coord") %>%
      dplyr::select(-all_of(cols_to_remove)) %>% 
      mutate(annotation_class = sub("\\s*\\(.*", "", annotation)) %>%
      dplyr::relocate("peak_coord") 
  
  return(df)
}

################################################################################################

anyhits <- function(df, n_required = 100, padj_thresh = 0.05, lfc_thresh = 0) {
  df <- as.data.frame(df)
  is_sig <- !is.na(df$padj) & df$padj < padj_thresh & abs(df$log2FoldChange) > lfc_thresh
  return(sum(is_sig, na.rm = TRUE) > n_required)
}


collect_input_info <- function(gr, ncores=1){
  
  #all_seqs <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, all_peaks)
  #fails because of require(SeqInfo) and R version isnt available on Sherlock. Go the long way. parallelization usually not needed either. consider removing in R-only script. 
  chunk_size = 10000
  chunks <- split(seq_along(gr), ceiling(seq_along(gr)/chunk_size))
  seq_list <- mclapply(chunks, function(idx) {
    sub_gr = gr[idx] 
    seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38,
                             names = as.character(seqnames(sub_gr)), start = start(sub_gr), end = end(sub_gr),
                             as.character = FALSE
    )
    names(seqs) <- sub_gr$peak_coord
    seqs
  }, mc.cores = ncores)
  
  failed <- !vapply(seq_list, is, logical(1), "DNAStringSet")
  if (any(failed)) warning(sprintf("%d/%d sequence chunks failed and will be dropped.", sum(failed), length(seq_list)))
  all_seqs <- Reduce(c, seq_list[!failed])
  
  pwm_list <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", 
                                                   tax_group = "vertebrates", 
                                                   all_versions = FALSE))
  
  cat("Number of sequences loaded:", length(all_seqs), "\n")
  cat("Number of motifs loaded:", length(pwm_list), "\n")
  return(list("seq"=all_seqs, "pwms"=pwm_list))
}


perform_motif_enrichment <- function(df, input_info, my_seed, padj_thresh=0.05, lfc_thresh=0,
                                      bg_multiplier=10, min_bg=20000){
  
  is_sig  <- !is.na(df$padj) & df$padj < padj_thresh & abs(df$log2FoldChange) > lfc_thresh
  fg_seqs <- input_info$seq[is_sig]
  #bg_seqs <- <- input_info$seq[!is_sig]
  
  #Memory issues pop up when the background set is too large. For now, use 10x the foreground set size OR 20K as the set to reduce the overhead. Can remove this later if needed.
  bg_all  <- input_info$seq[!is_sig]
  n_bg    <- min(max(length(fg_seqs) * bg_multiplier, min_bg), length(bg_all) - 1L)
  bg_seqs <- bg_all[sample(length(bg_all), max(n_bg, 1L))]

  cat("Foreground (DE) peaks:", length(fg_seqs), "\n")
  cat("Background peaks:     ", length(bg_seqs), "of", length(bg_all), "\n")
  
  set.seed(my_seed)
  enrich_res <- suppressWarnings(monaLisa::calcBinnedMotifEnrR(
    seqs     = c(fg_seqs, bg_seqs),
    bins     = factor(
      c(rep("foreground", length(fg_seqs)),
        rep("background", length(bg_seqs))),
      levels = c("background", "foreground")  # background is reference level
    ),
    pwmL     = toPWM(input_info$pwms)
  ))
  
  results_df <- cbind(
    data.frame(log2enr= assay(enrich_res, "log2enr")[, "foreground"],
               negLog10Padj = assay(enrich_res, "negLog10Padj")[, "foreground"]),
    rowData(enrich_res)[c("motif.name","motif.percentGC")]
  )
  
  results_df <- results_df[order(-results_df$negLog10Padj), ]
  
  return(list("monalisa_obj"=enrich_res, "enrich_df"=results_df))
}


output_pwms <- function(sel, outpath){
  if (!file.exists(outpath)) {
    dir.create(outpath, recursive = TRUE)
  }
  
  pwms <- rowData(sel)$motif.pwm
  names <- rowData(sel)$motif.name
  motif <- rowData(sel)$motif.id
  
  for (i in seq_along(pwms)) {
    png(sprintf("%s/%s_%s.png", outpath, names[i], motif[i]), width=800, height=400)
    mat <- exp(pwms[[i]]@profileMatrix)
    mat <- sweep(mat, 2, colSums(mat), "/")
    seqLogo::seqLogo(mat)  
    dev.off()
  }
}

################################################################################################
#Perform GSEA
run_gsea <- function(de_res_list){
  gsea_output = list(GO_BP = list(), DO = list())
  
  gene_rank_list = list()
  for (cond in names(de_res_list)){
    gene_rank_list[[cond]] <- de_res_list[[cond]] %>%
      as.data.frame() %>%
      filter(!is.na(padj), !is.na(geneId)) %>%
      mutate( weight = log2FoldChange * -log10(padj)) %>%
      group_by(geneId) %>%
      summarise(score = sum(weight, na.rm = TRUE)) %>%
      arrange(desc(score)) %>%
      { setNames(.$score, .$geneId) }
  }
  
  for (cond in names(gene_rank_list)){
    gene_rank = gene_rank_list[[cond]]
    
    ego <- gseGO(geneList     = gene_rank,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 verbose      = FALSE)
    ego_df <- as.data.frame(ego)
    if (nrow(ego_df) > 1){gsea_output[["GO_BP"]][[cond]] <- ego_df}
    
    ego <- gseDO(geneList     = gene_rank,
                 minGSSize    = 50,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 verbose      = FALSE)
    ego_df <- as.data.frame(ego)
    if (nrow(ego_df) > 1){gsea_output[["DO"]][[cond]] <- ego_df}
  }
  
  return(gsea_output)
}



output_gffs <- function(gene_name, consensus_gr, de_res_list, output_dir, padj_thresh = 0.05) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pretty_p <- function(p) {
    if (is.na(p) | p == ".") {return(NA)}
    if (p >= 1e-3) {
      sprintf("%.3f", p)
    } else {
      formatC(p, format = "e", digits = 2)
    }
  }

  get_de <- function(tp) {
    if (!tp %in% names(de_res_list)) {
      tmp <- allpeaks[, "peak_coord", drop=FALSE]
      tmp$padj <- NA
      tmp$log2FoldChange <- NA
    } else {
      tmp <- tibble::rownames_to_column(as.data.frame(de_res_list[[tp]]), "peak_coord")
      tmp$padj <- sapply(tmp$padj, pretty_p)
    }
    return(tmp[c("peak_coord", "padj", "log2FoldChange")])
  }

  gene     <- gene_name
  allpeaks <- as.data.frame(consensus_gr)

  de6  <- get_de("6h")
  de24 <- get_de("24h")

  df = full_join(de6, de24, by=c("peak_coord"), suffix=c("_6hr", "_24hr"))
  df = full_join(df, allpeaks, by=c("peak_coord"))

  min_idx <- apply(cbind(
    ifelse((!is.na(df$padj_6hr)  & df$padj_6hr  < padj_thresh), df$padj_6hr,  NA),
    ifelse((!is.na(df$padj_24hr) & df$padj_24hr < padj_thresh), df$padj_24hr, NA)
  ), 1, function(x) {
    if (all(is.na(x))) NA_integer_ else which.min(x)
  })
  min_idx <- as.data.frame(min_idx)
  min_idx$best_fc <- df[c("log2FoldChange_6hr", "log2FoldChange_24hr")][cbind(seq_len(nrow(df)), min_idx$min_idx)]
  df$DE_status <- ifelse(
    is.na(min_idx$best_fc), "not DE",
    ifelse(min_idx$best_fc > 0, "up", "down")
  )

  attr_cols <- c("SYMBOL", "annotation_class", "DE_status", "log2FoldChange_6hr", "padj_6hr", "log2FoldChange_24hr", "padj_24hr")
  attr_names <- c("GeneSymbol", "PeakAnnot", "DE_status", "log2FC_6hr", "padj_6hr", "log2FC_24hr", "padj_24hr")
  df$attributes <- apply(df[, attr_cols], 1, function(x) paste0(attr_names, "=", x, collapse = ";"))

  gff <- data.frame(
    seqname    = df$seqnames,
    source     = ".",
    feature    = "peak",
    start      = df$start + 1L,   # BED 0-based → GFF3 1-based
    end        = df$end,
    score      = ".",
    strand     = ".",
    frame      = ".",
    attributes = df$attributes,
    stringsAsFactors = FALSE
  )
  gff <- gff[order(gff$seqname, gff$start), ]

  out_bed <- file.path(output_dir, paste0(gene, ".gff3"))
  write.table(gff, out_bed, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  return(out_bed)
}



################################################################################################


################################################################################################
# Write result metadata sheet

make_dracc_metadata <- function(gene_name, deseq_object, has_hits, gsea_list) {
  
  result_spec <- data.frame(
    result_type = c("da_peaks", "da_peak_enriched_motifs", "da_peak_gene_set_enrichment_GO.BP", "da_peak_gene_set_enrichment_DO"),
    plot_type   = c("volcano",  "heatmap",   "barchart",  "barchart"),
    gsea_cat    = c(NA,          NA,          "GO_BP",     "DO"),
    stringsAsFactors = FALSE
  )
  
  base <- data.frame(gene = gene_name, 
                     assay = "ATAC", 
                     BAM_source = "internal_release::MorPhiC_Release_September_2025/NW_ATAC-Seq02/processed/", 
                     PEAK_source = "internal_release::MorPhiC_Release_May_2026/NW_ATACseq02_DRACC_processed_v3/processed",
                     stringsAsFactors = FALSE)
  
  rows <- lapply(c("6h", "24h"), function(tp) {
    tp_ran      <- tp %in% names(deseq_object$res_list)
    tp_has_hits <- tp_ran && isTRUE(has_hits[tp])
    
    lapply(seq_len(nrow(result_spec)), function(i) {
      rt  <- result_spec$result_type[i]
      pt  <- result_spec$plot_type[i]
      cat <- result_spec$gsea_cat[i]
      
      if (!tp_ran) {
        fn <- NA_character_;  st <- "not_run (missing timepoint)"
      } else if (rt == "da_peaks") {
        fn <- sprintf("%s_vs_0h.tsv.gz", tp);  st <- "generated"
      } else if (rt == "da_peak_enriched_motifs") {
        fn <- ifelse(tp_has_hits, sprintf("%s_vs_0h_DE_MotifEnrich.tsv.gz", tp), NA_character_)
        st <- ifelse(tp_has_hits, "generated", "not_generated (insufficient DE peaks)")
      } else {
        was <- tp %in% names(gsea_list[[cat]])
        fn  <- ifelse(was, sprintf("%s_vs_0h_DE_GSEA_%s.tsv.gz", tp, cat), NA_character_)
        st  <- ifelse(was, "generated", "not_generated (no significant enrichment)")
      }
      
      cbind(base, data.frame(timepoint = tp, 
                             comparison = paste0(tp, "_vs_0h"),
                             result_type = rt, 
                             plot_type = pt,
                             file_name = fn, 
                             status = st, stringsAsFactors = FALSE))
    })
  })
  
  return(do.call(rbind, unlist(rows, recursive = FALSE)))
}



