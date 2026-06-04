# NW_degronATAC_analysis

ATAC-seq analysis pipeline for the MorPhiC degron screen (NW_ATAC-Seq02) to evaluate chromatin accessibility changes following auxin-inducible degradation of chromatin regulatory proteins at 0h, 6h, and 24h post-auxin treatment in KOLF2.1J. 

## Repository structure

```
NW_degronATAC_analysis/
├── 00_Peaks/
│   ├── 00_inputs/                 # Input BAMs (see file_structure.txt; not included for size, can be downloaded from Globus)
│   ├── 01_outputs/                # Per-gene/timepoint peak files (see file_structure.txt; not included for size, can be made available)
│   └── scripts/
│       ├── 00_run_pipeline.sh     # merges BAMs, call peaks, filters by exclusion list
│       ├── 01_submit_jobs.sh      # SLURM submission wrapper
│       └── env.yaml               # conda env
│
├── 01_DApeaks/
│   ├── 00_inputs/                                 # Per-gene metadata TSVs linking BAMs and peaks (see file_structure.txt; not included; can be made available)
│   ├── 01_outputs/
│   │   └── {GENE}/
│   │       ├── {tp}_vs_0h.tsv.gz                  # Differential accessibility results (DESeq2)
│   │       ├── {tp}_vs_0h_DE_MotifEnrich*.tsv.gz  # Motif enrichment in DA peaks
│   │       ├── {tp}_vs_0h_DE_GSEA_GO_BP.tsv.gz    # GO biological process enrichment
│   │       ├── {tp}_vs_0h_DE_GSEA_DO.tsv.gz       # Disease ontology enrichment
│   │       ├── QC_report_{$GENE}.html             # HTML output report for each gene. It is self-contained, but you will need to download locally for it to render. 
│   │       ├── dracc_metadata.tsv                 # Metadata for Data Portal Plotting
│   │       ├── consensus_peaks.bed.gz             # not included; Annotated dataset consensus peaks
│   │       └── peak_counts.tsv.gz                 # not included; Peak x Sample count matrix. 
│   └── scripts/
│       ├── 00_make_metadata.sh             # Builds per-gene metadata TSVs from input dirs
│       ├── 01_atac_pipeline.R              # Main DA pipeline (DESeq2, motif, GSEA)
│       ├── 01_atac_pipeline_functions.R    # Functions for pipeline above. 
│       ├── 02_submitJobs.sh                # SLURM submission wrapper
│       └── QC_report.Rmd                   # Per-gene QC report template
│
└── 02_e2g/
    ├── 00_inputs/                # Input config and qsub files (see file_structure.txt; not included)
    ├── 01_outputs/               # TBD -- will update soon.
    ├── ENCODE_rE2G/              # External GitHub repository (dev branch; not included)
    └── scripts/
        └── 00_make_metadata.sh   # Write config, qsub, and creates neccessary symlinks to peak files. 
```

## Pipeline overview

### Step 0 — Peak calling (`00_Peaks/`)

Run `01_submit_jobs.sh` to submit one SLURM job per gene/timepoint. Each job (`00_run_pipeline.sh`) does the following:

1. **Merge & sort BAMs** — replicate BAMs merged with `samtools merge` and position-sorted
2. **Convert to TagAlign** — BAM → BEDPE → tagAlign (paired-end aware, mates separated)
3. **Call peaks** — MACS3 with ATAC-seq shift/extension parameters (`--shift -75 --extsize 150 --call-summits`)
4. **Filter exclusion list** — peaks overlapping the ENCODE hg38 exclusion list are removed

Outputs per gene/timepoint: `ATAC-{GENE}-{tp}_peaks.filtered.narrowPeak.gz`, `_peaks.narrowPeak.gz`, `_summits.bed.gz`

### Step 1 — Differential accessibility (`01_DApeaks/`)

1. **Build metadata** — run `00_make_metadata.sh` to generate `{GENE}.metadata.tsv` files linking BAM paths and peak paths
2. **Run DA pipeline** — `02_submitJobs.sh` submits one SLURM job per gene running `01_atac_pipeline.R`, which:
   - Merges per-timepoint peaks into one consensus peak set
   - Counts reads over peaks with `Rsubread::featureCounts`
   - Filters to peaks with >1 CPM in at least 2 samples (replicate size)
   - Annotates peaks with `ChIPseeker` (gene, distance to TSS)
   - Generates a 500bp fixed peak coordinate (peak_coord500bp) in addition to the default peak coordinate (peak_coord)
   - Runs DESeq2 for each timepoint vs 0h
   - Runs motif enrichment on DA peaks with `monaLisa` (JASPAR2020 motifs)
   - Runs GSEA on DA-peak-associated genes with `clusterProfiler`
   - Renders a per-gene QC HTML report

#### Output file formats

| File pattern | Contents |
|---|---|
| `{tp}_vs_0h.tsv.gz` | DESeq2 results: `peak_coord`, `log2FoldChange`, `padj`, peak annotation columns |
| `*_DE_MotifEnrich*.tsv.gz` | Motif enrichment in gained/lost DA peaks (monaLisa scores) |
| `*_DE_GSEA_GO_BP.tsv.gz` | GO biological process enrichment (clusterProfiler) |
| `*_DE_GSEA_DO.tsv.gz` | Disease ontology enrichment (clusterProfiler) |
| `dracc_metadata.tsv` | Links metadata to input files. This gives a reason code for missing files and should be useful for plotting. |

### Step 2 — Enhancer-to-Gene Links (`02_e2g/`)

1. **Build Config Files** — run `00_make_qsubs.sh` to generate `00_inputs/qsubs` and `00_inputs/configs/` files for input to [ENCODE E2G platform](https://github.com/EngreitzLab/ENCODE_rE2G)

2. **Apply E2G** — details TBD


## Other Info:

1. There are [IGV tracks vizualizations](https://mitra.stanford.edu/engreitz/oak/Users/aeveritt/morphic_degron/IGV_tracks_public/) of the data
   - The top table presents quick links to the top DA peaks. 
   - The per-sample pileups are shown
   - The consensus peak set w/ variable width peaks colored by DA status. 
   - The E2G links per condition are shown as links. 

For questions/comments please contact Jesse Engreitz @ Stanford DAV

