# AIV Dataset Processing

Scripts for cleaning avian influenza virus metadata, preparing alignments, subsampling, computing distances, and summarizing phylogenetic topology agreement. Bash helpers wrap MAFFT/FastTree, and R scripts handle preprocessing, grouping, and visualization.

## Repository Layout
- `coding_upload/` raw drop of sequence/metadata files (by year range).
- `function.R` shared helpers (I/O, header parsing, plotting utilities, package imports).
- `meta_clean_v2.R` merge/clean of GISAID + IRD metadata and sequences.
- `align.sh` MAFFT alignment runner for all `.fa` in a working directory.
- `run_fast_tree.sh` FastTree GTR nucleotide trees for `_sampling.fa` inputs.
- `topology_assignment.R` + `toplogy_assignment.sh` IQ-TREE topology agreement summaries across distance thresholds.
- `patristic_distance_analysis.R`, `distance.sh`, `inter_intra_GTR_gd.R` distance matrices and intra/inter-segment summaries.
- `sampling_as.sh`, `sampling_aa_Identicle_reflection.R`, `PB2_NS_sampling_result.R` amino-acid and segment subsampling utilities.
- `plots.R`, `agreement.sh`, `save_data.R` plot creation and result export.

## Data Expectations
- Metadata: directory of CSVs from GISAID/IRD (`-m`), columns include `Isolate_Id`, `Subtype`, `Clade`, `Location`, `Host`, `Collection_Date`, `Submission_Date`, and INSDC accessions.
- Sequences: combined FASTA from GISAID (`-gis`) and IRD (`-ird`). Headers must contain `EPI_ISL` IDs and segment info separated by pipes (see parsing in `meta_clean_v2.R` and helpers in `function.R`).
- Segment naming: `seg_sub_level` in `function.R` defines expected segments (`PB2`, `PB1`, `PA`, `H1-16`, `NP`, `N1-9`, `MP`, `NS`). Update if adding new targets.

## Environment Setup
- R (tested with 4.3). Install required packages:
```
Rscript -e "install.packages(c('argparse','dplyr','purrr','fs','seqinr','stringr','data.table','stringdist','ggplot2','ComplexHeatmap','circlize','ape','RColorBrewer','ggtree','tidytree','ggnewscale','TDbook','aplot','phytools','ggthemes','ggplotify','smplot2','scales','ggpubr','gridExtra','treemap','castor'))"
```
- Command-line tools: MAFFT, FastTree, IQ-TREE available on PATH.
- Paths in scripts are currently user-specific (`/home/eric/...`). Adjust `sdir`, `rdir`, and argument defaults to match your filesystem before running.

## Workflow Overview
1) **Clean metadata and sequences** (produces filtered FASTA + cleaned metadata):
```
Rscript meta_clean_v2.R \
  -m /path/to/meta/ \
  -gis /path/to/gisaid_all.fasta \
  -ird /path/to/IRD_Sequence.fa \
  -bc 10 \
  -prop 5 \
  -ns 4 \
  -p gisaid_IRD_merged_ \
  -o /path/to/output_dir/
```
2) **Alignment** (align all FASTA in a directory; edit `sdir` in `align.sh`):
```
bash align.sh
```
3) **Tree building** (FastTree on sampled FASTA; edit `sdir` in `run_fast_tree.sh`):
```
bash run_fast_tree.sh
```
4) **Topology agreement** (compare IQ-TREE groupings vs distance thresholds):
```
Rscript topology_assignment.R \
  -seg PA \
  -iqt /path/to/PA_5p_sampling.nexus \
  -iqg /path/to/PA_iq_MPD_groups.RData \
  -op _as_sampling \
  -o /path/to/output_dir/
```
5) **Distances and summaries**:
- `distance.sh` wraps distance calculations per segment.
- `inter_intra_GTR_gd.R` and `patristic_distance_analysis.R` compute intra/inter group distances from trees.
- `sampling_as.sh` and related R scripts perform amino-acid/segment subsampling before distance or topology steps.

## Script Cheat Sheet (inputs → outputs)
- `meta_clean_v2.R` merge/clean
  - Inputs: `-m` metadata dir, `-gis` GISAID FASTA, `-ird` IRD FASTA, `-bc` box cut, `-prop` max non-NT %, `-ns` min segment count, `-p` prefix, `-o` output dir.
  - Outputs: cleaned FASTA(s) and metadata with standardized headers.
- `align.sh` alignment
  - Inputs: `.fa` files in `sdir` (set inside script).
  - Outputs: `_aligned.fa` per input using `mafft --auto --anysymbol --thread -8`.
- `run_fast_tree.sh` FastTree
  - Inputs: `*_sampling.fa` files in `sdir`.
  - Outputs: `_sampling_ft.tree` (GTR, nucleotide).
- `distance.sh` distance + grouping + agreement
  - Inputs: IQ/FT trees and patristic matrices (paths hard-coded).
  - Steps: (commented) patristic distance generation; active parts partition distance matrices (`partistic_partition*.R`) and compute IQ vs FT agreement (`agreement_grouping_0924.R`).
  - Outputs: grouped RData files and agreement tables in `~/Analysis/aiv/merge/0307/distance/ns_sampling/`.
- `inter_intra_GTR_gd.R`, `patristic_distance_analysis.R`
  - Inputs: tree directory, file prefix, output suffix/dir.
  - Outputs: patristic distance matrices and intra/inter group summaries for each segment.
- `topology_assignment.R`
  - Inputs: `-seg`, `-iqt` (nexus), `-iqg` (grouping RData), `-op` suffix, `-o` out dir.
  - Outputs: agreement tables/plots across distance thresholds for a segment.
- `sampling_as.sh`
  - Inputs: segment list (`elements`), sampled IQ/FT results.
  - Outputs: calls `PB2_NS_sampling_result.R` to write agreement RData/plots for each segment.
- `PB2_NS_sampling_result.R`
  - Inputs: `--segment`, `--iqtree`, `-iqg`, `-ftg`, `-op`, `-o`.
  - Outputs: agreement metrics between IQ and FT groupings (RData + figures).
- `sampling_aa_Identicle_reflection.R`
  - Inputs: amino-acid alignment and sampling parameters inside the script.
  - Outputs: balanced/reflected amino-acid sampling sets for downstream analyses.
- `plots.R`, `save_data.R`, `agreement.sh`, `iq_sampling.sh`
  - Plot/export helpers and wrappers; edit paths/segment lists before running.

## Key Inputs/Outputs
- Inputs: raw metadata CSVs, combined FASTA files, optional subsampling lists.
- Intermediate: `_aligned.fa`, `_sampling.fa`, `*_ft.tree`, distance matrices (`*.RData`), IQ-TREE outputs (`*.nexus`).
- Outputs: cleaned FASTA and metadata, agreement tables, plots, and exported RData objects (`save_data.R`).

## Customization Tips
- Update hard-coded paths (`sdir`, `rdir`, file prefixes) near the top of each script.
- If header formats differ, adjust parsing in `function.R` (`read_as_list`, header cleaning helpers, `organ_subtype`).
- Ensure output directories exist before running scripts; most do not create parent paths.

## Troubleshooting
- Missing package errors: rerun the install command above.
- Empty or NA segment assignments: confirm header contains segment tokens and matches `seg_sub_level`.
- MAFFT/FastTree not found: verify PATH or install via your package manager.
- RData shape changes: some scripts assume specific column names (`Strain_number`, `Isolate_Id`, `segment`); keep these consistent across inputs.

## End-to-End Example (edit paths)
```
# 1) Clean metadata + sequences
Rscript meta_clean_v2.R \
  -m ~/Analysis/aiv/gisaid_all/meta/ \
  -gis ~/Analysis/aiv/gisaid_all/all.fasta \
  -ird ~/Analysis/aiv/ird/IRD_Sequence.fa \
  -bc 10 -prop 5 -ns 4 \
  -p gisaid_IRD_merged_ \
  -o ~/Analysis/aiv/merge/0307/

# 2) Align the cleaned FASTA
cd ~/Analysis/aiv/merge/0307/
bash /path/to/repo/align.sh        # set sdir in script to this folder

# 3) Build FastTree on sampled FASTA
bash /path/to/repo/run_fast_tree.sh  # set sdir to sampling directory

# 4) IQ-TREE runs (if needed) and distance partitions
bash /path/to/repo/iq_sampling.sh    # edit IQ settings/segments
bash /path/to/repo/distance.sh       # edit tree/patristic paths

# 5) Agreement summaries and plots
Rscript topology_assignment.R \
  -seg PA \
  -iqt ~/Analysis/aiv/merge/0307/sampling_subset/0924/iq/PA_5p_sampling.nexus \
  -iqg ~/Analysis/aiv/merge/0307/distance/sampling_0924/PA_iq_MPD_groups.RData \
  -op _as_sampling \
  -o ~/Analysis/aiv/merge/0307/distance/sampling_0924/
```

## File Naming Conventions
- Cleaned FASTA prefix: `gisaid_IRD_merged_` (set by `-p` in `meta_clean_v2.R`).
- Alignments: `<prefix>_aligned.fa` from `align.sh`.
- Sampled FASTA: `<seg>_5p_sampling.fa` (example; created upstream).
- Trees: `<seg>_5p_sampling.nexus` (IQ-TREE), `<seg>_sampling_ft.tree` (FastTree).
- Patristic matrices: `<seg>_*_patristic_dist_matrix.RData`.
- Groupings: `<seg>_*_patristic_group_bympd_percentile.RData`.
- Agreement outputs: `<seg>_*_agreement_result.RData` plus plots.

## Directory Layout Tips
- Keep raw inputs separate from derived outputs:
  - `raw/` for original metadata/FASTA.
  - `work/` for cleaned FASTA and alignments.
  - `trees/` for IQ-TREE and FastTree outputs.
  - `distance/` for patristic matrices and groupings.
  - `figures/` or `reports/` for plots and tables.
- Mirror these paths in the scripts (`sdir`, `rdir`, output flags) to avoid path confusion.

## Data Flow (short version)
- Metadata + raw FASTA → `meta_clean_v2.R` → cleaned FASTA + metadata tables.
- Cleaned FASTA → `align.sh` → `_aligned.fa`.
- Aligned/sampled FASTA → IQ-TREE / FastTree (`iq_sampling.sh`, `run_fast_tree.sh`) → tree files.
- Trees → `patristic_distance_analysis.R` / `inter_intra_GTR_gd.R` → patristic matrices + intra/inter summaries.
- Matrices → `partistic_partition*.R` (via `distance.sh`) → grouping RData by MPD percentiles.
- Groupings + trees → `PB2_NS_sampling_result.R` / `topology_assignment.R` → agreement metrics + plots.

## Helpful Commands
- Check FASTA header patterns:
```
grep -m 5 ">" /path/to/file.fasta
```
- Count sequences per segment (after cleaning):
```
grep ">" cleaned.fa | cut -d"|" -f4 | sort | uniq -c
```
- Validate MAFFT/FastTree availability:
```
mafft --version
fasttree -help | head -n 1
```

## Reproducibility Checklist
- Record tool versions: `mafft --version`, `fasttree -help | head -n 1`, `iqtree2 --version`, `R --version`.
- Snapshot R package versions with `sessionInfo()` when running key scripts.
- Keep the exact command lines (with arguments) used for each run; store logs alongside outputs.
- Avoid modifying raw inputs; write all derived data to new directories with timestamps (e.g., `merge/0307/`, `distance/0924/`).

## Assumptions
- FASTA headers contain `EPI_ISL_` IDs and segment tokens separated by `|`.
- Metadata CSVs share column names expected in `meta_clean_v2.R`.
- Segments conform to `seg_sub_level`; expand if new segments/hosts are added.

## FAQ
- **How do I switch to a new working directory?** Edit `sdir` (and `rdir` if needed) near the top of each `.sh` script to point to your run directory; keep paths consistent across scripts.
- **Can I add more segments?** Update `seg_sub_level` in `function.R`, ensure metadata/headers contain the new segment tokens, and include them in the `elements` arrays inside the shell scripts.
- **What if metadata has different column names?** Adjust column selections in `meta_clean_v2.R` where it selects `Isolate_Id`, `Subtype`, `Clade`, `Location`, `Host`, `Collection_Date`, `Submission_Date`, and INSDC columns.
- **Distance thresholds seem off.** Inspect the MPD percentile thresholds in `partistic_partition*.R`; retune percentiles or add additional cut points, then re-run `distance.sh`.
- **Plots look empty or sparse.** Check that grouping RData files are non-empty; confirm segments had sufficient sequences after filtering; verify header parsing didn’t drop segment labels.
- **I need to rerun only agreement plots.** Reuse existing grouping RData and trees; rerun `PB2_NS_sampling_result.R` or `topology_assignment.R` with updated `-op` suffix to avoid overwriting prior outputs.

## Where to Edit First
- Paths: `sdir`, `rdir`, and argument defaults in `*.sh` and top of R scripts.
- Segment lists: `elements=(...)` arrays inside shell wrappers.
- Prefixes/suffixes: file naming in `meta_clean_v2.R`, `align.sh`, and any distance/plotting scripts to keep outputs organized by run date.
