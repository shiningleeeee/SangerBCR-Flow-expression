# =============================================================================
# 01_pair_and_complete.R
#
# Step 1 of IgG BCR expression-ready selection.
#
# Input : one IgBLAST AIRR TSV per donor in data/, named after the donor
#         (data/sample01.tsv, data/sample02.tsv, ...).
# Output: output/*.xlsx
#
# Self-contained: settings, input discovery, reference acquisition and the
# analysis all live in this one file, in that order, so it can be read top to
# bottom. R/02_extract_simplified.R is a separate, thin script that turns this
# script's workbook into the tables a human reads.
#
# What it does, per donor:
#   * keeps productive, stop-free, in-frame rearrangements
#   * keeps heavy chains whose subclass is explicitly called as IGHG1-4
#   * pairs heavy and light chains by cell id
#   * completes the N terminus from the V germline and the C terminus from the
#     J germline, recording exactly how many residues were added
#   * runs whole-variable-domain QC (frame, translation, CIGAR, framework
#     profile, V-prefix identity) and grades every pair A/B/C/D
#
# Only A_priority and B_acceptable pairs are expression-ready. Every C/D pair
# keeps its evidence in manual_check / whole_variable_risk, and its
# final_VH_aa / final_VL_aa are set to NA so that an incomplete sequence can
# never be read from a `final_` column by mistake.
#
# The completion and grading rules are unchanged from the validated prototype.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(writexl)
})


# =============================================================================
# SETTINGS -- everything you might want to change lives in this block
# =============================================================================

# -----------------------------------------------------------------------------
# project layout
# -----------------------------------------------------------------------------

# Resolve the project root without depending on the caller's working directory.
# Rscript does not expose its own path, so we walk the call stack for the
# --file= argument; when that fails we fall back to the working directory.
bcr_script_path <- function() {
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)$ofile
    if (!is.null(f) && nzchar(f)) return(normalizePath(f, mustWork = FALSE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) {
    return(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE))
  }
  normalizePath(file.path(getwd(), "01_pair_and_complete.R"), mustWork = FALSE)
}

# This script sits in <root>/R/, so the project root is one level up.
bcr_root_dir <- dirname(dirname(bcr_script_path()))

CFG <- list(
  root_dir     = bcr_root_dir,
  data_dir     = file.path(bcr_root_dir, "data"),
  ref_dir      = file.path(bcr_root_dir, "reference"),
  output_dir   = file.path(bcr_root_dir, "output"),
  scripts_dir  = file.path(bcr_root_dir, "R")
)

# Rscript on Windows cannot read its script from %TEMP% on this machine
# (observed access violation), so warn early with a readable message instead of
# letting the process die without output.
if (grepl("[/\\\\](Temp|tmp)[/\\\\]", CFG$root_dir)) {
  warning(
    "The project appears to be running from a temporary directory (", CFG$root_dir, ").\n",
    "  On this machine Rscript may crash (0xC0000005) when reading its script from %TEMP%.\n",
    "  Run the scripts from inside the project folder instead.",
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# config.ini reader
#
# Deliberately dependency-free: adding ini to Imports for a handful of numbers
# is not worth the install burden. Flat [section] / key = value files only.
# -----------------------------------------------------------------------------

read_ini_flat <- function(path) {
  if (!file.exists(path)) return(list())
  lines <- readLines(path, warn = FALSE)
  lines <- sub("[#;].*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  out <- list()
  for (ln in lines) {
    if (grepl("^\\[.*\\]$", ln)) next          # ignore section headers
    pos <- regexpr("=", ln, fixed = TRUE)
    if (pos < 1) next
    key <- trimws(substr(ln, 1, pos - 1))
    val <- trimws(substr(ln, pos + 1, nchar(ln)))
    val <- sub('^"(.*)"$', "\\1", val)
    val <- sub("^'(.*)'$", "\\1", val)
    if (nzchar(key)) out[[key]] <- val
  }
  out
}

ini_values <- read_ini_flat(file.path(CFG$root_dir, "config.ini"))


# -----------------------------------------------------------------------------
# typed getters: env var > config.ini > default
# -----------------------------------------------------------------------------

cfg_lgl <- function(key, default) {
  env <- Sys.getenv(key, unset = NA_character_)
  raw <- if (!is.na(env)) env else ini_values[[key]]
  if (is.null(raw) || is.na(raw) || !nzchar(raw)) return(default)
  toupper(trimws(raw)) %in% c("TRUE", "T", "1", "YES", "Y")
}

cfg_chr <- function(key, default) {
  env <- Sys.getenv(key, unset = NA_character_)
  raw <- if (!is.na(env)) env else ini_values[[key]]
  if (is.null(raw) || is.na(raw) || !nzchar(raw)) return(default)
  trimws(raw)
}

cfg_num <- function(key, default, min = -Inf, max = Inf) {
  env <- Sys.getenv(key, unset = NA_character_)
  raw <- if (!is.na(env)) env else ini_values[[key]]
  if (is.null(raw) || is.na(raw) || !nzchar(raw)) return(default)
  v <- suppressWarnings(as.numeric(raw))
  if (is.na(v)) {
    warning("config: '", key, "' = '", raw, "' is not numeric, using default ",
            default, call. = FALSE)
    return(default)
  }
  if (v < min || v > max) {
    stop("config: '", key, "' = ", v, " is outside the allowed range [",
         min, ", ", max, "].", call. = FALSE)
  }
  v
}

cfg_int <- function(key, default, min = 0L, max = .Machine$integer.max) {
  as.integer(round(cfg_num(key, default, min = min, max = max)))
}


# -----------------------------------------------------------------------------
# run control
# -----------------------------------------------------------------------------

CFG$strict_igg_subclass_only <- cfg_lgl("BCR_STRICT_IGG_SUBCLASS_ONLY", TRUE)

# Minimum IGH records a file must contain to be accepted as a donor sample.
# Below this the file is reported rather than silently used: it usually means
# the wrong tsv was dropped into data/.
CFG$min_igh_rows <- cfg_int("BCR_MIN_IGH_ROWS", 20L, min = 1L)

# How to recover the cell/well identity from a sequence_id.
#
# A Sanger record id identifies the well and the chain, for example
# "A11H1.<run>.<...>" for the heavy chain of well A11 or "A11K.<run>.<...>"
# for its kappa chain. Heavy and light reads are paired on the well name, so the
# chain marker has to be stripped from the first dot-separated field.
#
# default_well_id() takes everything before the first "." and then removes a
# trailing chain marker: an optional separator followed by H with optional
# digits, or K, or L. A well genuinely called "H10" is left intact, because the
# H there is followed by more of the name rather than being a suffix.
#
# If your ids use another convention -- a different separator, "VH"/"VL"
# instead of H/K/L, the chain marker at the end rather than the beginning --
# override BCR_CHAIN_SUFFIX_PATTERN rather than editing the scripts.
CFG$chain_suffix_pattern <- cfg_chr("BCR_CHAIN_SUFFIX_PATTERN", "[-_ ]?(H\\d*|K|L)$")

# Separator between the well name and the chain marker when the two are split
# apart, used only to keep the pattern readable in error messages.
CFG$chain_suffix_examples <- cfg_chr("BCR_CHAIN_SUFFIX_EXAMPLES", "H, H1, K, L")

# Terminal completion limits.
CFG$max_added_N_clean <- cfg_int("BCR_MAX_ADDED_N_CLEAN", 10L, min = 0L)
CFG$max_added_N_hard  <- cfg_int("BCR_MAX_ADDED_N_HARD",  15L, min = 0L)
CFG$max_added_C_clean <- cfg_int("BCR_MAX_ADDED_C_CLEAN",  6L, min = 0L)
CFG$max_added_C_hard  <- cfg_int("BCR_MAX_ADDED_C_HARD",   8L, min = 0L)

# Alignment consistency thresholds.
CFG$min_anchor_clean    <- cfg_num("BCR_MIN_ANCHOR_CLEAN",    0.75, min = 0, max = 1)
CFG$min_anchor_hard     <- cfg_num("BCR_MIN_ANCHOR_HARD",     0.65, min = 0, max = 1)
CFG$min_j_overlap_clean <- cfg_num("BCR_MIN_J_OVERLAP_CLEAN", 0.75, min = 0, max = 1)
CFG$min_j_overlap_hard  <- cfg_num("BCR_MIN_J_OVERLAP_HARD",  0.65, min = 0, max = 1)

# Plausible variable-domain lengths (amino acids).
CFG$VH_len_min <- cfg_int("BCR_VH_LEN_MIN", 105L, min = 1L)
CFG$VH_len_max <- cfg_int("BCR_VH_LEN_MAX", 140L, min = 1L)
CFG$VL_len_min <- cfg_int("BCR_VL_LEN_MIN",  95L, min = 1L)
CFG$VL_len_max <- cfg_int("BCR_VL_LEN_MAX", 125L, min = 1L)
if (CFG$VH_len_min >= CFG$VH_len_max || CFG$VL_len_min >= CFG$VL_len_max) {
  stop("config: VH/VL length ranges must have min < max.", call. = FALSE)
}

# Whole-variable-domain frame / structure QC.
CFG$whole_v_prefix_len <- cfg_int("BCR_WHOLE_V_PREFIX_LEN", 80L, min = 40L)
CFG$min_whole_v_identity_hard <-
  cfg_num("BCR_MIN_WHOLE_V_IDENTITY_HARD", 0.60, min = 0, max = 1)
CFG$max_whole_v_mismatch_run_hard <-
  cfg_int("BCR_MAX_WHOLE_V_MISMATCH_RUN_HARD", 10L, min = 1L)

# Framework regions use deliberately loose ranges: they exist to spot gross
# anomalies, not to reproduce IMGT numbering.
CFG$FWR1_len_range <- c(cfg_int("BCR_FWR1_LEN_MIN", 15L, 1L),
                        cfg_int("BCR_FWR1_LEN_MAX", 35L, 1L))
CFG$FWR2_len_range <- c(cfg_int("BCR_FWR2_LEN_MIN", 10L, 1L),
                        cfg_int("BCR_FWR2_LEN_MAX", 25L, 1L))
CFG$FWR3_len_range <- c(cfg_int("BCR_FWR3_LEN_MIN", 25L, 1L),
                        cfg_int("BCR_FWR3_LEN_MAX", 55L, 1L))


# -----------------------------------------------------------------------------
# germline reference acquisition
# -----------------------------------------------------------------------------

# Only human is wired up. Another species needs its own germline FASTAs and a
# matching locus prefix table in R/01_pair_and_complete.R; the parameter exists
# so that call sites do not have to change later.
CFG$species <- cfg_chr("BCR_SPECIES", "human")
if (!identical(CFG$species, "human")) {
  stop(
    "config: unsupported species '", CFG$species, "'.\n",
    "  Only 'human' is implemented. For another species, place the IMGT\n",
    "  germline FASTAs in ", CFG$ref_dir, " manually.",
    call. = FALSE
  )
}

CFG$download_reference <- cfg_lgl("BCR_DOWNLOAD_REFERENCE", TRUE)

CFG$imgt_base_url <- cfg_chr(
  "BCR_IMGT_BASE_URL",
  paste0("https://www.imgt.org/download/V-QUEST/",
         "IMGT_V-QUEST_reference_directory/Homo_sapiens/IG")
)

# Per-attempt download timeout in seconds.
CFG$download_timeout <- cfg_int("BCR_DOWNLOAD_TIMEOUT", 120L, min = 5L)

# Gene families the pipeline reads. IGHD is fetched for completeness even
# though terminal completion only uses V and J.
CFG$reference_genes <- c("IGHV", "IGKV", "IGLV", "IGHJ", "IGKJ", "IGLJ", "IGHD")

# Sanity floors, measured against the IMGT human release shipped in reference/
# (IGHV 460, IGKV 132, IGLV 119, IGHJ 15, IGKJ 10, IGLJ 11, IGHD 47). They sit
# well below those counts so a legitimate IMGT update is not rejected, while a
# truncated or HTML-error-page download still is.
CFG$ref_min_records <- c(
  IGHV = cfg_int("BCR_REF_MIN_IGHV", 300L),
  IGKV = cfg_int("BCR_REF_MIN_IGKV",  90L),
  IGLV = cfg_int("BCR_REF_MIN_IGLV",  80L),
  IGHJ = cfg_int("BCR_REF_MIN_IGHJ",  10L),
  IGKJ = cfg_int("BCR_REF_MIN_IGKJ",   6L),
  IGLJ = cfg_int("BCR_REF_MIN_IGLJ",   6L),
  IGHD = cfg_int("BCR_REF_MIN_IGHD",  30L)
)


# -----------------------------------------------------------------------------
# priority labels shared by both scripts
#
# Kept here so that renaming a priority cannot silently desynchronise step 2
# from step 1 (step 2 would otherwise classify every row as "CHECK_manual").
# -----------------------------------------------------------------------------

CFG$priority_ready <- c("A_priority", "B_acceptable")

CFG$priority_whole_variable_risk <- c(
  "D_whole_variable_frame_risk_manual_check",
  "D_abnormal_variable_domain_manual_check",
  "D_conserved_cys_manual_check"
)


# -----------------------------------------------------------------------------
# tiny shared helpers
# -----------------------------------------------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Cannot create directory: ", path, call. = FALSE)
  invisible(path)
}

log_msg <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}

# -----------------------------------------------------------------------------
# short aliases for the settings used inside the graded QC expressions
#
# The values live in CFG, which stays the single source of truth. These aliases
# exist so that the deeply nested case_when() that grades each pair reads as a
# formula rather than a wall of CFG$ lookups.
# -----------------------------------------------------------------------------

strict_igg_subclass_only <- CFG$strict_igg_subclass_only

whole_v_prefix_len            <- CFG$whole_v_prefix_len
min_whole_v_identity_hard     <- CFG$min_whole_v_identity_hard
max_whole_v_mismatch_run_hard <- CFG$max_whole_v_mismatch_run_hard

FWR1_len_range <- CFG$FWR1_len_range
FWR2_len_range <- CFG$FWR2_len_range
FWR3_len_range <- CFG$FWR3_len_range


# =============================================================================
# FINDING THE INPUT FILES
# =============================================================================

# --- column contract ---------------------------------------------------------

# Absent -> the pipeline cannot function at all.
AIRR_COLUMNS_REQUIRED <- c(
  "sequence_id", "locus", "v_call", "j_call",
  "productive", "stop_codon", "vj_in_frame",
  "junction", "junction_aa", "cdr3", "cdr3_aa",
  "sequence_alignment_aa"
)

# Absent -> one QC dimension silently degrades to NA. Worth telling the user
# about, but not worth refusing the file: `igblastn -outfmt 19` versions differ
# in which of these they emit.
AIRR_COLUMNS_OPTIONAL <- c(
  "v_identity", "d_call", "c_call",
  "fwr1_aa", "fwr2_aa", "fwr3_aa", "fwr4_aa",
  "v_cigar", "d_cigar", "j_cigar", "complete_vdj"
)

# Files that are pipeline products rather than inputs. Catching them by name
# gives a much clearer message than "missing required columns" would.
KNOWN_NON_INPUT_STEMS <- c(
  "summary", "discovery_report", "all_paired_strict_igg",
  "expression_ready_ab", "manual_check", "whole_variable_risk",
  "excluded_inferred_igg", "qc_details", "all_key_info", "attention"
)


# --- helpers -----------------------------------------------------------------

# "sample01.tsv" -> "sample01". The stem is used verbatim: nothing is stripped
# or inferred beyond the extension, so the donor name in the output is always
# exactly what the file was called.
stem_of <- function(file_name) {
  tools::file_path_sans_ext(basename(file_name))
}

natural_order <- function(x) {
  if (length(x) == 0) return(character(0))
  keys <- strsplit(tolower(x), "(?<=[0-9])(?=[a-z])|(?<=[a-z])(?=[0-9])|[^a-z0-9]+",
                   perl = TRUE)
  ord <- order(
    vapply(keys, function(k) paste(sprintf("%08d", suppressWarnings(as.integer(k))),
                                   collapse = ""), character(1)),
    x
  )
  x[ord]
}

locus_of_v_call <- function(v_call) {
  top <- sub(",.*$", "", as.character(v_call))
  top <- sub("\\*.*$", "", top)
  out <- rep(NA_character_, length(top))
  out[grepl("^IGH", top)] <- "IGH"
  out[grepl("^IGK", top)] <- "IGK"
  out[grepl("^IGL", top)] <- "IGL"
  out
}

# Count IGH / IGK / IGL records without loading the whole 110-column table.
# Returns NA counts when the file cannot supply a usable locus signal.
count_loci <- function(path, columns) {
  empty <- c(IGH = NA_integer_, IGK = NA_integer_, IGL = NA_integer_)
  tryCatch({
    if ("locus" %in% columns) {
      d <- readr::read_tsv(path, col_select = dplyr::all_of("locus"),
                           show_col_types = FALSE, progress = FALSE)
      loc <- toupper(trimws(as.character(d$locus)))
    } else if ("v_call" %in% columns) {
      d <- readr::read_tsv(path, col_select = dplyr::all_of("v_call"),
                           show_col_types = FALSE, progress = FALSE)
      loc <- locus_of_v_call(d$v_call)
    } else {
      return(empty)
    }
    c(
      IGH = sum(loc == "IGH", na.rm = TRUE),
      IGK = sum(loc == "IGK", na.rm = TRUE),
      IGL = sum(loc == "IGL", na.rm = TRUE)
    )
  }, error = function(e) empty)
}

# Read the header row only. read_tsv with n_max = 0 keeps the column names.
read_tsv_header <- function(path) {
  d <- readr::read_tsv(path, n_max = 0, show_col_types = FALSE, progress = FALSE)
  colnames(d)
}


# --- inspection --------------------------------------------------------------

inspect_candidate <- function(path, cfg = CFG) {
  stem <- stem_of(path)

  res <- list(
    sample         = stem,
    file           = basename(path),
    path           = path,
    status         = "accepted",
    reason         = "",
    n_rows         = NA_integer_,
    igh_rows       = NA_integer_,
    igk_rows       = NA_integer_,
    igl_rows       = NA_integer_,
    missing_optional = ""
  )

  # 1. pipeline's own outputs are not inputs
  if (tolower(stem) %in% KNOWN_NON_INPUT_STEMS) {
    res$status <- "rejected"
    res$reason <- "looks like a pipeline output, not an IgBLAST input"
    return(res)
  }

  if (!nzchar(stem)) {
    res$status <- "rejected"
    res$reason <- "empty file stem, cannot derive a donor name"
    return(res)
  }

  # 2. header / column contract
  columns <- tryCatch(read_tsv_header(path), error = function(e) NULL)
  if (is.null(columns) || length(columns) == 0) {
    res$status <- "rejected"
    res$reason <- "unreadable as a tab-separated table"
    return(res)
  }

  missing_required <- setdiff(AIRR_COLUMNS_REQUIRED, columns)
  if (length(missing_required) > 0) {
    res$status <- "rejected"
    res$reason <- paste0("missing required AIRR column(s): ",
                         paste(head(missing_required, 6), collapse = ", "),
                         if (length(missing_required) > 6) ", ..." else "")
    return(res)
  }

  res$missing_optional <- paste(setdiff(AIRR_COLUMNS_OPTIONAL, columns), collapse = ",")

  # 3. does it hold immunoglobulin records at all?
  n <- tryCatch(
    nrow(readr::read_tsv(path, col_select = dplyr::all_of("sequence_id"),
                         show_col_types = FALSE, progress = FALSE)),
    error = function(e) NA_integer_
  )
  res$n_rows <- n

  loci <- count_loci(path, columns)
  res$igh_rows <- loci[["IGH"]]
  res$igk_rows <- loci[["IGK"]]
  res$igl_rows <- loci[["IGL"]]

  total_ig <- sum(loci, na.rm = TRUE)
  if (is.na(total_ig) || total_ig == 0) {
    res$status <- "rejected"
    res$reason <- "no IGH/IGK/IGL records found"
    return(res)
  }

  # 4. weak but usable: too few heavy chains to be the intended donor file
  if (!is.na(res$igh_rows) && res$igh_rows < cfg$min_igh_rows) {
    res$status <- "rejected"
    res$reason <- paste0("only ", res$igh_rows, " IGH record(s), below the ",
                         cfg$min_igh_rows, "-record floor (BCR_MIN_IGH_ROWS)")
    return(res)
  }

  res
}


# --- main entry point --------------------------------------------------------

#' Discover donor TSV files in a directory.
#'
#' @return list(samples = data.frame of accepted donors,
#'              report  = data.frame of every candidate with its status,
#'              rejected = subset of report with status == "rejected")
#' @export
discover_samples <- function(data_dir = CFG$data_dir, cfg = CFG) {
  if (!dir.exists(data_dir)) {
    stop(
      "Input directory not found: ", data_dir, "\n",
      "  Create it and drop one renamed IgBLAST TSV per donor into it,\n",
      "  for example data/sample01.tsv, data/sample02.tsv.",
      call. = FALSE
    )
  }

  files <- list.files(data_dir, pattern = "\\.tsv$", full.names = TRUE,
                      ignore.case = TRUE)

  if (length(files) == 0) {
    stop(
      "No .tsv files found in ", data_dir, "\n",
      "  Drop one renamed IgBLAST TSV per donor into that folder,\n",
      "  for example data/sample01.tsv, data/sample02.tsv.",
      call. = FALSE
    )
  }

  files <- natural_order(basename(files))
  files <- file.path(data_dir, files)

  # --- stem collisions are fatal -------------------------------------------
  # One file = one donor. Two files sharing a stem would be merged into a
  # single `sample` group and silently reported as one donor.
  stems <- vapply(files, stem_of, character(1))
  dup_stems <- unique(stems[duplicated(stems)])
  if (length(dup_stems) > 0) {
    detail <- vapply(dup_stems, function(s) {
      hits <- basename(files[stems == s])
      paste0("  '", s, "' <- ", paste(hits, collapse = ", "))
    }, character(1))
    stop(
      "Sample name collision in ", data_dir, ":\n",
      paste(detail, collapse = "\n"), "\n",
      "  The file name is the donor identity, so two files sharing a stem would\n",
      "  be merged into one donor without any error. Rename them so that every\n",
      "  donor has a unique name, for example 02_igblast_results.tsv -> sample01.tsv.",
      call. = FALSE
    )
  }

  inspected <- lapply(files, inspect_candidate, cfg = cfg)

  report <- do.call(rbind, lapply(inspected, function(x) {
    data.frame(
      sample           = x$sample,
      file             = x$file,
      status           = x$status,
      reason           = x$reason,
      n_rows           = x$n_rows,
      igh_rows         = x$igh_rows,
      igk_rows         = x$igk_rows,
      igl_rows         = x$igl_rows,
      missing_optional = x$missing_optional,
      stringsAsFactors = FALSE
    )
  }))

  accepted <- report[report$status == "accepted", , drop = FALSE]
  rejected <- report[report$status == "rejected", , drop = FALSE]

  if (nrow(accepted) == 0) {
    stop(
      "None of the ", nrow(report), " .tsv file(s) in ", data_dir,
      " passed the AIRR preflight check.\n",
      paste0("  - ", rejected$file, ": ", rejected$reason, collapse = "\n"),
      call. = FALSE
    )
  }

  # `report` already follows the natural order established above, because the
  # files were sorted before inspection. Re-sorting here would be worse than
  # redundant: natural_order() returns sorted *values*, not a permutation, so
  # order(natural_order(x)) compares the values as strings and turns the
  # numeric ordering back into a lexicographic one (sample10 before sample2).
  rownames(accepted) <- NULL
  rownames(report) <- NULL
  rownames(rejected) <- NULL

  list(
    samples  = data.frame(
      sample = accepted$sample,
      path   = file.path(data_dir, accepted$file),
      stringsAsFactors = FALSE
    ),
    report   = report,
    rejected = rejected
  )
}


#' Print the discovery outcome in a form a human can act on.
#' @export
report_discovery <- function(disc, cfg = CFG) {
  log_msg("Discovered ", nrow(disc$samples), " donor sample(s) in ", cfg$data_dir)

  for (i in seq_len(nrow(disc$samples))) {
    r <- disc$report[disc$report$sample == disc$samples$sample[i], , drop = FALSE][1, ]
    log_msg("  sample '", r$sample, "'  rows=", r$n_rows,
            "  IGH=", r$igh_rows, " IGK=", r$igk_rows, " IGL=", r$igl_rows,
            "  (", r$file, ")")
    if (nzchar(r$missing_optional)) {
      log_msg("    note: optional column(s) absent, related QC will read NA: ",
              r$missing_optional)
    }
  }

  if (nrow(disc$rejected) > 0) {
    log_msg("Skipped ", nrow(disc$rejected), " file(s) that are not usable donor input:")
    for (i in seq_len(nrow(disc$rejected))) {
      r <- disc$rejected[i, ]
      log_msg("  SKIP ", r$file, " -- ", r$reason,
              if (is.na(r$n_rows)) "" else paste0("  (rows=", r$n_rows, ")"))
    }
    log_msg("  These files are listed in output/discovery_report.tsv.")
  }

  invisible(disc)
}


#' Persist the discovery report so a run can be audited afterwards.
#' @export
write_discovery_report <- function(disc, output_dir = CFG$output_dir) {
  ensure_dir(output_dir)
  path <- file.path(output_dir, "discovery_report.tsv")
  readr::write_tsv(disc$report, path, na = "")
  path
}

# =============================================================================
# GERMLINE REFERENCES
# =============================================================================

#' Count FASTA records (lines starting with '>').
#' @export
count_fasta_records <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  lines <- readLines(path, warn = FALSE)
  as.integer(sum(startsWith(lines, ">")))
}


#' Read the header lines of a FASTA file.
#' @export
read_fasta_headers <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines[startsWith(lines, ">")]
}


#' Inspect one reference file. Never throws: returns a verdict.
#'
#' @return list(ok, records, headers, with_id, missing_id, digest, reason)
#' @export
inspect_reference_file <- function(path, gene, min_records) {
  verdict <- list(
    gene = gene, path = path, ok = FALSE, records = NA_integer_,
    headers = NA_integer_, with_id = NA_integer_, missing_id = NA_integer_,
    digest = NA_character_, reason = ""
  )

  if (!file.exists(path)) {
    verdict$reason <- "file not found"
    return(verdict)
  }
  if (file.size(path) == 0) {
    verdict$reason <- "file is empty"
    return(verdict)
  }

  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
  if (is.null(lines) || length(lines) == 0) {
    verdict$reason <- "file cannot be read as text"
    return(verdict)
  }

  headers <- lines[startsWith(lines, ">")]
  verdict$headers <- length(headers)

  if (length(headers) == 0) {
    verdict$reason <- "no FASTA header found (not a FASTA file?)"
    return(verdict)
  }

  # An HTML error page or a directory listing saved as .fasta lands here.
  if (grepl("<html|<!DOCTYPE|<head>", paste(head(headers, 3), collapse = " "),
            ignore.case = TRUE)) {
    verdict$reason <- "content looks like an HTML page, not FASTA"
    return(verdict)
  }

  # IMGT numbering gap characters are tolerated; digits/asterisks are not.
  body <- lines[!startsWith(lines, ">")]
  if (length(body) > 0) {
    body_clean <- toupper(gsub("[^A-Za-z.*-]", "", body))
    if (any(grepl("[^ACGTNacgtn.*-]", body_clean))) {
      verdict$reason <- "sequence lines contain non-nucleotide characters"
      return(verdict)
    }
  }

  verdict$records <- length(headers)

  # a-z is required, not decorative: IMGT names the OR genes with lowercase
  # letters (IGHD1/OR15-1a*01), so an uppercase-only pattern reports 10 of the
  # 47 IGHD headers as id-less and would reject a perfectly good file.
  verdict$with_id <- sum(grepl("IG[HKL][VDJ][A-Za-z0-9/-]+\\*[0-9]+", headers))
  verdict$missing_id <- verdict$headers - verdict$with_id

  # pick_ref() falls back to a gene-level match, so a handful of id-less
  # headers is survivable; a systematic loss is not.
  if (verdict$with_id < 0.9 * verdict$records) {
    verdict$reason <- paste0(
      "only ", verdict$with_id, " of ", verdict$records,
      " headers carry a standard *allele id"
    )
    return(verdict)
  }

  if (verdict$records < min_records) {
    verdict$reason <- paste0(
      "only ", verdict$records, " records, below the floor of ", min_records,
      " (truncated or outdated download)"
    )
    return(verdict)
  }

  verdict$digest <- tryCatch(
    as.character(tools::md5sum(path)),
    error = function(e) NA_character_
  )
  verdict$ok <- TRUE
  verdict$reason <- "ok"
  verdict
}


#' Download one file, trying several transports before giving up.
#'
#' Returns TRUE on success. The caller validates the result; a 200 response
#' carrying an error page must not count as success.
#' @export
download_with_fallback <- function(url, dest, timeout = 120L) {
  partial <- paste0(dest, ".partial")
  unlink(partial)

  attempt <- function(label, expr) {
    ok <- tryCatch({
      force(expr)
      TRUE
    }, error = function(e) {
      log_msg("    ", label, " failed: ", conditionMessage(e))
      FALSE
    }, warning = function(w) {
      log_msg("    ", label, " warning: ", conditionMessage(w))
      FALSE
    })
    ok
  }

  methods <- c("libcurl", "wininet", "curl")

  for (m in methods) {
    if (!capabilities("libcurl") && m == "libcurl") next
    log_msg("    trying download method '", m, "' ...")
    unlink(partial)
    ok <- attempt(m, suppressWarnings(utils::download.file(
      url, partial, method = m, quiet = TRUE, mode = "wb", timeout = timeout,
      headers = c("User-Agent" = "SangerBCR-Flow-expression/1.0")
    )))
    if (ok && file.exists(partial) && file.size(partial) > 0) {
      if (file.rename(partial, dest)) return(TRUE)
      # Cross-device rename can fail; fall back to copy.
      if (file.copy(partial, dest, overwrite = TRUE)) {
        unlink(partial)
        return(TRUE)
      }
    }
  }

  unlink(partial)
  FALSE
}


#' Build a copy-pasteable instruction block for obtaining a reference by hand.
#' @export
manual_reference_help <- function(genes, base_url, ref_dir) {
  lines <- c(
    "  Download the germline FASTA file(s) manually and put them in:",
    paste0("    ", ref_dir),
    "",
    "  Bash / WSL:",
    paste0("    mkdir -p '", ref_dir, "'"),
    paste0("    cd '", ref_dir, "'"),
    paste0("    for g in ", paste(genes, collapse = " "), "; do"),
    paste0("      curl -fL -o \"$g.fasta\" '", base_url, "/$g.fasta'"),
    "    done",
    "",
    "  PowerShell:",
    paste0("    New-Item -ItemType Directory -Force '", ref_dir, "' | Out-Null"),
    paste0("    foreach ($g in @('", paste(genes, collapse = "','"), "')) {"),
    paste0("      Invoke-WebRequest \"", base_url,
           "/$g.fasta\" -OutFile (Join-Path '", ref_dir, "' \"$g.fasta\")"),
    "    }",
    "",
    "  The 7 files are also downloadable from the IMGT V-QUEST reference",
    "  directory for Homo sapiens (IG). A machine with working TLS may be",
    "  needed: this pipeline cannot download when R, curl and PowerShell all",
    "  fail the TLS handshake."
  )
  paste(lines, collapse = "\n")
}


#' Ensure all needed references are present and valid.
#'
#' @param genes character vector of gene families
#' @param cfg   configuration list
#' @return list(tables = named list of file paths, provenance = data.frame)
#' @export
ensure_references <- function(genes = CFG$reference_genes, cfg = CFG) {
  ensure_dir(cfg$ref_dir)

  paths <- file.path(cfg$ref_dir, paste0(genes, ".fasta"))
  min_records <- cfg$ref_min_records[genes]

  verdicts <- lapply(seq_along(genes), function(i) {
    inspect_reference_file(paths[i], genes[i], min_records[i])
  })
  names(verdicts) <- genes

  missing <- genes[!vapply(verdicts, function(v) isTRUE(v$ok), logical(1))]

  # Say what is wrong before trying to fix it: an invalid local file is a
  # different situation from an absent one.
  for (g in missing) {
    v <- verdicts[[g]]
    log_msg("Reference ", g, " unusable: ", v$reason)
  }

  if (length(missing) > 0 && isTRUE(cfg$download_reference)) {
    log_msg("Fetching ", length(missing), " reference file(s) from IMGT ...")
    for (g in missing) {
      url <- paste0(cfg$imgt_base_url, "/", g, ".fasta")
      log_msg("  ", g, " <- ", url)
      ok <- download_with_fallback(url, file.path(cfg$ref_dir, paste0(g, ".fasta")),
                                  timeout = cfg$download_timeout)
      if (!ok) {
        log_msg("  ", g, ": every download method failed")
        next
      }
      verdicts[[g]] <- inspect_reference_file(
        file.path(cfg$ref_dir, paste0(g, ".fasta")), g, min_records[[g]]
      )
      if (isTRUE(verdicts[[g]]$ok)) {
        log_msg("  ", g, ": ok, ", verdicts[[g]]$records, " records")
      } else {
        log_msg("  ", g, ": downloaded file rejected -- ", verdicts[[g]]$reason)
      }
    }
  }

  still_bad <- genes[!vapply(verdicts, function(v) isTRUE(v$ok), logical(1))]
  if (length(still_bad) > 0) {
    detail <- vapply(still_bad, function(g) {
      paste0("  - ", g, ".fasta: ", verdicts[[g]]$reason)
    }, character(1))
    stop(
      "Cannot obtain ", length(still_bad), " germline reference file(s):\n",
      paste(detail, collapse = "\n"), "\n\n",
      manual_reference_help(still_bad, cfg$imgt_base_url, cfg$ref_dir),
      call. = FALSE
    )
  }

  provenance <- do.call(rbind, lapply(genes, function(g) {
    v <- verdicts[[g]]
    data.frame(
      gene     = g,
      records  = v$records,
      headers  = v$headers,
      with_id  = v$with_id,
      md5      = v$digest,
      file     = basename(v$path),
      stringsAsFactors = FALSE
    )
  }))
  rownames(provenance) <- NULL

  list(
    paths      = stats::setNames(paths, genes),
    provenance = provenance
  )
}


# ============================================================
# 1. basic helpers
# ============================================================

flag_true <- function(x) {
  toupper(as.character(x)) %in% c("T", "TRUE", "1", "YES")
}

flag_false <- function(x) {
  toupper(as.character(x)) %in% c("F", "FALSE", "0", "NO")
}

is_blank <- function(x) {
  is.na(x) | str_trim(as.character(x)) == "" | toupper(as.character(x)) == "NA"
}

top_call <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- str_split_fixed(x, ",", 2)[, 1]
  str_trim(x)
}

clean_gene <- function(x) {
  top_call(x) %>%
    str_remove("\\*.*$")
}

clean_aa <- function(x) {
  x <- toupper(as.character(x))
  x[is.na(x)] <- ""
  x <- str_replace_all(x, "[\\.\\-\\s_]", "")
  x <- str_replace_all(x, "[^A-Z\\*]", "")
  x
}

valid_aa <- function(x) {
  x <- clean_aa(x)
  x != "" & str_detect(x, "^[ACDEFGHIKLMNPQRSTVWY]+$")
}

clean_nt <- function(x) {
  x <- toupper(as.character(x))
  x[is.na(x)] <- ""
  str_replace_all(x, "[^ACGTN]", "")
}

nt_frame_len_ok <- function(x) {
  nt <- clean_nt(x)
  out <- rep(NA, length(nt))
  has <- nt != ""
  out[has] <- nchar(nt[has]) %% 3 == 0
  out
}

raw_aa_has_stop <- function(x) {
  x <- toupper(as.character(x))
  x[is.na(x)] <- ""
  str_detect(x, "\\*")
}

# An I/D run in an AIRR CIGAR whose length is not a multiple of 3 hints at a
# frameshift indel. A missing CIGAR returns NA and is not treated as a failure.
cigar_nontriplet_indel_one <- function(x) {
  if (is.na(x) || str_trim(as.character(x)) == "") return(NA)
  x <- toupper(as.character(x))
  m <- str_match_all(x, "([0-9]+)([ID])")[[1]]
  if (nrow(m) == 0) return(FALSE)
  lens <- suppressWarnings(as.integer(m[, 2]))
  any(!is.na(lens) & lens %% 3 != 0)
}

cigar_nontriplet_indel <- function(x) {
  vapply(x, cigar_nontriplet_indel_one, logical(1))
}

junction_translation_ok_one <- function(nt, aa) {
  if (is.na(nt) || is.na(aa)) return(NA)
  nt <- clean_nt(nt)
  aa <- clean_aa(aa)
  if (nt == "" || aa == "") return(NA)
  if (nchar(nt) %% 3 != 0) return(FALSE)
  if (str_detect(nt, "N")) return(NA)
  tr <- translate_nt_frame(nt, 1)
  if (tr == "" || str_detect(tr, "[X\\*]")) return(FALSE)
  identical(tr, aa)
}

junction_translation_ok <- function(nt_vec, aa_vec) {
  mapply(
    junction_translation_ok_one,
    nt_vec,
    aa_vec,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
}

framework_region_ok_one <- function(x, min_len, max_len) {
  if (is.na(x) || str_trim(as.character(x)) == "") return(NA)
  raw <- toupper(as.character(x))
  aa <- clean_aa(raw)
  if (str_detect(raw, "\\*|X")) return(FALSE)
  if (!valid_aa(aa)) return(FALSE)
  L <- nchar(aa)
  L >= min_len && L <= max_len
}

framework_profile_one <- function(fwr1, fwr2, fwr3) {
  q1 <- framework_region_ok_one(fwr1, FWR1_len_range[1], FWR1_len_range[2])
  q2 <- framework_region_ok_one(fwr2, FWR2_len_range[1], FWR2_len_range[2])
  q3 <- framework_region_ok_one(fwr3, FWR3_len_range[1], FWR3_len_range[2])
  q <- c(q1, q2, q3)
  if (all(is.na(q))) return("not_available")
  if (any(q == FALSE, na.rm = TRUE)) return("fail")
  "pass"
}

framework_profile <- function(fwr1, fwr2, fwr3) {
  mapply(
    framework_profile_one,
    fwr1,
    fwr2,
    fwr3,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
}

max_mismatch_run <- function(a, b) {
  a <- clean_aa(a)
  b <- clean_aa(b)
  L <- min(nchar(a), nchar(b))
  if (L < 1) return(NA_integer_)
  aa <- strsplit(substr(a, 1, L), "")[[1]]
  bb <- strsplit(substr(b, 1, L), "")[[1]]
  bad <- aa != bb
  r <- rle(bad)
  if (!any(r$values)) return(0L)
  as.integer(max(r$lengths[r$values]))
}

compare_v_prefix_one <- function(final_aa, v_ref_id, ref_tbl, prefix_len = 80) {
  final_aa <- clean_aa(final_aa)
  v_ref_id <- as.character(v_ref_id)
  if (final_aa == "" || is.na(v_ref_id) || v_ref_id == "") {
    return(list(identity = NA_real_, mismatch_run = NA_integer_, compared_len = NA_integer_))
  }

  hit <- ref_tbl %>% filter(.data$id == .env$v_ref_id)
  if (nrow(hit) == 0) {
    return(list(identity = NA_real_, mismatch_run = NA_integer_, compared_len = NA_integer_))
  }

  ref_aa <- clean_aa(hit$aa[1])
  L <- min(prefix_len, nchar(final_aa), nchar(ref_aa))
  if (L < 40) {
    return(list(identity = NA_real_, mismatch_run = NA_integer_, compared_len = L))
  }

  a <- substr(final_aa, 1, L)
  b <- substr(ref_aa, 1, L)
  aa <- strsplit(a, "")[[1]]
  bb <- strsplit(b, "")[[1]]

  list(
    identity = mean(aa == bb),
    mismatch_run = max_mismatch_run(a, b),
    compared_len = L
  )
}

compare_v_prefix_vec <- function(final_vec, ref_id_vec, ref_tbl, prefix) {
  z <- mapply(
    function(a, r) compare_v_prefix_one(a, r, ref_tbl, whole_v_prefix_len),
    final_vec,
    ref_id_vec,
    SIMPLIFY = FALSE
  )

  data.frame(
    whole_v_identity = sapply(z, function(x) x[["identity"]]),
    whole_v_max_mismatch_run = sapply(z, function(x) x[["mismatch_run"]]),
    whole_v_compared_len = sapply(z, function(x) x[["compared_len"]]),
    stringsAsFactors = FALSE
  ) %>%
    setNames(paste0(prefix, names(.)))
}

coalesce_text <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  a[is.na(a)] <- ""
  b[is.na(b)] <- ""
  idx <- is_blank(a)
  a[idx] <- b[idx]
  a
}

# AIRR defines cdr3_aa as excluding the conserved flanking residues, so every C
# in it is an internal CDR3 cysteine. junction_aa does include the 5' conserved
# Cys, so only that first position is discounted -- and only when it really is
# a C, to avoid undercounting on a truncated junction.
count_extra_cdr3_cys_one <- function(seq_aa, source) {
  aa <- clean_aa(seq_aa)
  source <- as.character(source)

  if (aa == "" || is.na(source) || source == "") {
    return(NA_integer_)
  }

  if (source == "cdr3_aa") {
    return(as.integer(str_count(aa, "C")))
  }

  if (source == "junction_aa") {
    aa_internal <- if (str_starts(aa, "C") && nchar(aa) > 1) {
      substr(aa, 2, nchar(aa))
    } else if (aa == "C") {
      ""
    } else {
      aa
    }

    return(as.integer(str_count(aa_internal, "C")))
  }

  NA_integer_
}

count_extra_cdr3_cys <- function(seq_vec, source_vec) {
  mapply(
    count_extra_cdr3_cys_one,
    seq_vec,
    source_vec,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
}

# Whether the extra-cysteine counter above says nothing about: can the two
# conserved cysteines that hold an immunoglobulin domain together actually be
# found? A chain that LOST one of them scores zero extras, so an extras-only
# check reads it as clean. junction_aa begins with the CDR3 cysteine whenever
# the boundary is right, so it answers the 3' half directly.
#
# NA means the input carried no junction_aa for that chain; that is left
# unknown rather than failed, and the caller falls back to the cysteine count.
junction_cys_ok <- function(x) {
  aa <- clean_aa(x)
  out <- rep(NA, length(aa))
  has <- aa != ""
  out[has] <- str_starts(aa[has], "C")
  out
}

# Recover the cell identity from a sequence_id.
#
# A Sanger record id identifies the well and the chain, for example
# "A11H1.<run>.<...>" for the heavy chain of well A11 and "A11K.<run>.<...>"
# for its kappa chain. Heavy and light reads are paired on the well name, so the
# chain marker is stripped from the first dot-separated field.
#
# The suffix pattern is a configuration value (BCR_CHAIN_SUFFIX_PATTERN),
# because the naming convention is a property of the sequencing run, not of this
# pipeline. See the SETTINGS block at the top of this file for the default and
# for how to override it.
extract_cell_id <- function(sequence_id) {
  x <- as.character(sequence_id)
  x[is.na(x)] <- ""

  x %>%
    str_split("\\.", simplify = TRUE) %>%
    .[, 1] %>%
    str_remove(CFG$chain_suffix_pattern) %>%
    str_remove("[-_]+$")
}

# The class includes a-z on purpose. IMGT names the OR genes with lowercase
# letters (IGHD1/OR15-1a*01, IGHD1/OR15-1b*01), so an uppercase-only pattern
# fails to match 10 of the 47 IGHD headers and those alleles fall back to a
# gene-level match.
IMGT_ALLELE_PATTERN <- "IG[HKL][VDJ][A-Za-z0-9/-]+\\*[0-9]+"

get_imgt_id <- function(x) {
  x <- top_call(x)
  id <- str_extract(x, IMGT_ALLELE_PATTERN)
  ifelse(is.na(id), x, id)
}

last_n <- function(x, n = 16) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  start <- pmax(1, nchar(x) - n + 1)
  substr(x, start, nchar(x))
}

classify_igg <- function(c_call, sequence_id) {
  c_top <- toupper(top_call(c_call))
  sid <- toupper(as.character(sequence_id))
  sid[is.na(sid)] <- ""

  subclass <- str_extract(c_top, "IGHG[1-4]")

  igg_source <- case_when(
    !is.na(subclass) ~ "IgG_subclass_called_by_c_call",
    str_detect(c_top, "^IGHG") ~ "IgG_called_by_c_call_no_subclass",
    str_detect(sid, "HIGG|IGHG") ~ "IgG_inferred_from_sequence_id",
    TRUE ~ "not_IgG"
  )

  data.frame(
    heavy_subclass = ifelse(is.na(subclass), "", subclass),
    igg_source = igg_source,
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 2. DNA translation
# ============================================================

genetic_code <- c(
  TTT="F", TTC="F", TTA="L", TTG="L",
  TCT="S", TCC="S", TCA="S", TCG="S",
  TAT="Y", TAC="Y", TAA="*", TAG="*",
  TGT="C", TGC="C", TGA="*", TGG="W",
  CTT="L", CTC="L", CTA="L", CTG="L",
  CCT="P", CCC="P", CCA="P", CCG="P",
  CAT="H", CAC="H", CAA="Q", CAG="Q",
  CGT="R", CGC="R", CGA="R", CGG="R",
  ATT="I", ATC="I", ATA="I", ATG="M",
  ACT="T", ACC="T", ACA="T", ACG="T",
  AAT="N", AAC="N", AAA="K", AAG="K",
  AGT="S", AGC="S", AGA="R", AGG="R",
  GTT="V", GTC="V", GTA="V", GTG="V",
  GCT="A", GCC="A", GCA="A", GCG="A",
  GAT="D", GAC="D", GAA="E", GAG="E",
  GGT="G", GGC="G", GGA="G", GGG="G"
)

translate_nt_frame <- function(nt, frame = 1) {
  nt <- toupper(as.character(nt))
  nt[is.na(nt)] <- ""
  nt <- str_replace_all(nt, "[^ACGTN]", "")

  if (nchar(nt) < frame + 2) return("")

  nt <- substr(nt, frame, nchar(nt))
  n <- nchar(nt)
  n <- n - n %% 3

  if (n < 3) return("")

  nt <- substr(nt, 1, n)
  codons <- substring(nt, seq(1, n, 3), seq(3, n, 3))
  aa <- genetic_code[codons]
  aa[is.na(aa)] <- "X"

  paste0(aa, collapse = "")
}

# Heuristic frame choice for short J references: the FR4 motif is a much
# stronger signal than raw length, hence the +100 reward. The weights are
# empirical, not model-based.
translate_best_frame <- function(nt) {
  aas <- sapply(1:3, function(fr) translate_nt_frame(nt, fr))

  scores <- sapply(aas, function(aa) {
    if (aa == "") return(-9999)

    score <- nchar(aa)
    score <- score - 50 * str_count(aa, "\\*")
    score <- score - 20 * str_count(aa, "X")

    if (str_detect(aa, "WG[QR]G|FGQG|FGGG|FGTG|FGPG|FGEG")) {
      score <- score + 100
    }

    score
  })

  aas[which.max(scores)]
}


# ============================================================
# 3. find the input files
#
# Done before the references are fetched: a missing or unusable input
# directory is the more fundamental problem, and it should be reported
# immediately rather than after a reference download.
# ============================================================

read_fasta_ref <- function(file, min_aa = 5, frame_mode = c("first", "best")) {
  frame_mode <- match.arg(frame_mode)

  if (!file.exists(file)) {
    stop(paste0("找不到参考文件: ", file))
  }

  lines <- readLines(file, warn = FALSE)
  head_idx <- which(str_starts(lines, ">"))

  if (length(head_idx) == 0) {
    stop(paste0("FASTA 文件没有检测到 header: ", file))
  }

  out <- lapply(seq_along(head_idx), function(i) {
    start <- head_idx[i]
    end <- ifelse(i == length(head_idx), length(lines), head_idx[i + 1] - 1)

    header <- lines[start]
    seq_lines <- if (start < end) lines[(start + 1):end] else character(0)

    nt <- paste0(seq_lines, collapse = "")
    nt <- toupper(str_replace_all(nt, "[^ACGTNacgtn]", ""))

    id <- str_extract(header, IMGT_ALLELE_PATTERN)
    gene <- str_remove(id, "\\*\\d+$")

    aa <- if (frame_mode == "first") {
      translate_nt_frame(nt, 1)
    } else {
      translate_best_frame(nt)
    }

    data.frame(
      id = id,
      gene = gene,
      nt = nt,
      aa = aa,
      header = header,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out) %>%
    filter(!is.na(id), aa != "") %>%
    filter(!str_detect(aa, "[X\\*]")) %>%
    filter(nchar(aa) >= min_aa)
}

# ============================================================
# 4. read the germline references
# ============================================================

ensure_dir(CFG$output_dir)

disc <- discover_samples()
report_discovery(disc)
disc_report_path <- write_discovery_report(disc)

refs <- ensure_references()
log_msg("References loaded from ", CFG$ref_dir)
ref_provenance <- refs$provenance

IGHV_ref <- read_fasta_ref(refs$paths[["IGHV"]], min_aa = 70, frame_mode = "first")
IGKV_ref <- read_fasta_ref(refs$paths[["IGKV"]], min_aa = 70, frame_mode = "first")
IGLV_ref <- read_fasta_ref(refs$paths[["IGLV"]], min_aa = 70, frame_mode = "first")

IGHJ_ref <- read_fasta_ref(refs$paths[["IGHJ"]], min_aa = 5, frame_mode = "best")
IGKJ_ref <- read_fasta_ref(refs$paths[["IGKJ"]], min_aa = 5, frame_mode = "best")
IGLJ_ref <- read_fasta_ref(refs$paths[["IGLJ"]], min_aa = 5, frame_mode = "best")

log_msg("Reference tables: IGHV ", nrow(IGHV_ref),
        ", IGKV ", nrow(IGKV_ref),
        ", IGLV ", nrow(IGLV_ref),
        ", IGHJ ", nrow(IGHJ_ref),
        ", IGKJ ", nrow(IGKJ_ref),
        ", IGLJ ", nrow(IGLJ_ref))


# ============================================================
# 5. reference matching
# ============================================================

pick_ref <- function(call, ref_tbl) {
  allele_query <- get_imgt_id(call)
  gene_query <- str_remove(allele_query, "\\*\\d+$")

  hit <- ref_tbl %>% filter(.data$id == .env$allele_query)
  if (nrow(hit) > 0) {
    return(list(ref = hit[1, ], match_type = "exact_allele"))
  }

  # An earlier version wrote filter(gene == gene), which compares the column
  # with itself and degenerates into "always TRUE", producing a wrong fallback.
  hit <- ref_tbl %>% filter(.data$gene == .env$gene_query)
  if (nrow(hit) > 0) {
    return(list(ref = hit[1, ], match_type = "gene_fallback"))
  }

  return(NULL)
}


# ============================================================
# 6. complete the V-domain N terminus
# ============================================================

find_start_in_v_ref <- function(
    obs_aa,
    ref_aa,
    anchor_len = 24,
    max_start_scan = 25
) {
  obs_aa <- clean_aa(obs_aa)
  ref_aa <- clean_aa(ref_aa)

  L <- min(nchar(obs_aa), anchor_len)

  if (L < 8) {
    return(list(start = NA, identity = NA))
  }

  query <- strsplit(substr(obs_aa, 1, L), "")[[1]]
  max_start_scan <- min(max_start_scan, nchar(ref_aa) - L + 1)

  if (max_start_scan < 1) {
    return(list(start = NA, identity = NA))
  }

  scores <- sapply(seq_len(max_start_scan), function(s) {
    ref_piece <- strsplit(substr(ref_aa, s, s + L - 1), "")[[1]]
    mean(query == ref_piece)
  })

  best_start <- which.max(scores)
  best_identity <- max(scores)

  list(start = best_start, identity = best_identity)
}

complete_N_terminal <- function(
    obs_aa,
    v_call,
    v_ref_tbl,
    min_anchor_hard = 0.65,
    max_added_N_hard = 15
) {
  obs_aa <- clean_aa(obs_aa)

  if (!valid_aa(obs_aa)) {
    return(list(
      final_aa = obs_aa,
      added_N_aa = "",
      added_N_len = NA,
      v_ref_id = NA,
      v_ref_gene = NA,
      v_ref_match_type = NA,
      start_in_v_ref = NA,
      v_anchor_identity = NA,
      N_note = "invalid_obs_aa_check_manually"
    ))
  }

  picked <- pick_ref(v_call, v_ref_tbl)

  if (is.null(picked)) {
    return(list(
      final_aa = obs_aa,
      added_N_aa = "",
      added_N_len = 0,
      v_ref_id = NA,
      v_ref_gene = NA,
      v_ref_match_type = NA,
      start_in_v_ref = NA,
      v_anchor_identity = NA,
      N_note = "no_V_reference_found_check_manually"
    ))
  }

  ref <- picked$ref
  pos <- find_start_in_v_ref(obs_aa, ref$aa)

  if (is.na(pos$start)) {
    return(list(
      final_aa = obs_aa,
      added_N_aa = "",
      added_N_len = 0,
      v_ref_id = ref$id,
      v_ref_gene = ref$gene,
      v_ref_match_type = picked$match_type,
      start_in_v_ref = NA,
      v_anchor_identity = NA,
      N_note = "cannot_align_to_V_reference_check_manually"
    ))
  }

  added_len <- pos$start - 1

  if (!is.na(pos$identity) && pos$identity < min_anchor_hard) {
    return(list(
      final_aa = obs_aa,
      added_N_aa = "",
      added_N_len = added_len,
      v_ref_id = ref$id,
      v_ref_gene = ref$gene,
      v_ref_match_type = picked$match_type,
      start_in_v_ref = pos$start,
      v_anchor_identity = pos$identity,
      N_note = "low_V_anchor_identity_check_manually"
    ))
  }

  if (added_len > max_added_N_hard) {
    return(list(
      final_aa = obs_aa,
      added_N_aa = "",
      added_N_len = added_len,
      v_ref_id = ref$id,
      v_ref_gene = ref$gene,
      v_ref_match_type = picked$match_type,
      start_in_v_ref = pos$start,
      v_anchor_identity = pos$identity,
      N_note = "too_much_N_terminal_missing_check_manually"
    ))
  }

  added <- if (added_len > 0) substr(ref$aa, 1, added_len) else ""
  final <- paste0(added, obs_aa)

  note <- if (added_len > 0) {
    "N_terminal_completed_by_V_germline"
  } else {
    "no_N_terminal_completion_needed"
  }

  list(
    final_aa = final,
    added_N_aa = added,
    added_N_len = added_len,
    v_ref_id = ref$id,
    v_ref_gene = ref$gene,
    v_ref_match_type = picked$match_type,
    start_in_v_ref = pos$start,
    v_anchor_identity = pos$identity,
    N_note = note
  )
}


# ============================================================
# 7. complete the J-domain C terminus
# ============================================================

find_end_in_j_ref <- function(
    obs_aa,
    j_ref_aa,
    min_overlap = 5,
    max_overlap = 24
) {
  obs_aa <- clean_aa(obs_aa)
  j_ref_aa <- clean_aa(j_ref_aa)

  obs_len <- nchar(obs_aa)
  ref_len <- nchar(j_ref_aa)
  max_L <- min(obs_len, ref_len, max_overlap)

  if (max_L < min_overlap) {
    return(list(
      start = NA,
      end = NA,
      identity = NA,
      overlap_len = NA
    ))
  }

  best_score <- -Inf
  best_identity <- NA
  best_L <- -1
  best_start <- NA
  best_end <- NA

  for (L in seq(min_overlap, max_L)) {
    obs_suffix <- substr(obs_aa, obs_len - L + 1, obs_len)
    obs_chars <- strsplit(obs_suffix, "")[[1]]

    for (s in seq_len(ref_len - L + 1)) {
      ref_piece <- substr(j_ref_aa, s, s + L - 1)
      ref_chars <- strsplit(ref_piece, "")[[1]]

      matches <- sum(obs_chars == ref_chars)
      identity <- matches / L
      e <- s + L - 1

      # Reward identity and overlap length together: ranking on match count
      # alone lets a long, low-identity overlap outrank a short, trustworthy
      # one.
      score <- identity * log2(L + 1)

      better <- FALSE
      if (score > best_score) {
        better <- TRUE
      } else if (isTRUE(all.equal(score, best_score)) && identity > best_identity) {
        better <- TRUE
      } else if (isTRUE(all.equal(score, best_score)) &&
                 isTRUE(all.equal(identity, best_identity)) && L > best_L) {
        better <- TRUE
      } else if (isTRUE(all.equal(score, best_score)) &&
                 isTRUE(all.equal(identity, best_identity)) && L == best_L && e > best_end) {
        better <- TRUE
      }

      if (better) {
        best_score <- score
        best_identity <- identity
        best_L <- L
        best_start <- s
        best_end <- e
      }
    }
  }

  list(
    start = best_start,
    end = best_end,
    identity = best_identity,
    overlap_len = best_L
  )
}

complete_C_terminal_by_J <- function(
    obs_aa,
    j_call,
    j_ref_tbl,
    min_j_overlap_hard = 0.65,
    max_added_C_hard = 8
) {
  obs_aa <- clean_aa(obs_aa)

  if (!valid_aa(obs_aa)) {
    return(list(
      final_aa = obs_aa,
      added_C_aa = "",
      added_C_len = NA,
      j_ref_id = NA,
      j_ref_gene = NA,
      j_ref_match_type = NA,
      j_overlap_identity = NA,
      j_overlap_len = NA,
      j_end_in_ref = NA,
      J_note = "invalid_obs_aa_for_J_check_manually"
    ))
  }

  picked <- pick_ref(j_call, j_ref_tbl)

  if (is.null(picked)) {
    return(list(
      final_aa = obs_aa,
      added_C_aa = "",
      added_C_len = 0,
      j_ref_id = NA,
      j_ref_gene = NA,
      j_ref_match_type = NA,
      j_overlap_identity = NA,
      j_overlap_len = NA,
      j_end_in_ref = NA,
      J_note = "no_J_reference_found_check_manually"
    ))
  }

  ref <- picked$ref
  pos <- find_end_in_j_ref(obs_aa, ref$aa)

  if (is.na(pos$end)) {
    return(list(
      final_aa = obs_aa,
      added_C_aa = "",
      added_C_len = 0,
      j_ref_id = ref$id,
      j_ref_gene = ref$gene,
      j_ref_match_type = picked$match_type,
      j_overlap_identity = NA,
      j_overlap_len = NA,
      j_end_in_ref = NA,
      J_note = "cannot_align_to_J_reference_check_manually"
    ))
  }

  ref_len <- nchar(ref$aa)
  added_len <- ref_len - pos$end

  if (!is.na(pos$identity) && pos$identity < min_j_overlap_hard) {
    return(list(
      final_aa = obs_aa,
      added_C_aa = "",
      added_C_len = added_len,
      j_ref_id = ref$id,
      j_ref_gene = ref$gene,
      j_ref_match_type = picked$match_type,
      j_overlap_identity = pos$identity,
      j_overlap_len = pos$overlap_len,
      j_end_in_ref = pos$end,
      J_note = "low_J_overlap_identity_check_manually"
    ))
  }

  if (added_len > max_added_C_hard) {
    return(list(
      final_aa = obs_aa,
      added_C_aa = "",
      added_C_len = added_len,
      j_ref_id = ref$id,
      j_ref_gene = ref$gene,
      j_ref_match_type = picked$match_type,
      j_overlap_identity = pos$identity,
      j_overlap_len = pos$overlap_len,
      j_end_in_ref = pos$end,
      J_note = "too_much_C_terminal_missing_check_manually"
    ))
  }

  added <- if (added_len > 0) {
    substr(ref$aa, pos$end + 1, ref_len)
  } else {
    ""
  }

  final <- paste0(obs_aa, added)

  note <- if (added_len > 0) {
    "C_terminal_completed_by_J_germline"
  } else {
    "no_C_terminal_completion_needed"
  }

  list(
    final_aa = final,
    added_C_aa = added,
    added_C_len = added_len,
    j_ref_id = ref$id,
    j_ref_gene = ref$gene,
    j_ref_match_type = picked$match_type,
    j_overlap_identity = pos$identity,
    j_overlap_len = pos$overlap_len,
    j_end_in_ref = pos$end,
    J_note = note
  )
}


# ============================================================
# 8. complete one VH/VL
# ============================================================

completion_error <- function(obs_aa, note) {
  list(
    final_aa = clean_aa(obs_aa),
    added_N_aa = "",
    added_N_len = NA,
    v_ref_id = NA,
    v_ref_gene = NA,
    v_ref_match_type = NA,
    start_in_v_ref = NA,
    v_anchor_identity = NA,
    N_note = note,
    added_C_aa = "",
    added_C_len = NA,
    j_ref_id = NA,
    j_ref_gene = NA,
    j_ref_match_type = NA,
    j_overlap_identity = NA,
    j_overlap_len = NA,
    j_end_in_ref = NA,
    J_note = note
  )
}

complete_chain <- function(obs_aa, v_call, j_call, v_ref_tbl, j_ref_tbl) {
  n <- complete_N_terminal(
    obs_aa = obs_aa,
    v_call = v_call,
    v_ref_tbl = v_ref_tbl,
    min_anchor_hard = CFG$min_anchor_hard,
    max_added_N_hard = CFG$max_added_N_hard
  )

  c <- complete_C_terminal_by_J(
    obs_aa = n$final_aa,
    j_call = j_call,
    j_ref_tbl = j_ref_tbl,
    min_j_overlap_hard = CFG$min_j_overlap_hard,
    max_added_C_hard = CFG$max_added_C_hard
  )

  list(
    final_aa = c$final_aa,

    added_N_aa = n$added_N_aa,
    added_N_len = n$added_N_len,
    v_ref_id = n$v_ref_id,
    v_ref_gene = n$v_ref_gene,
    v_ref_match_type = n$v_ref_match_type,
    start_in_v_ref = n$start_in_v_ref,
    v_anchor_identity = n$v_anchor_identity,
    N_note = n$N_note,

    added_C_aa = c$added_C_aa,
    added_C_len = c$added_C_len,
    j_ref_id = c$j_ref_id,
    j_ref_gene = c$j_ref_gene,
    j_ref_match_type = c$j_ref_match_type,
    j_overlap_identity = c$j_overlap_identity,
    j_overlap_len = c$j_overlap_len,
    j_end_in_ref = c$j_end_in_ref,
    J_note = c$J_note
  )
}

completion_fields <- c(
  "final_aa",
  "added_N_aa", "added_N_len", "v_ref_id", "v_ref_gene", "v_ref_match_type",
  "start_in_v_ref", "v_anchor_identity", "N_note",
  "added_C_aa", "added_C_len", "j_ref_id", "j_ref_gene", "j_ref_match_type",
  "j_overlap_identity", "j_overlap_len", "j_end_in_ref", "J_note"
)

completion_frame <- function(x, prefix) {
  out <- data.frame(
    lapply(completion_fields, function(f) sapply(x, function(z) z[[f]])),
    stringsAsFactors = FALSE
  )
  colnames(out) <- paste0(prefix, completion_fields)
  out
}

complete_chain_vec <- function(obs_vec, v_call_vec, j_call_vec, v_ref_tbl, j_ref_tbl, prefix) {
  x <- mapply(
    function(obs, v, j) {
      complete_chain(
        obs_aa = obs,
        v_call = v,
        j_call = j,
        v_ref_tbl = v_ref_tbl,
        j_ref_tbl = j_ref_tbl
      )
    },
    obs_vec,
    v_call_vec,
    j_call_vec,
    SIMPLIFY = FALSE
  )

  completion_frame(x, prefix)
}

complete_light_chain_vec <- function(obs_vec, v_call_vec, j_call_vec, locus_vec) {
  x <- mapply(
    function(obs, v, j, locus) {
      locus <- toupper(as.character(locus))

      if (locus == "IGK") {
        complete_chain(obs, v, j, IGKV_ref, IGKJ_ref)
      } else if (locus == "IGL") {
        complete_chain(obs, v, j, IGLV_ref, IGLJ_ref)
      } else {
        completion_error(obs, "invalid_light_locus_check_manually")
      }
    },
    obs_vec,
    v_call_vec,
    j_call_vec,
    locus_vec,
    SIMPLIFY = FALSE
  )

  completion_frame(x, "light_")
}


# ============================================================
# 9. process one donor
# ============================================================

process_one_sample <- function(sample_name, file) {

  log_msg("Processing sample '", sample_name, "'  (", basename(file), ")")

  df_raw <- read_tsv(file, show_col_types = FALSE)
  present_cols <- colnames(df_raw)

  needed_cols <- c(
    "productive", "stop_codon", "vj_in_frame",
    "sequence_id", "locus",
    "v_call", "d_call", "j_call", "c_call",
    "cdr3_aa", "junction_aa",
    "cdr3", "junction",
    "v_identity",
    "sequence_alignment", "sequence_alignment_aa",
    "v_cigar", "d_cigar", "j_cigar",
    "fwr1_aa", "fwr2_aa", "fwr3_aa", "fwr4_aa",
    "complete_vdj"
  )

  df <- df_raw

  for (cc in needed_cols) {
    if (!cc %in% colnames(df)) {
      df[[cc]] <- NA
    }
  }

  if ("productive" %in% present_cols && !all(is.na(df$productive))) {
    df <- df %>% filter(flag_true(productive))
  }

  if ("stop_codon" %in% present_cols && !all(is.na(df$stop_codon))) {
    df <- df %>% filter(flag_false(stop_codon))
  }

  if ("vj_in_frame" %in% present_cols && !all(is.na(df$vj_in_frame))) {
    df <- df %>% filter(flag_true(vj_in_frame))
  }

  df <- df %>%
    filter(!is.na(sequence_id), !is.na(v_call), !is.na(j_call)) %>%
    mutate(
      sample = sample_name,
      locus = toupper(as.character(locus)),
      cell_id = extract_cell_id(sequence_id),
      v_gene = clean_gene(v_call),
      j_gene = clean_gene(j_call),
      vdj_aa = clean_aa(sequence_alignment_aa),
      cdr3_aa_source = case_when(
        !is_blank(cdr3_aa) ~ "cdr3_aa",
        !is_blank(junction_aa) ~ "junction_aa",
        TRUE ~ "missing"
      ),
      cdr3_aa_use = coalesce_text(cdr3_aa, junction_aa),
      cdr3_aa_clean = clean_aa(cdr3_aa_use),
      junction_nt_clean = clean_nt(junction),
      cdr3_nt_clean = clean_nt(cdr3),
      junction_frame_ok = nt_frame_len_ok(junction),
      cdr3_frame_ok = nt_frame_len_ok(cdr3),
      junction_translation_ok = junction_translation_ok(junction, junction_aa),
      v_cigar_frameshift = cigar_nontriplet_indel(v_cigar),
      d_cigar_frameshift = cigar_nontriplet_indel(d_cigar),
      j_cigar_frameshift = cigar_nontriplet_indel(j_cigar),
      framework_profile = framework_profile(fwr1_aa, fwr2_aa, fwr3_aa),
      raw_alignment_aa_has_stop = raw_aa_has_stop(sequence_alignment_aa)
    ) %>%
    filter(valid_aa(cdr3_aa_clean)) %>%
    filter(valid_aa(vdj_aa))

  # -----------------------------
  # heavy IgG
  # -----------------------------

  heavy_all <- df %>%
    filter(locus == "IGH", str_detect(v_gene, "^IGHV"))

  igg_info <- classify_igg(heavy_all$c_call, heavy_all$sequence_id)
  heavy_all <- bind_cols(heavy_all, igg_info)

  if (strict_igg_subclass_only) {
    heavy_igg <- heavy_all %>%
      filter(heavy_subclass %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4"))
  } else {
    heavy_igg <- heavy_all %>%
      filter(igg_source != "not_IgG")
  }

  heavy_excluded_inferred <- heavy_all %>%
    filter(igg_source != "not_IgG") %>%
    filter(!(heavy_subclass %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4"))) %>%
    transmute(
      sample,
      cell_id,
      heavy_sequence_id = sequence_id,
      heavy_v_call = top_call(v_call),
      heavy_j_call = top_call(j_call),
      heavy_c_call = top_call(c_call),
      heavy_subclass,
      igg_source,
      heavy_cdr3_aa = cdr3_aa_clean,
      heavy_cdr3_aa_source = cdr3_aa_source,
      heavy_vdj_aa_raw = vdj_aa
    )

  heavy_igg <- heavy_igg %>%
    transmute(
      sample,
      cell_id,
      heavy_sequence_id = sequence_id,
      heavy_v_call = top_call(v_call),
      heavy_d_call = top_call(d_call),
      heavy_j_call = top_call(j_call),
      heavy_c_call = top_call(c_call),
      heavy_subclass,
      igg_source,
      heavy_v_gene = v_gene,
      heavy_j_gene = j_gene,
      heavy_cdr3_aa = cdr3_aa_clean,
      heavy_cdr3_aa_source = cdr3_aa_source,
      heavy_v_identity = v_identity,
      heavy_vdj_aa_raw = vdj_aa,
      heavy_junction_nt = junction_nt_clean,
      heavy_cdr3_nt = cdr3_nt_clean,
      heavy_junction_aa_raw = clean_aa(junction_aa),
      heavy_junction_frame_ok = junction_frame_ok,
      heavy_cdr3_frame_ok = cdr3_frame_ok,
      heavy_junction_translation_ok = junction_translation_ok,
      heavy_v_cigar_frameshift = v_cigar_frameshift,
      heavy_d_cigar_frameshift = d_cigar_frameshift,
      heavy_j_cigar_frameshift = j_cigar_frameshift,
      heavy_framework_profile = framework_profile,
      heavy_raw_alignment_aa_has_stop = raw_alignment_aa_has_stop,
      heavy_complete_vdj = complete_vdj,
      heavy_fwr1_aa = clean_aa(fwr1_aa),
      heavy_fwr2_aa = clean_aa(fwr2_aa),
      heavy_fwr3_aa = clean_aa(fwr3_aa)
    )

  # -----------------------------
  # light
  # -----------------------------

  light_productive <- df %>%
    filter(locus %in% c("IGK", "IGL")) %>%
    filter(str_detect(v_gene, "^(IGKV|IGLV)")) %>%
    transmute(
      sample,
      cell_id,
      light_sequence_id = sequence_id,
      light_locus = locus,
      light_v_call = top_call(v_call),
      light_j_call = top_call(j_call),
      light_c_call = top_call(c_call),
      light_v_gene = v_gene,
      light_j_gene = j_gene,
      light_cdr3_aa = cdr3_aa_clean,
      light_cdr3_aa_source = cdr3_aa_source,
      light_v_identity = v_identity,
      light_vj_aa_raw = vdj_aa,
      light_junction_nt = junction_nt_clean,
      light_cdr3_nt = cdr3_nt_clean,
      light_junction_aa_raw = clean_aa(junction_aa),
      light_junction_frame_ok = junction_frame_ok,
      light_cdr3_frame_ok = cdr3_frame_ok,
      light_junction_translation_ok = junction_translation_ok,
      light_v_cigar_frameshift = v_cigar_frameshift,
      light_j_cigar_frameshift = j_cigar_frameshift,
      light_framework_profile = framework_profile,
      light_raw_alignment_aa_has_stop = raw_alignment_aa_has_stop,
      light_complete_vdj = complete_vdj,
      light_fwr1_aa = clean_aa(fwr1_aa),
      light_fwr2_aa = clean_aa(fwr2_aa),
      light_fwr3_aa = clean_aa(fwr3_aa)
    )

  # -----------------------------
  # pairing
  # -----------------------------

  paired <- inner_join(
    heavy_igg,
    light_productive,
    by = c("sample", "cell_id")
  ) %>%
    arrange(sample, cell_id)

  empty_summary <- function(paired_n = 0) {
    data.frame(
      sample = sample_name,
      raw_rows = nrow(df_raw),
      after_basic_filter = nrow(df),
      strict_IgG_heavy = nrow(heavy_igg),
      inferred_or_unsubclassed_IgG_heavy_excluded = nrow(heavy_excluded_inferred),
      productive_light = nrow(light_productive),
      paired_strict_IgG_BCR = paired_n,
      A_priority = 0,
      B_acceptable = 0,
      manual_check = 0,
      whole_variable_frame_or_domain_risk = 0,
      multi_pair_cells = 0
    )
  }

  if (nrow(paired) == 0) {
    return(list(
      paired = paired,
      ready = paired,
      manual = paired,
      excluded = heavy_excluded_inferred,
      summary = empty_summary(0)
    ))
  }

  # -----------------------------
  # N-terminal and J-terminal completion
  # -----------------------------

  heavy_complete <- complete_chain_vec(
    obs_vec = paired$heavy_vdj_aa_raw,
    v_call_vec = paired$heavy_v_call,
    j_call_vec = paired$heavy_j_call,
    v_ref_tbl = IGHV_ref,
    j_ref_tbl = IGHJ_ref,
    prefix = "heavy_"
  )

  light_complete <- complete_light_chain_vec(
    obs_vec = paired$light_vj_aa_raw,
    v_call_vec = paired$light_v_call,
    j_call_vec = paired$light_j_call,
    locus_vec = paired$light_locus
  )

  heavy_whole_v <- compare_v_prefix_vec(
    heavy_complete$heavy_final_aa,
    heavy_complete$heavy_v_ref_id,
    IGHV_ref,
    prefix = "heavy_"
  )

  light_whole_v_list <- mapply(
    function(a, r, locus) {
      if (toupper(as.character(locus)) == "IGK") {
        compare_v_prefix_one(a, r, IGKV_ref, whole_v_prefix_len)
      } else if (toupper(as.character(locus)) == "IGL") {
        compare_v_prefix_one(a, r, IGLV_ref, whole_v_prefix_len)
      } else {
        list(identity = NA_real_, mismatch_run = NA_integer_, compared_len = NA_integer_)
      }
    },
    light_complete$light_final_aa,
    light_complete$light_v_ref_id,
    paired$light_locus,
    SIMPLIFY = FALSE
  )

  light_whole_v <- data.frame(
    light_whole_v_identity = sapply(light_whole_v_list, function(x) x[["identity"]]),
    light_whole_v_max_mismatch_run = sapply(light_whole_v_list, function(x) x[["mismatch_run"]]),
    light_whole_v_compared_len = sapply(light_whole_v_list, function(x) x[["compared_len"]]),
    stringsAsFactors = FALSE
  )

  paired <- bind_cols(
    paired,
    heavy_complete,
    light_complete,
    heavy_whole_v,
    light_whole_v
  ) %>%
    group_by(sample, cell_id) %>%
    mutate(pair_count = n()) %>%
    ungroup() %>%
    mutate(
      final_VH_aa = heavy_final_aa,
      final_VL_aa = light_final_aa,
      VH_len = nchar(final_VH_aa),
      VL_len = nchar(final_VL_aa),

      # extra-Cys must distinguish AIRR cdr3_aa from junction_aa: cdr3_aa
      # already excludes the 5' conserved Cys, junction_aa does not.
      heavy_extra_cys_n = count_extra_cdr3_cys(
        heavy_cdr3_aa,
        heavy_cdr3_aa_source
      ),
      light_extra_cys_n = count_extra_cdr3_cys(
        light_cdr3_aa,
        light_cdr3_aa_source
      ),
      extra_cys_total = heavy_extra_cys_n + light_extra_cys_n,

      cys_risk = case_when(
        extra_cys_total == 0 ~ "no_extra_Cys",
        extra_cys_total == 1 ~ "one_extra_Cys_warning",
        extra_cys_total == 2 ~ "two_extra_Cys_warning",
        TRUE ~ "multiple_extra_Cys_manual_check"
      ),

      # Presence of the two conserved cysteines, which the extra-cysteine
      # counting above cannot see. FR1's cysteine sits at 22 in the heavy chain
      # and at 22 or 23 in the light chain, and a completed sequence can start a
      # residue either side of it, so a window is searched rather than one
      # position. Without this, a chain broken between FR1 and FR2 keeps a
      # plausible length, an intact FR4 and a 0.9 whole-V identity, and reaches
      # expression unable to fold.
      heavy_n_cys = str_count(final_VH_aa, "C"),
      light_n_cys = str_count(final_VL_aa, "C"),
      heavy_fr1_cys_ok = str_detect(substr(final_VH_aa, 19, 26), "C"),
      light_fr1_cys_ok = str_detect(substr(final_VL_aa, 19, 26), "C"),
      heavy_junction_cys_ok = junction_cys_ok(heavy_junction_aa_raw),
      light_junction_cys_ok = junction_cys_ok(light_junction_aa_raw),

      cys_pair_ok =
        heavy_n_cys >= 2 & light_n_cys >= 2 &
        heavy_fr1_cys_ok & light_fr1_cys_ok &
        (is.na(heavy_junction_cys_ok) | heavy_junction_cys_ok) &
        (is.na(light_junction_cys_ok) | light_junction_cys_ok),

      VH_FR4_tail = last_n(final_VH_aa, 16),
      VL_FR4_tail = last_n(final_VL_aa, 16),

      VH_len_ok = !is.na(VH_len) & VH_len >= CFG$VH_len_min & VH_len <= CFG$VH_len_max,
      VL_len_ok = !is.na(VL_len) & VL_len >= CFG$VL_len_min & VL_len <= CFG$VL_len_max,

      VH_FR4_typical = str_detect(VH_FR4_tail, "WG[QR]G"),
      VL_FR4_typical = str_detect(VL_FR4_tail, "FGQG|FGGG|FGTG|FGPG|FGEG"),

      heavy_anchor_080 = !is.na(heavy_v_anchor_identity) & heavy_v_anchor_identity >= 0.80,
      light_anchor_080 = !is.na(light_v_anchor_identity) & light_v_anchor_identity >= 0.80,
      heavy_anchor_clean = !is.na(heavy_v_anchor_identity) & heavy_v_anchor_identity >= CFG$min_anchor_clean,
      light_anchor_clean = !is.na(light_v_anchor_identity) & light_v_anchor_identity >= CFG$min_anchor_clean,
      heavy_anchor_hard = !is.na(heavy_v_anchor_identity) & heavy_v_anchor_identity >= CFG$min_anchor_hard,
      light_anchor_hard = !is.na(light_v_anchor_identity) & light_v_anchor_identity >= CFG$min_anchor_hard,

      heavy_j_080 = !is.na(heavy_j_overlap_identity) & heavy_j_overlap_identity >= 0.80,
      light_j_080 = !is.na(light_j_overlap_identity) & light_j_overlap_identity >= 0.80,
      heavy_j_clean = !is.na(heavy_j_overlap_identity) & heavy_j_overlap_identity >= CFG$min_j_overlap_clean,
      light_j_clean = !is.na(light_j_overlap_identity) & light_j_overlap_identity >= CFG$min_j_overlap_clean,
      heavy_j_hard = !is.na(heavy_j_overlap_identity) & heavy_j_overlap_identity >= CFG$min_j_overlap_hard,
      light_j_hard = !is.na(light_j_overlap_identity) & light_j_overlap_identity >= CFG$min_j_overlap_hard,

      heavy_N6 = !is.na(heavy_added_N_len) & heavy_added_N_len <= 6,
      light_N6 = !is.na(light_added_N_len) & light_added_N_len <= 6,
      heavy_N10 = !is.na(heavy_added_N_len) & heavy_added_N_len <= CFG$max_added_N_clean,
      light_N10 = !is.na(light_added_N_len) & light_added_N_len <= CFG$max_added_N_clean,
      heavy_N_hard = !is.na(heavy_added_N_len) & heavy_added_N_len <= CFG$max_added_N_hard,
      light_N_hard = !is.na(light_added_N_len) & light_added_N_len <= CFG$max_added_N_hard,

      heavy_C3 = !is.na(heavy_added_C_len) & heavy_added_C_len <= 3,
      light_C3 = !is.na(light_added_C_len) & light_added_C_len <= 3,
      heavy_C6 = !is.na(heavy_added_C_len) & heavy_added_C_len <= CFG$max_added_C_clean,
      light_C6 = !is.na(light_added_C_len) & light_added_C_len <= CFG$max_added_C_clean,
      heavy_C_hard = !is.na(heavy_added_C_len) & heavy_added_C_len <= CFG$max_added_C_hard,
      light_C_hard = !is.na(light_added_C_len) & light_added_C_len <= CFG$max_added_C_hard,

      # -----------------------------
      # whole variable domain frame / structure QC
      # -----------------------------
      heavy_explicit_frame_risk =
        (!is.na(heavy_junction_frame_ok) & !heavy_junction_frame_ok) |
        (!is.na(heavy_cdr3_frame_ok) & !heavy_cdr3_frame_ok) |
        (!is.na(heavy_junction_translation_ok) & !heavy_junction_translation_ok) |
        (!is.na(heavy_v_cigar_frameshift) & heavy_v_cigar_frameshift) |
        (!is.na(heavy_d_cigar_frameshift) & heavy_d_cigar_frameshift) |
        (!is.na(heavy_j_cigar_frameshift) & heavy_j_cigar_frameshift) |
        heavy_raw_alignment_aa_has_stop,

      light_explicit_frame_risk =
        (!is.na(light_junction_frame_ok) & !light_junction_frame_ok) |
        (!is.na(light_cdr3_frame_ok) & !light_cdr3_frame_ok) |
        (!is.na(light_junction_translation_ok) & !light_junction_translation_ok) |
        (!is.na(light_v_cigar_frameshift) & light_v_cigar_frameshift) |
        (!is.na(light_j_cigar_frameshift) & light_j_cigar_frameshift) |
        light_raw_alignment_aa_has_stop,

      heavy_v_profile_bad =
        (!is.na(heavy_whole_v_identity) & heavy_whole_v_identity < min_whole_v_identity_hard) |
        (!is.na(heavy_whole_v_max_mismatch_run) &
           heavy_whole_v_max_mismatch_run > max_whole_v_mismatch_run_hard),

      light_v_profile_bad =
        (!is.na(light_whole_v_identity) & light_whole_v_identity < min_whole_v_identity_hard) |
        (!is.na(light_whole_v_max_mismatch_run) &
           light_whole_v_max_mismatch_run > max_whole_v_mismatch_run_hard),

      heavy_framework_bad = heavy_framework_profile == "fail",
      light_framework_bad = light_framework_profile == "fail",

      heavy_whole_variable_bad =
        heavy_explicit_frame_risk | heavy_v_profile_bad | heavy_framework_bad,
      light_whole_variable_bad =
        light_explicit_frame_risk | light_v_profile_bad | light_framework_bad,

      heavy_whole_variable_qc = case_when(
        heavy_explicit_frame_risk ~ "FAIL_explicit_frame_or_translation_risk",
        heavy_v_profile_bad ~ "FAIL_abnormal_whole_V_profile",
        heavy_framework_bad ~ "FAIL_framework_profile",
        TRUE ~ "PASS_no_negative_evidence"
      ),

      light_whole_variable_qc = case_when(
        light_explicit_frame_risk ~ "FAIL_explicit_frame_or_translation_risk",
        light_v_profile_bad ~ "FAIL_abnormal_whole_V_profile",
        light_framework_bad ~ "FAIL_framework_profile",
        TRUE ~ "PASS_no_negative_evidence"
      ),

      note_all = paste(
        heavy_N_note,
        heavy_J_note,
        light_N_note,
        light_J_note,
        sep = ";"
      ),

      note_bad = str_detect(
        note_all,
        "check_manually|invalid|empty|no_reference|cannot_align|too_much|low_"
      )
    ) %>%
    mutate(
      priority = case_when(
        pair_count > 1 ~ "D_multi_pair_manual_check",
        heavy_explicit_frame_risk | light_explicit_frame_risk ~
          "D_whole_variable_frame_risk_manual_check",
        heavy_v_profile_bad | light_v_profile_bad |
          heavy_framework_bad | light_framework_bad ~
          "D_abnormal_variable_domain_manual_check",
        note_bad ~ "D_sequence_completion_manual_check",
        !heavy_N_hard | !light_N_hard | !heavy_C_hard | !light_C_hard ~
          "D_large_terminal_gap_manual_check",
        !heavy_anchor_hard | !light_anchor_hard | !heavy_j_hard | !light_j_hard ~
          "D_low_alignment_manual_check",
        # A missing conserved cysteine is a whole-domain defect, so it is
        # graded D and sits above the C checks: length, FR4 tail and the
        # extra-cysteine count can all look acceptable on a chain that can no
        # longer fold.
        !cys_pair_ok ~ "D_conserved_cys_manual_check",
        !VH_len_ok | !VL_len_ok ~ "C_length_manual_check",
        !VH_FR4_typical | !VL_FR4_typical ~ "C_FR4_tail_manual_check",
        cys_risk == "multiple_extra_Cys_manual_check" ~
          "C_multiple_extra_Cys_manual_check",

        heavy_N6 & light_N6 &
          heavy_C3 & light_C3 &
          heavy_anchor_080 & light_anchor_080 &
          heavy_j_080 & light_j_080 &
          VH_len_ok & VL_len_ok &
          VH_FR4_typical & VL_FR4_typical &
          !heavy_whole_variable_bad & !light_whole_variable_bad ~ "A_priority",

        heavy_N10 & light_N10 &
          heavy_C6 & light_C6 &
          heavy_anchor_clean & light_anchor_clean &
          heavy_j_clean & light_j_clean &
          VH_len_ok & VL_len_ok &
          VH_FR4_typical & VL_FR4_typical &
          !heavy_whole_variable_bad & !light_whole_variable_bad ~ "B_acceptable",

        TRUE ~ "C_manual_check"
      )
    )

  # Retain the completed sequence for audit, then blank it out for anything
  # that is not expression-ready. Otherwise a truncated chain would still sit
  # in a column named final_*_aa, and any downstream reader that forgets to
  # filter on priority would silently pick it up.
  paired <- paired %>%
    mutate(
      final_VH_aa_completed = final_VH_aa,
      final_VL_aa_completed = final_VL_aa,
      final_VH_aa = if_else(priority %in% CFG$priority_ready, final_VH_aa, NA_character_),
      final_VL_aa = if_else(priority %in% CFG$priority_ready, final_VL_aa, NA_character_)
    ) %>%
    select(
      sample,
      cell_id,
      priority,
      pair_count,

      heavy_sequence_id,
      light_sequence_id,
      heavy_subclass,
      light_locus,

      heavy_v_call,
      heavy_d_call,
      heavy_j_call,
      heavy_c_call,
      light_v_call,
      light_j_call,
      light_c_call,

      heavy_cdr3_aa,
      light_cdr3_aa,
      heavy_cdr3_aa_source,
      light_cdr3_aa_source,
      heavy_v_identity,
      light_v_identity,

      heavy_whole_variable_qc,
      light_whole_variable_qc,
      heavy_whole_v_identity,
      light_whole_v_identity,
      heavy_whole_v_max_mismatch_run,
      light_whole_v_max_mismatch_run,
      heavy_whole_v_compared_len,
      light_whole_v_compared_len,
      heavy_junction_frame_ok,
      light_junction_frame_ok,
      heavy_cdr3_frame_ok,
      light_cdr3_frame_ok,
      heavy_junction_translation_ok,
      light_junction_translation_ok,
      heavy_v_cigar_frameshift,
      light_v_cigar_frameshift,
      heavy_d_cigar_frameshift,
      heavy_j_cigar_frameshift,
      light_j_cigar_frameshift,
      heavy_framework_profile,
      light_framework_profile,

      final_VH_aa,
      final_VL_aa,
      final_VH_aa_completed,
      final_VL_aa_completed,
      VH_len,
      VL_len,

      heavy_added_N_aa,
      light_added_N_aa,
      heavy_added_N_len,
      light_added_N_len,
      heavy_added_C_aa,
      light_added_C_aa,
      heavy_added_C_len,
      light_added_C_len,

      heavy_v_ref_id,
      light_v_ref_id,
      heavy_j_ref_id,
      light_j_ref_id,
      heavy_v_ref_match_type,
      light_v_ref_match_type,
      heavy_j_ref_match_type,
      light_j_ref_match_type,

      heavy_v_anchor_identity,
      light_v_anchor_identity,
      heavy_j_overlap_identity,
      light_j_overlap_identity,
      heavy_j_overlap_len,
      light_j_overlap_len,

      heavy_N_note,
      heavy_J_note,
      light_N_note,
      light_J_note,

      VH_FR4_tail,
      VL_FR4_tail,
      VH_FR4_typical,
      VL_FR4_typical,

      heavy_extra_cys_n,
      light_extra_cys_n,
      cys_risk,
      heavy_n_cys,
      light_n_cys,
      heavy_fr1_cys_ok,
      light_fr1_cys_ok,
      heavy_junction_cys_ok,
      light_junction_cys_ok,
      cys_pair_ok,

      heavy_junction_nt,
      light_junction_nt,
      heavy_cdr3_nt,
      light_cdr3_nt,
      heavy_junction_aa_raw,
      light_junction_aa_raw,
      heavy_fwr1_aa,
      heavy_fwr2_aa,
      heavy_fwr3_aa,
      light_fwr1_aa,
      light_fwr2_aa,
      light_fwr3_aa,
      heavy_vdj_aa_raw,
      light_vj_aa_raw
    )

  ready <- paired %>%
    filter(priority %in% CFG$priority_ready) %>%
    arrange(
      priority,
      heavy_added_N_len + light_added_N_len + heavy_added_C_len + light_added_C_len,
      desc(heavy_v_anchor_identity + light_v_anchor_identity +
             heavy_j_overlap_identity + light_j_overlap_identity)
    )

  manual <- paired %>%
    filter(!priority %in% CFG$priority_ready) %>%
    arrange(
      priority,
      heavy_added_N_len + light_added_N_len + heavy_added_C_len + light_added_C_len
    )

  summary <- data.frame(
    sample = sample_name,
    raw_rows = nrow(df_raw),
    after_basic_filter = nrow(df),
    strict_IgG_heavy = nrow(heavy_igg),
    inferred_or_unsubclassed_IgG_heavy_excluded = nrow(heavy_excluded_inferred),
    productive_light = nrow(light_productive),
    paired_strict_IgG_BCR = nrow(paired),
    A_priority = sum(paired$priority == "A_priority"),
    B_acceptable = sum(paired$priority == "B_acceptable"),
    manual_check = sum(!paired$priority %in% CFG$priority_ready),
    whole_variable_frame_or_domain_risk = sum(
      paired$priority %in% CFG$priority_whole_variable_risk
    ),
    multi_pair_cells = paired %>%
      count(sample, cell_id, name = "n_pair") %>%
      filter(n_pair > 1) %>%
      nrow()
  )

  list(
    paired = paired,
    ready = ready,
    manual = manual,
    excluded = heavy_excluded_inferred,
    summary = summary
  )
}


# ============================================================
# 10. run every donor, then export
# ============================================================

run_started <- Sys.time()

res <- lapply(seq_len(nrow(disc$samples)), function(i) {
  process_one_sample(
    sample_name = disc$samples$sample[i],
    file        = disc$samples$path[i]
  )
})
names(res) <- disc$samples$sample

summary_table <- bind_rows(lapply(res, `[[`, "summary"))

# Provenance for the run: which reference release produced these sequences.
ref_sheet <- ref_provenance
ref_sheet$downloaded_at <- format(run_started, "%Y-%m-%d %H:%M:%S")

summary_sheet <- bind_rows(
  summary_table,
  data.frame(
    sample = "TOTAL",
    raw_rows = sum(summary_table$raw_rows),
    after_basic_filter = sum(summary_table$after_basic_filter),
    strict_IgG_heavy = sum(summary_table$strict_IgG_heavy),
    inferred_or_unsubclassed_IgG_heavy_excluded =
      sum(summary_table$inferred_or_unsubclassed_IgG_heavy_excluded),
    productive_light = sum(summary_table$productive_light),
    paired_strict_IgG_BCR = sum(summary_table$paired_strict_IgG_BCR),
    A_priority = sum(summary_table$A_priority),
    B_acceptable = sum(summary_table$B_acceptable),
    manual_check = sum(summary_table$manual_check),
    whole_variable_frame_or_domain_risk =
      sum(summary_table$whole_variable_frame_or_domain_risk),
    multi_pair_cells = sum(summary_table$multi_pair_cells)
  )
)

sample_sheets <- function(r) {
  # Which priority labels this run actually used. Step 2 reads this sheet rather
  # than carrying its own copy of the strings, so the two scripts cannot drift
  # apart if a label is ever renamed.
  priority_labels <- data.frame(
    group = c(
      rep("ready", length(CFG$priority_ready)),
      rep("whole_variable_risk", length(CFG$priority_whole_variable_risk)),
      rep("manual", length(unique(r$paired$priority[
        !r$paired$priority %in% c(CFG$priority_ready,
                                  CFG$priority_whole_variable_risk)
      ])))
    ),
    priority = c(
      CFG$priority_ready,
      CFG$priority_whole_variable_risk,
      unique(r$paired$priority[
        !r$paired$priority %in% c(CFG$priority_ready,
                                  CFG$priority_whole_variable_risk)
      ])
    ),
    stringsAsFactors = FALSE
  )

  whole_variable_risk <- r$paired %>%
    filter(priority %in% CFG$priority_whole_variable_risk)

  list(
    expression_ready_AB = r$ready,
    all_paired_strict_IgG = r$paired,
    manual_check = r$manual,
    whole_variable_risk = whole_variable_risk,
    excluded_inferred_IgG = r$excluded,
    priority_labels = priority_labels
  )
}

written <- character(0)

for (s in names(res)) {
  path <- file.path(CFG$output_dir, paste0(s, "_IgG_BCR_expression_ready.xlsx"))
  write_xlsx(
    c(
      sample_sheets(res[[s]]),
      list(
        summary = summary_table[summary_table$sample == s, , drop = FALSE],
        reference_provenance = ref_sheet
      )
    ),
    path
  )
  written <- c(written, path)
}

# A combined workbook keeps the familiar single-file workflow when several
# donors are processed in one run.
all_paired <- bind_rows(lapply(res, function(r) r$paired))

all_priority_labels <- data.frame(
  group = c(
    rep("ready", length(CFG$priority_ready)),
    rep("whole_variable_risk", length(CFG$priority_whole_variable_risk)),
    rep("manual", length(unique(all_paired$priority[
      !all_paired$priority %in% c(CFG$priority_ready,
                                  CFG$priority_whole_variable_risk)
    ])))
  ),
  priority = c(
    CFG$priority_ready,
    CFG$priority_whole_variable_risk,
    unique(all_paired$priority[
      !all_paired$priority %in% c(CFG$priority_ready,
                                  CFG$priority_whole_variable_risk)
    ])
  ),
  stringsAsFactors = FALSE
)

combined_path <- file.path(CFG$output_dir, "combined_IgG_BCR_expression_ready.xlsx")
write_xlsx(
  list(
    expression_ready_AB = bind_rows(lapply(res, function(r) r$ready)),
    all_paired_strict_IgG = all_paired,
    manual_check = bind_rows(lapply(res, function(r) r$manual)),
    whole_variable_risk = all_paired %>%
      filter(priority %in% CFG$priority_whole_variable_risk),
    excluded_inferred_IgG = bind_rows(lapply(res, function(r) r$excluded)),
    summary = summary_sheet,
    reference_provenance = ref_sheet,
    priority_labels = all_priority_labels
  ),
  combined_path
)
written <- c(written, combined_path)

log_msg("Done.")
print(summary_sheet)
cat("\n")
log_msg("Files written:")
for (p in written) log_msg("  ", p)
log_msg("  ", disc_report_path)

cat("\n")
cat("Step 1 finished.\n")
cat("  output/<sample>_IgG_BCR_expression_ready.xlsx : per-donor result\n")
cat("  output/combined_IgG_BCR_expression_ready.xlsx : all donors merged\n")
cat("  sheet expression_ready_AB        : use these for expression (A/B only)\n")
cat("  sheet all_paired_strict_IgG      : complete result, includes C/D\n")
cat("  sheet manual_check               : needs a human decision\n")
cat("  sheet whole_variable_risk        : frame / domain anomalies, review first\n")
cat("  sheet excluded_inferred_IgG      : IgG not subclass-called, excluded\n")
cat("  sheet reference_provenance       : reference release and record counts\n")
cat("  output/discovery_report.tsv      : which files were used, which skipped\n")
cat("\n")
cat("Note: final_VH_aa / final_VL_aa are NA for anything that is not A/B.\n")
cat("      The completed sequence is kept in final_VH_aa_completed /\n")
cat("      final_VL_aa_completed for auditing.\n")
cat("\nNext: Rscript R/02_extract_simplified.R\n")
