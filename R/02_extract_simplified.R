# =============================================================================
# 02_extract_simplified.R
#
# Step 2 of IgG BCR expression-ready selection. A thin companion to
# 01_pair_and_complete.R: it reads step 1's workbook and turns it into the
# tables a human reads and a company receives. It does no analysis of its own.
#
# Input : output/combined_IgG_BCR_expression_ready.xlsx   (produced by step 1)
#         selected_cells.xlsx                             (optional)
# Output: output/IgG_BCR_simplified_for_expression.xlsx
#
# Only rows that are genuinely ready (A/B, single pairing, matched) reach the
# table handed out for expression.
#
# Self-contained: like step 1, its settings live at the top of this file. The
# two scripts share no file, so neither is a hidden dependency of the other.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(writexl)
})


# =============================================================================
# SETTINGS
# =============================================================================

# Project root: this file sits in <root>/R/.
project_root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) {
    dirname(dirname(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE)))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})

data_dir   <- file.path(project_root, "data")
output_dir <- file.path(project_root, "output")

# Optional request list: two columns, sample and cell_id. Absent -> every paired
# BCR is summarised instead of a requested subset.
selected_file <- Sys.getenv("BCR_SELECTED_CELLS",
                            unset = file.path(project_root, "selected_cells.xlsx"))

output_xlsx <- file.path(output_dir, "IgG_BCR_simplified_for_expression.xlsx")

# Fallback names for the two grading groups, used only if step 1's workbook
# predates the `priority_labels` sheet. Step 1 is the authority: see
# read_priority_labels() below.
ready_priorities_fallback <- c("A_priority", "B_acceptable")
whole_variable_risk_fallback <- c(
  "D_whole_variable_frame_risk_manual_check",
  "D_abnormal_variable_domain_manual_check"
)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Cannot create directory: ", path, call. = FALSE)
  invisible(path)
}


# ============================================================
# 1. basic helpers
# ============================================================

normalize_sample <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x[is.na(x)] <- ""
  str_replace_all(x, "\\s+", "")
}

# Normalisation used only when matching an optional request list against the
# results. It is deliberately conservative: case, surrounding whitespace and
# dash variety are ignored, and a single-digit numeric suffix is zero-padded so
# that a request written as "plateA-G1" still matches a result stored as
# "PLATEA-G01".
normalize_cell_id <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x[is.na(x)] <- ""
  x <- str_replace_all(x, "\u2013|\u2014", "-")
  x <- str_replace_all(x, "\\s+", "")

  # Zero-pad a one-character letter group followed by a single digit, written
  # with paste0 rather than a backreference so the intent is unambiguous.
  m <- str_match(x, "^(.*-)([A-Z])(\\d)$")
  idx <- !is.na(m[, 1])
  if (any(idx)) {
    x[idx] <- paste0(m[idx, 2], m[idx, 3], "0", m[idx, 4])
  }

  x
}

short_light_locus <- function(x) {
  x <- toupper(as.character(x))
  case_when(
    x %in% c("IGK", "K", "KAPPA") ~ "K",
    x %in% c("IGL", "L", "LAMBDA") ~ "L",
    TRUE ~ x
  )
}

fmt_num <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

fmt_int <- function(x) {
  x <- suppressWarnings(as.integer(x))
  ifelse(is.na(x), "", as.character(x))
}

safe_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

# Locate the step 1 workbook. The combined file is preferred because it holds
# every donor; a single per-donor file is accepted as well. Several per-donor
# files with no combined file is ambiguous, so the run stops rather than
# silently analysing one donor.
find_analysis_xlsx <- function(dir = output_dir) {
  combined <- file.path(dir, "combined_IgG_BCR_expression_ready.xlsx")
  if (file.exists(combined)) return(combined)

  cand <- list.files(
    path = dir,
    pattern = "_IgG_BCR_expression_ready\\.xlsx$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  cand <- cand[!str_detect(basename(cand), regex("simplified|combined", ignore_case = TRUE))]
  cand <- unique(cand)

  if (length(cand) == 1) {
    warning("No combined workbook found; using the only per-donor file: ",
            basename(cand), call. = FALSE)
    return(cand)
  }

  if (length(cand) > 1) {
    stop(
      "Found ", length(cand), " per-donor workbooks and no combined workbook in ",
      output_dir, ":\n",
      paste0("  - ", basename(cand), collapse = "\n"), "\n",
      "  Re-run step 1 to produce combined_IgG_BCR_expression_ready.xlsx, or\n",
      "  point BCR_SELECTED_CELLS at one donor if that is what you meant.",
      call. = FALSE
    )
  }

  stop(
    "No step 1 result found in ", output_dir, "\n",
    "  Run step 1 first:  Rscript R/01_pair_and_complete.R",
    call. = FALSE
  )
}


pick_sheet <- function(xlsx) {
  sheets <- excel_sheets(xlsx)

  # Read the full table rather than only the ready rows, so that a manually
  # requested C/D cell surfaces its risk instead of disappearing silently.
  if ("all_paired_strict_IgG" %in% sheets) {
    return("all_paired_strict_IgG")
  }

  if ("expression_ready_AB" %in% sheets) {
    warning("all_paired_strict_IgG is absent; reading expression_ready_AB, ",
            "so C/D risk rows cannot be shown.", call. = FALSE)
    return("expression_ready_AB")
  }

  stop("Neither all_paired_strict_IgG nor expression_ready_AB is present.", call. = FALSE)
}


read_selected_cells <- function(file) {
  # 真正兼容“有表头/无表头”：统一按无表头读取，再判断第一行是否是标题
  x <- read_excel(file, sheet = 1, col_names = FALSE)

  if (ncol(x) < 2) {
    stop("selected_cells.xlsx 至少需要两列：第一列 sample，第二列 cell_id。")
  }

  x <- x[, 1:2]
  colnames(x) <- c("sample_input", "cell_id_input")

  if (nrow(x) == 0) {
    return(data.frame(
      sample_input = character(0),
      cell_id_input = character(0),
      row_order = integer(0),
      sample_norm = character(0),
      cell_id_norm = character(0),
      stringsAsFactors = FALSE
    ))
  }

  first_sample <- normalize_sample(x$sample_input[1])
  first_cell <- normalize_sample(x$cell_id_input[1])

  sample_header <- first_sample %in% c("SAMPLE", "SAMPLE_ID", "样本", "样本ID")
  cell_header <- first_cell %in% c("CELL_ID", "CELLID", "CELL", "细胞", "细胞ID")

  if (sample_header && cell_header) {
    x <- x[-1, , drop = FALSE]
  }

  x %>%
    mutate(
      sample_input = str_trim(as.character(sample_input)),
      cell_id_input = str_trim(as.character(cell_id_input))
    ) %>%
    filter(
      !(is.na(sample_input) | sample_input == ""),
      !(is.na(cell_id_input) | cell_id_input == "")
    ) %>%
    mutate(
      row_order = row_number(),
      sample_norm = normalize_sample(sample_input),
      cell_id_norm = normalize_cell_id(cell_id_input)
    )
}


# ============================================================
# 2. read the step 1 result
# ============================================================

input_xlsx <- find_analysis_xlsx()

message("Step 1 result: ", input_xlsx)

source_sheet <- pick_sheet(input_xlsx)
df <- read_excel(input_xlsx, sheet = source_sheet)

hard_required <- c("sample", "cell_id", "final_VH_aa", "final_VL_aa",
                   "light_locus", "priority")
missing_hard <- setdiff(hard_required, colnames(df))
if (length(missing_hard) > 0) {
  stop("Step 1 result is missing required column(s): ",
       paste(missing_hard, collapse = ", "), call. = FALSE)
}

# The audit copy of the completed sequence. Older step 1 output has only the
# priority-gated final_*_aa columns, which are NA for every C/D row and so
# cannot show what completion produced.
if (!"final_VH_aa_completed" %in% colnames(df)) {
  warning(
    "Step 1 result has no final_VH_aa_completed column; this looks like output ",
    "from an older version of the pipeline. Non-ready rows will show an empty ",
    "sequence in QC_details.", call. = FALSE
  )
  df$final_VH_aa_completed <- NA_character_
  df$final_VL_aa_completed <- NA_character_
}

# Columns used for the readable QC summary; absent ones are filled with NA so
# an older input still works.
need_cols <- c(
  "pair_count",
  "heavy_subclass",
  "heavy_v_call", "heavy_d_call", "heavy_j_call",
  "light_v_call", "light_j_call",
  "heavy_cdr3_aa", "light_cdr3_aa",
  "heavy_cdr3_aa_source", "light_cdr3_aa_source",
  "VH_len", "VL_len",
  "heavy_added_N_len", "light_added_N_len",
  "heavy_added_C_len", "light_added_C_len",
  "heavy_v_anchor_identity", "light_v_anchor_identity",
  "heavy_j_overlap_identity", "light_j_overlap_identity",
  "VH_FR4_tail", "VL_FR4_tail",
  "VH_FR4_typical", "VL_FR4_typical",
  "heavy_whole_variable_qc", "light_whole_variable_qc",
  "heavy_whole_v_identity", "light_whole_v_identity",
  "heavy_whole_v_max_mismatch_run", "light_whole_v_max_mismatch_run",
  "heavy_junction_frame_ok", "light_junction_frame_ok",
  "heavy_cdr3_frame_ok", "light_cdr3_frame_ok",
  "heavy_junction_translation_ok", "light_junction_translation_ok",
  "heavy_v_cigar_frameshift", "light_v_cigar_frameshift",
  "heavy_d_cigar_frameshift", "heavy_j_cigar_frameshift", "light_j_cigar_frameshift",
  "heavy_framework_profile", "light_framework_profile",
  "cys_risk"
)

for (cc in need_cols) {
  if (!cc %in% colnames(df)) df[[cc]] <- NA
}

# Step 1 records which priority labels its run used; read them from the workbook
# instead of keeping a second copy of the strings in this script.
read_priority_labels <- function(xlsx) {
  sheets <- excel_sheets(xlsx)
  if (!"priority_labels" %in% sheets) {
    message("Workbook has no priority_labels sheet (older step 1 output); ",
            "falling back to the built-in label names.")
    return(list(ready = ready_priorities_fallback,
                whole_variable_risk = whole_variable_risk_fallback))
  }
  lab <- read_excel(xlsx, sheet = "priority_labels")
  if (!all(c("group", "priority") %in% colnames(lab)) || nrow(lab) == 0) {
    message("priority_labels sheet is empty or malformed; ",
            "falling back to the built-in label names.")
    return(list(ready = ready_priorities_fallback,
                whole_variable_risk = whole_variable_risk_fallback))
  }
  ready <- lab$priority[lab$group == "ready"]
  risk <- lab$priority[lab$group == "whole_variable_risk"]
  list(
    ready = if (length(ready) > 0) ready else ready_priorities_fallback,
    whole_variable_risk = if (length(risk) > 0) risk else whole_variable_risk_fallback
  )
}

labels <- read_priority_labels(input_xlsx)
ready_priorities <- labels$ready
whole_variable_risk_priorities <- labels$whole_variable_risk

message("Grading labels in use: ready = ",
        paste(ready_priorities, collapse = ", "))


# ============================================================
# 3. build the key-result overview
# ============================================================

all_key_info <- df %>%
  mutate(
    sample_norm = normalize_sample(sample),
    cell_id_norm = normalize_cell_id(cell_id),
    light_locus = short_light_locus(light_locus),

    pair_count = suppressWarnings(as.integer(pair_count)),
    pair_count = ifelse(is.na(pair_count), 1L, pair_count),

    whole_variable_QC = case_when(
      str_starts(safe_text(heavy_whole_variable_qc), "FAIL") |
        str_starts(safe_text(light_whole_variable_qc), "FAIL") ~ "FAIL_check",
      str_starts(safe_text(heavy_whole_variable_qc), "PASS") &
        str_starts(safe_text(light_whole_variable_qc), "PASS") ~ "PASS_no_negative_evidence",
      TRUE ~ "not_available"
    ),

    expression_status = case_when(
      pair_count > 1 ~ "CHECK_multi_pair",
      priority %in% whole_variable_risk_priorities ~ "CHECK_whole_variable_risk",
      priority %in% ready_priorities ~ "READY",
      is.na(priority) | priority == "" ~ "CHECK_missing_priority",
      TRUE ~ "CHECK_manual"
    ),

    terminal_completion = paste0(
      "H:N+", fmt_int(heavy_added_N_len), "/C+", fmt_int(heavy_added_C_len),
      " | L:N+", fmt_int(light_added_N_len), "/C+", fmt_int(light_added_C_len)
    ),

    alignment_QC = paste0(
      "V H/L=", fmt_num(heavy_v_anchor_identity), "/", fmt_num(light_v_anchor_identity),
      " | J H/L=", fmt_num(heavy_j_overlap_identity), "/", fmt_num(light_j_overlap_identity)
    ),

    whole_V_QC = paste0(
      "H:", safe_text(heavy_whole_variable_qc),
      " | L:", safe_text(light_whole_variable_qc)
    )
  ) %>%
  transmute(
    sample,
    cell_id,
    expression_status,
    priority,
    pair_count,
    heavy_subclass,
    light_locus,

    final_VH_aa,
    final_VL_aa,
    # Carried through from step 1 so that a non-ready pair, whose final_*_aa is
    # NA by design, still has an auditable completed sequence here.
    final_VH_aa_completed,
    final_VL_aa_completed,
    VH_len,
    VL_len,

    heavy_v_call,
    heavy_d_call,
    heavy_j_call,
    light_v_call,
    light_j_call,

    heavy_cdr3_aa,
    light_cdr3_aa,
    heavy_cdr3_aa_source,
    light_cdr3_aa_source,

    whole_variable_QC,
    whole_V_QC,
    heavy_whole_v_identity,
    light_whole_v_identity,
    heavy_whole_v_max_mismatch_run,
    light_whole_v_max_mismatch_run,

    terminal_completion,
    alignment_QC,
    VH_FR4_tail,
    VL_FR4_tail,
    VH_FR4_typical,
    VL_FR4_typical,
    cys_risk,

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

    sample_norm,
    cell_id_norm
  ) %>%
  arrange(sample, cell_id, priority)


# ============================================================
# 4. 根据 selected_cells.xlsx 提取；没有该文件则直接展示全部
# ============================================================

if (file.exists(selected_file)) {
  selected <- read_selected_cells(selected_file)

  # 为防止同一 cell 多配对导致 left_join 静默扩行，显式记录每个选择项匹配了几行
  selected_joined <- selected %>%
    left_join(
      all_key_info %>%
        rename(
          sample_result = sample,
          cell_id_result = cell_id
        ),
      by = c("sample_norm", "cell_id_norm")
    ) %>%
    group_by(row_order) %>%
    mutate(
      matched_rows = sum(!is.na(sample_result)),
      match_status = case_when(
        matched_rows == 0 ~ "not_found",
        matched_rows > 1 ~ "multiple_rows_check",
        is.na(final_VH_aa) | is.na(final_VL_aa) ~ "partial_match_check",
        TRUE ~ "matched"
      )
    ) %>%
    ungroup() %>%
    mutate(
      sample = sample_input,
      cell_id = cell_id_input,
      review_status = case_when(
        match_status == "not_found" ~ "NOT_FOUND",
        match_status == "partial_match_check" ~ "CHECK_partial_match",
        match_status == "multiple_rows_check" ~ "CHECK_multi_pair_or_duplicate",
        expression_status == "READY" ~ "READY",
        TRUE ~ expression_status
      )
    ) %>%
    arrange(row_order, priority)

  simplified_with_QC <- selected_joined %>%
    select(
      sample,
      cell_id,
      review_status,
      match_status,
      matched_rows,
      expression_status,
      priority,
      pair_count,
      heavy_subclass,
      light_locus,
      final_VH_aa,
      final_VL_aa,
      final_VH_aa_completed,
      final_VL_aa_completed,
      VH_len,
      VL_len,
      heavy_v_call,
      heavy_d_call,
      heavy_j_call,
      light_v_call,
      light_j_call,
      heavy_cdr3_aa,
      light_cdr3_aa,
      heavy_cdr3_aa_source,
      light_cdr3_aa_source,
      whole_variable_QC,
      whole_V_QC,
      terminal_completion,
      alignment_QC,
      VH_FR4_tail,
      VL_FR4_tail,
      cys_risk,
      heavy_whole_v_identity,
      light_whole_v_identity,
      heavy_whole_v_max_mismatch_run,
      light_whole_v_max_mismatch_run,
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
      row_order
    )

} else {
  message("未检测到 selected_cells.xlsx：将输出全部配对 BCR 的简化总览。")

  simplified_with_QC <- all_key_info %>%
    mutate(
      review_status = expression_status,
      match_status = "all_records",
      matched_rows = 1L,
      row_order = row_number()
    ) %>%
    select(
      sample,
      cell_id,
      review_status,
      match_status,
      matched_rows,
      expression_status,
      priority,
      pair_count,
      heavy_subclass,
      light_locus,
      final_VH_aa,
      final_VL_aa,
      final_VH_aa_completed,
      final_VL_aa_completed,
      VH_len,
      VL_len,
      heavy_v_call,
      heavy_d_call,
      heavy_j_call,
      light_v_call,
      light_j_call,
      heavy_cdr3_aa,
      light_cdr3_aa,
      heavy_cdr3_aa_source,
      light_cdr3_aa_source,
      whole_variable_QC,
      whole_V_QC,
      terminal_completion,
      alignment_QC,
      VH_FR4_tail,
      VL_FR4_tail,
      cys_risk,
      heavy_whole_v_identity,
      light_whole_v_identity,
      heavy_whole_v_max_mismatch_run,
      light_whole_v_max_mismatch_run,
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
      row_order
    )
}


# ============================================================
# 5. export the simplified tables
# ============================================================

# Keep the full technical detail for troubleshooting; not for daily reading.
# final_*_aa_completed is the authoritative sequence for a non-ready row, since
# step 1 blanks final_*_aa for anything that is not A/B.
QC_details <- simplified_with_QC

# Main reading table: only what actually drives a decision.
simplified_view <- QC_details %>%
  transmute(
    sample,
    cell_id,
    review_status,
    priority,
    pair_count,
    heavy_subclass,
    light_locus,

    heavy_v_call,
    heavy_d_call,
    heavy_j_call,
    light_v_call,
    light_j_call,

    heavy_cdr3_aa,
    light_cdr3_aa,
    VH_len,
    VL_len,

    whole_variable_QC,
    terminal_completion,
    alignment_QC,
    cys_risk,

    final_VH_aa,
    final_VL_aa,

    match_status,
    matched_rows,
    row_order
  ) %>%
  arrange(row_order, priority)

# The table actually handed out for expression: A/B only, single pairing,
# successfully matched. Unlike an earlier version, C/D and multi-pair rows are
# never silently included.
simplified_for_company <- QC_details %>%
  filter(
    review_status == "READY",
    priority %in% ready_priorities,
    pair_count <= 1,
    !is.na(final_VH_aa), final_VH_aa != "",
    !is.na(final_VL_aa), final_VL_aa != ""
  ) %>%
  transmute(
    sample,
    cell_id,
    final_VH_aa,
    final_VL_aa,
    light_locus
  )

# Rows a human must look at: C/D, multi-pair, partial match, not found.
attention <- simplified_view %>%
  filter(review_status != "READY") %>%
  arrange(row_order, review_status)

not_found <- simplified_view %>%
  filter(review_status == "NOT_FOUND")

# Compact overview of every paired BCR, for scanning the whole run quickly.
all_key_info_export <- all_key_info %>%
  transmute(
    sample,
    cell_id,
    expression_status,
    priority,
    pair_count,
    heavy_subclass,
    light_locus,
    heavy_v_call,
    heavy_d_call,
    heavy_j_call,
    light_v_call,
    light_j_call,
    heavy_cdr3_aa,
    light_cdr3_aa,
    VH_len,
    VL_len,
    whole_variable_QC,
    terminal_completion,
    alignment_QC,
    cys_risk,
    final_VH_aa,
    final_VL_aa
  )

summary_table <- data.frame(
  metric = c(
    "input_file",
    "input_sheet",
    "selected_file_used",
    "cdr3_source_fields_available",
    "rows_in_simplified_view",
    "unique_requested_cells",
    "ready_for_expression_rows",
    "attention_rows",
    "not_found_rows",
    "multi_pair_or_duplicate_rows",
    "whole_variable_risk_rows"
  ),
  value = c(
    input_xlsx,
    source_sheet,
    ifelse(file.exists(selected_file), "YES", "NO"),
    ifelse(
      all(c("heavy_cdr3_aa_source", "light_cdr3_aa_source") %in% colnames(df)),
      "YES",
      "NO_legacy_input"
    ),
    as.character(nrow(simplified_view)),
    as.character(n_distinct(paste(simplified_view$sample, simplified_view$cell_id, sep = "||"))),
    as.character(nrow(simplified_for_company)),
    as.character(nrow(attention)),
    as.character(nrow(not_found)),
    as.character(sum(simplified_view$review_status == "CHECK_multi_pair_or_duplicate", na.rm = TRUE)),
    as.character(sum(QC_details$priority %in% whole_variable_risk_priorities, na.rm = TRUE))
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 6. export
# ============================================================

ensure_dir(output_dir)

write_xlsx(
  list(
    simplified_for_company = simplified_for_company,
    simplified_with_QC = simplified_view %>% select(-row_order),
    attention = attention %>% select(-row_order),
    not_found = not_found %>% select(-row_order),
    all_key_info = all_key_info_export,
    QC_details = QC_details %>% select(-row_order),
    summary = summary_table
  ),
  output_xlsx
)

cat("Done!\n")
cat("Input file:", input_xlsx, "\n")
cat("Input sheet:", source_sheet, "\n")
cat("Selected file used:", ifelse(file.exists(selected_file), "YES", "NO"), "\n")
cat("Output:", output_xlsx, "\n")
cat("Rows in simplified_with_QC:", nrow(simplified_view), "\n")
cat("Ready for expression:", nrow(simplified_for_company), "\n")
cat("Need attention:", nrow(attention), "\n")
cat("Not found:", nrow(not_found), "\n")
