library(tidyverse)
library(data.table)
library(optparse)

option_list <- list(
    optparse::make_option(c("--AnnotationData"), type = "character", default = NULL,
                        help = "Annotation data for variants containing ancestry specific AC/AN and cohort level AC/AN", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type = "character", default = NULL,
                        help = "Output file name prefix", metavar = "type"),
    optparse::make_option(c("--PipThreshold"), type = "double", default = 0.9,
                        help = "Minimum PIP threshold used to select variants [default %default]", metavar = "number"),
    optparse::make_option(c("--AdmixedSubpops"), type = "character", default = "oth",
                        help = "Comma-separated GVS subpopulation labels to remove for the no-admixed skew calculation [default %default]", metavar = "string")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (is.null(opt$AnnotationData)) {
    stop("--AnnotationData is required", call. = FALSE)
}
if (is.null(opt$OutputPrefix)) {
    stop("--OutputPrefix is required", call. = FALSE)
}

AnnotationPath <- opt$AnnotationData
OutputName <- paste0(opt$OutputPrefix, ".AncestrySkew.tsv.gz")
PipThreshold <- opt$PipThreshold
AdmixedSubpops <- opt$AdmixedSubpops %>%
    str_split(",", simplify = FALSE) %>%
    unlist() %>%
    str_trim() %>%
    str_to_lower()
AdmixedSubpops <- AdmixedSubpops[AdmixedSubpops != ""]

required_columns <- c("variant", "pip", "gvs_all_ac", "gvs_all_an")

####### HELPERS ########
maf <- function(af) {
    pmin(af, 1 - af)
}

missing_columns <- function(df, cols) {
    setdiff(cols, names(df))
}

is_gzip_file <- function(path) {
    con <- file(path, "rb")
    on.exit(close(con))
    magic <- readBin(con, what = "raw", n = 2)
    length(magic) == 2 && identical(as.integer(magic), c(0x1f, 0x8b))
}

read_annotation <- function(path) {
    if (is_gzip_file(path)) {
        return(fread(cmd = paste("gzip -dc", shQuote(path))))
    }
    fread(path)
}

require_columns <- function(df, cols, context) {
    missing <- missing_columns(df, cols)
    if (length(missing) > 0) {
        stop(
            paste0(context, " requires missing column(s): ", paste(missing, collapse = ", ")),
            call. = FALSE
        )
    }
}

get_subpop_from_af_col <- function(col) {
    col %>%
        str_remove("^gvs_") %>%
        str_remove("_af$")
}

choose_max_subpop <- function(df, subpops, prefix) {
    maf_cols <- paste0("gvs_", subpops, "_maf")
    max_maf_col <- paste0(prefix, "_max_maf")
    max_subpop_col <- paste0(prefix, "_max_subpop")

    maf_matrix <- as.matrix(df[, maf_cols, drop = FALSE])
    max_idx <- apply(maf_matrix, 1, function(row) {
        if (all(is.na(row))) {
            return(NA_integer_)
        }
        which.max(replace(row, is.na(row), -Inf))
    })

    df[[max_maf_col]] <- ifelse(
        is.na(max_idx),
        NA_real_,
        maf_matrix[cbind(seq_len(nrow(maf_matrix)), max_idx)]
    )
    df[[max_subpop_col]] <- ifelse(is.na(max_idx), NA_character_, subpops[max_idx])
    df
}

add_max_counts <- function(df, prefix) {
    max_subpop_col <- paste0(prefix, "_max_subpop")
    max_ac_col <- paste0(prefix, "_max_ac")
    max_an_col <- paste0(prefix, "_max_an")

    df[[max_ac_col]] <- NA_real_
    df[[max_an_col]] <- NA_real_

    for (subpop in unique(na.omit(df[[max_subpop_col]]))) {
        ac_col <- paste0("gvs_", subpop, "_ac")
        an_col <- paste0("gvs_", subpop, "_an")
        rows <- df[[max_subpop_col]] == subpop & !is.na(df[[max_subpop_col]])
        df[[max_ac_col]][rows] <- df[[ac_col]][rows]
        df[[max_an_col]][rows] <- df[[an_col]][rows]
    }

    df
}

add_background_counts <- function(df, prefix, all_ac_col, all_an_col) {
    max_ac_col <- paste0(prefix, "_max_ac")
    max_an_col <- paste0(prefix, "_max_an")
    bg_ac_col <- paste0(prefix, "_background_ac")
    bg_an_col <- paste0(prefix, "_background_an")

    df[[bg_ac_col]] <- df[[all_ac_col]] - df[[max_ac_col]]
    df[[bg_an_col]] <- df[[all_an_col]] - df[[max_an_col]]
    df
}

run_fisher <- function(df, prefix) {
    max_ac_col <- paste0(prefix, "_max_ac")
    max_an_col <- paste0(prefix, "_max_an")
    bg_ac_col <- paste0(prefix, "_background_ac")
    bg_an_col <- paste0(prefix, "_background_an")
    odds_col <- paste0(prefix, "_odds_ratio")
    p_col <- paste0(prefix, "_p_value")

    res <- apply(df[, c(max_ac_col, max_an_col, bg_ac_col, bg_an_col), drop = FALSE], 1, function(x) {
        if (any(is.na(x)) || any(x < 0)) {
            return(c(odds_ratio = NA_real_, p_value = NA_real_))
        }
        test <- fisher.test(matrix(x, nrow = 2, byrow = TRUE))
        c(odds_ratio = unname(test$estimate), p_value = test$p.value)
    })

    df[[odds_col]] <- as.numeric(res["odds_ratio", ])
    df[[p_col]] <- as.numeric(res["p_value", ])
    df
}

####### LOAD DATA ########
AnnotationDf <- read_annotation(AnnotationPath)

require_columns(AnnotationDf, required_columns, "Ancestry skew")

af_cols <- names(AnnotationDf) %>%
    keep(~ str_detect(.x, "^gvs_[^_]+_af$") && .x != "gvs_all_af" && .x != "gvs_max_af")

if (length(af_cols) == 0) {
    stop("No GVS subpopulation allele frequency columns found. Expected columns like gvs_afr_af.", call. = FALSE)
}

subpops <- map_chr(af_cols, get_subpop_from_af_col)
admixed_subpops <- intersect(AdmixedSubpops, subpops)
nonadmixed_subpops <- setdiff(subpops, admixed_subpops)

if (length(admixed_subpops) == 0) {
    stop(
        paste0(
            "None of the requested admixed subpops were found in GVS AF columns: ",
            paste(AdmixedSubpops, collapse = ", ")
        ),
        call. = FALSE
    )
}
if (length(nonadmixed_subpops) == 0) {
    stop("No non-admixed subpopulations remain after applying --AdmixedSubpops.", call. = FALSE)
}

count_cols <- c(paste0("gvs_", subpops, "_ac"), paste0("gvs_", subpops, "_an"))
require_columns(
    AnnotationDf,
    count_cols,
    paste0(
        "MAF-based max subpopulation AC/AN lookup. The attached annotation contains AF columns, ",
        "but exact skew needs matching per-subpopulation AC/AN columns"
    )
)

numeric_cols <- names(AnnotationDf) %>%
    keep(~ str_detect(.x, "^gvs_.*_(ac|an|af)$") || .x %in% c("pip"))
numeric_cols <- unique(numeric_cols)

removed_count_cols <- unlist(map(
    admixed_subpops,
    ~ c(paste0("gvs_no_admixed_removed_", .x, "_ac"), paste0("gvs_no_admixed_removed_", .x, "_an"))
))

OutputColumns <- c(
    "variant", "pip",
    "gvs_max_subpop", "gvs_max_maf", "gvs_max_ac", "gvs_max_an",
    "gvs_background_ac", "gvs_background_an", "gvs_odds_ratio", "gvs_p_value",
    "gvs_no_admixed_all_ac", "gvs_no_admixed_all_an",
    "gvs_no_admixed_max_subpop", "gvs_no_admixed_max_maf",
    "gvs_no_admixed_max_ac", "gvs_no_admixed_max_an",
    "gvs_no_admixed_background_ac", "gvs_no_admixed_background_an",
    "gvs_no_admixed_odds_ratio", "gvs_no_admixed_p_value",
    removed_count_cols
)

SkewInput <- AnnotationDf %>%
    select(variant, all_of(numeric_cols)) %>%
    mutate(across(all_of(numeric_cols), ~ as.numeric(.))) %>%
    filter(pip >= PipThreshold) %>%
    mutate(across(all_of(af_cols), maf, .names = "{.col}_maf"))

maf_col_names <- paste0(af_cols, "_maf")
names(SkewInput)[match(maf_col_names, names(SkewInput))] <- paste0("gvs_", subpops, "_maf")

if (nrow(SkewInput) == 0) {
    warning("No variants passed the PIP threshold; writing an empty output.")
    empty_output <- as_tibble(setNames(replicate(length(OutputColumns), logical(0), simplify = FALSE), OutputColumns))
    empty_output %>% write_tsv(OutputName)
    quit(save = "no", status = 0)
}

SkewInput <- SkewInput %>%
    choose_max_subpop(subpops, "gvs") %>%
    add_max_counts("gvs") %>%
    add_background_counts("gvs", "gvs_all_ac", "gvs_all_an") %>%
    run_fisher("gvs")

for (subpop in admixed_subpops) {
    SkewInput[[paste0("gvs_no_admixed_removed_", subpop, "_ac")]] <- SkewInput[[paste0("gvs_", subpop, "_ac")]]
    SkewInput[[paste0("gvs_no_admixed_removed_", subpop, "_an")]] <- SkewInput[[paste0("gvs_", subpop, "_an")]]
}

SkewInput <- SkewInput %>%
    mutate(
        gvs_no_admixed_all_ac = gvs_all_ac - rowSums(across(all_of(paste0("gvs_", admixed_subpops, "_ac"))), na.rm = FALSE),
        gvs_no_admixed_all_an = gvs_all_an - rowSums(across(all_of(paste0("gvs_", admixed_subpops, "_an"))), na.rm = FALSE)
    ) %>%
    choose_max_subpop(nonadmixed_subpops, "gvs_no_admixed") %>%
    add_max_counts("gvs_no_admixed") %>%
    add_background_counts("gvs_no_admixed", "gvs_no_admixed_all_ac", "gvs_no_admixed_all_an") %>%
    run_fisher("gvs_no_admixed")

SkewInput %>%
    select(all_of(OutputColumns)) %>%
    write_tsv(OutputName)
