################################################################################
# GENOMIC HALLMARKS OF DEPOT MEDROXYPROGESTERONE ACETATE-ASSOCIATED MENINGIOMAS

# SUMMARY
# DNA methylation analysis of n=10 female patients with depot medroxyprogesterone
# acetate (DMPA)-associated meningiomas, compared against two external reference
# cohorts. Patients are labeled DMPA-1 through DMPA-10 throughout, matching the
# manuscript.
#
# Reference cohorts:
#   Baylor (Bayley et al., Science Advances 2022): 110 meningiomas (n=90 WHO I), clinically annotated
#   Heidelberg (Capper et al., Nature 2018): 90 meningiomas (n=66 WHO I), clinically annotated; non-annotated cases excluded
#
# Arrays:  Baylor: EPIC v1 | DMPA: EPIC v2 | Heidelberg: 450k
# pd:      Baylor: pd1     | DMPA: pd2     | Heidelberg: pd3
#
# STRUCTURE
#   Sections 1-9 build the combined, batch-corrected beta matrix and phenotype (pd) table.
#   Sections 10+ are downstream analyses and figures. The script can be run end-to-end,
#   or restarted at Section 10 by reading the saved beta/pd objects from disk.
#
# REPRODUCIBILITY
#   Set `base_dir` in Section 1 to the project root that contains the data folders.
#   All input/output paths are derived from `base_dir`; no other paths require editing.
#   Raw arrays are deposited on GEO (DMPA: this study; Baylor: GSE189521; Heidelberg: GSE109381).
#
# PIPELINE OUTLINE
#    1. Load environment and set paths
#    2. Build phenotype table per cohort (pd1/pd2/pd3)
#    3. Baylor data: clinical annotation and filtering
#    4. Heidelberg data: clinical annotation and filtering
#    5. Read IDATs (minfi) and extract beta matrices
#    6. ChAMP probe filtering
#    7. SVD diagnostics (pre-ComBat)
#    8. ComBat cross-batch correction
#    9. SVD diagnostics (post-ComBat)
#   10. Restrict to top median absolute deviation (MAD) probes
#   11. Consensus clustering
#   12. Concordance plots
#   13. t-SNE plot
#   14. Principal component analysis (PCA) plot
#   15. Assign DMPA samples to Baylor methylation groups (random forest)
#   16. SeSAMe copy number analysis
#   17. PGR gene differential methylation
#   18. 11q22.1 cytoband differential methylation
#   19. Lollipop plots (GlioSeq targeted panel)
#   20. Oncoplot (GlioSeq mutation analysis)
#   21. Master oncoprint
#   22. PGR & 11q22.1 methylation in restricted Baylor reference subgroups
#   23. Progesterone signaling pathway differential methylation (GO:0032570, GO:0050847)

################################################################################

##==============================================================================
###-----1. LOAD ENVIRONMENT-----###
##==============================================================================

###-----Install packages and load libraries-----###
pkgs <- c(
  "ComplexHeatmap","circlize","dendsort","randomForest","ChAMP","DESeq2","biomaRt","Rtsne",
  "caret","pls","sesame","dendextend","cluster", "gdata","survival","survminer",
  "minfi","RColorBrewer","ConsensusClusterPlus",  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICv2manifest", "lumi", "ChAMPdata", "CopyNeutralIMA", 
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38", "Illumina450ProbeVariants.db", "DMRcate", 
  "AnnotationDbi", "org.Hs.eg.db", "limma", "matrixStats", "readr", "statmod", "NMF",
  "ggplot2", "tidyverse", "dplyr", "tidyr", "plyr", "plotly", "ineq", "DBI", "RSQLite", "missMethyl",
  "scales", "magick", "grid", "ggforce", "umap", "GenomicRanges", "ggalluvial", "stringr", 
  "scales", "sesame", "sesameData", "dplyr", "readxl")

bioc_avail <- tryCatch(BiocManager::available(), error = function(e) character())
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% bioc_avail) {
      BiocManager::install(p, ask = FALSE, update = FALSE)
    } else {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
  suppressPackageStartupMessages(
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  )
}

#Make sure dplyr verbs are used, not plyr
select     <- dplyr::select
filter     <- dplyr::filter
distinct   <- dplyr::distinct
mutate     <- dplyr::mutate
left_join  <- dplyr::left_join
arrange    <- dplyr::arrange
group_by   <- dplyr::group_by
ungroup    <- dplyr::ungroup
count      <- dplyr::count
summarise  <- dplyr::summarise

###------Set paths (EDIT base_dir ONLY)------###
# base_dir is the project root containing the cohort data folders referenced below.
# Every other path in this script is derived from it.
base_dir <- path.expand("~/Depo_Meningiomas")   # <-- EDIT THIS to your project root

# Results/working directory: all script outputs are written here.
working_directory <- file.path(base_dir, "Combined Analysis_WHO I_Final")
ResultsDir        <- working_directory
dir.create(working_directory, recursive = TRUE, showWarnings = FALSE)
setwd(working_directory)
message("Working directory: ", normalizePath(working_directory))

###-----IDAT input directories (one per cohort/array)-----###
epic2_dir <- file.path(base_dir, "Beta values/DMPA Idat files")     # DMPA       (EPIC v2)
epic_dir  <- file.path(base_dir, "Baylor/GSE189521_RAW")            # Baylor     (EPIC v1)
heid_dir  <- file.path(base_dir, "GSE109381_MNG_2018_Nature_paper") # Heidelberg (450k)

#Clear existing plots
while (!is.null(dev.list())) dev.off()


###-----Helpers-----###
#Standardize Sentrix positions
normalize_pos <- function(x){
  x <- toupper(as.character(x))
  x <- sub("^R([0-9])C", "R0\\1C", x)
  x <- sub("C([0-9])$", "C0\\1", x)
  x
}

#Enumerate paired IDATs from a directory (robust to .gz)
idat_index_from_dir <- function(idat_dir) {
  f <- list.files(idat_dir, pattern="\\.idat(\\.gz)?$", recursive=TRUE, full.names=TRUE, ignore.case=TRUE)
  if (!length(f)) stop("No IDAT files under: ", idat_dir)
  is_red <- grepl("(_Red\\.idat(\\.gz)?)$", f, ignore.case=TRUE)
  is_grn <- grepl("(_Grn\\.idat(\\.gz)?)$", f, ignore.case=TRUE)
  typ <- ifelse(is_red, "Red", ifelse(is_grn, "Grn", NA_character_))
  keep <- !is.na(typ); f <- f[keep]; typ <- typ[keep]
  basenames <- sub("(_Red|_Grn)\\.idat(\\.gz)?$", "", f, ignore.case=TRUE)
  tab <- as.data.frame.matrix(table(basenames, typ)); tab[is.na(tab)] <- 0
  have_red <- if ("Red" %in% colnames(tab)) tab[,"Red"] > 0 else FALSE
  have_grn <- if ("Grn" %in% colnames(tab)) tab[,"Grn"] > 0 else FALSE
  paired <- rownames(tab)[have_red & have_grn]
  if (!length(paired)) stop("Found IDATs in ", idat_dir, " but no Red/Green pairs.")
  leaf <- basename(paired)
  sid  <- sub("_.*$", "", leaf)
  pos  <- normalize_pos(sub("^[^_]*_", "", leaf))
  parent <- basename(dirname(paired))
  gsm_like <- grepl("^GSM\\d+$", parent, ignore.case=TRUE)
  data.frame(
    Sample_Name      = ifelse(gsm_like, parent, paste0(sid, "_", pos)),
    Sentrix_ID       = sid,
    Sentrix_Position = pos,
    Basename         = paired,
    stringsAsFactors = FALSE
  )
}

#Build PD from files
pd_from_sheet_or_files <- function(idat_dir, sheet_path=NULL) {
  idx <- idat_index_from_dir(idat_dir)
  if (is.null(sheet_path)) {
    pd <- idx
    return(pd)
  }
  delim <- if (grepl("\\.tsv$", sheet_path, ignore.case=TRUE)) "\t" else ","
  ss <- suppressMessages(readr::read_delim(sheet_path, delim=delim, show_col_types=FALSE))
  cn <- names(ss)
  pick <- function(cn, exact, rx){
    hit <- intersect(exact, cn); if (length(hit)) return(hit[1])
    for (r in rx) { m <- grep(r, cn, value=TRUE, ignore.case=TRUE, perl=TRUE); if (length(m)) return(m[1]) }
    NA_character_
  }
  slide_col <- pick(cn,
                    exact=c("Sentrix_ID","Slide","Array_ID","Chip_ID","Barcode","BeadChip","BeadChip_ID"),
                    rx=c("sentrix.*id","^slide(_id)?$","^array(_id)?$","chip(_id)?","barcode","beadchip"))
  pos_col <- pick(cn,
                  exact=c("Sentrix_Position","Array","Array_Position","Position","Well","Sentrix_Array","Sentrix_Array_Position"),
                  rx=c("sentrix.*pos","^array(_position)?$","^position$","^well$","^r\\d{1,2}c\\d{1,2}$"))
  if (is.na(slide_col) || is.na(pos_col)) {
    # fall back to files only
    return(idx)
  }
  t_sid <- as.character(ss[[slide_col]])
  t_pos <- normalize_pos(ss[[pos_col]])
  keys <- data.frame(Sentrix_ID=t_sid, Sentrix_Position=t_pos, .row=seq_len(nrow(ss)))
  joined <- merge(keys, idx, by=c("Sentrix_ID","Sentrix_Position"), all.x=TRUE)
  joined <- joined[order(joined$.row),]
  pd <- ss
  pd$Sentrix_ID <- t_sid
  pd$Sentrix_Position <- t_pos
  pd$Basename <- joined$Basename
  if (anyNA(pd$Basename)) {
    stop("Some rows in the sheet did not match actual IDAT pairs in ", idat_dir,
         ". Example:\n", capture.output(print(head(pd[is.na(pd$Basename), c("Sentrix_ID","Sentrix_Position")], 8))))
  }
  pd
}

#Uppercase & strip spaces, for robust matching
norm <- function(x) toupper(gsub("\\s+", "", as.character(x)))

#Basename leaf from file path
basename_from_path <- function(p) {
  leaf <- basename(as.character(p))
  sub("(_Red|_Grn)?\\.idat(\\.gz)?$", "", leaf, ignore.case = TRUE)
}

#Split basename leaf into Sentrix_ID / Sentrix_Position
split_sid_pos <- function(bn) {
  sid <- sub("_.*$", "", bn)
  pos_raw <- sub("^[^_]*_", "", bn)
  data.frame(Sentrix_ID = sid, Sentrix_Position = normalize_pos(pos_raw), stringsAsFactors = FALSE)
}

#Cast selected columns to character (safe merges)
to_char <- function(df, cols) { for (nm in intersect(cols, names(df))) df[[nm]] <- as.character(df[[nm]]); df }

#Align data.frame to target set of columns
align_cols <- function(x, target_cols) {missing <- setdiff(target_cols, names(x)); for (nm in missing) x[[nm]] <- NA_character_; x[, target_cols, drop=FALSE]}

#Safe getter
pick_col <- function(df, nm, default = NA_character_) {
  if (!is.na(nm) && nzchar(nm) && nm %in% names(df)) df[[nm]] else default
}

#Normalize WHO grade
normalize_who <- function(x){
  x0 <- norm(x)
  x0 <- sub("\\(2016\\)|\\(2021\\)", "", x0)  # strip year tags
  x0 <- sub("^WHO", "", x0)                    # "WHO I" -> " I"
  x0 <- gsub("[^IV0-9]", "", x0)               # keep numerals
  dplyr::case_when(
    x0 %in% c("I","1")   ~ "WHO I",
    x0 %in% c("II","2")  ~ "WHO II",
    x0 %in% c("III","3") ~ "WHO III",
    TRUE ~ NA_character_
  )
}

#Normalize sex
normalize_sex <- function(x){
  x0 <- toupper(trimws(collapse_spaces(x)))
  dplyr::case_when(
    grepl("^F", x0) ~ "F",
    grepl("^M", x0) ~ "M",
    TRUE            ~ NA_character_
  )
}

pull_col <- function(df, colname) {
  if (is.na(colname)) return(rep(NA_character_, nrow(df)))
  as.character(df[[colname]])
}

###-----GLOBAL COHORT HARMONIZATION (ONE-STOP)-----###
PRIMARY_COHORT <- "DMPA"  # the cohort of interest in this study

harmonize_cohort <- function(x) {
  x <- gsub("\u00A0"," ", as.character(x))
  x <- trimws(gsub("[[:space:]]+"," ", x))
  
  dplyr::case_when(
    tolower(x) %in% c("dmpa","dmpa cohort") ~ "DMPA",
    tolower(x) %in% c("baylor","baylor cohort") ~ "Baylor",
    tolower(x) %in% c("heidelberg","heidelberg cohort","hd") ~ "Heidelberg",
    TRUE ~ x
  )
}

apply_pd_harmonization <- function(pd) {
  stopifnot("Cohort" %in% names(pd))
  pd$Cohort <- factor(harmonize_cohort(pd$Cohort))
  pd
}

COHORT_COLORS <- c(
  DMPA       = "#582C83",  # purple
  Baylor     = "grey95",
  Heidelberg = "grey85"
)

##==============================================================================
###-----2. BUILD PD PER COHORT (pd1/pd2/pd3)-----###
##==============================================================================

#---Build separate pd file for each batch/cohort---#
pd1 <- pd_from_sheet_or_files(epic_dir);   pd2 <- pd_from_sheet_or_files(epic2_dir);   pd3 <- pd_from_sheet_or_files(heid_dir)

# Tag cohorts (add a Cohort column)
pd1$ArrayType <- "EPIC";   pd1$Batch <- "EPIC_v1";          pd1$Cohort <- "Baylor"
pd2$ArrayType <- "EPICv2"; pd2$Batch <- "EPIC_v2";          pd2$Cohort <- "DMPA"
pd3$ArrayType <- "450K";   pd3$Batch <- "Heidelberg_MNG";   pd3$Cohort <- "Heidelberg"

message(sprintf("Total per cohort: EPIC (Baylor): %d | EPICv2 (DMPA): %d | 450k (Heidelberg): %d", nrow(pd1), nrow(pd2), nrow(pd3)))

char_cols <- c("Sample_Name","Sentrix_ID_(.idat)","Sentrix_Position","Basename",
               "Batch","ArrayType","Cohort","WHO_Grade","Gender", "Sex", "WHO Grade")

pd1 <- to_char(pd1, char_cols)
pd2 <- to_char(pd2, char_cols)
pd3 <- to_char(pd3, char_cols)

all_cols <- Reduce(union, list(names(pd1), names(pd2), names(pd3)))
pd1 <- align_cols(pd1, all_cols)
pd2 <- align_cols(pd2, all_cols)
pd3 <- align_cols(pd3, all_cols)

stopifnot(identical(names(pd1), names(pd2)), identical(names(pd2), names(pd3)))


##==============================================================================
###-----3. BAYLOR DATA-----###
#Attach clinical annotations and choose filter
##==============================================================================

###-----Load Baylor clinical annotations and join to pd1-----###
#Baylor clinical annotations
meta_path  <- file.path(base_dir, "Baylor/GSE189521_Clinical_data (Bayley et al).xlsx")
meta_sheet <- "2. Clinical and genomic dataSH"                            
meta <- readxl::read_excel(meta_path, sheet = meta_sheet, skip = 2) #Read Excel sheet with headers on row 3

#Identify the required columns in the Excel sheet
cn <- names(meta)
col_idat <- cn[match(tolower("Idat file"), tolower(cn))]
col_who  <- cn[match(tolower("WHO grade"), tolower(cn))]
col_gender<- cn[match(tolower("Gender"),    tolower(cn))]
if (any(is.na(c(col_idat, col_who)))) {
  stop("Couldn't find columns 'Idat file' (B) and 'WHO grade' (G) in header row 3 of the Excel tab.")
}

#Build keys from Excel
meta$BasenameLeaf <- basename_from_path(meta[[col_idat]])
sidpos <- split_sid_pos(meta$BasenameLeaf)
meta_keys <- cbind(
  data.frame(BasenameLeaf = meta$BasenameLeaf, stringsAsFactors = FALSE),
  sidpos,
  WHO_Grade = meta[[col_who]],
  Gender    = meta[[col_gender]],
  stringsAsFactors = FALSE
)

#Prepare pd1 join keys (pd1$Basename comes from idat discovery)
pd1$BasenameLeaf <- basename_from_path(pd1$Basename)

#Join by BasenameLeaf (exact filename match)
pd1 <- merge(pd1,
             meta_keys[, c("BasenameLeaf","WHO_Grade","Gender","Sentrix_ID","Sentrix_Position")],
             by = "BasenameLeaf", all.x = TRUE, sort = FALSE)


#Fallback join if WHO_Grade or Gender are still NA
if (sum(!is.na(pd1$WHO_Grade)) < nrow(pd1) || sum(!is.na(pd1$Gender)) < nrow(pd1)) {
  pd1$WHO_Grade <- NULL; pd1$Gender <- NULL
  pd1 <- merge(pd1,
               unique(meta_keys[, c("Sentrix_ID","Sentrix_Position","WHO_Grade","Gender")]),
               by = c("Sentrix_ID","Sentrix_Position"),
               all.x = TRUE, sort = FALSE)
}

###-----Filter Baylor data (EPIC v1, pd1)-----###
#WHO I, both sexes 
allowed_grades <- c("WHO I")                        
allowed_genders <- c("M", "F")                     

#Apply filters set above to pd1 sheet
pd1 <- subset(pd1,
              norm(WHO_Grade) %in% norm(allowed_grades) &
              norm(Gender)    %in% norm(allowed_genders))

message(sprintf("EPIC (Baylor) retained %d samples after WHO Grade (%s) + Gender (%s) filter.",
                nrow(pd1),
                paste(allowed_grades, collapse = ", "),
                paste(allowed_genders, collapse = ", ")))

if (nrow(pd1) == 0) stop("Filter removed all EPIC samples. Check Excel values or adjust 'allowed_grades'/'allowed_genders'.")

#Clean up helper column
pd1$BasenameLeaf <- NULL


##==============================================================================
###-----4. HEIDELBERG DATA-----###
#Attach clinical annotations and choose filter
##==============================================================================

###-----Load Heidelberg clinical annotations and join to pd3-----###
#Heidelberg clinical annotations
meta_heid_path  <- file.path(base_dir, "DMPABaylorHeidelbergCombined/HeidelbergClinicalAnnotations_MNG.xlsx")
meta_heid_sheet <- "Heidelberg_MNG_CLEAN"
skip_rows_heid  <- 1        # change to 2 if header starts on row-3

metaH <- readxl::read_excel(meta_heid_path, sheet = "Heidelberg_MNG_CLEAN", skip = 1)

###-----Detect header columns in Excel sheet (robust matching)-----###
collapse_spaces <- function(x){
  x <- gsub("\u00A0", " ", x)       # NBSP -> space
  x <- gsub("[[:space:]]+", " ", x) # collapse spaces
  trimws(x)
}

hdr_raw   <- names(metaH)
hdr_norm  <- tolower(collapse_spaces(hdr_raw))
hdr_simpl <- gsub("[^a-z0-9]+", "", hdr_norm)

choose_exact <- function(cands){
  # prefer exact header match
  hit <- intersect(cands, hdr_raw)
  if (length(hit)) return(hit[1])
  # then match in simplified space (strip punctuation/spacing)
  simp_cands <- gsub("[^a-z0-9]+","", tolower(cands))
  idx <- which(hdr_simpl %in% simp_cands)
  if (length(idx)) return(hdr_raw[idx[1]])
  NA_character_
}

pick_regex <- function(rx){
  idx <- grep(rx, hdr_norm, perl = TRUE)
  if (length(idx)) hdr_raw[idx[1]] else NA_character_
}

col_basename <- choose_exact(c(
  "BasenameLeaf", "Basename", "IDAT", "IDAT file"
))
if (is.na(col_basename)) col_basename <- pick_regex("basename|idat")

col_sid <- choose_exact(c(
  "Sentrix_ID", "Sentrix_ID_(.idat)", "Array ID", "Slide", "Array_ID", "Chip_ID", "Barcode", "BeadChip"
))
if (is.na(col_sid)) col_sid <- pick_regex("sentrix.*id|^slide(_id)?$|array(_id)?|chip(_id)?|barcode|beadchip")

col_pos <- choose_exact(c(
  "Sentrix_Position", "Array", "Array_Position", "Position", "Well", "Sentrix_Array_Position"
))
if (is.na(col_pos)) col_pos <- pick_regex("sentrix.*pos|array.*pos|^position$|^well$|r\\d{1,2}c\\d{1,2}")

col_who <- choose_exact(c(
  "WHO_Grade", "WHO Grade (2016)", "Pathological Diagnosis (WHO 2016)", "WHO grade", "Grade"
))
if (is.na(col_who)) col_who <- pick_regex("\\bwho\\b.*\\bgrade\\b|\\bgrade\\b.*\\bwho\\b|^grade$")

col_sex <- choose_exact(c(
  "Sex", "Gender"
))
if (is.na(col_sex)) col_sex <- pick_regex("^sex$|^gender$|\\bsex\\b|\\bgender\\b")

message("Detected columns: ",
        paste(c(
          paste0("basename=", col_basename),
          paste0("sid=", col_sid),
          paste0("pos=", col_pos),
          paste0("who=", col_who),
          paste0("sex=", col_sex)
        ), collapse=" | "))

###-----Build clinical key table-----###
#From metaH using  detected headers 
BasenameLeaf_vec <- if (!is.na(col_basename)) basename_from_path(pull_col(metaH, col_basename)) else rep(NA_character_, nrow(metaH))
Sentrix_ID_vec   <- pull_col(metaH, col_sid)
Sentrix_Pos_vec  <- normalize_pos(pull_col(metaH, col_pos))
WHO_raw_vec      <- pull_col(metaH, col_who)
Sex_raw_vec      <- pull_col(metaH, col_sex)

meta_heid_keys <- data.frame(
  BasenameLeaf     = BasenameLeaf_vec,
  Sentrix_ID       = Sentrix_ID_vec,
  Sentrix_Position = Sentrix_Pos_vec,
  WHO_Grade        = normalize_who(WHO_raw_vec),
  Gender           = normalize_sex(Sex_raw_vec),
  stringsAsFactors = FALSE
)

cat("Heidelberg clinical rows in key table: ", nrow(meta_heid_keys), "\n")
if (nrow(meta_heid_keys) < 5) {
  warning("Unusually few clinical rows detected; confirm your header detection and `skip_rows_heid`.")
}

###-----Join clinical key table to pd3-----###
req_cols <- c("Basename","Sentrix_ID","Sentrix_Position")
if (!all(req_cols %in% names(pd3))) {
  stop("pd3 is missing required columns: ", paste(setdiff(req_cols, names(pd3)), collapse=", "),
       ". Ensure pd3 was created via pd_from_sheet_or_files(heid_dir) before this block.")
}
pd3$BasenameLeaf     <- basename_from_path(pd3$Basename)
pd3$Sentrix_ID       <- as.character(pd3$Sentrix_ID)
pd3$Sentrix_Position <- normalize_pos(pd3$Sentrix_Position)

#Primary merge by BasenameLeaf
pd3_j <- merge(
  pd3,
  meta_heid_keys[, c("BasenameLeaf","WHO_Grade","Gender")],
  by = "BasenameLeaf", all.x = TRUE, sort = FALSE
)

#Fallback merge by (Sentrix_ID, Sentrix_Position) where needed
need_fallback <- which(is.na(pd3_j$WHO_Grade) | is.na(pd3_j$Gender))
if (length(need_fallback)) {
  add_sidpos <- unique(meta_heid_keys[, c("Sentrix_ID","Sentrix_Position","WHO_Grade","Gender")])
  pd3_j <- merge(
    pd3_j, add_sidpos,
    by = c("Sentrix_ID","Sentrix_Position"),
    all.x = TRUE, sort = FALSE, suffixes = c("", ".sidpos")
  )
  # Coalesce (prefer BasenameLeaf match)
  pd3_j$WHO_Grade <- ifelse(is.na(pd3_j$WHO_Grade), pd3_j$WHO_Grade.sidpos, pd3_j$WHO_Grade)
  pd3_j$Gender    <- ifelse(is.na(pd3_j$Gender),    pd3_j$Gender.sidpos,    pd3_j$Gender)
  pd3_j$WHO_Grade.sidpos <- NULL
  pd3_j$Gender.sidpos    <- NULL
}

#Diagnostics
message("Attached WHO_Grade for ", sum(!is.na(pd3_j$WHO_Grade)), " / ", nrow(pd3_j), " samples.")
message("Attached Gender for ", sum(!is.na(pd3_j$Gender)),
        " / ", nrow(pd3_j), " samples.")


###------Filter Heidelberg data-----###
#Option 1 (Default): WHO I, both sexes
allowed_heid_grades  <- c("WHO I")
allowed_heid_genders <- c("M","F")

grade_ok  <- !is.na(pd3_j$WHO_Grade) & normalize_who(pd3_j$WHO_Grade) %in% normalize_who(allowed_heid_grades)
gender_ok <- !is.na(pd3_j$Gender)    & norm(pd3_j$Gender)        %in% norm(allowed_heid_genders)
keep_mask <- grade_ok & gender_ok
pd3_keep  <- pd3_j[keep_mask, , drop = FALSE]

#Diagnostics after filtering
message(sprintf(
  "Heidelberg (450k) retained %d of %d samples after WHO Grade (%s) + Gender (%s).",
  nrow(pd3_keep), nrow(pd3_j),
  paste(allowed_heid_grades, collapse = ", "),
  paste(allowed_heid_genders, collapse = ", ")
))

#Filter pd3
pd3_keep$BasenameLeaf <- NULL
pd3 <- pd3_keep

#pd3 diagnostics 
message("Heidelberg WHO x Gender table:")
print(with(pd3, table(WHO_Grade, Gender)))

###-----Make pd1, pd2, pd3 have identical columns and types-----###
#Adds WHO Grade and Gender to pd2 (DMPA data) so the datasets combine together cleanly
pd1 <- to_char(pd1, char_cols)  # Baylor (EPIC)
pd2 <- to_char(pd2, char_cols)  # DMPA  (EPICv2)
pd3 <- to_char(pd3, char_cols)  # Heidelberg (450k)

#Add any missing columns to each with NA, then order columns identically
align_cols <- function(x, target_cols) {
  missing <- setdiff(target_cols, names(x))
  for (nm in missing) x[[nm]] <- NA_character_
  x[, target_cols, drop = FALSE]
}

all_cols <- Reduce(union, list(names(pd1), names(pd2), names(pd3)))
pd1 <- align_cols(pd1, all_cols)
pd2 <- align_cols(pd2, all_cols)
pd3 <- align_cols(pd3, all_cols)

# Sanity check - make sure all three have identical columns and order
stopifnot(identical(names(pd1), names(pd2)), identical(names(pd2), names(pd3)))

pd1$Sentrix_Position <- normalize_pos(pd1$Sentrix_Position)
pd2$Sentrix_Position <- normalize_pos(pd2$Sentrix_Position)
pd3$Sentrix_Position <- normalize_pos(pd3$Sentrix_Position)

message("Aligned PDs -> Baylor: ", nrow(pd1), " | DMPA: ", nrow(pd2), " | Heidelberg: ", nrow(pd3))
message("You now have a pd1 file containing filtered patients from Baylor, a pd2 file containing all DMPA patients, and a pd3 file containing filtered patients from Heidelberg")
 
               
##==============================================================================
###-----5. READ INTO MINFI AND EXTRACT BETA MATRICES-----###
#Starts by treating pd1, pd2, and pd3 (because they were run on different arrays). 
#Will combine them later
##==============================================================================

#Read raw data (RGChannelSet) with minfi#
rg1 <- minfi::read.metharray.exp(targets = pd1, extended = TRUE, force = TRUE)
rg2 <- minfi::read.metharray.exp(targets = pd2, extended = TRUE, force = TRUE)
rg3 <- minfi::read.metharray.exp(targets = pd3, extended = TRUE, force = TRUE)

#Detection P (from raw rg)#
detP1 <- minfi::detectionP(rg1)
detP2 <- minfi::detectionP(rg2)
detP3 <- minfi::detectionP(rg3)

#Noob normalize EACH batch (single-sample; safe to do separately)#
gr1 <- minfi::preprocessNoob(rg1)
gr2 <- minfi::preprocessNoob(rg2)
gr3 <- minfi::preprocessNoob(rg3)

#Extract matrices & INTERSECT probes#
beta1 <- minfi::getBeta(gr1); beta2 <- minfi::getBeta(gr2); beta3 <- minfi::getBeta(gr3)
M1    <- minfi::getM(gr1);    M2    <- minfi::getM(gr2); M3    <- minfi::getM(gr3)

#Strip EPIC v2 suffixes#
# Keep your EPICv2 suffix stripping; safe to apply to all (no effect on 450k)
strip_suffix <- function(v) sub("_(BC|TC)[0-9]*$", "", v, perl = TRUE)
for (nm in c("beta1","beta2","beta3","M1","M2","M3","detP1","detP2","detP3")) {
  m <- get(nm); rownames(m) <- strip_suffix(rownames(m)); assign(nm, m)
}

collapse_dups <- function(mat){
  rn <- rownames(mat); sp <- split(seq_len(nrow(mat)), rn, drop=TRUE)
  out <- vapply(sp, function(ix) colMeans(mat[ix, , drop=FALSE], na.rm=TRUE), FUN.VALUE = numeric(ncol(mat)))
  out <- t(out); rownames(out) <- names(sp); colnames(out) <- colnames(mat); out[order(rownames(out)), , drop=FALSE]
}
beta1 <- collapse_dups(beta1); beta2 <- collapse_dups(beta2); beta3 <- collapse_dups(beta3)
M1    <- collapse_dups(M1);    M2    <- collapse_dups(M2);    M3    <- collapse_dups(M3)
detP1 <- collapse_dups(detP1); detP2 <- collapse_dups(detP2); detP3 <- collapse_dups(detP3)

#Three-way intersection of datasets 
#Build combined matrices with a single columns order reference
common <- Reduce(intersect, list(rownames(beta1), rownames(beta2), rownames(beta3)))
cat("Common CpGs across EPIC, EPICv2, and 450k:", length(common), "\n")

#Define the exact combined column order once
cols_beta <- c(colnames(beta1), colnames(beta2), colnames(beta3))

beta <- cbind(beta1[common, , drop = FALSE], beta2[common, , drop = FALSE], beta3[common, , drop = FALSE])

M <- cbind(M1[common, , drop = FALSE], M2[common, , drop = FALSE], M3[common, , drop = FALSE])

#Build detP in the SAME order as 'cols_beta'
detP <- cbind(
  detP1[common, colnames(beta1), drop = FALSE],
  detP2[common, colnames(beta2), drop = FALSE],
  detP3[common, colnames(beta3), drop = FALSE]
)

stopifnot(identical(colnames(beta), cols_beta),
          identical(colnames(M),    cols_beta),
          identical(colnames(detP), cols_beta))

#Build pd in the same order
pd <- rbind(
  as.data.frame(pData(rg1), stringsAsFactors = FALSE),
  as.data.frame(pData(rg2), stringsAsFactors = FALSE),
  as.data.frame(pData(rg3), stringsAsFactors = FALSE)
)

#Ensure rownames are sample names, then reorder to 'cols_beta'
if (is.null(rownames(pd)) || !all(cols_beta %in% rownames(pd))) {
  # fall back to sampleNames from RG sets
  rn <- c(colnames(beta1), colnames(beta2), colnames(beta3))
  rownames(pd) <- rn
}
pd <- pd[cols_beta, , drop = FALSE]

#Final invariants
stopifnot(identical(colnames(beta), rownames(pd)),
          identical(colnames(detP), colnames(beta)),
          identical(colnames(M),    colnames(beta)))

cat("Combined beta/M/detP/pd built with perfectly matching columns/rows.\n")

#Save pd and beta files
write.csv(pd, "pd_Combined_WHO I.csv")
write.csv(beta, "beta_Combined_WHO I.csv")


##==============================================================================
###-----6. ChAMP FILTER-----###
##==============================================================================

###-----Sanitize data before ChAMP Filter-----###
#Keep everything as basic R objects
beta <- as.matrix(beta); storage.mode(beta) <- "double"
detP <- as.matrix(detP); storage.mode(detP) <- "double"
pd   <- as.data.frame(pd, stringsAsFactors = FALSE)

#Ensure we have Sample_Name and set it to the beta colnames
pd$Sample_Name <- rownames(pd)

#Re-order pd so Sample_Name order == beta columns
pd <- pd[match(colnames(beta), pd$Sample_Name), , drop = FALSE]

#Final invariants
stopifnot(
  identical(pd$Sample_Name,   colnames(beta)),
  identical(colnames(detP),   colnames(beta)),
  !any(is.na(pd$Sample_Name)),
  !any(duplicated(pd$Sample_Name))
)

###-----ChAMP Filter-----###
#Install and load ChAMP dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ChAMPdata", quietly = TRUE)) {
  BiocManager::install("ChAMPdata", ask = FALSE, update = FALSE)
}
library(ChAMPdata) # Load the package for access

data_list <- data(package = "ChAMPdata")$results[, "Item"]

for (data_name in data_list) {
  suppressMessages(suppressWarnings(utils::data(data_name, package = "ChAMPdata")))
}

#Run ChAMP Filter
myFilt <- ChAMP::champ.filter(
  beta         = beta,
  pd           = pd,
  detP         = detP,
  filterDetP   = TRUE,
  detPcut      = 0.01,
  ProbeCutoff  = 0.01,   # Drop probes failing in >1% of samples
  SampleCutoff = 0.01,   #Drop samples with >1% probes failing
  autoimpute   = TRUE,  # TRUE to impute missing values in ChAMP
  filterNoCG   = TRUE,
  filterSNPs   = TRUE,
  filterXY = TRUE, 
  filterMultiHit = TRUE,
  arraytype    = "450k"  #Common probe universe
)

#Initial data sanitizing#
beta_filt <- as.matrix(myFilt$beta)
pd_filt   <- as.data.frame(myFilt$pd)

#Quality check on filtered data#
cat("Probes:", nrow(beta), "Samples:", ncol(beta), "\n")
cat("DetP NA count:", sum(is.na(detP)), "\n")
apply(detP < 0.01, 2, mean) |> summary()  # fraction of passing probes per sample
apply(detP < 0.01, 1, mean) |> summary()  # fraction of passing samples per probe


##==============================================================================
###-----7. RUN SVD in ChAMP (Pre-ComBat)-----###
##==============================================================================

###-----Sanitize data prior to SVD - make beta/pd exactly what champ.SVD expects-----###
#Force beta to a plain *data.frame* with class length 1
#Ensure pd is also a plain data.frame with rownames = sample IDs
#Align rownames to beta_df cols
beta_df <- as.data.frame(beta_filt, stringsAsFactors = FALSE, check.names = FALSE); class(beta_df) <- "data.frame"
pd_df   <- pd_filt; if (is.null(rownames(pd_df))) rownames(pd_df) <- pd_df$Sample_Name
pd_df   <- pd_df[colnames(beta_df), , drop = FALSE]
stopifnot(identical(colnames(beta_df), rownames(pd_df)))

#Keep only covariates with ≥2 unique non-NA levels
nuniq <- function(x) length(unique(x[!is.na(x)]))
vary_cols <- names(pd_df)[vapply(pd_df, nuniq, integer(1)) >= 2]
pd_svd <- pd_df[, vary_cols, drop = FALSE]

##Run SVD on the sanitized, imputed data##
svd_res <- ChAMP::champ.SVD(
  beta       = beta_df,
  pd         = pd_svd,
  resultsDir = file.path(getwd(), "SVD_imputed_3way"),  #Adjust this file name based on how many datasets included
  PDFplot    = TRUE
)


##==============================================================================
###-----8. RUN ComBat IN ChAMP-----###
#Cross batch normalization - if deemed necessary by SVD#
##==============================================================================

#Sanity check 
beta <- as.matrix(myFilt$beta); storage.mode(beta) <- "double"
pd   <- as.data.frame(pd_svd, stringsAsFactors = FALSE)
stopifnot(identical(colnames(beta), rownames(pd)))

###-----Batch factor-----###
#Option 1: for 3 datasheets (Default. Batch alone carries three cohorts - Epic_v1, Epic_v2, Heidelberg)
batch <- factor(pd$Batch)     # <- use Batch only
if (anyNA(batch) || nlevels(batch) < 2) 
  stop("Batch factor needs ≥2 levels and no NA")

#Option 2: for 2 datasheets - create composite batch (Batch × ArrayType)
#By default this is not run. This would be if we were only comparing DMPA to Baylor OR Heidelberg, not both)
#batch <- interaction(pd$Batch, pd$ArrayType, drop = TRUE)
#if (anyNA(batch) || nlevels(batch) < 2) stop("Batch/ArrayType need ≥ 2 levels and no NA.")

###-----ComBat Pipeline----###
#Convert β → M 
eps <- 1e-6
beta_clipped <- pmin(pmax(beta, eps), 1 - eps)
M <- log2(beta_clipped / (1 - beta_clipped))

#Run ComBat on M-values
M_combat <- sva::ComBat(
  dat         = M,
  batch       = batch,
  par.prior   = TRUE,
  prior.plots = FALSE
)

#Back-transform to β
beta_combat <- 2^M_combat / (1 + 2^M_combat)


##==============================================================================
###-----9. RUN SVD in ChAMP (Post-ComBat)-----###
##==============================================================================

###-----SVD after ComBat-----###
svd_res <- ChAMP::champ.SVD(
  beta       = as.data.frame(beta_combat),          # pass the data.frame version
  pd         = as.data.frame(pd),
  resultsDir = file.path(getwd(), "SVD_after_ComBat_3way"),
  PDFplot    = TRUE
)

#OPTIONAL - Save Beta Combat and pd files to csv 
write.csv(beta_combat, "beta_combat_Combined Analysis_WHO I.csv", quote = FALSE)     
write.csv(pd_df,       "pd_Combined Analysis_WHO I.csv",           row.names = TRUE)


##==============================================================================
###-----DOWNSTREAM ANALYSIS-----###
##==============================================================================

###-----Optional - read in saved beta combat sheet (to start here, if beta_combat is not in memory)-----###
#If doing this, need to also load packages and directories at the top of this script
#Can alternatively skip this section if running the whole pipeline from the top

#Read saved beta sheet
beta_df <- read.csv(file.path(base_dir, "Combined Analysis_WHO I_Final/beta_combat_Combined Analysis_WHO I.csv"), row.names = 1, check.names = FALSE)
num_cols <- vapply(beta_df, is.numeric, logical(1))
beta_combat <- as.matrix(beta_df[, num_cols, drop = FALSE]) #"beta_combat" downstream refers to whichever beta matrix was read in above

#Read pd (phenotype table)
pd <- readr::read_csv(
  file.path(base_dir, "Combined Analysis_WHO I_Final/pd_Combined Analysis_WHO I.csv"),
  show_col_types = FALSE)

pd <- as.data.frame(pd)
pd <- apply_pd_harmonization(pd)

###-----Harmonize beta_combat and pd (when starting from saved files)-----###
# Small helper to clean IDs
norm_id <- function(x) {
  x <- gsub("\u00A0", " ", as.character(x))  # replace NBSP
  trimws(x)                                  # trim outer spaces
}

#Make sure pd has Sample_Name and clean it
stopifnot("Sample_Name" %in% names(pd))
pd$Sample_Name <- norm_id(pd$Sample_Name)
rownames(pd)   <- pd$Sample_Name

#Clean beta_combat column names
colnames(beta_combat) <- norm_id(colnames(beta_combat))

#Restrict to the intersection of sample IDs, in the same order
common_ids <- intersect(colnames(beta_combat), rownames(pd))
if (length(common_ids) == 0L) {
  stop("No overlapping sample IDs between beta_combat columns and pd rows.")
}

if (length(common_ids) < ncol(beta_combat)) {
  warning("Dropping ", ncol(beta_combat) - length(common_ids),
          " samples from beta_combat with no matching row in pd.")
}
if (length(common_ids) < nrow(pd)) {
  warning("Dropping ", nrow(pd) - length(common_ids),
          " rows from pd with no matching column in beta_combat.")
}

beta_combat <- beta_combat[, common_ids, drop = FALSE]
pd          <- pd[common_ids, , drop = FALSE]

stopifnot(identical(colnames(beta_combat), rownames(pd)))
cat("Harmonized beta_combat and pd — samples:", length(common_ids), "\n")

##==============================================================================
###-----10. FILTER BETAS TO X MEDIAN ABSOLUTE DEVIATION (MAD) PROBES-----###
###-----Create top 2k, 5k, 10k matrices-----###
##==============================================================================

#Calculate the Median Absolute Deviation (MAD) for each probe
mads <- apply(beta_combat, 1, mad, na.rm = TRUE)  # The '1' tells the 'apply' function to operate on each row

#Order the probes by MAD and select the top 50,000
ordered_probes <- order(mads, decreasing = TRUE)
top_50k_probes <- beta_combat[ordered_probes[1:50000], ]
cat("Final dimensions of the data matrix (Probes x Columns):\n")
print(dim(top_50k_probes))

#Create matrices for top 2k, top 5k, top 10k most variable probes
top_2k_betas <- top_50k_probes[1:2000, ]
top_2k_matrix <- as.matrix(top_2k_betas)

top_5k_betas <- top_50k_probes[1:5000, ]
top_5k_matrix <- as.matrix(top_5k_betas)

top_10k_betas <- top_50k_probes[1:10000, ]
top_10k_matrix <- as.matrix(top_10k_betas)

top_50k_betas <- top_50k_probes[1:50000, ]
top_50k_matrix <- as.matrix(top_50k_betas)

#(Optional): Save  2K, 5k, 10k data matrices
#2k
df_top2k <- as.data.frame(top_2k_matrix)
df_top2k <- cbind(ProbeID = rownames(df_top2k), df_top2k)
write.csv(
  df_top2k,
  file = "top_2k_matrix_Combined_WHO I.csv",
  row.names = FALSE
)
cat("Saved top_2k_matrix to: top_2k_matrix_Combined_WHO I.csv\n")

#5k
df_top5k <- as.data.frame(top_5k_matrix)
df_top5k <- cbind(ProbeID = rownames(df_top5k), df_top5k)
write.csv(
  df_top5k,
  file = "top_5k_matrix_Combined_WHO I.csv",
  row.names = FALSE
)
cat("Saved top_5k_matrix to: top_5k_matrix_Combined_WHO I.csv\n")

#10k
df_top10k <- as.data.frame(top_10k_matrix)
df_top10k <- cbind(ProbeID = rownames(df_top10k), df_top10k)
write.csv(
  df_top10k,
  file = "top_10k_matrix_Combined_WHO I.csv",
  row.names = FALSE
)
cat("Saved top_10k_matrix to: top_10k_matrix_Combined_WHO I.csv\n")


##==============================================================================
####-----11. CONSENSUS CLUSTERING-----####
#Runs for 2k, 5k, 10k MAD automatically#
#Paste and run this whole block together. It is memory-intensive
##==============================================================================

mvp_list <- list(`5k_WHOI` = top_5k_matrix)
chosen_k_value <- 2

###-----COHORT FACTOR + COLORS-----###
pd <- pd %>% mutate(Cohort = factor(as.character(Cohort)))

COHORT_COLORS <- c(
  DMPA       = "#582C83",
  Baylor     = "grey95",
  Heidelberg = "grey85"
)

#Helper: run one CC block and save outputs
run_cc_block <- function(
    mat, label, base_dir, pd,
    chosen_k,                        # K specified manually
    maxK        = 5,
    reps        = 1000,
    pItem       = 0.8,
    pFeature    = 1,
    clusterAlg  = "hc",
    distance    = "pearson",
    seed        = 12345,
    save_pdf    = FALSE,
    
    #Color palettes
    heatmap_colors = c("navy", "white", "firebrick3"),   # β-value gradient
    cluster_colors = c("red1", "#1F77B4", "forestgreen",
                       "#E69F00", "#CC79A7", "#56B4E9",
                       "#F0E442", "#0072B2"),            # cluster hues
    cohort_colors = c(DMPA="#582C83", Baylor="grey95", Heidelberg="grey85")


) {
  
message("\n===== Running Consensus Clustering: ", label, " =====")
  
#Sanity checks
stopifnot(is.matrix(mat) || is.data.frame(mat))
stopifnot(is.data.frame(pd), "Cohort" %in% colnames(pd))
stopifnot(all(colnames(mat) %in% rownames(pd)))

#Create results folder
run_dir <- file.path(base_dir, paste0("CC_", label))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

#Save plots as high resolution
#Keep colors true for Illustrator, more robust raster device
grDevices::pdf.options(colormodel = "srgb", useDingbats = FALSE)
#options(bitmapType = "cairo")

###-----Run ConsensusClusterPlus-----###
#Run CC
set.seed(seed)
  res <- ConsensusClusterPlus(
    as.matrix(mat),
    maxK       = maxK,
    reps       = reps,
    pItem      = pItem,
    pFeature   = pFeature,
    title      = run_dir,
    clusterAlg = clusterAlg,
    distance   = distance,
    seed       = seed,
    plot       = "png"
  )
  
  
###-----Extract chosen cluster assignments-----###
if (is.null(res[[chosen_k]])) {
  stop(sprintf("chosen_k=%d not available (maxK=%d).", chosen_k, maxK))
}

# Force the cluster factor to have levels 1..K in order (preserve names!)
clu_vec <- res[[chosen_k]]$consensusClass
nm <- names(clu_vec)                         # keep names
clu <- factor(as.integer(clu_vec),
              levels = seq_len(chosen_k),
              labels = as.character(seq_len(chosen_k)))
if (!is.null(nm)) {
  names(clu) <- nm                           # restore names
} else {
  names(clu) <- colnames(mat)                # fallback if CCP returned no names
}
  
###-----Save cluster assignments-----###
  write.csv(
    data.frame(Sample_Name = names(clu), Cluster = clu),
    file = file.path(run_dir, sprintf("CC_Assignments_%s_k%d.csv", label, chosen_k)),
    row.names = FALSE
  )
  
###-----CC heatmap setup-----###

#Heatmap continuous colors
hm_cols <- colorRampPalette(heatmap_colors)(255)
  
#Cluster colors: extend if K > length provided
  K <- nlevels(clu)
  if (length(cluster_colors) < K) {
    extra <- RColorBrewer::brewer.pal(max(3, K), "Set3")
    cluster_colors <- c(cluster_colors, setdiff(extra, cluster_colors))
    if (length(cluster_colors) < K) {
      # still short — recycle safely (last resort)
      cluster_colors <- rep(cluster_colors, length.out = K)
    }
    message("Note: extended cluster palette to cover K = ", K)
  }
  clu_cols <- setNames(cluster_colors[seq_len(K)], levels(clu))

  
  sg <- factor(pd[colnames(mat), "Cohort"])
  sg <- droplevels(sg)
  
  cohort_lvls <- levels(droplevels(sg))
  if (!all(cohort_lvls %in% names(cohort_colors))) {
    missing <- setdiff(cohort_lvls, names(cohort_colors))
    fallback <- setNames(
      RColorBrewer::brewer.pal(max(3, length(cohort_lvls)), "Set2"),
      cohort_lvls
    )
  #Start with fallback, then overwrite with user-specified where present
    cohort_cols <- fallback
    overlap <- intersect(names(cohort_cols), names(cohort_colors))
    cohort_cols[overlap] <- cohort_colors[overlap]
    if (length(missing)) {
      warning("No explicit color set for cohort(s): ",
              paste(missing, collapse=", "),
              ". Using fallback colors.")
    }
  } else {
    cohort_cols <- cohort_colors[cohort_lvls]
    names(cohort_cols) <- cohort_lvls
  }
  
  #Column annotations
  ha <- HeatmapAnnotation(
    df = data.frame(
      Cohort  = sg,
      Cluster = clu[colnames(mat)]
    ),
    col = list(Cohort = cohort_cols, Cluster = clu_cols),
    annotation_name_gp = gpar(fontsize = 9),
    simple_anno_size   = unit(3.5, "mm")
  )
  
###-----Create consensus cluster heatmap-----###
  ht <- Heatmap(
    as.matrix(mat),
    name              = "Beta",
    col               = hm_cols,
    column_split      = clu,
    cluster_column_slices = FALSE,
    top_annotation    = ha,
    show_row_names    = FALSE,
    show_column_names = FALSE,
    cluster_rows      = TRUE,
    cluster_columns   = TRUE,
    use_raster        = TRUE,
    raster_quality    = 0.8,
    raster_device     = "png",
    heatmap_legend_param = list(
      title = "Beta",
      title_gp   = gpar(fontsize = 9),
      labels_gp  = gpar(fontsize = 8)
    )
  )
  
  # ----- Save heatmap (robust) -----
  png(
    file.path(run_dir, sprintf("Heatmap_%s_k%d.png", label, chosen_k)),
    width = 2200, height = 1600, res = 300
  )
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  }

###-----Loop consensus cluster script for 2k / 5k / 10k MVPs-----###
#Creates and saves all clustering figures and heatmaps for 2k, 5k, 10k probes

#Select matrices to run
mvp_list <- list(
  `2k_WHOI`  = top_2k_matrix,
  `5k_WHOI`  = top_5k_matrix,
  `10k_WHOI` = top_10k_matrix
)

#Select K (same K for all)
chosen_k_value <- 2

#Run Conensus Clustering
cc_runs <- lapply(names(mvp_list), function(lbl) {
  tryCatch({
    run_cc_block(
      mat       = mvp_list[[lbl]],
      label     = lbl,
      base_dir  = working_directory,
      pd        = pd,
      chosen_k  = chosen_k_value
    )
  }, error = function(e) {
    message(sprintf("❌ %s failed: %s", lbl, e$message))
    NULL
  })
})
names(cc_runs) <- names(mvp_list)

message("\n🎯 All clustering runs attempted.")

# ---- REPLAY CONSENSUS CLUSTER PLOTS IN PLOTS PANE ----
replay_cc_pngs <- function(label, base_dir = working_directory) {
  
  # Re-create the folder name used inside run_cc_block()
  run_dir <- file.path(base_dir, paste0("CC_", label))
  
  if (!dir.exists(run_dir)) {
    stop("run_dir does not exist: ", run_dir)
  }
  
  if (!requireNamespace("png", quietly = TRUE)) {
    install.packages("png")
  }
  library(png)
  library(grid)
  
  png_files <- list.files(run_dir, pattern = "\\.png$", full.names = TRUE)
  
  if (!length(png_files)) {
    message("No PNG files found in: ", run_dir)
    return(invisible(NULL))
  }
  
  message("Replaying all PNGs from: ", run_dir)
  
  for (f in png_files) {
    message("Showing: ", basename(f))
    img <- png::readPNG(f)
    grid::grid.newpage()
    grid::grid.raster(img)
  }
  
  invisible(png_files)
}

replay_cc_pngs("2k_WHOI")
replay_cc_pngs("5k_WHOI")
replay_cc_pngs("10k_WHOI")

##==============================================================================
###-----12. CONCORDANCE PLOTS-----###
#Creates Alluvial plots showing cluster crossover when running across 2k/5k/10k
#Uses chosen k set in prior block 
##==============================================================================

  library(readr)
  library(dplyr)       # brings in %>% pipe and verbs
  library(tibble)
  library(ggplot2)     # ggplot()
  library(ggalluvial)
  
  #Set output directory for Alluvial plots
  concord_dir <- file.path(working_directory, "Concordance_K2")
  dir.create(concord_dir, showWarnings = FALSE, recursive = TRUE)
  
  ###-----Helpers-----###
  #Align K=2 labels (preserve names)
  align_k2_to_ref <- function(ref, other) {
    common <- intersect(names(ref), names(other))
    r <- as.integer(ref[common])
    o <- as.integer(other[common])
    
    agree_identity <- sum(r == o, na.rm = TRUE)
    agree_swap     <- sum(r == (3 - o), na.rm = TRUE)
    
    if (agree_swap > agree_identity) {
      nm <- names(other)  # preserve sample names
      other <- factor(3 - as.integer(other), levels = c(1,2), labels = levels(ref))
      names(other) <- nm
    }
    return(other)
  }
  
  ###-----Read assignment CSVs-----###
  #Robust assignment reader 
  read_assign3 <- function(label, k = 2) {
    path <- file.path(working_directory, paste0("CC_", label),
                      sprintf("CC_Assignments_%s_k%d.csv", label, k))
    if (!file.exists(path)) stop("Missing assignments file: ", path)
    df <- readr::read_csv(path, show_col_types = FALSE)
    
    #Normalized column names
    cn <- names(df)
    cn_norm <- tolower(gsub("[^a-z0-9]+", "_", trimws(cn)))
    
    #Locate CLUSTER column
    #By name…
    cl_cands <- which(cn_norm %in% c("cluster","consensusclass","consensus_class","consensus_cluster","k2_cluster"))
    #Or by values being only {1,2} (ignoring NA)
    if (!length(cl_cands)) {
      cl_cands <- which(vapply(df, function(x) {
        ux <- unique(na.omit(as.character(x)))
        length(ux) > 0 && all(ux %in% c("1","2"))
      }, logical(1)))
    }
    stopifnot(length(cl_cands) >= 1)
    cl_idx <- cl_cands[1]
    
    #Locate SAMPLE column
    #Heuristics by header name first
    sn_name_cands <- which(cn_norm %in% c(
      "sample_name","samplename","sample","sampleid","sample_id",
      "id","name","gsm","basename","basenameleaf","colnames"
    ))
    
    #Heuristics by values (look like GSM… or Sentrix position tail …_R##C##)
    looks_like_sample <- function(v) {
      v <- as.character(v)
      frac_non_na <- mean(!is.na(v))
      uniq        <- length(unique(na.omit(v)))
      has_sentrix <- any(grepl("_R\\d\\dC\\d\\d$", v))            # …_R05C01
      has_gsm     <- any(grepl("^GSM\\d+", v, ignore.case = TRUE)) # GSM5702860…
      (frac_non_na > 0.9 && uniq > 10 && (has_sentrix || has_gsm))
    }
    sn_value_cands <- which(vapply(df, looks_like_sample, logical(1)))
    
    #Typical rowname-dump columns (empty header, …1 / X1 / ï..Sample_Name)
    typical_rowname_cols <- which(cn %in% c("", "...1", "X1", "\u00EF..\u00BB\u00BFSample_Name", "ï..Sample_Name"))
    
    #Combine candidates, exclude the cluster column
    sn_cands <- setdiff(unique(c(sn_name_cands, sn_value_cands, typical_rowname_cols)), cl_idx)
    if (!length(sn_cands)) {
      # Last resort: first non-cluster texty column
      others <- setdiff(seq_along(cn), cl_idx)
      is_texty <- vapply(df[others], function(x) any(grepl("[A-Za-z0-9]", as.character(x))), logical(1))
      sn_cands <- others[which(is_texty)]
    }
    stopifnot(length(sn_cands) >= 1)
    
    #Choose the best sample column: highest uniqueness
    uniq_counts <- vapply(df[sn_cands], function(x) length(unique(na.omit(as.character(x)))), integer(1))
    sn_idx <- sn_cands[ which.max(uniq_counts) ]
    
    out <- data.frame(
      Sample_Name = as.character(df[[sn_idx]]),
      Cluster     = as.factor(df[[cl_idx]]),
      stringsAsFactors = FALSE
    )
    out <- subset(out, !is.na(Sample_Name) & nzchar(Sample_Name))
    
    #Guardrail: ensure we truly have many unique sample names
    if (length(unique(out$Sample_Name)) < nrow(out) / 2) {
      stop(sprintf("Sample column detection failed for %s: only %d unique of %d rows.\nColumns were: %s",
                   basename(path), length(unique(out$Sample_Name)), nrow(out), paste(names(df), collapse=", ")))
    }
    
    #Coerce cluster levels to 1/2
    out$Cluster <- factor(as.integer(as.character(out$Cluster)), levels = c(1,2))
    out
  }
  
  ###-----Re-import and verify cluster assignment outputs (2k/5k/10k MVP)-----###
  a2  <- read_assign3("2k_WHOI", 2)
  a5  <- read_assign3("5k_WHOI", 2)
  a10 <- read_assign3("10k_WHOI", 2)
  
  message("Rows: 2k=", nrow(a2), " | 5k=", nrow(a5), " | 10k=", nrow(a10))
  message("Unique samples: 2k=", length(unique(a2$Sample_Name)),
          " | 5k=", length(unique(a5$Sample_Name)),
          " | 10k=", length(unique(a10$Sample_Name)))
  
  ###-----Build named vectors + align labels-----###
  cl2  <- stats::setNames(a2$Cluster,  a2$Sample_Name)
  cl5  <- stats::setNames(a5$Cluster,  a5$Sample_Name)
  cl10 <- stats::setNames(a10$Cluster, a10$Sample_Name)
  
  cl5_aligned  <- align_k2_to_ref(cl2,  cl5)
  cl10_aligned <- align_k2_to_ref(cl2, cl10)
  
  ###-----Rebuild runs_long-----###
  runs_long <- dplyr::bind_rows(
    tibble::tibble(Sample_Name = names(cl2),         Run = "2k",  Cluster = as.integer(cl2)),
    tibble::tibble(Sample_Name = names(cl5_aligned),  Run = "5k",  Cluster = as.integer(cl5_aligned)),
    tibble::tibble(Sample_Name = names(cl10_aligned), Run = "10k", Cluster = as.integer(cl10_aligned))
  ) |>
    dplyr::mutate(Run = factor(Run, levels = c("2k","5k","10k")),
                  Cluster = factor(Cluster, levels = c(1,2)))
  
  #Sanity check: should include all samples
  print(table(runs_long$Run))
  print(tapply(runs_long$Sample_Name, runs_long$Run, function(x) length(unique(x))))
  
  
  ###-----Attach Cohort information to flow and cluster data-----###
  #Includes helpers
  pd_csv_path <- file.path(working_directory, "pd_Combined_WHO I.csv")
  pd_map <- tryCatch(readr::read_csv(pd_csv_path, show_col_types = FALSE), error = function(e) NULL)
  
  if (!is.null(pd_map)) {
    nm <- names(pd_map)
    if (!"Sample_Name" %in% nm) {
      if (!anyDuplicated(pd_map[[1]]) && !grepl("^V\\d+$", names(pd_map)[1])) {
        names(pd_map)[1] <- "Sample_Name"
      } else {
        pd_map$Sample_Name <- pd_map[[1]]
      }
    }
    stopifnot("Cohort" %in% names(pd_map))
    pd_map <- dplyr::transmute(pd_map,
                               Sample_Name      = as.character(.data$Sample_Name),
                               Cohort           = as.character(.data$Cohort),
                               Sentrix_ID       = dplyr::coalesce(as.character(.data$Sentrix_ID), NA_character_),
                               Sentrix_Position = dplyr::coalesce(as.character(.data$Sentrix_Position), NA_character_))
  } else {
    stopifnot(exists("pd"))
    tmp <- as.data.frame(pd, stringsAsFactors = FALSE)
    if (!"Sample_Name" %in% names(tmp)) tmp$Sample_Name <- rownames(tmp)
    stopifnot("Cohort" %in% names(tmp))
    pd_map <- dplyr::transmute(tmp,
                               Sample_Name      = as.character(.data$Sample_Name),
                               Cohort           = as.character(.data$Cohort),
                               Sentrix_ID       = dplyr::coalesce(as.character(.data$Sentrix_ID), NA_character_),
                               Sentrix_Position = dplyr::coalesce(as.character(.data$Sentrix_Position), NA_character_))
  }
  
  normalize_pos <- function(x){
    x <- toupper(as.character(x))
    x <- sub("^R([0-9])C", "R0\\1C", x)
    x <- sub("C([0-9])$", "C0\\1", x)
    x
  }
  extract_sidpos <- function(x){
    sid <- sub("_.*$", "", x)
    pos <- sub("^[^_]*_", "", x)
    tibble::tibble(Sentrix_ID = sid, Sentrix_Position = normalize_pos(pos))
  }
  
  runs_annot <- dplyr::left_join(
    runs_long,
    dplyr::select(pd_map, Sample_Name, Cohort),
    by = "Sample_Name"
  )
  
  if (anyNA(runs_annot$Cohort)) {
    need <- which(is.na(runs_annot$Cohort))
    sidpos <- extract_sidpos(runs_annot$Sample_Name[need])
    sidpos$ix <- need
    pd_sidpos <- dplyr::filter(pd_map, !is.na(Sentrix_ID), !is.na(Sentrix_Position)) |>
      dplyr::mutate(Sentrix_Position = normalize_pos(Sentrix_Position)) |>
      dplyr::select(Sentrix_ID, Sentrix_Position, Cohort)
    joined <- dplyr::left_join(sidpos, pd_sidpos, by = c("Sentrix_ID","Sentrix_Position"))
    if (any(!is.na(joined$Cohort))) {
      runs_annot$Cohort[joined$ix] <- joined$Cohort
    }
  }
  
  if (anyNA(runs_annot$Cohort)) {
    bad <- unique(runs_annot$Sample_Name[is.na(runs_annot$Cohort)])
    warning("Could not attach Cohort for ", length(bad), " sample(s). Examples: ",
            paste(utils::head(bad, 5), collapse = ", "))
    runs_annot$Cohort[is.na(runs_annot$Cohort)] <- "Unknown"
  }
  
  ###-----Generate concordance summary tables-----###
  concord_tbl <- runs_annot |>
    dplyr::count(Sample_Name, Cluster, name = "count_in_cluster") |>
    dplyr::group_by(Sample_Name) |>
    dplyr::mutate(
      majority_count   = max(count_in_cluster, na.rm = TRUE),
      majority_cluster = Cluster[which.max(count_in_cluster)]
    ) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      runs_total = 3,
      percent_concordance = 100 * majority_count / runs_total
    ) |>
    dplyr::select(Sample_Name, majority_cluster, percent_concordance)
  
  cross_flags <- runs_annot |>
    dplyr::group_by(Sample_Name) |>
    dplyr::summarize(num_unique = dplyr::n_distinct(Cluster), .groups = "drop") |>
    dplyr::mutate(crossover = num_unique > 1) |>
    dplyr::select(Sample_Name, crossover)
  
  sample_summary <- runs_annot |>
    dplyr::distinct(Sample_Name, Cohort) |>
    dplyr::left_join(concord_tbl, by = "Sample_Name") |>
    dplyr::left_join(cross_flags, by = "Sample_Name") |>
    dplyr::arrange(Cohort, dplyr::desc(percent_concordance), Sample_Name)
  
  readr::write_csv(sample_summary, file.path(concord_dir, "K2_concordance_per_sample.csv"))
  
  cohort_stability <- runs_annot |>
    dplyr::group_by(Sample_Name, Cohort) |>
    dplyr::summarize(num_unique = dplyr::n_distinct(Cluster), .groups = "drop_last") |>
    dplyr::mutate(crossover = num_unique > 1) |>
    dplyr::count(Cohort, crossover, name = "n") |>
    dplyr::group_by(Cohort) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup()
  
  #Save concordance summary stables
  readr::write_csv(cohort_stability, file.path(concord_dir, "Cohort_stability_summary.csv"))
  
  ###-----Generate Alluvium plot for concordance between MVP runs-----###
  #User knobs - EDIT THESE
  working_directory <- file.path(base_dir, "Combined Analysis_WHO I_Final")
  cluster_cols <- c("1"="red1","2"="#1F77B4")  #Colors to match heatmap in prior code 
  legend_xy     <- c(0.845, 0.53)                 #Legend in white space in 3rd column
  bar_width     <- 0.90
  gap_expand    <- c(0.006, 0)
  
  concord_dir <- file.path(working_directory, "Concordance_K2")
  dir.create(concord_dir, recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf.options(colormodel = "srgb", useDingbats = FALSE)
  options(bitmapType = "cairo")
  
  #Helpers
  norm_ws <- function(x) { x <- gsub("\u00A0"," ", x); x <- gsub("[[:space:]]+"," ", x); trimws(x) }
  norm_pos <- function(x){
    x <- toupper(as.character(x))
    x <- sub("^R([0-9])C", "R0\\1C", x)
    x <- sub("C([0-9])$", "C0\\1", x)
    x
  }
  sid_from_name <- function(s) sub("_.*$", "", s)
  pos_from_name <- function(s) norm_pos(sub("^[^_]*_", "", s))
  sidpos_key    <- function(sid, pos) paste0(sid, "_", pos)
  
  align_k2_to_ref <- function(ref, other){
    common <- intersect(names(ref), names(other))
    r <- as.integer(ref[common]); o <- as.integer(other[common])
    if (sum(r == (3 - o), na.rm = TRUE) > sum(r == o, na.rm = TRUE)) {
      nm <- names(other)
      other <- factor(3 - as.integer(other), levels = c(1,2), labels = levels(ref))
      names(other) <- nm
    }
    other
  }
  
  ###-----Build master_pd with a reliable join_key-----###
  pd_raw <- readr::read_csv(file.path(working_directory, "pd_Combined_WHO I.csv"),
                            show_col_types = FALSE)
  if (!"Sample_Name" %in% names(pd_raw)) names(pd_raw)[1] <- "Sample_Name"
  
  pd <- pd_raw %>%
    mutate(
      Sample_Name = as.character(Sample_Name),
      Cohort = norm_ws(Cohort),
      Cohort = dplyr::case_when(
        Cohort %in% c("baylor","BAYLEY","Baylor cohort") ~ "Baylor",
        Cohort %in% c("heidelberg","Heidelberg cohort","HD") ~ "Heidelberg",
        Cohort %in% c("dmpa","DMPA cohort") ~ "DMPA",
        TRUE ~ Cohort
      ),
      Sentrix_Position = norm_pos(Sentrix_Position),
      # primary key: if a canonical Sample_Name exists, use it;
      # otherwise (Baylor/Heidelberg rows that only have SID/POS), synthesize one
      join_key = dplyr::if_else(
        !is.na(Sample_Name) & nzchar(Sample_Name),
        Sample_Name,
        sidpos_key(Sentrix_ID, Sentrix_Position)
      ),
      sidpos_key = sidpos_key(Sentrix_ID, Sentrix_Position)
    ) %>%
    filter(!is.na(Cohort)) %>%
    # keep *one* row per join_key; prefer rows that actually have a Sample_Name
    arrange(desc(nzchar(Sample_Name))) %>%
    distinct(join_key, .keep_all = TRUE)
  
  
  #Robust assignment reader using same keys
  read_assign_robust <- function(label, k = 2){
    path <- file.path(working_directory, paste0("CC_", label),
                      sprintf("CC_Assignments_%s_k%d.csv", label, k))
    if (!file.exists(path)) stop("Missing assignments file: ", path)
    df <- readr::read_csv(path, show_col_types = FALSE)
    cn <- names(df); cn_norm <- tolower(gsub("[^a-z0-9]+","_", trimws(cn)))
    
    # detect cluster column
    cl_idx <- which(cn_norm %in% c("cluster","consensusclass","consensus_class","consensus_cluster","k2_cluster"))
    if (!length(cl_idx)) {
      cl_idx <- which(vapply(df, function(x){
        ux <- unique(na.omit(as.character(x))); length(ux)>0 && all(ux %in% as.character(1:k))
      }, logical(1)))
    }
    stopifnot(length(cl_idx)>=1); cl_idx <- cl_idx[1]
    
    # detect sample-name column
    sn_idx <- which(cn_norm %in% c("sample_name","samplename","sample","sample_id","id","name","gsm","basename","basenameleaf","colnames"))
    if (!length(sn_idx)) sn_idx <- 1  # fall back to first col if needed
    stopifnot(length(sn_idx)>=1); sn_idx <- sn_idx[1]
    
    raw_names <- as.character(df[[sn_idx]])
    cluster   <- factor(as.integer(as.character(df[[cl_idx]])), levels = 1:k)
    
    tibble::tibble(
      file_name  = raw_names,
      join_key   = ifelse(grepl("_R\\d\\dC\\d\\d$", raw_names), raw_names, NA_character_),
      sidpos_key = sidpos_key(sid_from_name(raw_names), pos_from_name(raw_names)),
      Cluster    = cluster
    ) %>%
      # primary join by join_key (exact Sample_Name or SID_POS string)
      left_join(pd %>% select(join_key, Sample_Name, Cohort), by = "join_key") %>%
      # fill via SID/POS if Cohort still NA
      left_join(pd %>% select(sidpos_key, Sample_Name2 = Sample_Name, Cohort2 = Cohort),
                by = "sidpos_key") %>%
      mutate(
        Sample_Name = dplyr::coalesce(Sample_Name, Sample_Name2, join_key, file_name),
        Cohort      = dplyr::coalesce(Cohort, Cohort2)
      ) %>%
      select(Sample_Name, Cohort, Cluster) %>%
      filter(!is.na(Cohort)) %>%                    # drop rows we cannot place
      distinct(Sample_Name, .keep_all = TRUE)
  }
  
  #Read & align assignments (k=2)
  a2  <- read_assign_robust("2k_WHOI",  2)
  a5  <- read_assign_robust("5k_WHOI",  2)
  a10 <- read_assign_robust("10k_WHOI", 2)
  
  #Sanity check
  print(table(a2$Cohort)); print(table(a5$Cohort)); print(table(a10$Cohort))
  
  cl2  <- stats::setNames(a2$Cluster,  a2$Sample_Name)
  cl5  <- stats::setNames(a5$Cluster,  a5$Sample_Name)
  cl10 <- stats::setNames(a10$Cluster, a10$Sample_Name)
  
  #Align labels of 5k/10k to 2k
  cl5  <- align_k2_to_ref(cl2,  cl5)
  cl10 <- align_k2_to_ref(cl2, cl10)
  
  #Long format, one row per Sample × Run
  runs_long <- bind_rows(
    tibble::tibble(Sample_Name = names(cl2),  Run="2k",  Cluster = as.integer(cl2)),
    tibble::tibble(Sample_Name = names(cl5),  Run="5k",  Cluster = as.integer(cl5)),
    tibble::tibble(Sample_Name = names(cl10), Run="10k", Cluster = as.integer(cl10))
  ) %>%
    distinct(Sample_Name, Run, .keep_all = TRUE) %>%         # <- remove dup rows (root of double bars)
    mutate(Run = factor(Run, levels = c("2k","5k","10k")),
           Cluster = factor(Cluster, levels = c(1,2)))
  
  #Attach cohort from pd (using same join_key logic)
  runs_annot <- runs_annot %>%
    mutate(
      Cohort = factor(as.character(Cohort), levels = c("Baylor","Heidelberg","DMPA"))
    )
  
  #Keep only samples present in all 3 runs (needed for continuous ribbons)
  flows_df <- runs_annot %>%
    group_by(Sample_Name, Cohort) %>%
    filter(dplyr::n_distinct(Run) == 3) %>%
    ungroup() %>%
    distinct(Sample_Name, Cohort, Run, .keep_all = TRUE)
  
  #Same y axis across groups 
  ymax_global <- flows_df %>%
    count(Cohort, Run) %>%
    summarise(max_n = max(n), .groups = "drop") %>%
    pull(max_n)
  
  ###-----Generate Alluvial plot-----###
  bar_width  <- 0.70          #Edit bar width
  moat       <- 0.02        #Moat - extra width for the white mask over flows (tune 0.04–0.08)
  gap_expand <- c(0.02, 0.02) #Extra horizontal space 
  
  p_alluvial <- ggplot(
    flows_df,
    aes(x = Run, stratum = Cluster, alluvium = Sample_Name, y = 1, fill = Cluster)
  ) +
    ggalluvial::geom_flow(knot.pos = 0.33, size = 0.10, alpha = 0.55, colour = "grey80") +
    ggalluvial::geom_stratum(width = bar_width + moat, fill = "white", colour = NA, alpha = 1) +
    ggalluvial::geom_stratum(width = bar_width, colour = "grey30", size = 0.28, alpha = 0.90) +
    geom_text(stat = "stratum", aes(label = after_stat(n)), vjust = 0.5, size = 3,
              colour = "white", fontface = "bold") +
    
    # ⬇️ move Baylor / Heidelberg / DMPA to the bottom
    facet_wrap(~ Cohort, nrow = 1, drop = FALSE, strip.position = "bottom") +
    
    scale_x_discrete(expand = gap_expand) +
    scale_y_continuous(limits = c(0, ymax_global),
                       breaks = seq(0, ymax_global, 10),
                       expand = c(0, 0)) +
    scale_fill_manual(values = cluster_cols, name = "Cluster") +
    labs(
      title    = "Cluster Flow Across MVP Sizes (K = 2)",
      subtitle = "2k → 5k → 10k; counts shown inside bars",
      y = "Sample count", x = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      
      # ⬇️ put strips outside so they sit below the axis tick labels
      strip.placement  = "outside",
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", margin = ggplot2::margin(t = 6)),
      
      panel.spacing.x  = grid::unit(10, "pt"),
      axis.line.x      = element_line(color = "black", linewidth = 0.4),
      axis.line.y      = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      axis.text.x      = element_text(margin = ggplot2::margin(t = 4)),
      
      legend.position  = c(0.845, 0.53),
      legend.direction = "horizontal",
      legend.justification = c("center","center"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key       = element_rect(fill = "transparent", colour = NA)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  
  
  #Save Alluvial plot
  ggsave(file.path(concord_dir, "Alluvial_K2_ribbons.png"),
         p_alluvial, width = 6, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(concord_dir, "Alluvial_K2_ribbons.pdf"),
         p_alluvial, width = 6, height = 6, device = cairo_pdf, bg = "white")
  
  #Print Alluvial plot
  print(p_alluvial)
  
  
##========================================================
###--DIMENSIONALITY REDUCTION PLOTS - SET GLOBAL THEME---###
##========================================================
library(ggplot2)

###-----User inputs-----###

X             <- as.matrix(top_5k_matrix)   # choose: top_2k_matrix / top_5k_matrix / top_10k_matrix
chosen_k      <- 2                          # consensus cluster K
n_pcs_for_tsne <- 30                        # number of PCs to use for t-SNE on PCA
seed_global   <- 42

#make pd rownames match X colnames (robust) ----
pd <- readr::read_csv(
  file.path(base_dir, "Combined Analysis_WHO I_Final/pd_Combined Analysis_All Grades.csv"),
  show_col_types = FALSE
)
pd <- as.data.frame(pd)

norm <- function(s) { s <- gsub("\u00A0"," ", as.character(s)); trimws(s) }
stopifnot("Sample_Name" %in% names(pd))
pd$Sample_Name <- norm(pd$Sample_Name)
rownames(pd)   <- pd$Sample_Name

# normalize X colnames and align
colnames(X) <- norm(colnames(X))

missing_in_pd  <- setdiff(colnames(X), rownames(pd))
missing_in_X   <- setdiff(rownames(pd), colnames(X))

if (length(missing_in_pd)) {
  warning("In X but not in pd: ", length(missing_in_pd),
          " (e.g., ", paste(utils::head(missing_in_pd, 5), collapse=", "), ")")
}
if (length(missing_in_X)) {
  message("In pd but not in X: ", length(missing_in_X),
          " (e.g., ", paste(utils::head(missing_in_X, 5), collapse=", "), ")")
}

keep <- intersect(colnames(X), rownames(pd))
X    <- X[, keep, drop = FALSE]
pd   <- pd[keep, , drop = FALSE]

stopifnot(all(colnames(X) %in% rownames(pd)))

#Theme settings
theme_pub <- function(base_size = 14, legend_x = 0.92, legend_y = 0.92) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = c(legend_x, legend_y),
      legend.justification = c("right", "top"),
      legend.background     = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.key            = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.text   = ggplot2::element_text(size = base_size * 0.9),
      legend.spacing.y = grid::unit(2, "pt"),
      legend.margin = ggplot2::margin(t = 2, b = 2, unit = "pt"),
      axis.text.x   = ggplot2::element_blank(),
      axis.text.y   = ggplot2::element_blank(),
      axis.ticks.x  = ggplot2::element_blank(),
      axis.ticks.y  = ggplot2::element_blank(),
      plot.title    = ggplot2::element_text(size = base_size, face = "bold", hjust = 0.5),
      axis.title.x  = ggplot2::element_text(size = base_size * 0.9),
      axis.title.y  = ggplot2::element_text(size = base_size * 0.9),
      plot.margin   = ggplot2::margin(5, 5, 5, 5, unit = "pt")
    )
}
theme_set(theme_pub())
update_geom_defaults("point", list(size = 1.9, alpha = 0.85, shape = 16))


#Flexible Naming and Helper Functions
infer_matrix_label <- function(X) {
  r <- nrow(X)
  if (r >= 9000) "10k" else if (r >= 4000) "5k" else if (r >= 1500) "2k" else paste0(r, "rows")
}
matrix_label <- infer_matrix_label(X)
plots_root   <- file.path(working_directory, "Plots")
out_dir      <- file.path(plots_root, matrix_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, prefix, matrix_label, chosen_k) {
  outfile <- file.path(out_dir, sprintf("%s_%s_k%d.png", prefix, matrix_label, chosen_k))
  ggsave(outfile, plot = p, width = 7.5, height = 5.5, dpi = 300)
  message("✅ Saved: ", normalizePath(outfile))
}


#Load Cluster Assignments
cc_csv_path <- file.path(
  working_directory,
  paste0("CC_", matrix_label, "_WHOI"),
  sprintf("CC_Assignments_%s_WHOI_k%d.csv", matrix_label, chosen_k)
)
if (!file.exists(cc_csv_path))
  stop("Could not find cluster assignments at: ", cc_csv_path)

clu_tbl <- read_csv(cc_csv_path, show_col_types = FALSE)
clu <- setNames(as.factor(clu_tbl$Cluster), clu_tbl$Sample_Name)

stopifnot(all(colnames(X) %in% rownames(pd)))
pd_plot <- pd[colnames(X), , drop = FALSE]
clu <- droplevels(clu[colnames(X)])

###-----SHARED PLOT LAYERS (PCA/tSNE)-----###
# Cohort labels (DMPA/Baylor/Heidelberg) are already harmonized upstream.
pd_plot$Cohort <- factor(pd_plot$Cohort)

HILITE_COHORT <- "DMPA"     # cohort shown in the legend and given the convex hull
purple        <- "#582C83"  # DMPA signature purple
hilite_all    <- rownames(pd_plot)[pd_plot$Cohort == HILITE_COHORT]

add_shared_layers <- function(df, x_aes, y_aes, title_text, x_lab, y_lab) {
  
  # make sure LegendGroup is a character for comparisons, then factor for legend order
  df$LegendGroup <- as.character(df$LegendGroup)
  df$LegendGroup <- factor(df$LegendGroup, levels = c("Cluster 1", "Cluster 2", HILITE_COHORT))
  
  ggplot(df, aes_string(x = x_aes, y = y_aes, color = "LegendGroup")) +
    
    # non-highlight points  (NOTE: no quotes around HILITE_COHORT)
    geom_point(
      data = subset(df, LegendGroup != HILITE_COHORT),
      size = 2.5, alpha = 0.9
    ) +
    
    # highlight points
    geom_point(
      data = subset(df, LegendGroup == HILITE_COHORT),
      size = 3.0, alpha = 0.95
    ) +
    
    # highlight hull
    ggforce::geom_mark_hull(
      data = subset(df, LegendGroup == HILITE_COHORT),
      aes_string(x = x_aes, y = y_aes),
      concavity = 5, expand = grid::unit(2, "mm"),
      colour = purple, linetype = "dashed",
      size = 0.7, fill = NA, inherit.aes = FALSE
    ) +
    
    scale_color_manual(
      values = c(
        "Cluster 1" = "red1",
        "Cluster 2" = "#1F77B4",
        setNames(purple, HILITE_COHORT)   # <- ensures DMPA is mapped
      ),
      breaks = c("Cluster 1", "Cluster 2", HILITE_COHORT),
      drop = FALSE
    ) +
    
    guides(color = guide_legend(title = NULL, override.aes = list(size = 3))) +
    labs(title = title_text, x = x_lab, y = y_lab)
}

##========================================================
###-----13. t-SNE PLOT-----###
##========================================================

set.seed(seed_global)
R <- suppressWarnings(cor(X, method = "spearman", use = "pairwise.complete.obs"))
D <- as.dist(1 - R)
perp <- max(5, min(30, floor((ncol(X) - 1) / 3)))
tsne_direct <- Rtsne(D, is_distance = TRUE, perplexity = perp,
                     theta = 0.5, max_iter = 1000, verbose = TRUE, pca = FALSE)

tsne_df <- data.frame(
  TSNE1   = tsne_direct$Y[,1],
  TSNE2   = tsne_direct$Y[,2],
  Sample  = colnames(X),
  Cohort  = factor(pd_plot$Cohort),
  Cluster = factor(clu)
) %>%
  mutate(
    IsHILITE    = Sample %in% hilite_all,
    LegendGroup = ifelse(IsHILITE, HILITE_COHORT, paste0("Cluster ", Cluster)),
    LegendGroup = factor(LegendGroup, levels = c("Cluster 1", "Cluster 2", HILITE_COHORT))
  )

p_tsne <- add_shared_layers(tsne_df, "TSNE1", "TSNE2", "t-SNE", "t-SNE 1", "t-SNE 2")
print(p_tsne)
save_plot(p_tsne, "tsne", matrix_label, chosen_k)


##==============================================================================
###-----14. PRINCIPAL COMPONENT ANALYSIS (PCA) PLOT-----###
##==============================================================================

set.seed(seed_global)
mat_pca <- t(X)
pca_res <- prcomp(mat_pca, center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  PC1    = pca_res$x[,1],
  PC2    = pca_res$x[,2],
  Sample = rownames(pca_res$x),
  Cohort = factor(pd_plot$Cohort[rownames(pca_res$x)]),
  Cluster = factor(clu[rownames(pca_res$x)])
) %>%
  mutate(
    IsHILITE    = Sample %in% hilite_all,
    LegendGroup = ifelse(IsHILITE, HILITE_COHORT, paste0("Cluster ", Cluster)),
    LegendGroup = factor(LegendGroup, levels = c("Cluster 1", "Cluster 2", HILITE_COHORT))
  )

p_pca <- add_shared_layers(pca_df, "PC1", "PC2", "PCA", "PC1", "PC2")
print(p_pca)
save_plot(p_pca, "pca", matrix_label, chosen_k)


##==============================================================================
###-----15. ASSIGN DMPA SAMPLES TO BAYLOR GROUPS-----###
#Trains a random forest classifier based on all 110 Baylor samples (all grades)
##==============================================================================

###-----User inputs-----###
set.seed(20251013)

#Baylor IDA file location (your earlier epic_dir)
epic_dir <- file.path(base_dir, "Baylor/GSE189521_RAW")

#Baylor clinical workbook + sheet location
meta_path  <- file.path(base_dir, "Baylor/GSE189521_Clinical_data (Bayley et al).xlsx")
meta_sheet <- "2. Clinical and genomic dataSH"

#Previously saved combined CSVs (from Block 5) earlier in script
working_directory <- file.path(base_dir, "Combined Analysis_WHO I_Final")
combined_pd_csv   <- file.path(working_directory,  "pd_Combined_WHO I.csv")
combined_beta_csv <- file.path(working_directory,  "beta_Combined_WHO I.csv")

#Output folder
ResultsDir <- working_directory
out_dir <- file.path(ResultsDir, "MenG_RF_fromClinical_ALLBaylor")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

###-----Helpers-----###
normalize_pos <- function(x){
  x <- toupper(as.character(x))
  x <- sub("^R([0-9])C", "R0\\1C", x)
  x <- sub("C([0-9])$", "C0\\1", x)
  x
}
basename_from_path <- function(p) {
  leaf <- basename(as.character(p))
  sub("(_Red|_Grn)?\\.idat(\\.gz)?$", "", leaf, ignore.case = TRUE)
}
split_sid_pos <- function(bn) {
  sid <- sub("_.*$", "", bn)
  pos <- normalize_pos(sub("^[^_]*_", "", bn))
  data.frame(Sentrix_ID = sid, Sentrix_Position = pos, stringsAsFactors = FALSE)
}

###-----Read ALL Baylor IDATs-----###
#Note that this runs independent of the code at beginning of script, 
#which reads in Baylor IDATs filtered by WHO grade

message("📂 Reading Baylor IDATs (all WHO grades) from: ", epic_dir)
rg_bay <- minfi::read.metharray.exp(epic_dir, recursive = TRUE, extended = TRUE, force = TRUE)
gr_bay <- minfi::preprocessNoob(rg_bay)
beta_bay_all <- minfi::getBeta(gr_bay)                      # CpGs x samples
pd_bay_all   <- as.data.frame(pData(rg_bay), stringsAsFactors = FALSE)
pd_bay_all$BetaCol <- colnames(beta_bay_all)

#Ensure pd has Sentrix columns (derive from BetaCol if necessary)
if (!("Sentrix_ID" %in% names(pd_bay_all)) || !("Sentrix_Position" %in% names(pd_bay_all))) {
  sidpos_pd <- split_sid_pos(pd_bay_all$BetaCol)
  pd_bay_all$Sentrix_ID       <- sidpos_pd$Sentrix_ID
  pd_bay_all$Sentrix_Position <- sidpos_pd$Sentrix_Position
} else {
  pd_bay_all$Sentrix_Position <- normalize_pos(pd_bay_all$Sentrix_Position)
}

###-----Read Baylor clinical Excel sheet-----###
stopifnot(file.exists(meta_path))
meta0 <- readxl::read_excel(meta_path, sheet = meta_sheet, skip = 2)
meta0 <- as.data.frame(meta0, stringsAsFactors = FALSE)

#Robust MenG extractor (handles 'MenG A', 'A', etc.)
extract_meng_ABC <- function(x){
  x0 <- toupper(trimws(as.character(x)))
  # keep only letters and spaces to simplify search
  x1 <- gsub("[^A-Z]", " ", x0)
  ifelse(grepl("\\bA\\b", x1), "A",
         ifelse(grepl("\\bB\\b", x1), "B",
                ifelse(grepl("\\bC\\b", x1), "C", NA_character_)))
}

if (ncol(meta0) < 3) stop("Clinical sheet has <3 columns; need the third column for MenG.")
MenG_col <- extract_meng_ABC(meta0[[3]])

#Helpers re-used here 
normalize_pos <- function(x){
  x <- toupper(as.character(x))
  x <- sub("^R([0-9])C", "R0\\1C", x)
  x <- sub("C([0-9])$", "C0\\1", x)
  x
}
basename_from_path <- function(p) {
  leaf <- basename(as.character(p))
  sub("(_Red|_Grn)?\\.idat(\\.gz)?$", "", leaf, ignore.case = TRUE)
}
split_sid_pos <- function(bn) {
  sid <- sub("_.*$", "", bn)
  pos <- normalize_pos(sub("^[^_]*_", "", bn))
  data.frame(Sentrix_ID = sid, Sentrix_Position = pos, stringsAsFactors = FALSE)
}

#Try to find IDAT/basename or Sentrix columns in clinical sheet
cn     <- names(meta0)
cn_low <- tolower(gsub("[^a-z0-9]+", " ", cn))

idx_idat <- which(grepl("\\bidat\\b|basename|file|filename|path", cn_low))

#Heuristic scanners for Sentrix
looks_like_sid <- function(v){
  v <- as.character(v)
  sum(grepl("^[0-9]{8,}$", v))  # many digits (8+)
}
looks_like_pos <- function(v){
  v <- toupper(as.character(v))
  sum(grepl("^R0?\\dC0?\\d$", v))
}

#Score every column
sid_scores <- sapply(meta0, looks_like_sid)
pos_scores <- sapply(meta0, looks_like_pos)
sid_idx    <- if (max(sid_scores, na.rm=TRUE) > 5) which.max(sid_scores) else integer(0)
pos_idx    <- if (max(pos_scores, na.rm=TRUE) > 5) which.max(pos_scores) else integer(0)

#Build clinical mapping keys (without filtering rows yet)
if (length(idx_idat)) {
  # Prefer basename / IDAT column
  baseleaf <- basename_from_path(meta0[[ idx_idat[1] ]])
  sidpos   <- split_sid_pos(baseleaf)
  meta_key <- data.frame(
    Sentrix_ID       = sidpos$Sentrix_ID,
    Sentrix_Position = sidpos$Sentrix_Position,
    stringsAsFactors = FALSE
  )
} else if (length(sid_idx) && length(pos_idx)) {
  meta_key <- data.frame(
    Sentrix_ID       = as.character(meta0[[ sid_idx[1] ]]),
    Sentrix_Position = normalize_pos(meta0[[ pos_idx[1] ]]),
    stringsAsFactors = FALSE
  )
} else {
  #last-ditch: try any pair of columns that *together* look like keys
  cand_sid <- which(sid_scores > 0)
  cand_pos <- which(pos_scores > 0)
  if (length(cand_sid) && length(cand_pos)) {
    meta_key <- data.frame(
      Sentrix_ID       = as.character(meta0[[ cand_sid[1] ]]),
      Sentrix_Position = normalize_pos(meta0[[ cand_pos[1] ]]),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Could not find an IDAT/basename column nor Sentrix_ID/Sentrix_Position in the clinical sheet.\n",
         "Columns seen:\n - ", paste(cn, collapse = " | "), "\n",
         "Tip: ensure the sheet has either a file path to the IDATs or explicit Sentrix ID and Position.")
  }
}

#Attach MenG (but do not filter yet. Keep size alignment)
meta_key$MenG <- MenG_col

###-----Diagnostics on clinical keys-----###
cat("Clinical keys detected: ", sum(!is.na(meta_key$Sentrix_ID) & !is.na(meta_key$Sentrix_Position)), "\n")
cat("Clinical MenG A/B/C counts: ",
    "A=", sum(meta_key$MenG=="A", na.rm=TRUE),
    " B=", sum(meta_key$MenG=="B", na.rm=TRUE),
    " C=", sum(meta_key$MenG=="C", na.rm=TRUE), "\n")

###-----Build Baylor pd keys from the IDATs we just read-----###
pd_bay_all$Sentrix_Position <- normalize_pos(pd_bay_all$Sentrix_Position)
key_meta <- paste(meta_key$Sentrix_ID, meta_key$Sentrix_Position, sep = "_")
key_pd   <- paste(pd_bay_all$Sentrix_ID, pd_bay_all$Sentrix_Position, sep = "_")

#Map MenG → Baylor pd
m <- match(key_pd, key_meta)
pd_bay_all$MenG_fromClinical <- meta_key$MenG[m]

#Filter to labeled A/B/C (drop Unknown/NA)
keep_labeled <- !is.na(pd_bay_all$MenG_fromClinical) & pd_bay_all$MenG_fromClinical %in% c("A","B","C")

cat("Mapped MenG onto Baylor betas: ", sum(keep_labeled), " labeled samples.\n")

keep_train <- keep_labeled

if (sum(keep_labeled) < 50) {
  # show a few example keys to debug quickly
  cat("Example Baylor pd keys (first 5):\n")
  print(head(data.frame(key_pd = key_pd, SID = pd_bay_all$Sentrix_ID,
                        POS = pd_bay_all$Sentrix_Position), 5))
  cat("Example clinical keys (first 5):\n")
  print(head(data.frame(key_meta = key_meta, SID = meta_key$Sentrix_ID,
                        POS = meta_key$Sentrix_Position, MenG = meta_key$MenG), 5))
  stop("Too few labeled Baylor samples matched. The clinical sheet likely uses different columns/format for Sentrix or IDAT. See diagnostics above.")
}

#build Baylor training matrix/labels
B_bay_full <- beta_bay_all[, keep_train, drop = FALSE]
y_bay      <- factor(pd_bay_all$MenG_fromClinical[keep_train], levels = c("A","B","C"))

message("Baylor training samples: ", ncol(B_bay_full),
        " | A=", sum(y_bay=="A"), " B=", sum(y_bay=="B"), " C=", sum(y_bay=="C"))

###-----Load combined pd-----###
#Harmonize pd <-> beta names
pd <- read.csv(combined_pd_csv, check.names = FALSE, stringsAsFactors = FALSE)
if (!"Sample_Name" %in% names(pd)) names(pd)[1] <- "Sample_Name"

beta <- tryCatch(
  as.matrix(read.csv(combined_beta_csv, check.names = FALSE, row.names = 1)),
  error = function(e) {
    if (file.exists(paste0(combined_beta_csv, ".gz"))) {
      as.matrix(read.csv(gzfile(paste0(combined_beta_csv, ".gz")), check.names = FALSE, row.names = 1))
    } else stop(e)
  }
)
storage.mode(beta) <- "double"

#Harmonize pd to beta colnames 
looks_like_sidpos <- function(v) grepl("^[0-9A-Z]+_R\\d\\dC\\d\\d$", v)
extract_sidpos <- function(v){
  sid <- sub("_.*$", "", v)
  pos <- normalize_pos(sub("^[^_]*_", "", v))
  data.frame(Sentrix_ID = sid, Sentrix_Position = pos, stringsAsFactors = FALSE)
}
for (nm in c("Sample_Name","Sentrix_ID","Sentrix_Position","Cohort")) {
  if (!nm %in% names(pd)) pd[[nm]] <- NA_character_
  pd[[nm]] <- as.character(pd[[nm]])
}
pd$Sentrix_Position <- normalize_pos(pd$Sentrix_Position)

#Try direct sample name match; else via Sentrix_ID/Position
beta_cols <- colnames(beta)
if (sum(pd$Sample_Name %in% beta_cols) < length(beta_cols) * 0.6 && any(looks_like_sidpos(beta_cols))) {
  sidpos_df <- cbind(data.frame(beta_col = beta_cols, stringsAsFactors = FALSE),
                     extract_sidpos(beta_cols))
  key_pd   <- paste(pd$Sentrix_ID, pd$Sentrix_Position, sep="_")
  key_beta <- paste(sidpos_df$Sentrix_ID, sidpos_df$Sentrix_Position, sep="_")
  m <- match(beta_cols, pd$Sample_Name)  # direct
  pd2 <- pd[m, , drop = FALSE]
  need <- which(is.na(pd2$Sample_Name))
  if (length(need)) {
    m2 <- match(paste(sidpos_df$Sentrix_ID, sidpos_df$Sentrix_Position, sep="_"), key_pd)
    ok <- !is.na(m2)
    pd2[ok, ] <- pd[m2[ok], , drop = FALSE]
    rownames(pd2) <- beta_cols
    pd <- pd2
  }
} else {
  rownames(pd) <- pd$Sample_Name
  pd <- pd[match(beta_cols, rownames(pd)), , drop = FALSE]
}
#Force alignment
rownames(pd) <- colnames(beta)
if (!"Sample_Name_Original" %in% names(pd)) pd$Sample_Name_Original <- pd$Sample_Name
pd$Sample_Name <- colnames(beta)

#DMPA columns for prediction
stopifnot("Cohort" %in% names(pd))
dmpa_cols <- rownames(pd)[pd$Cohort == "DMPA"]
dmpa_cols <- intersect(dmpa_cols, colnames(beta))
if (!length(dmpa_cols)) stop("No DMPA columns found in combined beta.")
B_dmpa_full <- beta[, dmpa_cols, drop = FALSE]

###-----Harmonize loci+ MVP selection-----###
message("Baylor probes: ", nrow(B_bay_full), " | DMPA probes: ", nrow(B_dmpa_full))

#Intersect EPIC v1 (Baylor) with EPIC v2 (DMPA) 
common0 <- intersect(rownames(B_bay_full), rownames(B_dmpa_full))
if (!length(common0)) stop("No probe overlap between Baylor and DMPA matrices. Check preprocessing and rownames (CpG IDs).")
B_bay_i  <- B_bay_full [common0, , drop = FALSE]
B_dmpa_i <- B_dmpa_full[common0, , drop = FALSE]
message("Common probes before sex-filter: ", length(common0))

#Build sex-probe set from BOTH manifests (v1 and v2)
get_sex_ids_safe <- function() {
  sex_ids <- character()
  # EPIC v1 (hg19/38 manifest)
  sex_ids <- tryCatch({
    suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICmanifest))
    ann1 <- minfi::getAnnotation(IlluminaHumanMethylationEPICmanifest)
    c(sex_ids, rownames(ann1)[which(ann1$chr %in% c("chrX","chrY"))])
  }, error = function(e) sex_ids)
  # EPIC v2 (hg38 manifest)
  sex_ids <- tryCatch({
    suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICv2anno.20a1.hg38))
    ann2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
    ann2 <- as.data.frame(ann2)
    c(sex_ids, rownames(ann2)[which(ann2$chr %in% c("chrX","chrY"))])
  }, error = function(e) sex_ids)
  unique(sex_ids)
}
sex_ids_union <- get_sex_ids_safe()

#Apply sex filter ONLY to probes we actually have
auto_mask <- !(rownames(B_bay_i) %in% sex_ids_union)
if (sum(auto_mask) == 0) {
  warning("Sex-chrom filter removed all probes; continuing WITHOUT sex filtering.")
  auto_mask <- rep(TRUE, nrow(B_bay_i))
}
B_bay  <- B_bay_i [auto_mask, , drop = FALSE]
B_dmpa <- B_dmpa_i[auto_mask, , drop = FALSE]
message("Common autosomal probes after sex-filter: ", nrow(B_bay))

#Drop CpGs with lots of NA in Baylor (looser threshold to be safe)
na_frac <- rowMeans(is.na(B_bay))
keep_na <- na_frac <= 0.10
if (sum(keep_na) < 1000) {
  warning("Strict NA filter would leave <1000 probes; relaxing to <=20% NA.")
  keep_na <- rowMeans(is.na(B_bay)) <= 0.20
}
B_bay  <- B_bay [keep_na, , drop = FALSE]
B_dmpa <- B_dmpa[keep_na, , drop = FALSE]
message("Probes after NA-filter: ", nrow(B_bay))

#Most-variable probes (MVPs) computed on Baylor, then z-scale using Baylor stats
if (nrow(B_bay) < 200) stop("Too few probes after filtering (", nrow(B_bay), "). Review preprocessing or relax filters.")

n_mvp <- min(10000, max(2000, nrow(B_bay)))  # aim for up to 10k, but keep at least 2k if possible
vars  <- matrixStats::rowVars(B_bay, na.rm = TRUE)
sel   <- order(vars, decreasing = TRUE)[seq_len(n_mvp)]
B_bay  <- B_bay [sel, , drop = FALSE]
B_dmpa <- B_dmpa[sel, , drop = FALSE]

mu <- rowMeans(B_bay, na.rm = TRUE)
sd <- matrixStats::rowSds(B_bay, na.rm = TRUE); sd[!is.finite(sd) | sd == 0] <- 1

Z_bay  <- t( sweep( sweep(B_bay,  1, mu, "-"), 1, sd, "/") )
Z_dmpa <- t( sweep( sweep(B_dmpa, 1, mu, "-"), 1, sd, "/") )

message("Training features: ", ncol(Z_bay),
        " | Baylor n=", nrow(Z_bay),
        " | DMPA n=", nrow(Z_dmpa))

###-----Train Random Forest model, predict DMPA classification into Baylor groups-----###
tbl <- table(y_bay)
classwt <- as.numeric(median(tbl) / tbl); names(classwt) <- names(tbl)

p <- ncol(Z_bay)
mtry_grid <- unique(pmax(1, floor(c(sqrt(p), p/10, p/5))))

rf_fits <- lapply(mtry_grid, function(mtry) {
  message("Training RF (mtry=", mtry, ", ntree=1000)…")
  randomForest::randomForest(
    x = Z_bay, y = y_bay,
    ntree = 1000, mtry = mtry,
    importance = TRUE, proximity = FALSE, na.action = na.omit,
    classwt = classwt, keep.forest = TRUE
  )
})
oob_err  <- sapply(rf_fits, function(f) tail(f$err.rate[,"OOB"], 1))
best_ix  <- which.min(oob_err); rf_fit <- rf_fits[[best_ix]]
best_mtry <- mtry_grid[best_ix]
message(sprintf("Best mtry = %d | OOB error = %.3f", best_mtry, oob_err[best_ix]))

#Predict DMPA
pred_lab <- predict(rf_fit, Z_dmpa, type = "response")
pred_prb <- predict(rf_fit, Z_dmpa, type = "prob")

# --- your crosswalk (paste here if it isn't defined above this block) ---
id_crosswalk <- data.frame(
  Idat = c(
    "208271520022_R05C01",
    "209429010079_R08C01",
    "209547720159_R02C01",
    "209547720159_R03C01",
    "209547720159_R04C01",
    "209547720159_R05C01",
    "209547720159_R06C01",
    "209547720159_R07C01",
    "209547720159_R08C01",
    "209514950092_R08C01"
  ),
  DMPA_ID = c(1,10,8,9,6,7,4,2,3,5),
  stringsAsFactors = FALSE
)

# helper: normalize Sentrix positions & strip EPICv2 suffixes
.norm_sidpos <- function(v) {
  v <- toupper(as.character(v))
  v <- sub("_(BC|TC)\\d+$", "", v, perl = TRUE)  # EPICv2 suffixes
  v <- sub("^R([0-9])C", "R0\\1C", v)
  v <- sub("C([0-9])$", "C0\\1", v)
  v
}

# match predicted sample names to crosswalk
samples_pred <- rownames(Z_dmpa)                 # these are the DMPA sample names (SID_POS)
key_pred     <- .norm_sidpos(samples_pred)
key_cw       <- .norm_sidpos(id_crosswalk$Idat)

mm        <- match(key_pred, key_cw)
dmpa_num  <- id_crosswalk$DMPA_ID[mm]
dmpa_lab  <- ifelse(is.na(dmpa_num), NA_character_, paste0("DMPA-", dmpa_num))

# confidence = max class prob across A/B/C
Confidence <- apply(pred_prb, 1, max)

# final table (keep "DMPA-ID" exactly)
pred_df <- data.frame(
  Sample         = samples_pred,
  `DMPA-ID`      = dmpa_lab,
  Predicted_MenG = as.character(pred_lab),
  Prob_A         = pred_prb[, "A", drop = TRUE],
  Prob_B         = pred_prb[, "B", drop = TRUE],
  Prob_C         = pred_prb[, "C", drop = TRUE],
  Confidence     = Confidence,
  check.names    = FALSE,
  stringsAsFactors = FALSE
)

# sort by confidence (optional) and save
pred_df <- pred_df[order(-pred_df$Confidence), ]
write.csv(pred_df, file.path(out_dir, "DMPA_MenG_predictions.csv"), row.names = FALSE)
message("✅ Saved with DMPA-ID: ", normalizePath(file.path(out_dir, "DMPA_MenG_predictions.csv")))


##==============================================================================
###-----16. SeSAMe COPY NUMBER ANALYSIS-----###
##==============================================================================

#Idat directory and results directories
idat_dir <- file.path(base_dir, "Beta values/DMPA Idat files")

out_root <- file.path(getwd(), "SeSAMe_CNV_FINAL")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

fig_dir         <- file.path(out_root, "figures")
diag_dir        <- file.path(out_root, "diagnostics")
dir.create(fig_dir,        showWarnings = FALSE, recursive = TRUE)
dir.create(diag_dir,       showWarnings = FALSE, recursive = TRUE)

###-----User inputs/knobs-----###
# Preprocessing toggles
USE_POOBAH          <- TRUE
USE_DYE_BIAS_TINORM <- TRUE

# Calling knobs 
LOSS_MIN_ABS_THR <- 0.25   # log2 threshold for losses (absolute value)
GAIN_MIN_ABS_THR <- 0.30   # log2 threshold for gains (absolute value)
LOSS_MAD_K       <- 0   # adaptive add-on for losses
GAIN_MAD_K       <- 0    # adaptive add-on for gains (threshold is max of ABS_THr and  K x MAD)

# Coverage gate (single knobs; applied to both gain and loss)
SEG_MAG_THR    <- 0.25
COVER_FRAC_ARM <- 0.50
COVER_FRAC_CHR <- 0.50

# Heatmap display-only cutoffs (separate from calling logic - this is graphed separately)
DISPLAY_LOSS_THR <- -0.25
DISPLAY_GAIN_THR <- +0.25

# DMPA crosswalk (IDAT basename → DMPA_ID) for ordered arm heatmap
id_crosswalk <- data.frame(
  Idat = c(
    "208271520022_R05C01",
    "209429010079_R08C01",
    "209547720159_R02C01",
    "209547720159_R03C01",
    "209547720159_R04C01",
    "209547720159_R05C01",
    "209547720159_R06C01",
    "209547720159_R07C01",
    "209547720159_R08C01",
    "209514950092_R08C01"
  ),
  DMPA_ID = c(1,10,8,9,6,7,4,2,3,5),
  stringsAsFactors = FALSE
)

#Install/load
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c(
  "sesame","readr","GenomicRanges","IRanges","GenomeInfoDb","S4Vectors",
  "ComplexHeatmap","circlize","grid","methods","ggplot2","dplyr","tidyr","tibble"
)
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) {
  BiocManager::install(p, ask = FALSE, update = FALSE, quiet = TRUE)
}
suppressPackageStartupMessages({
  library(sesame); library(readr); library(GenomicRanges); library(IRanges)
  library(GenomeInfoDb); library(S4Vectors); library(ComplexHeatmap)
  library(circlize); library(grid); library(methods)
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
})

###-----Helpers-----###
normalize_pos <- function(x){
  x <- toupper(as.character(x))
  x <- sub("^R([0-9])C","R0\\1C",x); x <- sub("C([0-9])$","C0\\1",x); x
}
idat_index_from_dir <- function(idat_dir) {
  f <- list.files(idat_dir, pattern="\\.idat(\\.gz)?$", recursive=TRUE, full.names=TRUE, ignore.case=TRUE)
  stopifnot(length(f) > 0)
  is_red <- grepl("(_Red\\.idat(\\.gz)?)$", f, ignore.case=TRUE)
  is_grn <- grepl("(_Grn\\.idat(\\.gz)?)$", f, ignore.case=TRUE)
  keep <- is_red | is_grn; f <- f[keep]
  basenames <- sub("(_Red|_Grn)\\.idat(\\.gz)?$", "", f, ignore.case=TRUE)
  tab <- as.data.frame.matrix(table(basenames, ifelse(is_red[keep],"Red","Grn"))); tab[is.na(tab)] <- 0
  paired <- rownames(tab)[(tab$Red > 0) & (tab$Grn > 0)]
  leaf <- basename(paired); sid <- sub("_.*$", "", leaf); pos <- normalize_pos(sub("^[^_]*_", "", leaf))
  data.frame(Sample_Name=paste0(sid,"_",pos), Basename=paired, stringsAsFactors = FALSE)
}
safe_write_csv <- function(df, out_path) {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  tmp <- tempfile(tmpdir = dirname(out_path))
  readr::write_csv(df, tmp)
  ok <- suppressWarnings(file.rename(tmp, out_path))
  if (!ok) { ok2 <- file.copy(tmp, out_path, overwrite = TRUE); unlink(tmp); if (!ok2) stop("Failed to write: ", out_path) }
}

add_dmpa_cols <- function(df, crosswalk = id_crosswalk) {
  # df must have a column named "Sample"
  up <- crosswalk$DMPA_ID[match(df$Sample, crosswalk$Idat)]
  df$DMPA_ID    <- up
  df$DMPA_Label <- ifelse(is.na(up), NA_character_, paste0("DMPA-", up))
  # Put the new columns right after Sample
  keep_order <- c("Sample","DMPA_ID","DMPA_Label",
                  setdiff(names(df), c("Sample","DMPA_ID","DMPA_Label")))
  df[, keep_order, drop = FALSE]
}

#Robust seg_to_df (diagnostic-friendly)
seg_to_df <- function(seg, sample_id, diag_dir = NULL) {
  if (!is.null(diag_dir)) {
    dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)
    txt_path <- file.path(diag_dir, paste0(sample_id, ".CNSegment.txt"))
    if (!file.exists(txt_path)) {
      capture.output({
        cat("class(seg): ", paste(class(seg), collapse=" / "), "\n", sep = "")
        if (is.list(seg)) {
          cat("list names: ", paste(names(seg), collapse=", "), "\n", sep = "")
          for (nm in names(seg)) {
            cat("  $", nm, " : ", paste(class(seg[[nm]]), collapse="/"), "\n", sep = "")
            if (is.data.frame(seg[[nm]])) cat("    cols: ", paste(colnames(seg[[nm]]), collapse=", "), "\n", sep = "")
          }
        }
        st <- try(getFromNamespace("segTable","sesame"), silent = TRUE)
        if (!inherits(st,"try-error")) {
          cat("\nsegTable() available: yes\n")
          tt <- try(st(seg), silent = TRUE)
          cat("segTable() class: ", paste(class(tt), collapse="/"), "\n", sep = "")
          if (is.data.frame(tt)) cat("segTable() cols: ", paste(colnames(tt), collapse=", "), "\n", sep = "")
        } else cat("\nsegTable() available: no\n")
      }, file = txt_path)
    }
  }
  .normalize_df <- function(df) {
    nn <- tolower(gsub("[^a-z0-9]+","", colnames(df)))
    pick <- function(cands) { i <- which(nn %in% cands)[1]; if (length(i)) i else NA_integer_ }
    i_chr <- pick(c("chrom","chr","chromosome","seqnames"))
    i_sta <- pick(c("locstart","start","startpos","startposition"))
    i_end <- pick(c("locend","end","endpos","endposition"))
    i_val <- pick(c("segmean","value","log2","mean","log2ratio","cnlrmedian","signal"))
    if (any(is.na(c(i_chr,i_sta,i_end,i_val)))) return(NULL)
    data.frame(
      Sample     = sample_id,
      Chromosome = as.character(df[[i_chr]]),
      Start      = as.integer(df[[i_sta]]),
      End        = as.integer(df[[i_end]]),
      Log2       = as.numeric(df[[i_val]]),
      stringsAsFactors = FALSE
    )
  }
  # segTable preferred
  st <- try(getFromNamespace("segTable","sesame"), silent = TRUE)
  if (!inherits(st, "try-error")) {
    tt <- try(st(seg), silent = TRUE)
    if (!inherits(tt, "try-error") && !is.null(tt)) {
      df <- as.data.frame(tt, stringsAsFactors = FALSE)
      out <- .normalize_df(df); if (!is.null(out)) return(out)
    }
  }
  # Common cases
  if (is.list(seg)) {
    for (nm in c("seg.signals","segments","table","df"))
      if (!is.null(seg[[nm]]) && is.data.frame(seg[[nm]])) {
        out <- .normalize_df(seg[[nm]]); if (!is.null(out)) return(out)
      }
    for (nm in c("seg.signals.gr","segments.gr","gr","ranges"))
      if (!is.null(seg[[nm]]) && methods::is(seg[[nm]],"GRanges")) {
        gr <- seg[[nm]]; mm <- S4Vectors::mcols(gr); vv <- NULL
        for (cand in c("value","seg.mean","signal","log2","cnlr.median","mean","log2ratio"))
          if (cand %in% colnames(mm)) { vv <- mm[[cand]]; break }
        if (is.null(vv) && ncol(mm)>=1) vv <- mm[[1]]
        return(data.frame(
          Sample=sample_id,
          Chromosome=as.character(GenomeInfoDb::seqnames(gr)),
          Start=as.integer(IRanges::start(gr)), End=as.integer(IRanges::end(gr)),
          Log2=as.numeric(vv), stringsAsFactors = FALSE))
      }
    for (nm in names(seg)) if (is.data.frame(seg[[nm]])) {
      out <- .normalize_df(seg[[nm]]); if (!is.null(out)) return(out)
    }
  }
  gr <- try(methods::as(seg,"GRanges"), silent = TRUE)
  if (!inherits(gr,"try-error")) {
    mm <- S4Vectors::mcols(gr); vv <- NULL
    for (cand in c("value","seg.mean","signal","log2","cnlr.median","mean","log2ratio"))
      if (cand %in% colnames(mm)) { vv <- mm[[cand]]; break }
    if (is.null(vv) && ncol(mm)>=1) vv <- mm[[1]]
    return(data.frame(
      Sample=sample_id,
      Chromosome=as.character(GenomeInfoDb::seqnames(gr)),
      Start=as.integer(IRanges::start(gr)), End=as.integer(IRanges::end(gr)),
      Log2=as.numeric(vv), stringsAsFactors = FALSE))
  }
  # last resort: 0-row df
  data.frame(Sample=character(0), Chromosome=character(0),
             Start=integer(0), End=integer(0), Log2=numeric(0), stringsAsFactors = FALSE)
}

infer_build_from_seg <- function(seg_csv) {
  df <- try(suppressMessages(readr::read_csv(seg_csv, show_col_types = FALSE, n_max = 20000)), silent = TRUE)
  if (inherits(df,"try-error") || !ncol(df)) return("hg38")
  chr <- gsub("^chr","", as.character(df$Chromosome)); chr <- suppressWarnings(as.integer(chr))
  df1 <- df[chr == 1 & is.finite(df$End), , drop = FALSE]
  if (!nrow(df1)) return("hg38")
  m1 <- max(as.numeric(df1$End), na.rm = TRUE)
  if (abs(m1 - 248956422) < abs(m1 - 249250621)) "hg38" else "hg19"
}
make_cytoband_arms <- function(genome = c("hg38","hg19")) {
  genome <- match.arg(genome)
  params <- switch(genome,
                   hg38 = list(
                     len=c(`1`=248956422,`2`=242193529,`3`=198295559,`4`=190214555,`5`=181538259,`6`=170805979,
                           `7`=159345973,`8`=145138636,`9`=138394717,`10`=133797422,`11`=135086622,`12`=133275309,
                           `13`=114364328,`14`=107043718,`15`=101991189,`16`=90338345,`17`=83257441,`18`=80373285,
                           `19`=58617616,`20`=64444167,`21`=46709983,`22`=50818468),
                     cen=c(`1`=123400000,`2`=93800000,`3`=90900000,`4`=50400000,`5`=48400000,`6`=61000000,
                           `7`=59900000,`8`=45600000,`9`=49200000,`10`=40200000,`11`=53700000,`12`=35800000,
                           `13`=17900000,`14`=17600000,`15`=19000000,`16`=36600000,`17`=24000000,`18`=17200000,
                           `19`=26500000,`20`=27500000,`21`=13200000,`22`=14700000)),
                   hg19 = list(
                     len=c(`1`=249250621,`2`=243199373,`3`=198022430,`4`=191154276,`5`=180915260,`6`=171115067,
                           `7`=159138663,`8`=146364022,`9`=141213431,`10`=135534747,`11`=135006516,`12`=133851895,
                           `13`=115169878,`14`=107349540,`15`=102531392,`16`=90354753,`17`=81195210,`18`=78077248,
                           `19`=59128983,`20`=63025520,`21`=48129895,`22`=51304566),
                     cen=c(`1`=121535434,`2`=92326171,`3`=90504854,`4`=49660117,`5`=46405641,`6`=58830166,
                           `7`=58054331,`8`=43838887,`9`=47367679,`10`=39254935,`11`=51644205,`12`=34856694,
                           `13`=16000000,`14`=16000000,`15`=17000000,`16`=35335801,`17`=22263006,`18`=15460898,
                           `19`=24681782,`20`=26369569,`21`=11288129,`22`=13000000))
  )
  ch <- as.integer(names(params$len))
  gr_p <- GRanges(seqnames=paste0("chr", ch), ranges=IRanges(1L, pmax(1L, as.integer(params$cen[ch])-1L)))
  gr_q <- GRanges(seqnames=paste0("chr", ch), ranges=IRanges(as.integer(params$cen[ch]), as.integer(params$len[ch])))
  mcols(gr_p)$arm <- "p"; mcols(gr_q)$arm <- "q"
  out <- c(gr_p, gr_q); GenomeInfoDb::seqlevelsStyle(out) <- "UCSC"; out
}

#Calling logic (mode-centering, asymmetric thresholds)
center_segments <- function(seg_df) {
  df <- seg_df
  chr <- gsub("^chr","", as.character(df$Chromosome))
  keep <- chr %in% as.character(1:22) & is.finite(df$Log2)
  if (!any(keep)) return(df)
  w <- pmax(1, as.numeric(df$End) - as.numeric(df$Start) + 1)
  x <- df$Log2[keep]; w <- w[keep]
  mode0 <- NA_real_
  try({
    w2 <- pmax(1, round(1e5 * (w / sum(w))))
    xx <- rep(x, times = w2)
    d  <- stats::density(xx, n = 2048, bw = "nrd0")
    mode0 <- d$x[which.max(d$y)]
  }, silent = TRUE)
  if (!is.finite(mode0)) {
    d <- stats::density(x, n = 2048, bw = "nrd0")
    mode0 <- d$x[which.max(d$y)]
  }
  df$Log2 <- df$Log2 - mode0
  df
}
sample_thresholds <- function(seg_df,
                              loss_min = LOSS_MIN_ABS_THR, gain_min = GAIN_MIN_ABS_THR,
                              k_loss = LOSS_MAD_K,       k_gain = GAIN_MAD_K) {
  x <- seg_df$Log2[is.finite(seg_df$Log2)]
  mad0 <- stats::mad(x, center = 0, constant = 1, na.rm = TRUE)
  thr_loss <- max(loss_min, k_loss * mad0)
  thr_gain <- max(gain_min, k_gain * mad0)
  c(loss = -thr_loss, gain = +thr_gain)
}
chrom_wm <- function(seg_df, autosomes = as.character(1:22)) {
  out <- setNames(rep(NA_real_, length(autosomes)), autosomes)
  if (!nrow(seg_df)) return(out)
  chr <- gsub("^chr","", as.character(seg_df$Chromosome))
  chr <- ifelse(chr %in% c("X","Y","M","MT"), chr, suppressWarnings(as.character(as.integer(chr))))
  width <- pmax(1, as.numeric(seg_df$End) - as.numeric(seg_df$Start) + 1)
  keep  <- chr %in% c(autosomes, "X","Y")
  if (!any(keep)) return(out)
  num <- tapply(seg_df$Log2[keep] * width[keep], chr[keep], sum, na.rm = TRUE)
  den <- tapply(width[keep],                     chr[keep], sum, na.rm = TRUE)
  wm  <- num/den; hit <- intersect(names(wm), autosomes); out[hit] <- wm[hit]; out
}
arm_wm <- function(seg_df, cyto) {
  if (!nrow(seg_df)) return(NULL)
  q <- GRanges(seqnames = seg_df$Chromosome,
               ranges   = IRanges(start = seg_df$Start, end = seg_df$End),
               value    = seg_df$Log2)
  GenomeInfoDb::seqlevelsStyle(q) <- GenomeInfoDb::seqlevelsStyle(cyto)[1]
  hits <- GenomicRanges::findOverlaps(q, cyto, ignore.strand = TRUE)
  if (!length(hits)) return(NULL)
  qi <- S4Vectors::queryHits(hits); si <- S4Vectors::subjectHits(hits)
  inter <- GenomicRanges::pintersect(q[qi], cyto[si], ignore.strand = TRUE)
  w     <- as.numeric(BiocGenerics::width(inter))
  labs  <- paste0(gsub("^chr","", as.character(GenomeInfoDb::seqnames(cyto)[si])),
                  S4Vectors::mcols(cyto)$arm[si])
  vals  <- S4Vectors::mcols(q)$value[qi]
  num <- tapply(vals * w, labs, sum, na.rm = TRUE)
  den <- tapply(w,        labs, sum, na.rm = TRUE)
  out <- num/den; out[order(names(out))]
}
frac_over_threshold <- function(seg_df, target_gr, mag_thr = SEG_MAG_THR) {
  q <- GRanges(seqnames = seg_df$Chromosome,
               ranges   = IRanges(start = seg_df$Start, end = seg_df$End),
               value    = seg_df$Log2)
  GenomeInfoDb::seqlevelsStyle(q) <- GenomeInfoDb::seqlevelsStyle(target_gr)[1]
  hits <- GenomicRanges::findOverlaps(q, target_gr, ignore.strand = TRUE)
  if (!length(hits)) return(rep(0, length(target_gr)))
  qi <- S4Vectors::queryHits(hits); si <- S4Vectors::subjectHits(hits)
  inter <- GenomicRanges::pintersect(q[qi], target_gr[si], ignore.strand = TRUE)
  w     <- as.numeric(BiocGenerics::width(inter))
  vals  <- S4Vectors::mcols(q)$value[qi]
  by_tgt <- split(seq_along(w), si)
  out <- numeric(length(target_gr)); tgt_len <- as.numeric(BiocGenerics::width(target_gr))
  for (k in names(by_tgt)) {
    idx <- by_tgt[[k]]
    w_keep <- w[idx][abs(vals[idx]) >= mag_thr]
    out[as.integer(k)] <- sum(w_keep) / tgt_len[as.integer(k)]
  }
  out
}
to_calls_disp <- function(M) ifelse(M <= DISPLAY_LOSS_THR, -1L, ifelse(M >= DISPLAY_GAIN_THR, +1L, 0L))

#Segment all samples (NO external ref)
idx <- idat_index_from_dir(idat_dir)
if (!nrow(idx)) stop("No paired IDATs found in: ", idat_dir)
message("SeSAMe CNV: found ", nrow(idx), " paired IDATs.")

seg_paths <- character(nrow(idx))
for (i in seq_len(nrow(idx))) {
  bs <- idx$Basename[i]; sm <- idx$Sample_Name[i]
  message(sprintf("[%d/%d] %s", i, nrow(idx), sm))
  tryCatch({
    sdf <- sesame::readIDATpair(bs)
    sdf <- sesame::noob(sdf)
    if (USE_POOBAH)          sdf <- sesame::pOOBAH(sdf)
    if (USE_DYE_BIAS_TINORM) sdf <- sesame::dyeBiasCorrTypeINorm(sdf)
    seg <- sesame::cnSegmentation(sdf)                # internal reference
    saveRDS(seg, file.path(diag_dir, paste0(sm, ".CNSegment.rds")))
    seg_df <- seg_to_df(seg, sm, diag_dir = diag_dir)
    seg_df <- center_segments(seg_df)                 # robust re-centering
    out_seg <- file.path(out_root, paste0(sm, ".seg.csv"))
    safe_write_csv(seg_df, out_seg)
    seg_paths[i] <- out_seg
  }, error = function(e) {
    message("  ⚠️  Skipping ", sm, " (", conditionMessage(e), "). See diagnostics.")
    seg_paths[i] <<- NA_character_
  })
}
seg_paths <- seg_paths[!is.na(seg_paths)]
if (!length(seg_paths)) stop("No segment files were produced. See ", normalizePath(diag_dir))

#Per-arm definitions
build <- infer_build_from_seg(seg_paths[1])
cyto  <- make_cytoband_arms(build)
message("Using arm definitions for genome build: ", build)

#Summaries & calls
load_seg <- function(pth) suppressMessages(readr::read_csv(pth, show_col_types = FALSE))
seg_list <- lapply(seg_paths, load_seg)
names(seg_list) <- sub("\\.seg\\.csv$", "", basename(seg_paths))

autosomes <- as.character(1:22)
arm_names <- as.vector(outer(autosomes, c("p","q"), paste0))

chr_mat   <- matrix(NA_real_, nrow = length(seg_list), ncol = length(autosomes),
                    dimnames = list(names(seg_list), autosomes))
arm_mat   <- matrix(NA_real_, nrow = length(seg_list), ncol = length(arm_names),
                    dimnames = list(names(seg_list), arm_names))
chr_calls <- matrix(NA_integer_, nrow = length(seg_list), ncol = length(autosomes),
                    dimnames = list(names(seg_list), autosomes))
arm_calls <- matrix(NA_integer_, nrow = length(seg_list), ncol = length(arm_names),
                    dimnames = list(names(seg_list), arm_names))

###-----Compute per-sample MAD and effective thresholds, including DMPA IDs-----###
mad_table <- lapply(names(seg_list), function(nm) {
  s <- seg_list[[nm]]
  x <- s$Log2[is.finite(s$Log2)]
  mad0 <- stats::mad(x, center = 0, constant = 1, na.rm = TRUE)
  thr_loss <- max(LOSS_MIN_ABS_THR, LOSS_MAD_K * mad0)
  thr_gain <- max(GAIN_MIN_ABS_THR, GAIN_MAD_K * mad0)
  
  #Lookup DMPA ID (if present)
  dmpa_id <- id_crosswalk$DMPA_ID[id_crosswalk$Idat == nm]
  if (length(dmpa_id) == 0) dmpa_id <- NA  # fallback if no match
  
  data.frame(
    Sample = nm,
    DMPA_ID = dmpa_id,
    MAD = mad0,
    thr_loss = thr_loss,
    thr_gain = thr_gain
  )
}) |> dplyr::bind_rows()

# Save to Excel-style CSV
safe_write_csv(mad_table, file.path(out_root, "PerSample_MAD_and_Thresholds_with_DMPA.csv"))
message("✅ Saved MAD table with DMPA IDs: ",
        normalizePath(file.path(out_root, "PerSample_MAD_and_Thresholds_with_DMPA.csv")))

#Mask acrocentric p-arms (for calling matrices)
ACRO_P <- c("13p","14p","15p","21p","22p")
if (length(intersect(ACRO_P, colnames(arm_mat)))) {
  arm_mat[,   intersect(ACRO_P, colnames(arm_mat))]   <- NA_real_
  arm_calls[, intersect(ACRO_P, colnames(arm_calls))] <- NA_integer_
}

#Chromosome lengths for coverage (hg38)
chr_len <- c(`1`=248956422,`2`=242193529,`3`=198295559,`4`=190214555,`5`=181538259,`6`=170805979,
             `7`=159345973,`8`=145138636,`9`=138394717,`10`=133797422,`11`=135086622,`12`=133275309,
             `13`=114364328,`14`=107043718,`15`=101991189,`16`=90338345,`17`=83257441,`18`=80373285,
             `19`=58617616,`20`=64444167,`21`=46709983,`22`=50818468)
chr_gr <- GRanges(seqnames = paste0("chr", autosomes),
                  ranges   = IRanges(1L, as.integer(chr_len[autosomes])))

for (nm in names(seg_list)) {
  s <- seg_list[[nm]]
  if (!nrow(s)) next
  
  thr <- sample_thresholds(s,
                           loss_min = LOSS_MIN_ABS_THR, gain_min = GAIN_MIN_ABS_THR,
                           k_loss = LOSS_MAD_K,          k_gain = GAIN_MAD_K)
  
  chr_mean <- chrom_wm(s, autosomes = autosomes)
  a_mean   <- arm_wm(s, cyto)
  
  chr_mat[nm, ] <- chr_mean
  if (!is.null(a_mean)) {
    keep <- intersect(names(a_mean), colnames(arm_mat))
    arm_mat[nm, keep] <- a_mean[keep]
  }
  
  chr_frac <- frac_over_threshold(s, chr_gr, mag_thr = SEG_MAG_THR); names(chr_frac) <- autosomes
  arm_frac_vec <- frac_over_threshold(s, cyto,    mag_thr = SEG_MAG_THR)
  names(arm_frac_vec) <- paste0(gsub("^chr","", as.character(seqnames(cyto))), mcols(cyto)$arm)
  arm_frac <- tapply(arm_frac_vec, names(arm_frac_vec), mean, na.rm = TRUE)
  
  chr_calls[nm, ] <- 0L
  chr_calls[nm, chr_mean <= thr["loss"] & chr_frac >= COVER_FRAC_CHR] <- -1L
  chr_calls[nm, chr_mean >= thr["gain"] & chr_frac >= COVER_FRAC_CHR] <- +1L
  
  if (!is.null(a_mean)) {
    a_col <- intersect(names(a_mean), colnames(arm_calls))
    arm_calls[nm, a_col] <- 0L
    arm_calls[nm, a_col[a_mean[a_col] <= thr["loss"] & arm_frac[a_col] >= COVER_FRAC_ARM]] <- -1L
    arm_calls[nm, a_col[a_mean[a_col] >= thr["gain"] & arm_frac[a_col] >= COVER_FRAC_ARM]] <- +1L
  }
}

###-----Three-platform CNV cross-validation concordance-----###
# Arm-level CNV calls from SeSAMe were cross-validated against two independent
# platforms: (1) the Heidelberg Epignostix classifier (Conumee 2.0) and
# (2) targeted next-generation sequencing (GlioSeq for 9 patients, Oncomine
# for 1 patient). The concordance table below summarizes results for all 10
# patients and identifies DMPA-10 as the sole discordant case.

cnv_concordance <- data.frame(
  DMPA_ID         = 1:10,
  SeSAMe_CNV      = c("Stable",                            # DMPA-1
                      "Stable",                            # DMPA-2
                      "1p/q,2p/q,7p/q,13q,18p/q,22q",     # DMPA-3
                      "Stable",                            # DMPA-4
                      "1p,2p",                             # DMPA-5
                      "Stable",                            # DMPA-6
                      "Stable",                            # DMPA-7
                      "Stable",                            # DMPA-8
                      "1p,2p",                             # DMPA-9
                      "20p/q,22q"),                        # DMPA-10
  Heidelberg_CNV  = c("Stable",                            # DMPA-1
                      "Stable",                            # DMPA-2
                      "Confirmed",                         # DMPA-3
                      "Stable",                            # DMPA-4
                      "Confirmed",                         # DMPA-5
                      "Stable",                            # DMPA-6
                      "Stable",                            # DMPA-7
                      "Stable",                            # DMPA-8
                      "Confirmed",                         # DMPA-9
                      "Flat"),                             # DMPA-10
  NGS_Platform    = c("GlioSeq (clinical)",                # DMPA-1
                      "GlioSeq",                           # DMPA-2
                      "GlioSeq",                           # DMPA-3
                      "GlioSeq",                           # DMPA-4
                      "Oncomine (clinical)",               # DMPA-5
                      "GlioSeq",                           # DMPA-6
                      "GlioSeq",                           # DMPA-7
                      "GlioSeq",                           # DMPA-8
                      "GlioSeq",                           # DMPA-9
                      "GlioSeq"),                          # DMPA-10
  NGS_CNV         = c("Negative",                          # DMPA-1
                      "No CNV (TRAF7 only)",               # DMPA-2
                      "1p/2p/2q/7p/7q/13q/22q",           # DMPA-3
                      "No CNV (TRAF7 only)",               # DMPA-4
                      "1p,2p",                             # DMPA-5
                      "Negative (limited quality)",        # DMPA-6
                      "No CNV (TRAF7 only)",               # DMPA-7
                      "No CNV (PIK3CA + FGFR1 only)",     # DMPA-8
                      "1p/2p",                             # DMPA-9
                      "No CNV (TRAF7 only)"),              # DMPA-10
  Concordance     = c("Full (3/3 stable)",                 # DMPA-1
                      "Full (3/3 stable)",                 # DMPA-2
                      "Full (3/3 confirmed)",              # DMPA-3
                      "Full (3/3 stable)",                 # DMPA-4
                      "Full (3/3 confirmed)",              # DMPA-5
                      "Full (3/3 stable)",                 # DMPA-6
                      "Full (3/3 stable)",                 # DMPA-7
                      "Full (3/3 stable)",                 # DMPA-8
                      "Full (3/3 confirmed)",              # DMPA-9
                      "Discordant (1/3 — SeSAMe only)"),  # DMPA-10
  stringsAsFactors = FALSE
)

safe_write_csv(cnv_concordance,
               file.path(out_root, "CNV_three_platform_concordance.csv"))
message("✅ Saved three-platform CNV concordance table: ",
        normalizePath(file.path(out_root, "CNV_three_platform_concordance.csv")))

# Print summary
message("\n=== Three-platform CNV cross-validation summary ===")
message("Patients with CNV events confirmed by ≥2 platforms: DMPA-3, DMPA-5, DMPA-9")
message("Patients CNV-stable across all 3 platforms:         DMPA-1, 2, 4, 6, 7, 8")
message("Discordant (SeSAMe-only, not confirmed):            DMPA-10")
message("Action: DMPA-10 arm-level calls overridden to CNV-stable\n")


###-----Cross-validation override: DMPA-10-----###
# Three-platform cross-validation (SeSAMe, Heidelberg Epignostix/Conumee 2.0,
# and GlioSeq targeted sequencing) identified DMPA-10 arm-level calls (20p, 20q,
# 22q) as false positives: Heidelberg profile is flat, GlioSeq reports no CNV
# events (TRAF7 K389N only), and SeSAMe weighted means are marginal (-0.30 to
# -0.39). All other patients with SeSAMe calls (DMPA-3, DMPA-5, DMPA-9) are
# confirmed by at least two independent platforms. Override DMPA-10 to CNV-stable.
dmpa10_idat <- id_crosswalk$Idat[id_crosswalk$DMPA_ID == 10]
if (dmpa10_idat %in% rownames(arm_calls)) {
  arm_calls[dmpa10_idat, ]  <- ifelse(
    is.na(arm_calls[dmpa10_idat, ]), NA_integer_, 0L
  )
  chr_calls[dmpa10_idat, ]  <- 0L
  message("DMPA-10 (", dmpa10_idat, "): arm + chr calls overridden to ",
          "CNV-stable (three-platform cross-validation)")
}

#Save tables
safe_write_csv(
  add_dmpa_cols(data.frame(Sample = rownames(chr_mat),   chr_mat,   check.names = FALSE)),
  file.path(out_root, "CNV_chromosome_weighted_means.csv")
)

safe_write_csv(
  add_dmpa_cols(data.frame(Sample = rownames(arm_mat),   arm_mat,   check.names = FALSE)),
  file.path(out_root, "CNV_arms_weighted_means.csv")
)

safe_write_csv(
  add_dmpa_cols(data.frame(Sample = rownames(chr_calls), chr_calls, check.names = FALSE)),
  file.path(out_root, "CNV_chromosome_calls_discrete.csv")
)

safe_write_csv(
  add_dmpa_cols(data.frame(Sample = rownames(arm_calls), arm_calls, check.names = FALSE)),
  file.path(out_root, "CNV_arms_calls_discrete.csv")
)


###-----CNV Publication figures - preparation-----###
# Two versions:
# 1) Display thresholds on continuous arm means 
# 2) Formal calls (uses arm_calls)


# --- shared aesthetics ---
col_disp  <- c("-1"="#08306B","0"="#f2f2f2","+1"="#A50F15")   # Loss / Neutral / Gain
na_tile   <- "#F2f2f2"                                        # for NA tiles (formal calls)
COMPACT_DROP_P <- c(13, 14, 15, 21, 22)                              # compact layout from your v2
autosomes <- as.character(1:22)

# helper: build plot from a -1/0/+1 (or NA) matrix, rows = samples, cols = arms "1p","1q",...
.build_pub_plot <- function(Mplot, title_text, outfile_png,
                            neutralize_acro_p_for_display = FALSE,
                            na_value = na_tile) {
  
  #Optionally neutralize acrocentric p arms for display
  ACRO_P <- c("13p","14p","15p","21p","22p")
  if (neutralize_acro_p_for_display) {
    hit <- intersect(ACRO_P, colnames(Mplot))
    if (length(hit)) Mplot[, hit] <- 0L
  }
  
  #Order rows by DMPA_ID using your crosswalk
  common_ids <- intersect(rownames(Mplot), id_crosswalk$Idat)
  stopifnot(length(common_ids) > 0)
  Mplot <- Mplot[match(id_crosswalk$Idat, rownames(Mplot)), , drop = FALSE]
  ord   <- order(id_crosswalk$DMPA_ID)
  Mplot <- Mplot[ord, , drop = FALSE]
  rownames(Mplot) <- paste0("DMPA-", id_crosswalk$DMPA_ID[ord])
  
  #Ensure arm order is 1p,1q,2p,2q,... (and drop any missing)
  autosomes <- as.character(1:22)
  arm_order <- as.vector(rbind(paste0(autosomes,"p"), paste0(autosomes,"q")))
  arm_order <- intersect(arm_order, colnames(Mplot))
  Mplot <- Mplot[, arm_order, drop = FALSE]
  
  #Long-format data
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(grid)
  })
  
  df_long <- as.data.frame(Mplot) |>
    rownames_to_column("Sample") |>
    pivot_longer(-Sample, names_to = "ArmFull", values_to = "call") |>
    separate(ArmFull, into = c("Chr","Arm"), sep = "(?<=\\d)(?=[pq])") |>
    mutate(
      Chr    = as.integer(Chr),
      call   = as.integer(call),
      Sample = factor(Sample, levels = rownames(Mplot)),
      LOSS   = (call == -1L)
    ) |>
    filter(!(Chr %in% COMPACT_DROP_P & Arm == "p"))
  
  #Map arms to x positions for compact layout
  arm_map <- lapply(1:22, function(ch) {
    aa <- if (ch %in% COMPACT_DROP_P) "q" else c("p","q")
    data.frame(Chr = ch, Arm = aa, stringsAsFactors = FALSE)
  }) |> bind_rows() |> mutate(x = dplyr::row_number())
  
  df_long <- dplyr::left_join(df_long, arm_map, by = c("Chr","Arm"))
  
  #Axis ticks, boxes, grid lines
  xticks <- arm_map |> dplyr::group_by(Chr) |> dplyr::summarise(x = mean(x), .groups = "drop")
  nS  <- nrow(Mplot); gap <- 0.10
  
  box <- arm_map |>
    dplyr::group_by(Chr) |>
    dplyr::summarise(xmin = min(x) - 0.5 + gap,
                     xmax = max(x) + 0.5 - gap,
                     .groups = "drop") |>
    dplyr::mutate(ymin = 0.5, ymax = nS + 0.5)
  
  edges <- arm_map |>
    dplyr::group_by(Chr) |>
    dplyr::summarise(xmin_chr = min(x), xmax_chr = max(x), .groups = "drop")
  
  df_tile <- dplyr::left_join(df_long, edges, by = "Chr") |>
    dplyr::mutate(
      pad_left  = ifelse(x == xmin_chr, gap, 0),
      pad_right = ifelse(x == xmax_chr, gap, 0),
      row_nat   = match(as.character(Sample), rownames(Mplot)),
      row_idx   = nS - row_nat + 1,
      xmin = x - 0.5 + pad_left,
      xmax = x + 0.5 - pad_right,
      ymin = row_idx - 0.5,
      ymax = row_idx + 0.5
    )
  
  grid_h <- box |>
    tidyr::expand_grid(y = seq(1, nS, by = 1)) |>
    dplyr::mutate(y = y + 0.5)
  
  #Vertical split between p and q within a chromosome box
  vmap <- arm_map |>
    dplyr::group_by(Chr) |>
    dplyr::arrange(x, .by_group = TRUE) |>
    dplyr::mutate(xnext = dplyr::lead(x)) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(xnext)) |>
    dplyr::left_join(box |> dplyr::select(Chr, ymin, ymax), by = "Chr") |>
    dplyr::transmute(
      xmin = (x + xnext)/2,
      xmax = (x + xnext)/2,
      ymin = ymin,
      ymax = ymax
    )
  
  guide_col <- "grey85"
  guide_lwd <- 0.35
  
  #Plot
  p <- ggplot() +
    # neutral background for each tile
    geom_rect(
      data = df_tile,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey95",
      color = NA,
      linewidth = 0
    ) +
    # losses overlay
    geom_rect(
      data = df_tile,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = LOSS),
      color = NA, linewidth = 0
    ) +
    # horizontal grid lines
    geom_segment(
      data = grid_h, inherit.aes = FALSE,
      aes(x = xmin, xend = xmax, y = y, yend = y),
      color = guide_col, linewidth = guide_lwd
    ) +
    # vertical split lines within a chromosome box
    geom_segment(
      data = vmap, inherit.aes = FALSE,
      aes(x = xmin, xend = xmax, y = ymin, yend = ymax),
      color = guide_col, linewidth = guide_lwd, lineend = "butt"
    ) +
    # chromosome box borders
    geom_rect(
      data = box, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      color = "#000000", fill = NA, linewidth = 0.3
    ) +
    scale_fill_manual(
      values = c(`TRUE` = "#08306B", `FALSE` = "#F2F2F2"),
      breaks = "TRUE",
      labels = "Copy number loss",
      name   = NULL,
      drop = FALSE, na.translate = FALSE
    ) +
    guides(fill = guide_legend(
      direction = "horizontal", nrow = 1, byrow = TRUE,
      keywidth = unit(16, "pt"), keyheight = unit(10, "pt")
    )) +
    scale_x_continuous(
      breaks = xticks$x,
      labels = xticks$Chr,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(
      breaks = 1:nS,
      labels = rev(rownames(Mplot)),
      expand = c(0, 0)
    ) +
    coord_fixed(ratio = 0.9, clip = "on", ylim = c(0.5, nS + 0.5)) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(
        size = 13, face = "plain",
        margin = ggplot2::margin(t = 6, unit = "pt")
      ),
      axis.text.y = element_text(size = 13, face = "plain", hjust = 1),
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text  = element_text(size = 14, face = "bold"),
      legend.box   = "horizontal"
    )
  
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(fig_dir, outfile_png), p,
         width = 16, height = 5, units = "in", dpi = 300, device = "png")
  
  message("✅ Saved: ", normalizePath(file.path(fig_dir, outfile_png)))
  invisible(Mplot)
}
  
  
###-----Plot figure with display thresholds (based on user inputs/knobs)-----###
M_disp <- ifelse(arm_mat <= DISPLAY_LOSS_THR, -1L,
                 ifelse(arm_mat >= DISPLAY_GAIN_THR, +1L, 0L))
storage.mode(M_disp) <- "integer"

.build_pub_plot(
  Mplot = M_disp,
  title_text = "Chromosome arm-level CNV (p/q paired, DMPA ordered) — display thresholds",
  outfile_png = "CNV_arms_publication_DISPLAY.png",
  neutralize_acro_p_for_display = TRUE,  # keep the old look
  na_value = na_tile
)

###-----Plot figure with formal calls (based on user inputs/knobs)-----###
M_calls <- arm_calls
storage.mode(M_calls) <- "integer"

#Figure from calls
.build_pub_plot(
  Mplot = M_calls,
  title_text = "Chromosome arm-level CNV (p/q paired, DMPA ordered) — FORMAL CALLS",
  outfile_png = "CNV_arms_publication_CALLS.png",
  neutralize_acro_p_for_display = FALSE, # show NA where appropriate
  na_value = na_tile
)

#Recreate the arranged matrix used inside .build_pub_plot()
arrange_calls_for_export <- function(M) {
  common_ids <- intersect(rownames(M), id_crosswalk$Idat)
  M <- M[match(id_crosswalk$Idat, rownames(M)), , drop = FALSE]
  ord <- order(id_crosswalk$DMPA_ID)
  M <- M[ord, , drop = FALSE]
  rownames(M) <- paste0("DMPA-", id_crosswalk$DMPA_ID[ord])
  arm_order <- as.vector(rbind(paste0(autosomes,"p"), paste0(autosomes,"q")))
  arm_order <- intersect(arm_order, colnames(M))
  M <- M[, arm_order, drop = FALSE]
  # drop p for 13/21/22 to match the compact figure
  keep_cols <- !(sub("p$","", colnames(M)) %in% as.character(COMPACT_DROP_P) & grepl("p$", colnames(M)))
  M[, keep_cols, drop = FALSE]
}

M_calls_figmatrix <- arrange_calls_for_export(M_calls)
safe_write_csv(
  data.frame(Sample = rownames(M_calls_figmatrix), M_calls_figmatrix, check.names = FALSE),
  file.path(fig_dir, "CNV_arms_calls_matrix_DMPA_ordered_forFigure.csv")
)
message("✅ Saved calls matrix (DMPA-ordered, p/q paired, compact): ",
        normalizePath(file.path(fig_dir, "CNV_arms_calls_matrix_DMPA_ordered_forFigure.csv")))



###-----SeSAMe-native genome-wide CNV plots (green/red dots + segments)-----###
suppressPackageStartupMessages({
  library(sesame); library(sesameData)
  library(GenomicRanges); library(GenomeInfoDb); library(IRanges)
  library(ggplot2); library(dplyr); library(tibble); library(grid)
})

#Output paths
out_dir    <- file.path(base_dir, "Combined Analysis_WHO I_Final/SeSAMe_CNV_Final")
native_dir <- file.path(out_dir, "_sesame_native_plots")
dir.create(native_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(exists("idx"), all(c("Sample_Name","Basename") %in% names(idx)))
stopifnot(exists("id_crosswalk"), all(c("Idat","DMPA_ID") %in% names(id_crosswalk)))

#Figure geometry
w_in <- 2000/150; h_in <- 700/150; dpi <- 150

#Helpers
drop_sex_chr <- function(seg){
  sex <- c("chrX","chrY")
  bc  <- seg$bin.coords
  keep_bc <- !(as.character(GenomeInfoDb::seqnames(bc)) %in% sex)
  bc <- bc[keep_bc]
  bc <- GenomeInfoDb::keepSeqlevels(bc, setdiff(GenomeInfoDb::seqlevels(bc), sex),
                                    pruning.mode = "coarse")
  seg$bin.coords  <- bc
  seg$bin.signals <- seg$bin.signals[names(bc)]
  if (!is.null(seg$seg.signals) && nrow(seg$seg.signals)) {
    seg$seg.signals <- seg$seg.signals[!(seg$seg.signals$chrom %in% sex), , drop = FALSE]
  }
  seg
}

make_dmpa_title <- local({
  cw <- setNames(paste0("DMPA-", id_crosswalk$DMPA_ID), id_crosswalk$Idat)
  function(sample_name, basename_path){
    leaf <- basename(as.character(basename_path))
    if (!is.na(leaf) && leaf %in% names(cw)) return(cw[[leaf]])
    if (!is.na(sample_name) && sample_name %in% names(cw)) return(cw[[sample_name]])
    sample_name
  }
})

chr_vlines_from_seg <- function(seg){
  si <- GenomeInfoDb::seqinfo(seg$bin.coords)
  nm <- GenomeInfoDb::seqlevels(si); len <- GenomeInfoDb::seqlengths(si)
  keep <- nm[nm %in% paste0("chr",1:22) & !is.na(len)]
  if (!length(keep)) return(list(
    x_left=0, x_right=1,
    chr_internal=data.frame(x=numeric(0)),
    arm_split=data.frame(x=numeric(0)),
    chr_mid=data.frame(x=numeric(0), lab=character(0))
  ))
  si <- si[keep]
  
  seqname <- GenomeInfoDb::seqlevels(si)
  seqlen  <- as.numeric(GenomeInfoDb::seqlengths(si))
  ord <- order(as.integer(sub("^chr","", seqname)))
  seqname <- seqname[ord]; seqlen <- seqlen[ord]
  
  totlen   <- sum(seqlen)
  seqcum   <- cumsum(seqlen)
  seqstart <- c(0, seqcum[-length(seqcum)]); names(seqstart) <- seqname
  
  x_left  <- (seqstart[1])    / totlen
  x_right <- (tail(seqcum,1)) / totlen
  chr_int <- data.frame(x = (seqstart[-1]) / totlen)
  chr_mid <- data.frame(x = ((seqstart + seqcum)/2)/totlen, lab = seqname)
  
  #Arm_split kept for completeness (optional draw)
  cen <- c(`chr1`=123400000,`chr2`=93800000,`chr3`=90900000,`chr4`=50400000,`chr5`=48400000,
           `chr6`=61000000,`chr7`=59900000,`chr8`=45600000,`chr9`=49200000,`chr10`=40200000,
           `chr11`=53700000,`chr12`=35800000,`chr13`=17900000,`chr14`=17600000,`chr15`=19000000,
           `chr16`=36600000,`chr17`=24000000,`chr18`=17200000,`chr19`=26500000,`chr20`=27500000,
           `chr21`=13200000,`chr22`=14700000)
  hit <- intersect(names(seqstart), names(cen))
  arm_split <- data.frame(x = (seqstart[hit] + cen[hit]) / totlen)
  
  list(x_left=x_left, x_right=x_right,
       chr_internal=chr_int, arm_split=arm_split, chr_mid=chr_mid)
}

###-----Plotting loop - SeSAMe native plots-----###
for (i in seq_len(nrow(idx))) {
  sm <- idx$Sample_Name[i]
  bs <- idx$Basename[i]
  title_id <- make_dmpa_title(sm, bs)
  
  if (!file.exists(paste0(bs,"_Grn.idat")) || !file.exists(paste0(bs,"_Red.idat"))) {
    message("Skipping ", sm, " (missing IDAT pair)"); next
  }
  message("▶ Native CNV plot for ", sm, " → title: ", title_id)
  
  #Read + preprocess
  sdf <- sesame::readIDATpair(bs)
  sdf <- sesame::noob(sdf)
  sdf <- sesame::dyeBiasCorrTypeINorm(sdf)
  if ("pCutoff" %in% names(formals(sesame::pOOBAH))) {
    sdf <- sesame::pOOBAH(sdf, pCutoff = 0.005)  # stricter than default
  } else {
    sdf <- sesame::pOOBAH(sdf)
  }
  
  #Drop low-intensity probes (~bottom 2%) BEFORE CNV calling
  if ("totalIntensities" %in% getNamespaceExports("sesame")) {
    tot <- sesame::totalIntensities(sdf)
  } else if ("meanIntensity" %in% getNamespaceExports("sesame")) {
    tot <- sesame::meanIntensity(sdf)
  } else {
    mu  <- sesame::signalMU(sdf); tot <- mu$M + mu$U; names(tot) <- rownames(mu)
  }
  q_lo <- stats::quantile(tot, 0.02, na.rm = TRUE)
  keep_probes <- which(tot >= q_lo)
  sdf <- sdf[keep_probes, ]
  
  #Segmentation
  seg <- tryCatch(sesame::cnSegmentation(sdf), error = function(e) {
    message("  Segmentation failed: ", e$message); NULL
  })
  if (is.null(seg)) next
  
  seg_auto <- drop_sex_chr(seg)
  VL <- chr_vlines_from_seg(seg_auto)
  
  #Base SeSAMe plot 
  p <- sesame::visualizeSegments(seg_auto) +
    ggplot2::ggtitle(title_id) +
    ggplot2::scale_x_continuous(
      limits = c(VL$x_left, VL$x_right),
      breaks = VL$chr_mid$x,
      labels = sub("^chr","", VL$chr_mid$lab),   # 1..22, horizontal
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-1.35, 1.35),
      breaks  = c(-1.2,-0.8,-0.4,0,0.4,0.8,1.2),
      expand  = c(0, 0)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position  = "none",
      axis.title       = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "black", linewidth = 0.75, fill = NA),
      plot.title       = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x      = element_text(angle = 0, vjust = 1, hjust = 0.5, size = 12, face = "bold"),
      axis.text.y      = element_text(size = 12, face = "bold"),
      axis.ticks.x     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.ticks.y.left     = element_line(colour = "grey55", linewidth = 0.28),
      axis.ticks.length.y.left = grid::unit(2.2, "mm")
    )
  
  #Dot tweak 
  pt_ix <- which(vapply(p$layers, function(L) inherits(L$geom, "GeomPoint"), logical(1)))[1]
  if (!is.na(pt_ix)) {
    p$layers[[pt_ix]]$aes_params$size   <- 0.7
    p$layers[[pt_ix]]$aes_params$alpha  <- 0.7
    p$layers[[pt_ix]]$aes_params$shape  <- 16     # solid circle
    p$layers[[pt_ix]]$aes_params$stroke <- 0      # no outline
  }
  
  #Draw verticals LAST so they sit on top of points
  p$layers <- Filter(function(L) !inherits(L$geom, "GeomVline"), p$layers)
  
  #Solid chromosome boundaries
  if (nrow(VL$chr_internal)) {
    p <- p + ggplot2::geom_vline(
      data = VL$chr_internal, ggplot2::aes(xintercept = x),
      linewidth = 0.30, colour = "grey75", lineend = "square", inherit.aes = FALSE
    )
  }
  
  #dashed p/q arm splits (optional; kept aligned off borders)
  if (nrow(VL$arm_split)) {
    boundary_x <- c(VL$chr_internal$x, VL$x_left, VL$x_right)
    keep <- vapply(VL$arm_split$x, function(xx) all(abs(xx - boundary_x) >= 0.002), logical(1))
    VL$arm_split <- VL$arm_split[keep, , drop = FALSE]
    if (nrow(VL$arm_split)) {
      p <- p + ggplot2::geom_vline(
        data = VL$arm_split, ggplot2::aes(xintercept = x),
        linewidth = 0.26, colour = "grey85", linetype = "longdash",
        lineend = "square", inherit.aes = FALSE
      )
    }
  }
  
  #Save
  outfile <- file.path(native_dir, paste0(sm, "_CNV_sesameNative.png"))
  ggsave(outfile, plot = p, width = w_in, height = h_in, dpi = dpi, units = "in", bg = "white")
  message("✅ Saved: ", normalizePath(outfile))
}

print(mad_table[order(mad_table$DMPA_ID), c("DMPA_ID", "MAD", "thr_loss")])


##==============================================================================
###-----17. PGR GENE DIFFERENTIAL METHYLATION ANALYSIS-----###
##==============================================================================

#Load beta combat and pd csv files
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(openxlsx)     # install once outside the script
})

beta_path <- file.path(base_dir, "Combined Analysis_WHO I_Final/beta_combat_Combined Analysis_All Grades.csv")
pd_path   <- file.path(base_dir, "Combined Analysis_WHO I_Final/pd_Combined Analysis_All Grades.csv")

ResultsDir <- file.path(base_dir, "Combined Analysis_WHO I_Final")
out_root <- file.path(ResultsDir, "PGR and 11q22 DMR")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

load_beta_table <- function(path) {
  stopifnot(file.exists(path))
  info <- file.info(path); stopifnot(info$size > 0)
  sig <- readBin(path, "raw", n = 4)
  is_gzip  <- function(b) length(b)>=2 && b[1]==as.raw(0x1f) && b[2]==as.raw(0x8b)
  is_xlsx  <- function(b) length(b)>=4 && all(b[1:2] == charToRaw("PK"))
  looks_tsv <- function(p){ con <- file(p,"rt"); on.exit(close(con), add=TRUE)
  ln <- ""; while(nzchar(ln)==0 && length(ln)==1) ln <- readLines(con, n=1); grepl("\t", ln) }
  if (is_xlsx(sig)) {
    suppressPackageStartupMessages(require(readxl))
    df <- readxl::read_excel(path) |> as.data.frame(check.names = FALSE)
  } else if (is_gzip(sig)) {
    df <- read.csv(gzfile(path), check.names = FALSE)
  } else if (looks_tsv(path)) {
    df <- readr::read_tsv(path, show_col_types = FALSE) |> as.data.frame(check.names = FALSE)
  } else {
    df <- tryCatch(readr::read_csv(path, show_col_types = FALSE) |> as.data.frame(check.names = FALSE),
                   error = function(e) read.csv(path, check.names = FALSE))
  }
  
  if (nrow(df)>0) {
    cpg_col <- intersect(c("probe","Probe","CpG","cg","rowname","ID","IlmnID"), names(df))
    if (!length(cpg_col)) {
      first <- names(df)[1]
      if (grepl("^cg\\d+", as.character(df[[first]])[1])) cpg_col <- first
    }
    if (length(cpg_col)) { rownames(df) <- as.character(df[[cpg_col[1]]]); df[[cpg_col[1]]] <- NULL }
  }
  stopifnot(is.data.frame(df), nrow(df)>0, ncol(df)>0)
  df
}

beta_combat <- load_beta_table(beta_path) |> as.matrix()
storage.mode(beta_combat) <- "double"

pd_df <- read.csv(pd_path, row.names = 1, check.names = FALSE)
stopifnot(identical(colnames(beta_combat), rownames(pd_df)))

grp <- factor(ifelse(pd_df$Cohort == "DMPA", "DMPA", "Reference"),
              levels = c("DMPA","Reference"))
names(grp) <- rownames(pd_df)

#Add annotations to beta/pd sheet. Include all CpG probes for PGR gene from 450k (hg19) 
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_all <- ann450k[ intersect(rownames(ann450k), rownames(beta_combat)), , drop = FALSE ]

`%||%` <- function(a,b) if (is.null(a)) b else a
has_PGR <- function(x) grepl("(^|;)PGR(;|$)", x %||% "", ignore.case = FALSE)

pgr_cpgs <- rownames(ann_all)[ has_PGR(ann_all$UCSC_RefGene_Name) ]
stopifnot(length(pgr_cpgs) > 0)

keep_cols <- intersect(
  c("Name","chr","pos","strand","UCSC_RefGene_Name","UCSC_RefGene_Group",
    "Relation_to_Island","UCSC_CpG_Islands_Name","Regulatory_Feature_Group",
    "Methyl27_Loci","Probe_rs"),
  colnames(ann_all)
)
pgr_ann <- ann_all[pgr_cpgs, keep_cols, drop = FALSE] |>
  as.data.frame() |>
  rownames_to_column("CpG")

#Per-CpG statistics (Welch + Wilcoxon)
test_one <- function(cg) {
  x <- beta_combat[cg, grp=="DMPA"]
  y <- beta_combat[cg, grp=="Reference"]
  tt <- t.test(x, y, var.equal = FALSE)
  ww <- wilcox.test(x, y, exact = FALSE)
  tibble(
    CpG = cg,
    mean_DMPA  = mean(x, na.rm=TRUE),
    mean_Ref   = mean(y, na.rm=TRUE),
    delta_beta = mean(x, na.rm=TRUE) - mean(y, na.rm=TRUE),
    t_stat     = unname(tt$statistic), pval_t = tt$p.value,
    W_stat     = unname(ww$statistic), pval_wil = ww$p.value
  )
}
res_pgr <- bind_rows(lapply(pgr_cpgs, test_one)) |>
  mutate(FDR_t = p.adjust(pval_t, method="BH"),
         FDR_wil = p.adjust(pval_wil, method="BH")) |>
  left_join(select(pgr_ann, CpG, any_of(c("chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group"))), by="CpG") |>
  arrange(FDR_t, desc(abs(delta_beta)))

#Gene-level summary (mean β across PGR CpGs)
pgr_mean_beta <- colMeans(beta_combat[pgr_cpgs, , drop=FALSE], na.rm=TRUE)
gene_t <- t.test(pgr_mean_beta[grp=="DMPA"], pgr_mean_beta[grp=="Reference"], var.equal=FALSE)
gene_w <- wilcox.test(pgr_mean_beta[grp=="DMPA"], pgr_mean_beta[grp=="Reference"], exact=FALSE)

gene_summary <- tibble(
  group     = c("Reference","DMPA"),
  mean_beta = c(mean(pgr_mean_beta[grp=="Reference"]), mean(pgr_mean_beta[grp=="DMPA"])),
  sd_beta   = c(sd(pgr_mean_beta[grp=="Reference"]),   sd(pgr_mean_beta[grp=="DMPA"])),
  n         = c(sum(grp=="Reference"),                 sum(grp=="DMPA"))
)
gene_tests <- tibble(
  metric     = c("Welch t-test","Wilcoxon rank-sum"),
  statistic  = c(unname(gene_t$statistic), unname(gene_w$statistic)),
  p_value    = c(gene_t$p.value, gene_w$p.value),
  delta_beta = mean(pgr_mean_beta[grp=="DMPA"]) - mean(pgr_mean_beta[grp=="Reference"])
)

#Save data into Excel workbook
wb <- createWorkbook()
addWorksheet(wb, "Per_CpG_Beta_Stats"); writeData(wb, "Per_CpG_Beta_Stats", res_pgr); freezePane(wb, "Per_CpG_Beta_Stats", 2)
addWorksheet(wb, "PGR_Gene_Summary");    writeData(wb, "PGR_Gene_Summary", gene_summary, startRow=1, startCol=1)
writeData(wb, "PGR_Gene_Summary", gene_tests,  startRow=nrow(gene_summary)+3, startCol=1)
for (sh in names(wb$worksheets)) setColWidths(wb, sh, cols=1:50, widths="auto")
saveWorkbook(wb, file.path(out_root,"PGR_beta_perCpG_and_stats.xlsx"), overwrite=TRUE)


###-----Butterfly plot of PGR CpG β values - DMPA vs. reference cohort-----###
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2)
})

stopifnot(exists("beta_combat"), exists("pd_df"))
beta_combat <- as.matrix(beta_combat); storage.mode(beta_combat) <- "double"

grp <- factor(ifelse(pd_df$Cohort == "DMPA", "DMPA", "Reference"),
              levels = c("DMPA","Reference"))

names(grp) <- rownames(pd_df)
stopifnot(identical(colnames(beta_combat), names(grp)))

#Helpers - safely get annotations
.safe_ann <- function(pkg_name, getter = "getAnnotation") {
  if (!requireNamespace(pkg_name, quietly = TRUE)) return(NULL)
  suppressPackageStartupMessages(require(pkg_name, character.only = TRUE))
  ga <- get(getter, asNamespace(pkg_name))
  out <- tryCatch(ga(get(pkg_name)), error = function(e) NULL)
  if (is.null(out)) {
    obj <- tryCatch(get(pkg_name), error = function(e) NULL)
    if (is.null(obj)) return(NULL)
    out <- tryCatch(ga(obj), error = function(e) NULL)
  }
  out
}

#Normalize various annotation objects to a common data.frame schema
norm_ann <- function(ann, source_tag) {
  if (is.null(ann)) return(NULL)
  df <- as.data.frame(ann)
  name_col <- c("Name","IlmnID","TargetID","probe_id")
  chr_col  <- c("chr","CHR","Chromosome","seqnames")
  pos_col  <- c("pos","MAPINFO","Start","start","Position")
  gene_col <- c("UCSC_RefGene_Name","gene_name","Gene_Name","GeneSymbol","Symbol","Genes")
  group_col<- c("UCSC_RefGene_Group","gene_group","Gene_Group","Group")
  pick <- function(x, choices) { hit <- intersect(choices, colnames(x)); if (length(hit)) hit[1] else NA_character_ }
  cn_name <- pick(df, name_col); if (is.na(cn_name)) return(NULL)
  cn_chr  <- pick(df, chr_col)
  cn_pos  <- pick(df, pos_col)
  cn_gene <- pick(df, gene_col)
  cn_grp  <- pick(df, group_col)
  out <- tibble(
    CpG   = df[[cn_name]],
    chr   = if (!is.na(cn_chr))  as.character(df[[cn_chr]])  else NA_character_,
    pos   = if (!is.na(cn_pos))  suppressWarnings(as.numeric(df[[cn_pos]])) else NA_real_,
    gene  = if (!is.na(cn_gene)) as.character(df[[cn_gene]]) else NA_character_,
    group = if (!is.na(cn_grp))  as.character(df[[cn_grp]])  else NA_character_,
    source = source_tag
  )
  out <- out[!is.na(out$CpG) & nzchar(out$CpG), , drop = FALSE]
  out
}

#Try multiple annotation packages (use whichever are installed)
cand_pkgs <- c(
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
)
anns <- list()
for (pkg in cand_pkgs) {
  ann_obj <- .safe_ann(pkg)
  if (!is.null(ann_obj)) anns[[pkg]] <- norm_ann(ann_obj, pkg)
}
anns <- Filter(Negate(is.null), anns)
if (!length(anns)) stop("No annotation packages available. Install at least one of: ",
                        paste(cand_pkgs, collapse = ", "))

#Find PGR CpGs across available annotations and intersect with our data
has_PGR <- function(x) grepl("(^|;)PGR(;|$)", ifelse(is.na(x), "", x), perl = TRUE)
hits_all <- dplyr::bind_rows(lapply(anns, \(df) dplyr::filter(df, has_PGR(gene)))) %>%
  distinct(CpG, .keep_all = TRUE)

pgr_cpgs <- intersect(hits_all$CpG, rownames(beta_combat))
stopifnot(length(pgr_cpgs) > 0)

pgr_ann <- hits_all %>% filter(CpG %in% pgr_cpgs)

#Build long table & genomic order for y-axis
per_cpg_beta_long <- beta_combat[pgr_cpgs, , drop = FALSE] %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("CpG") %>%
  pivot_longer(-CpG, names_to = "Sample", values_to = "Beta") %>%
  mutate(Group = grp[Sample]) %>%
  left_join(select(pgr_ann, CpG, chr, pos, gene, group), by = "CpG") %>%
  filter(!is.na(Group))

pgr_ann2 <- per_cpg_beta_long %>%
  distinct(CpG, chr, pos) %>%
  mutate(chr_clean = sub("^chr","", chr),
         pos_num   = suppressWarnings(as.numeric(pos)))
chr_levels <- c(as.character(1:22), "X", "Y", "MT", "M")
pgr_ann2 <- pgr_ann2 %>%
  arrange(factor(chr_clean, levels = chr_levels, ordered = TRUE), pos_num)
cpg_order <- pgr_ann2$CpG
per_cpg_beta_long$CpG <- factor(per_cpg_beta_long$CpG, levels = cpg_order)

###-----One-stop labels for butterfly plots-----###
HILITE_LABEL <- "DMPA"                 # what you want printed on the plot
REF_LABEL    <- "Baylor/Heidelberg"

COL_HILITE   <- "#582C83"              # DMPA signature purple
COL_REF      <- "#555555"

###-----Plot butterfly graph - mean β per CpG site by cohort-----###
summ <- per_cpg_beta_long %>%
  group_by(CpG, Group) %>%
  summarise(n = sum(!is.na(Beta)),
            mean_beta = mean(Beta, na.rm = TRUE),
            .groups = "drop")

y_key  <- tibble(CpG = cpg_order, y = seq_along(cpg_order))
summ   <- left_join(summ, y_key, by = "CpG")
leftDF  <- filter(summ, Group == "DMPA")
rightDF <- filter(summ, Group == "Reference")

library(grid)  # for unit()

cols <- c(DMPA = "#582C83", Reference = "#555555")

p <- ggplot() +
  geom_segment(aes(x = 0, xend = 0, y = 0.5, yend = nrow(y_key) + 0.5),
               color = "grey65", linewidth = 0.6) +
  geom_segment(data = y_key,
               aes(x = -0.02, xend = 0.02, y = y, yend = y),
               color = "grey70", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", linewidth = 0.6, color = "grey45") +
  
  # DMPA (left)
  geom_segment(data = leftDF,
               aes(x = 0, xend = -mean_beta, y = y, yend = y),
               color = cols["DMPA"], linewidth = 2.2, lineend = "round") +
  geom_point(data = leftDF,
             aes(x = -mean_beta, y = y),
             color = cols["DMPA"], size = 2.2) +
  
  # Reference (right)
  geom_segment(data = rightDF,
               aes(x = 0, xend = mean_beta, y = y, yend = y),
               color = cols["Reference"], linewidth = 2.2, lineend = "round") +
  geom_point(data = rightDF,
             aes(x = mean_beta, y = y),
             color = cols["Reference"], size = 2.2) +
  
  annotate("text", x = -0.65, y = nrow(y_key) + 1.3, label = HILITE_LABEL,
           color = COL_HILITE, fontface = "bold", size = 5.2, hjust = 0.5) +
  annotate("text", x =  0.65, y = nrow(y_key) + 1.3, label = REF_LABEL,
           color = COL_REF, fontface = "bold", size = 5.2, hjust = 0.5) +
  
  scale_y_continuous(breaks = y_key$y, labels = y_key$CpG,
                     expand = expansion(mult = c(0.02, 0.10))) +
  scale_x_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, by = 0.25),
                     labels = c("1.00","0.75","0.50","0.25","0",
                                "0.25","0.50","0.75","1.00")) +
  labs(
    title = "Progesterone receptor methylation",
    x = "Mean β magnitude",
    y = "CpG sites along PGR gene"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line.x = element_line(linewidth = 1.1, color = "black"),
    axis.line.y = element_line(linewidth = 1.1, color = "black"),
    axis.ticks  = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(5, "pt"),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 15, margin = ggplot2::margin(t = 10)),
    axis.title.y = element_text(size = 15, margin = ggplot2::margin(r = 16)),
    panel.grid = element_blank()
  )

print(p)

#Save
out_dir <- if (exists("ResultsDir")) ResultsDir else getwd()
png_path <- file.path(out_root, "Fig_PGR_butterfly_mean_beta_with_axes.png")
pdf_path <- file.path(out_root, "Fig_PGR_butterfly_mean_beta_with_axes.pdf")
ggsave(png_path, p, width = 9, height = 7.5, dpi = 300)
ggsave(pdf_path, p, width = 9, height = 7.5)
message("Saved: ", normalizePath(png_path))
message("Saved: ", normalizePath(pdf_path))


##==============================================================================
###-----18. 11q22.1 CYTOBAND DIFFERENTIAL METHYLATION ANALYSIS-----###
##==============================================================================

###-----Create butterfly plot figure11q22.1 cytoband- same styling as PGR CpG site figure-----###
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

stopifnot(exists("beta_combat"), exists("pd_df"))
beta_combat <- as.matrix(beta_combat); storage.mode(beta_combat) <- "double"

#Cohorts (same mapping as before)
grp <- factor(ifelse(pd_df$Cohort == "DMPA", "DMPA", "Reference"),
              levels = c("DMPA","Reference"))
names(grp) <- rownames(pd_df)
stopifnot(identical(colnames(beta_combat), names(grp)))

#Get 11q22.1 probes from 450k annotation (hg19), fallback window 100–102 Mb 
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_df  <- as.data.frame(ann450k) |>
  transmute(CpG = Name, chr = as.character(chr), pos = as.numeric(pos),
            Genes = as.character(UCSC_RefGene_Name))

band_rng <- c(start = 100000000, end = 102000000)  # robust fallback around PGR (hg19)
cpg_band <- ann_df |>
  filter(chr == "chr11", pos >= band_rng["start"], pos <= band_rng["end"]) |>
  pull(CpG) |> intersect(rownames(beta_combat))
stopifnot(length(cpg_band) > 0)

#Build gene-level matrix: mean β per gene per sample (only genes with probes) 
probe_gene_tbl <- ann_df |>
  filter(CpG %in% cpg_band, !is.na(Genes), nzchar(Genes)) |>
  select(CpG, Genes, pos) |>
  separate_rows(Genes, sep = ";") |>
  mutate(Genes = trimws(Genes)) |>
  filter(nzchar(Genes))

stopifnot(nrow(probe_gene_tbl) > 0)

#per-gene genomic position (median CpG position per gene) to order y-axis
gene_pos <- probe_gene_tbl |>
  group_by(Genes) |>
  summarise(gene_pos = median(pos, na.rm = TRUE), .groups = "drop")

#aggregate sample-level mean β per gene
mat_sub <- beta_combat[unique(probe_gene_tbl$CpG), , drop = FALSE]
idx_by_gene <- split(probe_gene_tbl$CpG, probe_gene_tbl$Genes)

gene_mat <- vapply(idx_by_gene, function(cpgs) {
  colMeans(mat_sub[intersect(cpgs, rownames(mat_sub)), , drop = FALSE], na.rm = TRUE)
}, FUN.VALUE = numeric(ncol(mat_sub)))
gene_mat <- t(gene_mat)           # genes × samples
colnames(gene_mat) <- colnames(mat_sub)

#Build long table & summaries (mirrors the PGR section)
gene_long <- gene_mat |>
  as.data.frame(check.names = FALSE) |>
  rownames_to_column("Gene") |>
  pivot_longer(-Gene, names_to = "Sample", values_to = "Beta") |>
  mutate(Group = grp[Sample]) |>
  filter(!is.na(Group)) |>
  left_join(gene_pos, by = c("Gene" = "Genes"))

#Y order by genomic position
gene_order <- gene_long |>
  distinct(Gene, gene_pos) |>
  arrange(gene_pos) |>
  pull(Gene)

gene_long$Gene <- factor(gene_long$Gene, levels = gene_order)

#Per-group means (for the butterfly “wings”)
summ <- gene_long |>
  group_by(Gene, Group) |>
  summarise(n = sum(!is.na(Beta)),
            mean_beta = mean(Beta, na.rm = TRUE),
            .groups = "drop")

y_key <- tibble(Gene = gene_order, y = seq_along(gene_order))
summ  <- left_join(summ, y_key, by = "Gene")
leftDF  <- filter(summ, Group == "DMPA")
rightDF <- filter(summ, Group == "Reference")

##lot — EXACT styling as  PGR figure (colors, axes, dashed 0.3, top labels) 
col_DMPA <- "#582C83"
col_REF  <- "#555555"  # darker grey

#top label x-positions: centered between dashed line (±0.3) and axis ends (±1)
left_label_x  <- (-1 + -0.3) / 2   # -0.65
right_label_x <- ( 1 +  0.3) / 2   #  0.65

p_cyto <- ggplot() +
  # center vertical spine and y ticks
  geom_segment(aes(x = 0, xend = 0, y = 0.5, yend = nrow(y_key) + 0.5),
               color = "grey65", linewidth = 0.6) +
  geom_segment(data = y_key,
               aes(x = -0.02, xend = 0.02, y = y, yend = y),
               color = "grey70", linewidth = 0.5) +
  # dashed methylation thresholds
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed",
             linewidth = 0.6, color = "grey45") +
  # DMPA wing (left)
  geom_segment(data = leftDF,
               aes(x = 0, xend = -mean_beta, y = y, yend = y),
               color = col_DMPA, linewidth = 2.2, lineend = "round") +
  geom_point(data = leftDF,
             aes(x = -mean_beta, y = y),
             color = col_DMPA, size = 2.2) +
  # Reference wing (right)
  geom_segment(data = rightDF,
               aes(x = 0, xend = mean_beta, y = y, yend = y),
               color = col_REF, linewidth = 2.2, lineend = "round") +
  geom_point(data = rightDF,
             aes(x = mean_beta, y = y),
             color = col_REF, size = 2.2) +
  # Top cohort labels
  annotate("text", x = -0.65, y = nrow(y_key) + 1.3, label = HILITE_LABEL,
           color = COL_HILITE, fontface = "bold", size = 5.2, hjust = 0.5) +
  annotate("text", x =  0.65, y = nrow(y_key) + 1.3, label = REF_LABEL,
           color = COL_REF, fontface = "bold", size = 5.2, hjust = 0.5) +
  # axes & scales (same as PGR)
  scale_y_continuous(breaks = y_key$y, labels = y_key$Gene,
                     expand = expansion(mult = c(0.02, 0.10))) +
  scale_x_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, by = 0.25),
                     labels = c("1.00","0.75","0.50","0.25","0",
                                "0.25","0.50","0.75","1.00")) +
  labs(
    title = "11q22.1 cytoband methylation",
    x = "Mean β magnitude",
    y = "Genes within 11q22.1 cytoband"
  ) +
  # EXACT same theme details as the PGR plot
  theme_classic(base_size = 16) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.line.x  = element_line(linewidth = 1.1, color = "black"),
    axis.line.y  = element_line(linewidth = 1.1, color = "black"),
    axis.ticks   = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(5, "pt"),
    axis.text.x  = element_text(size = 13),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 15, margin = ggplot2::margin(t = 10)),
    axis.title.y = element_text(size = 15, margin = ggplot2::margin(r = 16)),
  )

print(p_cyto)

#Save
out_root <- file.path(getwd(), "PGR and 11q22 DMR")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)


ggsave(file.path(out_root, "Fig_11q22_1_cytoband_butterfly_mean_beta_by_cohort.png"),
       p_cyto, width = 9, height = 7.5, dpi = 300)
ggsave(file.path(out_root, "Fig_11q22_1_cytoband_butterfly_mean_beta_by_cohort.pdf"),
       p_cyto, width = 9, height = 7.5)
message("Saved cytoband butterfly to: ", normalizePath(out_dir))


##==============================================================================
###-----19. LOLLIPOP PLOTS FOR GLIOSEQ-----###
##==============================================================================

###-----Build enhanced mini-MAF from DepoMeningiomaGlioSeqResults.xlsx-----###
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

#Read the Excel sheet (new-computer path)
in_xlsx <- file.path(base_dir, "DepoMeningiomaGlioSeqResults.xlsx")
raw <- readxl::read_excel(in_xlsx)

#Tidy TRAF7 / PIK3CA / FGFR1 into long format
mut_long <- raw %>%
  dplyr::select(
    DMPA_ID,
    TRAF7_mutation, TRAF7_percent,
    PIK3CA_mutation, PIK3CA_percent,
    FGFR1_mutation, FGFR1_percent
  ) %>%
  tidyr::pivot_longer(
    cols            = -DMPA_ID,
    names_to        = c("Hugo_Symbol", "field"),
    names_pattern   = "(TRAF7|PIK3CA|FGFR1)_(mutation|percent)",
    values_to       = "value",
    values_transform = list(value = as.character)  # fix mixed type issue
  ) %>%
  tidyr::pivot_wider(
    names_from  = field,   # -> 'mutation' and 'percent'
    values_from = value
  ) %>%
  dplyr::mutate(
    percent = as.numeric(percent)
  ) %>%
  dplyr::filter(!is.na(mutation), mutation != "")

#Build the "mini-MAF" style table
mini_maf <- mut_long %>%
  dplyr::mutate(
    Tumor_Sample_Barcode = paste0("DMPA-", DMPA_ID),
    HGVSc       = stringr::str_extract(mutation, "c\\.[0-9]+[ACGT]>[ACGT]"),
    HGVSp_Short = stringr::str_extract(mutation, "p\\.[A-Z][0-9]+[A-Z]"),
    VAF         = percent / 100
  ) %>%
  dplyr::select(
    Hugo_Symbol,
    Tumor_Sample_Barcode,
    DMPA_ID,
    HGVSc,
    HGVSp_Short,
    VAF
  )

#Extract amino-acid position for lollipop x-axis
mini_maf <- mini_maf %>%
  dplyr::mutate(
    AA_position = as.integer(stringr::str_extract(HGVSp_Short, "(?<=p\\.[A-Z])[0-9]+"))
  )

mini_maf

#Save to disk for reproducibility
out_maf <- file.path(base_dir, "DepoMini_miniMAF_withAA.tsv")
readr::write_tsv(mini_maf, out_maf)
cat("Saved enhanced mini-MAF to:\n", out_maf, "\n")


###-----Build cBioPortal-compatible MAF for Mutation Mapper-----###
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

#Load the Excel file
in_xlsx <- file.path(base_dir, "DepoMeningiomaGlioSeqResults.xlsx")
raw <- readxl::read_excel(in_xlsx)

#Tidy mutation data
mut_long <- raw %>%
  select(
    DMPA_ID,
    TRAF7_mutation, TRAF7_percent,
    PIK3CA_mutation, PIK3CA_percent,
    FGFR1_mutation, FGFR1_percent
  ) %>%
  pivot_longer(
    cols = -DMPA_ID,
    names_to = c("Hugo_Symbol", "field"),
    names_pattern = "(TRAF7|PIK3CA|FGFR1)_(mutation|percent)",
    values_to = "value",
    values_transform = list(value = as.character)
  ) %>%
  pivot_wider(names_from = field, values_from = value) %>%
  filter(!is.na(mutation), mutation != "") %>%
  mutate(percent = as.numeric(percent))

#Extract protein change & amino acid position
mini_maf <- mut_long %>%
  mutate(
    Tumor_Sample_Barcode = paste0("DMPA-", DMPA_ID),
    HGVSc       = str_extract(mutation, "c\\.[0-9]+[ACGT]>[ACGT]"),
    Protein_Change = str_extract(mutation, "p\\.[A-Z][0-9]+[A-Z]"),
    AA_position  = as.integer(str_extract(Protein_Change, "(?<=p\\.[A-Z])[0-9]+")),
    VAF          = percent / 100
  )

#Build REQUIRED Mutation Mapper format
cbio_maf <- mini_maf %>%
  transmute(
    Hugo_Symbol,
    Tumor_Sample_Barcode,
    Variant_Classification = "Missense_Mutation",
    Variant_Type           = "SNP",
    Protein_Change,
    HGVSc,
    Reference_Allele       = "",
    Tumor_Seq_Allele2      = "",
    Mutation_Type          = "Missense_Mutation"
  )

#Save in tab-delimited format for easy copy/paste
out_file <- file.path(base_dir, "Depo_CBioPortal_MutationMapper_MAF.tsv")
write_tsv(cbio_maf, out_file)

cat("Saved cBioPortal-compatible MAF to:\n", out_file, "\n\n")
cat("You may now open this file and copy-paste the entire contents into Mutation Mapper.\n")

###-----Lollipop figures: TRAF7, NF2, PIK3CA, FGFR1 (using Track Viewer)-----###
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

###-----Read mini-MAF and summarise mutations-----###

maf_path <- file.path(base_dir, "DepoMini_miniMAF_withAA.tsv")

maf <- readr::read_tsv(maf_path, show_col_types = FALSE)

needed_cols <- c("Hugo_Symbol", "HGVSp_Short", "AA_position")
stopifnot(all(needed_cols %in% names(maf)))

###-----Summarize mutations and build display labels-----###

raw_mut_summary <- maf %>%
  dplyr::filter(Hugo_Symbol %in% c("TRAF7", "NF2", "PIK3CA", "FGFR1")) %>%
  dplyr::mutate(
    AA_Pos = as.numeric(AA_position)
  ) %>%
  dplyr::group_by(Hugo_Symbol, AA_Pos, Protein_Change = HGVSp_Short) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  )

mut_summary <- raw_mut_summary %>%
  dplyr::mutate(
    # final text that will be drawn above each lollipop
    Label = gsub("\\s*\\(n=2\\)", "", Protein_Change)
  )

###-----Domain definitions with different color palate options###

# Muted, but a bit more TCGA-like:
# warm peaches/oranges for functional domains, blues/teals for Ig/kinase, soft violet/pink for WD/TM

col_ring    <- "#F2B38C"  # RING / ABD-ish
col_cc      <- "#E57E5A"  # coiled-coil / RBD / FERM_C
col_wd      <- "#E7C2E8"  # WD repeats
col_fermN   <- "#F7D4A7"  # FERM_N
col_fermM   <- "#A7CFB5"  # FERM_M
col_fermC   <- "#E39A7C"  # FERM_C (warmer)
col_erm     <- "#F7CFAF"  # ERM
col_abd     <- "#F6D2A2"  # ABD
col_rbd     <- "#E9A77F"  # RBD
col_c2      <- "#F9E2B4"  # C2
col_helical <- "#B6D8B6"  # Helical (pale green)
col_kinase  <- "#8FB3E2"  # Kinase (mid blue)
col_ig      <- "#B9D4EF"  # Ig domains (light blue)
col_tm      <- "#C4A6E8"  # TM (soft violet)

# TRAF7 (UniProt Q6Q0C0; ~670 aa)
domains_TRAF7 <- tibble::tribble(
  ~start, ~end, ~label,         ~fill,
  60,  120, "RING",          col_ring,
  220,  320, "Coiled-coil",   col_cc,
  370,  400, "WD1",           col_wd,
  401,  430, "WD2",           col_wd,
  431,  460, "WD3",           col_wd,
  461,  490, "WD4",           col_wd,
  491,  520, "WD5",           col_wd,
  521,  550, "WD6",           col_wd,
  551,  640, "WD7",           col_wd
)
len_TRAF7 <- 670

# NF2 (~595 aa)
domains_NF2 <- tibble::tribble(
  ~start, ~end, ~label,   ~fill,
  20,  120, "FERM_N",  col_fermN,
  121,  310, "FERM_M",  col_fermM,
  311,  412, "FERM_C",  col_fermC,
  430,  595, "ERM",     col_erm
)
len_NF2 <- 595

# PIK3CA (~1068 aa)
domains_PIK3CA <- tibble::tribble(
  ~start, ~end, ~label,  ~fill,
  1,  108, "ABD",    col_abd,
  109,  314, "RBD",    col_rbd,
  335,  526, "C2",     col_c2,
  545,  720, "Helical",col_helical,
  721, 1068, "Kinase", col_kinase
)
len_PIK3CA <- 1068

# FGFR1 (~822 aa)
domains_FGFR1 <- tibble::tribble(
  ~start, ~end, ~label, ~fill,
  30,  110, "Ig1",   col_ig,
  111,  200, "Ig2",   col_ig,
  201,  300, "Ig3",   col_ig,
  370,  410, "TM",    col_tm,
  411,  822, "Kinase",col_kinase
)
len_FGFR1 <- 822

###-----Panel-building function (matched axis weights, stems at bar top)-----###

make_gene_panel <- function(gene, protein_len, domains_df) {
  
  gene_mut <- mut_summary %>% dplyr::filter(Hugo_Symbol == gene)
  
  #layout: x-axis at 0, bar mid a bit higher
  x_axis_y   <- 0.05
  bar_mid_y  <- 0.35
  bar_h      <- 0.18
  bar_ymin   <- bar_mid_y - bar_h / 2
  bar_ymax   <- bar_mid_y + bar_h / 2
  
  #domains slightly taller than grey bar
  dom_ymin   <- bar_mid_y - bar_h * 0.9
  dom_ymax   <- bar_mid_y + bar_h * 0.9
  
  #where the "No mutations detected" text goes
  no_mut_y <- dom_ymax + 0.28
  
  bar_bg <- data.frame(
    xmin = 0,
    xmax = protein_len,
    ymin = bar_ymin,
    ymax = bar_ymax
  )
  
  #bracket positions: 0 at top of bar, 1 and 2 above it
  y0    <- bar_ymax
  ystep <- 0.45
  y1    <- y0 + ystep
  y2    <- y0 + 2 * ystep
  
  #map counts to y positions, plus label positions
  if (nrow(gene_mut)) {
    gene_mut2 <- gene_mut %>%
      dplyr::mutate(
        y_point = dplyr::case_when(
          n == 1 ~ y1,
          n >= 2 ~ y2
        ),
        # default label positions
        label_y = y_point + 0.10,
        label_x = AA_Pos
      )
    
    #TRAF7 tweaks: stack WD1 doublet and nudge p.Y538C label
    if (gene == "TRAF7") {
      gene_mut2 <- gene_mut2 %>%
        dplyr::mutate(
          label_y = dplyr::case_when(
            AA_Pos == 389 ~ y_point + 0.10,
            AA_Pos == 390 ~ y_point + 0.30,
            TRUE          ~ label_y
          ),
          label_x = dplyr::case_when(
            Label == "p.Y538C" ~ AA_Pos + 12,
            TRUE               ~ label_x
          )
        )
    }
  } else {
    gene_mut2 <- gene_mut
  }
  
  #bracket geometry in x-units
  bracket_x    <- -protein_len * 0.035
  tick_len     <-  protein_len * 0.0125
  label_x_num  <-  bracket_x + tick_len + protein_len * 0.01
  
  #positions for left-side text
  x_no_mut <- bracket_x - protein_len * 0.02
  x_gene   <- bracket_x - protein_len * 0.05
  
  #build plot first
  p <- ggplot() +
    # x-axis line
    geom_segment(
      aes(x = 0, xend = protein_len, y = x_axis_y, yend = x_axis_y),
      linewidth = 0.8
    ) +
    
    #grey bar (gene backbone)
    geom_rect(
      data = bar_bg,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill    = "grey95",
      colour  = "black",
      linewidth = 0.6
    ) +
    
    #domains
    geom_rect(
      data = domains_df,
      aes(xmin = start, xmax = end, ymin = dom_ymin, ymax = dom_ymax, fill = fill),
      colour   = "black",
      linewidth = 0.6
    ) +
    geom_text(
      data = domains_df,
      aes(x = (start + end) / 2, y = bar_mid_y, label = label),
      size = 4, fontface = "bold"
    ) +
    
    #lollipops
    {
      if (nrow(gene_mut2)) {
        list(
          geom_segment(
            data = gene_mut2,
            aes(x = AA_Pos, xend = AA_Pos, y = dom_ymax, yend = y_point),
            colour   = "black",
            linewidth = 0.8
          ),
          geom_point(
            data = gene_mut2,
            aes(x = AA_Pos, y = y_point),
            colour = "#d7191c",
            fill   = "#d7191c",
            shape  = 21,
            stroke = 0.3,
            size   = 3.5
          ),
          geom_text(
            data = gene_mut2,
            aes(x = label_x, y = label_y, label = Label),
            size     = 4,
            fontface = "bold",
            vjust    = 0,
            nudge_y  = 0.03
          )
        )
      }
    } +
    
    #y-axis bracket
    geom_segment(aes(x = bracket_x, xend = bracket_x, y = y0, yend = y2), linewidth = 0.5) +
    geom_segment(aes(x = bracket_x, xend = bracket_x + tick_len, y = y0, yend = y0), linewidth = 0.5) +
    geom_segment(aes(x = bracket_x, xend = bracket_x + tick_len, y = y1, yend = y1), linewidth = 0.5) +
    geom_segment(aes(x = bracket_x, xend = bracket_x + tick_len, y = y2, yend = y2), linewidth = 0.5) +
    geom_text(
      data = data.frame(y = c(y0, y1, y2), lab = c("0","1","2")),
      aes(x = label_x_num, y = y, label = lab),
      hjust = 0, vjust = 0.4,
      size  = 4, fontface = "bold"
    ) +
    
    #"No. mutations"
    geom_text(
      aes(x = x_no_mut, y = (y0 + y2) / 2, label = "No. mutations"),
      angle = 90, vjust = 0.5, size = 4, fontface = "bold"
    ) +
    
    #gene name
    geom_text(
      aes(x = x_gene, y = (y0 + y2) / 2, label = gene),
      hjust = 1, vjust = 0.5,
      size = 4.5, fontface = "bold.italic"
    ) +
    
    coord_cartesian(
      xlim = c(bracket_x - protein_len * 0.18, protein_len * 1.04),
      ylim = c(x_axis_y, y2 + 0.5),
      expand = FALSE
    ) +
    
    theme_bw(base_size = 13) +
    theme(
      plot.title      = element_blank(),
      axis.title      = element_blank(),
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      panel.grid      = element_blank(),
      panel.border    = element_blank(),
      axis.line       = element_blank(),
      legend.position = "none",
      plot.margin     = margin(t = 5, r = 15, b = 5, l = 60)
    )
  
  #NOW add "No mutations detected" (after p exists) ---
  if (nrow(gene_mut) == 0) {
    p <- p +
      geom_text(
        aes(x = protein_len / 2, y = no_mut_y, label = "No mutations detected"),
        size = 4,
        colour = "grey45",
        fontface = "italic"
      )
  }
  
  p
}

###-----Build panels & save combined figure-----###

p_traf  <- make_gene_panel("TRAF7",  len_TRAF7,  domains_TRAF7)
p_nf2   <- make_gene_panel("NF2",    len_NF2,    domains_NF2)
p_pik3  <- make_gene_panel("PIK3CA", len_PIK3CA, domains_PIK3CA)
p_fgfr1 <- make_gene_panel("FGFR1",  len_FGFR1,  domains_FGFR1)

combined <- gridExtra::grid.arrange(p_traf, p_nf2, p_pik3, p_fgfr1, ncol = 1)

out_dir <- file.path(base_dir, "Combined Analysis_WHO I_Final/Lollipops")

ggsave(file.path(out_dir, "Lollipop_TRAF7_NF2_PIK3CA_FGFR1.pdf"),
       combined, width = 13, height = 7, device = cairo_pdf)

ggsave(file.path(out_dir, "Lollipop_TRAF7_NF2_PIK3CA_FGFR1.png"),
       combined, width = 13, height = 7, dpi = 300, bg = "white")


##==============================================================================
###-----20. ONCOPLOT FOR GLIOSEQ MUTATION ANALYSIS -----###
##==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(grid)   # for textGrob / annotation_custom
})

###-----Input + colors-----###
in_xlsx <- file.path(base_dir, "DepoMeningiomaGlioSeqResults.xlsx")

# Match your copy-number palette
col_nomut   <- "#EFEFEF"  # light gray tiles
col_mut     <- "#08306B"   # deep navy tiles
tile_border <- "#C5C5C5"  # slightly darker gray outline


###-----Read and reshape GlioSeq mutation data-----###

# Genes (will be reversed so TRAF7 is on top, AKT1 on bottom)
genes_alpha <- c(
  "TRAF7", "NF2", "AKT1", "PIK3CA", "SMO",
  "POLR2A","TERT", "CDKN2A", "SMARCB1", "FGFR1", "KLF4", "SUFU"
)

gene_pattern <- paste(genes_alpha, collapse = "|")

raw <- read_excel(in_xlsx)

mut_long <- raw %>%
  select(
    DMPA_ID,
    matches(paste0("^(", gene_pattern, ")_(mutation|percent)$"))
  ) %>%
  pivot_longer(
    cols          = -DMPA_ID,
    names_to      = c("Gene", "field"),
    names_pattern = paste0("(", gene_pattern, ")_(mutation|percent)"),
    values_to     = "value",
    values_transform = list(value = as.character)
  ) %>%
  pivot_wider(names_from = field, values_from = value) %>%
  mutate(
    Sample = as.character(DMPA_ID),
    Gene   = as.character(Gene)
  )

# Rows with an actual mutation
mut_any <- mut_long %>%
  filter(!is.na(mutation), mutation != "") %>%
  transmute(
    Gene,
    Sample,
    status = "Mutation"
  )

###-----Full gene × sample grid (fill in "No mutation")-----###
# Sort by numeric ID, not alphabetically — alphabetic sort puts DMPA-10 between DMPA-1 and DMPA-2
all_samples <- mut_long %>%
  dplyr::distinct(Sample) %>%
  dplyr::mutate(num = as.numeric(stringr::str_extract(Sample, "\\d+"))) %>%
  dplyr::arrange(num) %>%
  dplyr::pull(Sample)

genes_rev <- rev(genes_alpha)

grid_df <- tidyr::expand_grid(
  Gene   = genes_rev,
  Sample = all_samples
) %>%
  left_join(mut_any, by = c("Gene","Sample")) %>%
  mutate(
    status = if_else(is.na(status), "No mutation", status),
    Gene   = factor(Gene,   levels = genes_rev),
    Sample = factor(Sample, levels = all_samples)
  )

n_genes   <- length(genes_rev)
n_samples <- length(all_samples)

###-----Base heatmap-style oncoprint-----###

p <- ggplot(grid_df, aes(x = Sample, y = Gene, fill = status)) +
  geom_tile(
    color  = tile_border,
    width  = 0.70,  # rectangles taller than wide
    height = 0.90
  ) +
  scale_fill_manual(
    values = c("No mutation" = col_nomut,
               "Mutation"    = col_mut)
  ) +
  scale_x_discrete(
    position = "top",
    expand   = expansion(mult = c(0.001, 0.001)),
    labels = function(x) stringr::str_extract(x, "\\d+")  # extract the numeric ID from each sample name
  ) +
  scale_y_discrete(
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_equal(clip = "off") +  # allow label/legend below and to the left
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid       = element_blank(),
    axis.ticks       = element_blank(),
    axis.text.y      = element_text(face = "italic", size = 11, hjust = 1),
    axis.text.x      = element_text(size = 11, face = "bold", colour = "black"),
    plot.title       = element_blank(),
    legend.position  = "none",              
    plot.margin = ggplot2::margin(t = 40, r = 40, b = 20, l = 80),
    panel.border     = element_blank()
  )


###-----Horizontal rule under numbers (unchanged)-----###
p <- p +
  annotate(
    "segment",
    x    = 0.5,
    xend = n_samples + 0.5,
    y    = n_genes + 0.60,
    yend = n_genes + 0.60,
    linewidth = 0.6
  )


###-----“DMPA patient ID” label as a grob (same as before)-----###

label_grob <- textGrob(
  "Patient ID",
  x = unit(-0.02, "npc"),  # move right, closer to numbers
  y = unit(1.02, "npc"),  # align vertically with x-axis labels
  just = "right",
  gp = gpar(fontface = "bold", fontsize = 12)
)

p <- p +
  annotation_custom(
    grob = label_grob,
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf
  )

###-----Legend-----###

legend_label <- "Missense mutation"

# Visual sizing
box_w <- unit(6.0 * 1.0, "mm")   # 20% smaller than 6 mm
box_h <- unit(7.5 * 1.0, "mm")   # 20% smaller than 7.5 mm
gap   <- unit(2.4, "mm")         # small gap between box and text

# Move legend closer to plot:
# More negative = lower; less negative = closer to plot
legend_y_npc <- unit(-0.055, "npc")

# Text grob (match your ggplot styling)
tg <- textGrob(
  legend_label,
  x = unit(0, "npc"),
  just = "left",
  gp = gpar(fontface = "bold", fontsize = 12)
)

# Compute true total width of (box + gap + text)
total_w <- box_w + gap + grobWidth(tg)

# Left edge of legend block so that its CENTER is at 0.5 npc
x_left <- unit(0.5, "npc") - total_w / 2

legend_grob <- grobTree(
  # square
  rectGrob(
    x = x_left + box_w/2,
    y = legend_y_npc,
    width  = box_w,
    height = box_h,
    gp = gpar(fill = col_mut, col = tile_border, lwd = 1)
  ),
  # text (tight to square)
  textGrob(
    legend_label,
    x = x_left + box_w + gap,
    y = legend_y_npc,
    just = "left",
    gp = gpar(fontface = "bold", fontsize = 12)
  )
)

p <- p +
  annotation_custom(
    grob = legend_grob,
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf
  )

###-----Print and/or save-----###
p  # show in RStudio

ggsave(
  file.path(base_dir, "Manuscript/Figures/Mutation_Oncoprint.png"),
  p,
  width = 7,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(base_dir, "Manuscript/Figures/Mutation_Oncoprint.pdf"),
  p,
  width  = 7,
  height = 8,
  device = cairo_pdf
)


##==============================================================================
###-----21. MASTER ONCOPRINT-----###
##==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

###-----GLOBAL VISUAL SETTINGS-----###

tile_h_mm <- 6
tile_w_mm <- 6   # <<< square tiles (set to 5 for taller-than-wide)
anno_height <- unit(tile_h_mm, "mm")

bg_no <- "grey96"

#Inside-cell shrink factors (keep; these control tiny internal padding)
tile_w_frac <- 0.92
tile_h_frac <- 0.95
col_gap     <- unit(0, "mm")

#Gap between heatmap and legend (smaller = closer)
heatmap_right_pad_mm <- 1.0

#Legend area width (smaller = legend closer, but may wrap)
legend_width_mm <- 74

#FIXED BODY WIDTH (this is the square-tile magic)
hm_body_width <- function(n_cols) unit(n_cols * tile_w_mm, "mm")


###-----LOAD DATA-----###
in_xlsx <- file.path(base_dir, "26.4.14 Depo Meningioma Data.xlsx")
dat_raw <- read_excel(in_xlsx, sheet = "1. DataSH", skip = 1)


###-----RESTRICT TO DMPA 1–10-----###

dat <- dat_raw %>%
  filter(DeidentifiedName %in% paste0("DMPA-", 1:10)) %>%
  mutate(
    SampleID   = factor(DeidentifiedName, levels = paste0("DMPA-", 1:10)),
    PatientNum = as.numeric(str_extract(DeidentifiedName, "\\d+"))
  ) %>%
  arrange(SampleID)

sample_ids     <- as.character(dat$SampleID)
patient_labels <- dat$PatientNum

stopifnot(length(sample_ids) > 0)
hm_w <- hm_body_width(length(sample_ids))


###-----MUTATION MATRIX-----###

mut_long <- dat %>%
  select(SampleID, TRAF7, PIK3CA, NF2, FGFR1) %>%
  pivot_longer(cols = -SampleID, names_to = "Gene", values_to = "mut") %>%
  mutate(status = if_else(!is.na(mut) & mut != "" & mut != "No mutations", "Mut", ""))

mat_mut <- mut_long %>%
  select(SampleID, Gene, status) %>%
  pivot_wider(names_from = SampleID, values_from = status) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

mat_mut <- mat_mut[, sample_ids, drop = FALSE]

col_mut <- c("Mut" = "#d7191c")

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w * tile_w_frac, h * tile_h_frac,
              gp = gpar(fill = bg_no, col = NA))
  },
  Mut = function(x, y, w, h) {
    grid.rect(x, y, w * tile_w_frac, h * tile_h_frac,
              gp = gpar(fill = col_mut["Mut"], col = NA))
  }
)

###-----CNV MATRIX-----###

cnv_cols <- c("1p","1q","2p","2q","7p","7q","13q","18p","18q","22q")

cnv_mat <- dat %>%
  select(SampleID, all_of(cnv_cols)) %>%
  mutate(across(all_of(cnv_cols),
                ~ dplyr::case_when(
                  !is.na(.) & . == "loss" ~ "CNV loss",
                  !is.na(.) & . == "gain" ~ "CNV gain",
                  TRUE                    ~ "CNV stable"
                ))) %>%
  pivot_longer(cols = all_of(cnv_cols), names_to = "arm", values_to = "status") %>%
  pivot_wider(names_from = SampleID, values_from = status) %>%
  column_to_rownames("arm") %>%
  as.matrix()

cnv_mat <- cnv_mat[, sample_ids, drop = FALSE]

cnv_col <- c(
  "CNV stable" = bg_no,
  "CNV loss"   = "#b2182b",
  "CNV gain"   = "#2166AC"
)


###-----CLINICAL & MOLECULAR ANNOTATION DATA-----###
bin_to_yesno <- function(x) {
  out <- ifelse(x == 1, "Yes", "No")
  out[is.na(x)] <- "No"
  factor(out, levels = c("No", "Yes"))
}

#Robustly pull the Histopathology(Dan) column even if Excel duplicates rename it
histo_col <- c("Histopathology", "Histology")
histo_col <- histo_col[histo_col %in% names(dat)][1]
stopifnot(!is.na(histo_col))

histo_vec <- dat[[histo_col]]          # <<< REQUIRED
histo_vec <- trimws(histo_vec)

#Normalize non-breaking spaces (Excel loves these)
histo_vec <- gsub("\u00A0", " ", histo_vec)

#Normalize back to older wording
histo_vec <- gsub("\\(\\+myxoid features\\)", "with myxoid features", histo_vec)
histo_vec <- trimws(gsub("\\s+", " ", histo_vec))

anno_df <- data.frame(
  Age         = as.numeric(dat$Age),
  Sex         = factor(dat$Sex),
  WHO         = factor(dat$`WHO grade`),
  Histology   = factor(histo_vec),
    DMPACluster = factor(dat$DMPACluster),
  MenG        = factor(dat$MenG_Prediction),
  DKFZ2       = as.numeric(dat$DKFZ_2_Score),
  DKFZ3       = as.numeric(dat$DKFZ_3_Score),
  MultiMM     = bin_to_yesno(dat$Multiple_meningiomas),
  Regression  = bin_to_yesno(dat$Regression),
  SB_Present  = bin_to_yesno(dat$Skull_base_tumor_present),
  SB_Sequenced = bin_to_yesno(dat$Skull_base_tumor_sequenced),
  row.names   = dat$SampleID
)

###-----COLORS-----###
PALETTE <- "jci"  
# options: "nature", "nejm", "jci"

options(ComplexHeatmap.use_raster = FALSE)

bg_no <- "#F5F5F5"

if (PALETTE == "nature") {
  # Nature / Cell / Science (muted, modern, balanced)
  age_fun <- circlize::colorRamp2(
    range(anno_df$Age, na.rm = TRUE),
    c("#E8EEF4", "#274C77")
  )
  
  dkfz_fun <- circlize::colorRamp2(
    c(0, 0.5, 1),
    c("#FAFAFA", "#E9A44C", "#8E1B1B")
  )
  
  sex_col <- c(F = "#B71C1C", M = "#1E3A8A")
  
  who_col <- c(
    "1" = "#77C4B2",
    "2" = "#F2A272",
    "3" = "#7E8BC2"
  )
  
  cluster_base <- c("#4E79A7", "#F28E2B", "#59A14F", "#B07AA1")
  
  yesno_col <- c(No = bg_no, Yes = "#1B5E20")
  
  col_mut <- c(Mut = "#B71C1C")
  
  cnv_col <- c(
    "CNV stable" = bg_no,
    "CNV loss"   = "#8E1B1B",
    "CNV gain"   = "#274C77"
  )
  
  histology_levels <- levels(anno_df$Histology)
  
  
} else if (PALETTE == "nejm") {
  
  #NEJM / JAMA (high contrast, conservative, clinical)
  age_fun <- circlize::colorRamp2(
    range(anno_df$Age, na.rm = TRUE),
    c("#EDF2F7", "#08306B")
  )
  
  dkfz_fun <- circlize::colorRamp2(
    c(0, 0.5, 1),
    c("#FFFFFF", "#FDAE61", "#67000D")
  )
  
  sex_col <- c(F = "#9E1B32", M = "#003A8F")
  
  who_col <- c(
    "1" = "#66C2A5",
    "2" = "#FC8D62",
    "3" = "#8DA0CB"
  )
  
  cluster_base <- c("#08519C", "#CB181D", "#238B45")
  
  yesno_col <- c(No = bg_no, Yes = "#00441B")
  
  col_mut <- c(Mut = "#9E1B32")
  
  cnv_col <- c(
    "CNV stable" = bg_no,
    "CNV loss"   = "#67000D",
    "CNV gain"   = "#084594"
  )
  
} else if (PALETTE == "jci") {
  # JCI / Science Translational (clean, readable, slightly warmer)
  age_fun <- circlize::colorRamp2(
    range(anno_df$Age, na.rm = TRUE),
    c("#EFF3F8", "#336699")
  )
  
  dkfz_fun <- circlize::colorRamp2(
    c(0, 0.5, 1),
    c("#FFFFFF", "#FDB863", "#B2182B")
  )
  
  sex_col <- c(F = "#C0392B", M = "#2471A3")
  
  who_col <- c(
    "1" = "#73C6B6",
    "2" = "#F0B27A",
    "3" = "#85929E"
  )
  
  cluster_base <- c("#2874A6", "#D68910", "#229954")
  
  yesno_col <- c(No = bg_no, Yes = "#1E8449")
  
  col_mut <- c(Mut = "#C0392B")
  
  cnv_col <- c(
    "CNV stable" = bg_no,
    "CNV loss"   = "#922B21",
    "CNV gain"   = "#1F618D"
  )
}

histology_col <- c(
  "Metaplastic"                               = "#F28E2B",
  "Meningothelial"                            = "#4E79A7",
  "Meningothelial with myxoid features"     = "#A0CBE8",
  "Transitional"                              = "#59A14F",
  "Transitional with myxoid features"       = "#8CD17D",
  "Secretory"                                 = "#C39A6B"
)

#Derived mappings
cluster_levels <- levels(anno_df$DMPACluster)
cluster_col <- setNames(cluster_base[seq_along(cluster_levels)], cluster_levels)

meng_levels <- levels(anno_df$MenG)
meng_col <- setNames(rep(who_col["1"], length(meng_levels)), meng_levels)

###-----ANNOTATION MATRICES (1 x 10)-----###

m_age   <- matrix(anno_df$Age,         nrow = 1, dimnames = list("Age", sample_ids))
m_sex   <- matrix(anno_df$Sex,         nrow = 1, dimnames = list("Sex", sample_ids))
m_who   <- matrix(anno_df$WHO,         nrow = 1, dimnames = list("WHO Grade", sample_ids))
m_sbs   <- matrix(anno_df$SB_Sequenced, nrow = 1, dimnames = list("SB Tumor Sequenced", sample_ids))
m_sbp   <- matrix(anno_df$SB_Present,  nrow = 1, dimnames = list("SB Tumor Present", sample_ids))
m_mult  <- matrix(anno_df$MultiMM,     nrow = 1, dimnames = list("Multiple Meningiomas", sample_ids))
m_reg   <- matrix(anno_df$Regression,  nrow = 1, dimnames = list("Regression After Stopping DMPA", sample_ids))
m_hist <- matrix(anno_df$Histology, nrow = 1, dimnames = list("Histology", sample_ids))
m_clust <- matrix(anno_df$DMPACluster, nrow = 1, dimnames = list("Methylation Cluster", sample_ids))
m_meng  <- matrix(anno_df$MenG,        nrow = 1, dimnames = list("MenG Group", sample_ids))
m_dkfz2 <- matrix(anno_df$DKFZ2,       nrow = 1, dimnames = list("DKFZ Ben-2", sample_ids))
m_dkfz3 <- matrix(anno_df$DKFZ3,       nrow = 1, dimnames = list("DKFZ Ben-3", sample_ids))

###-----cell_fun HELPERS-----###

cell_fun_continuous <- function(mat, col_fun) {
  function(j, i, x, y, w, h, col) {
    val  <- mat[i, j]
    fill <- if (is.na(val)) bg_no else col_fun(val)
    grid.rect(x, y, w * tile_w_frac, h * tile_h_frac, gp = gpar(fill = fill, col = NA))
  }
}

cell_fun_yesno <- function(mat, col_map) {
  function(j, i, x, y, w, h, col) {
    val  <- as.character(mat[i, j])
    fill <- if (is.na(val) || val == "No") bg_no else col_map["Yes"]
    grid.rect(x, y, w * tile_w_frac, h * tile_h_frac, gp = gpar(fill = fill, col = NA))
  }
}

cell_fun_cat <- function(mat, col_map) {
  function(j, i, x, y, w, h, col) {
    val  <- as.character(mat[i, j])
    fill <- col_map[val]
    if (is.na(val) || is.na(fill)) fill <- bg_no
    grid.rect(x, y, w * tile_w_frac, h * tile_h_frac,
              gp = gpar(fill = fill, col = NA))
  }
}

###-----TOP ANNOTATION (Patient ID)-----###

top_anno_pid <- HeatmapAnnotation(
  `Patient ID` = anno_text(
    patient_labels,
    rot  = 0,
    gp   = gpar(fontsize = 10, fontface = "bold"),
    just = "center"
  ),
  annotation_name_side   = "left",
  annotation_name_rot    = 0,
  annotation_name_gp     = gpar(fontsize = 10, fontface = "bold"),
  annotation_name_offset = unit(2, "mm"),
  height                 = unit(tile_h_mm, "mm")
)

###-----HEATMAP BUILDERS (ALL WITH FIXED WIDTH)-----###

make_hm_cat <- function(mat, title, col_map) {
  Heatmap(
    mat,
    name = title,
    col  = col_map,
    cluster_rows = FALSE, cluster_columns = FALSE,
    rect_gp = gpar(fill = NA, col = NA),
    row_names_gp = gpar(fontsize = 10),
    show_row_names = TRUE, row_names_side = "left",
    show_column_names = FALSE,
    column_gap = col_gap,
    height = anno_height,
    width  = hm_w,                         # <<< FIXED COLUMN WIDTH
    cell_fun = cell_fun_cat(mat, col_map),
    show_heatmap_legend = FALSE
  )
}

make_hm_yesno <- function(mat, title) {
  Heatmap(
    mat,
    name = title,
    col  = yesno_col,
    cluster_rows = FALSE, cluster_columns = FALSE,
    rect_gp = gpar(fill = NA, col = NA),
    row_names_gp = gpar(fontsize = 10),
    show_row_names = TRUE, row_names_side = "left",
    show_column_names = FALSE,
    column_gap = col_gap,
    height = anno_height,
    width  = hm_w,                         # <<< FIXED COLUMN WIDTH
    cell_fun = cell_fun_yesno(mat, yesno_col),
    show_heatmap_legend = FALSE
  )
}

make_hm_cont <- function(mat, title, col_fun) {
  Heatmap(
    mat,
    name = title,
    col  = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    rect_gp = gpar(fill = NA, col = NA),
    row_names_gp = gpar(fontsize = 10),
    show_row_names = TRUE, row_names_side = "left",
    show_column_names = FALSE,
    column_gap = col_gap,
    height = anno_height,
    width  = hm_w,                         # <<< FIXED COLUMN WIDTH
    cell_fun = cell_fun_continuous(mat, col_fun),
    show_heatmap_legend = FALSE
  )
}

m_pid <- matrix(patient_labels, nrow = 1,
                dimnames = list("Patient ID", sample_ids))

ht_pid <- Heatmap(
  m_pid,
  name = "Patient ID",
  col  = c("1" = "white"),               # dummy
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  rect_gp = gpar(fill = NA, col = NA),
  row_names_side = "left",
  row_names_gp   = gpar(fontsize = 10, fontface = "bold"),
  height = anno_height,
  width  = hm_w,
  cell_fun = function(j, i, x, y, w, h, col) {
    # optional: very light background to match your style
    grid.rect(x, y, w * tile_w_frac, h * tile_h_frac,
              gp = gpar(fill = "white", col = NA))
    grid.text(m_pid[i, j], x, y, gp = gpar(fontsize = 10, fontface = "bold"))
  }
)

ht_age <- Heatmap(
  m_age,
  name = "Age",
  col  = age_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  rect_gp = gpar(fill = NA, col = NA),
  row_names_gp = gpar(fontsize = 10),
  show_row_names = TRUE, row_names_side = "left",
  show_column_names = FALSE,
  column_gap = col_gap,
  height = anno_height,
  width  = hm_w,                           # <<< FIXED COLUMN WIDTH
  cell_fun = cell_fun_continuous(m_age, age_fun),
  show_heatmap_legend = FALSE
)

ht_sex   <- make_hm_cat(m_sex, "Sex", sex_col)
ht_who   <- make_hm_cat(m_who, "WHO Grade", who_col)
ht_sbs   <- make_hm_yesno(m_sbs, "SB Tumor Sequenced")
ht_sbp   <- make_hm_yesno(m_sbp, "SB Tumor Present")
ht_mult  <- make_hm_yesno(m_mult, "Multiple Meningiomas")
ht_reg   <- make_hm_yesno(m_reg, "Regression")
ht_hist <- make_hm_cat(m_hist, "Histopathology", histology_col)
ht_clust <- make_hm_cat(m_clust, "Methylation Cluster", cluster_col)
ht_meng  <- make_hm_cat(m_meng, "MenG Group", meng_col)
ht_dkfz2 <- make_hm_cont(m_dkfz2, "DKFZ Ben-2", dkfz_fun)
ht_dkfz3 <- make_hm_cont(m_dkfz3, "DKFZ Ben-3", dkfz_fun)

ht_mut <- oncoPrint(
  mat_mut,
  alter_fun = alter_fun,
  alter_fun_is_vectorized = FALSE,
  col = col_mut,
  
  #italic gene labels
  row_names_side = "left",
  row_names_gp   = gpar(fontsize = 10, fontface = "italic"),
  
  remove_empty_columns = FALSE,
  show_pct = FALSE,
  
  #remove barplots cleanly
  top_annotation   = NULL,
  right_annotation = NULL,
  
  column_order = sample_ids,
  show_column_names = FALSE,
  column_gap = col_gap,
  
  #FORCE tile geometry to match others
  height  = unit(nrow(mat_mut) * tile_h_mm, "mm"),
  row_gap = unit(0, "mm"),   # prevents extra vertical slack
  width   = hm_w,
  
  show_heatmap_legend = FALSE
)

ht_cnv <- Heatmap(
  cnv_mat,
  name  = "Copy number alteration",
  col   = cnv_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp   = gpar(fontsize = 10),
  show_column_names = FALSE,
  rect_gp = gpar(fill = NA, col = NA),
  column_gap = col_gap,
  height = unit(nrow(cnv_mat) * tile_h_mm, "mm"),
  width  = hm_w,                           # <<< FIXED COLUMN WIDTH
  cell_fun = function(j, i, x, y, w, h, col) {
    val  <- cnv_mat[i, j]
    fill <- cnv_col[val]
    grid.rect(x, y, w * tile_w_frac, h * tile_h_frac, gp = gpar(fill = fill, col = NA))
  },
  show_heatmap_legend = FALSE
)

###-----SEPARATORS-----###

sep_mat <- matrix(1, nrow = 1, ncol = length(sample_ids), dimnames = list("", sample_ids))

sep_ht <- function(name) {
  Heatmap(
    sep_mat,
    name = name,
    col  = c("1" = "grey90"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    rect_gp = gpar(col = NA),
    height = unit(1.2, "mm"),
    width  = hm_w,                         # <<< keep aligned
    column_gap = col_gap,
    show_heatmap_legend = FALSE
  )
}
ht_sep_reg_hist <- sep_ht("sep_reg_hist")   # <<< ADD
ht_sep_top_mid  <- sep_ht("sep_top_mid")
ht_sep_clin_mut <- sep_ht("sep_clin_mut")
ht_sep_mut_cnv  <- sep_ht("sep_mut_cnv")

###-----STACK-----###

ht_all <- 
  ht_pid %v% ht_age %v%
  ht_sex %v% ht_who %v%
  ht_sbs %v% ht_sbp %v% ht_mult %v% ht_reg %v%
  ht_sep_reg_hist %v%          # <<< NEW GREY LINE
  ht_hist %v%                  # <<< NEW HISTOPATHOLOGY ROW
  ht_sep_top_mid %v%           # <<< EXISTING GREY LINE (now between histology and methylation)
  ht_clust %v% ht_meng %v% ht_dkfz2 %v% ht_dkfz3 %v%
  ht_sep_clin_mut %v%
  ht_mut %v%
  ht_sep_mut_cnv %v%
  ht_cnv

###-----LEGENDS (with grey separator bars + corrected spacing)-----###

lgd_title_gp <- gpar(fontface = "bold", fontsize = 8.5, lineheight = 0.95)
lgd_label_gp <- gpar(fontsize = 7.7)

#Thin grey horizontal separator bar (like your oncoprint separators)
lgd_sep_bar <- function(w_mm = 40, h_mm = 1.4, fill = "grey90") {
  Legend(
    title  = "",
    labels = "",
    legend_gp   = gpar(fill = fill, col = NA),
    grid_width  = unit(w_mm, "mm"),
    grid_height = unit(h_mm, "mm"),
    labels_gp   = gpar(col = "#FFFFFF00", fontsize = 1),
    title_gp    = gpar(col = "#FFFFFF00", fontsize = 1)
  )
}

lgd_age <- Legend(
  title = "Age", col_fun = age_fun, direction = "horizontal",
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp,
  legend_width = unit(20, "mm")
)
lgd_sex <- Legend(
  title = "Sex", labels = c("F","M"),
  legend_gp = gpar(fill = sex_col[c("F","M")]),
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp
)
lgd_who <- Legend(
  title = "WHO Grade", labels = names(who_col),
  legend_gp = gpar(fill = who_col),
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp
)
lgd_sbs <- Legend(
  title = "SB Tumor\nSequenced", labels = c("No","Yes"),
  legend_gp = gpar(fill = yesno_col[c("No","Yes")]),
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp
)
lgd_sbp <- Legend(
  title = "SB Tumor\nPresent", labels = c("No","Yes"),
  legend_gp = gpar(fill = yesno_col[c("No","Yes")]),
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp
)
lgd_mult <- Legend(
  title = "Multiple\nmeningiomas", labels = c("No","Yes"),
  legend_gp = gpar(fill = yesno_col[c("No","Yes")]),
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp
)
lgd_reg <- Legend(
  title = "Regression After\nStopping DMPA", labels = c("No","Yes"),
  legend_gp = gpar(fill = yesno_col[c("No","Yes")]),
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp
)

legend_left <- packLegend(
  lgd_age, lgd_sex, lgd_who, lgd_sbs, lgd_sbp, lgd_mult, lgd_reg,
  direction = "vertical",
  gap = unit(1.8, "mm")
)

lgd_hist <- Legend(
  title = "Histopathology",
  labels = names(histology_col),
  legend_gp = gpar(fill = histology_col),
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp
)

lgd_clust <- Legend(
  title = "Methylation\nCluster", labels = names(cluster_col),
  legend_gp = gpar(fill = cluster_col),
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp
)
lgd_meng <- Legend(
  title = "MenG Group", labels = names(meng_col),
  legend_gp = gpar(fill = meng_col),
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp
)
lgd_dkfz2 <- Legend(
  title = "DKFZ Ben-2", col_fun = dkfz_fun, at = c(0,0.5,1),
  direction = "horizontal",
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp,
  legend_width = unit(20, "mm")
)
lgd_dkfz3 <- Legend(
  title = "DKFZ Ben-3", col_fun = dkfz_fun, at = c(0,0.5,1),
  direction = "horizontal",
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp,
  legend_width = unit(20, "mm")
)
lgd_mut <- Legend(
  title = "Mutation Present", labels = "Missense",
  legend_gp = gpar(fill = col_mut["Mut"]),
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp
)
lgd_cnv <- Legend(
  title = "Copy Number\nAlterations",
  labels = c("CNV stable","CNV loss","CNV gain"),
  legend_gp = gpar(fill = cnv_col[c("CNV stable","CNV loss","CNV gain")]),
  title_gp = lgd_title_gp, labels_gp = lgd_label_gp
)


###-----LEGEND PACKING WITH SECTION DIVIDER BARS + CUSTOM SPACING-----###
#ht_mult display title capitalization (matches rowname) 
ht_mult  <- make_hm_yesno(m_mult, "Multiple Meningiomas")

#A “fake legend” that draws as a horizontal light-gray divider bar
lgd_bar <- function(w_mm = 30, h_mm = 1.2, col = "grey90") {
  Legend(
    title  = "",
    labels = "",
    legend_gp = gpar(fill = col, col = NA),
    grid_width  = unit(w_mm, "mm"),
    grid_height = unit(h_mm, "mm"),
    title_gp  = gpar(col = "#FFFFFF00", fontsize = 1),
    labels_gp = gpar(col = "#FFFFFF00", fontsize = 1),
    direction = "horizontal",
    legend_width  = unit(w_mm, "mm"),
    legend_height = unit(h_mm, "mm")
  )
}

#Left legend column (clinical)
legend_left <- packLegend(
  lgd_age, lgd_sex, lgd_who, lgd_sbr, lgd_sbp, lgd_mult, lgd_reg,
  direction = "vertical",
  gap = unit(1.8, "mm")
)

#DKFZ block:
dkfz_block <- packLegend(
  lgd_dkfz2,
  lgd_dkfz3,
  direction = "vertical",
  gap = unit(4, "mm")    # <<< space BETWEEN DKFZ Ben-2 and DKFZ Ben-3
)

#Methylation/molecular block: Cluster → MenG → DKFZ 
meth_block <- packLegend(
  lgd_clust,
  lgd_meng,
  dkfz_block,
  direction = "vertical",
  gap = unit(1.8, "mm")
)

#Mutation block (single legend)
mut_block <- packLegend(
  lgd_mut,
  direction = "vertical",
  gap = unit(1.8, "mm")
)

#CNV block
cnv_block <- packLegend(
  lgd_cnv,
  direction = "vertical",
  gap = unit(1.8, "mm")
)

#RIGHT legend column:
#methylation block → gray bar → mutation → gray bar → CNV
#Histopathology (TOP) → gray bar → methylation block → gray bar → mutation → gray bar → CNV
legend_right <- packLegend(
  lgd_hist,
  lgd_bar(w_mm = 30, h_mm = 0.9, col = "grey90"),
  
  meth_block,
  lgd_bar(w_mm = 30, h_mm = 0.9, col = "grey90"),
  
  mut_block,
  lgd_bar(w_mm = 30, h_mm = 0.9, col = "grey90"),
  
  cnv_block,
  direction = "vertical",
  gap = unit(2, "mm")
)

#Two-column legend (left clinical | right histology+molecular+mutation+CNV)
legend_2col <- packLegend(
  legend_left,
  legend_right,
  direction = "horizontal",
  gap = unit(6, "mm")
)

###-----DRAW IN RSTUDIO-----###
gap_mm_preview <- -60

grid::grid.newpage()

lg_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(legend_2col))
lg_w_mm <- grid::convertWidth(grid::grobWidth(lg_grob), "mm", valueOnly = TRUE)
if (is.na(lg_w_mm) || lg_w_mm < 1) lg_w_mm <- 120

lay <- grid::grid.layout(
  nrow = 1, ncol = 3,
  widths = grid::unit.c(
    grid::unit(1, "npc") - grid::unit(lg_w_mm + gap_mm_preview, "mm"),
    grid::unit(gap_mm_preview, "mm"),
    grid::unit(lg_w_mm, "mm")
  )
)

grid::pushViewport(grid::viewport(layout = lay))

grid::pushViewport(grid::viewport(layout.pos.col = 1))
ComplexHeatmap::draw(
  ht_all,
  newpage = FALSE,
  merge_legends = FALSE,
  show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE,
  padding = grid::unit(c(4, 0, 4, 2), "mm")
)
grid::popViewport()

grid::pushViewport(grid::viewport(layout.pos.col = 3))
grid::pushViewport(grid::viewport(x = unit(0, "npc"), y = unit(0.5, "npc"),
                                  just = c("left", "center")))
grid::grid.draw(lg_grob)
grid::popViewport(2)

grid::popViewport()


###-----EXPORT TO PDF-----###
out_pdf <- file.path(base_dir, "Combined Analysis_WHO I_Final/Oncoprint_Master.pdf")

#Vector-safe spacing controls
gap_mm_pdf      <- 0.8            # actual gap column (keep >= 0)
legend_shift_mm <- 18             # pull legend LEFT inside its column (try 10–30)
margins_mm      <- c(6, 6, 6, 6)  # top, right, bottom, left (mm)

#Must match how YOU draw the heatmap
left_padding_mm <- 2              # matches draw(... padding = unit(c(t,r,b,l),"mm")) LEFT
label_pad_mm    <- 0.8            # cushion between row labels and first tile (mm)

#helper: measure grob sizes in mm on a PDF device
.measure_mm <- function(expr) {
  g <- grid::grid.grabExpr(expr)  # only to measure; do not reuse grob for drawing
  w <- grid::convertWidth(grid::grobWidth(g),  "mm", valueOnly = TRUE)
  h <- grid::convertHeight(grid::grobHeight(g), "mm", valueOnly = TRUE)
  list(w = w, h = h)
}

#Force ComplexHeatmap to stay vector
options(ComplexHeatmap.use_raster = FALSE)

#Row labels exactly as shown (exclude separator rows)
row_labels <- c(
  "Patient ID",
  "Age","Sex","WHO Grade",
  "SB Tumor Sequenced","SB Tumor Present","Multiple Meningiomas","Regression After Stopping DMPA",
  "Histology",
  "Methylation Cluster","MenG Group","DKFZ Ben-2","DKFZ Ben-3",
  "TRAF7","PIK3CA","FGFR1","NF2",
  "1p","1q","2p","2q","7p","7q","13q","18p","18q","22q"
)

#Measure on a temporary PDF device
tmp <- tempfile(fileext = ".pdf")
pdf(tmp, width = 10, height = 10, useDingbats = FALSE)
grid::grid.newpage()

#Legend size
m_lg <- .measure_mm(ComplexHeatmap::draw(legend_2col))

#Heatmap BODY width (tiles only)
hm_body_mm <- length(sample_ids) * tile_w_mm

#Max row-label width using same font settings as row names
row_gp <- grid::gpar(fontsize = 10, fontface = "plain")  # match your row_names_gp
grid::pushViewport(grid::viewport(gp = row_gp))
label_w_mm <- grid::convertWidth(max(grid::stringWidth(row_labels)), "mm", valueOnly = TRUE)
grid::popViewport()

#Heatmap height (for page height)
m_ht <- .measure_mm(
  ComplexHeatmap::draw(
    ht_all,
    newpage = FALSE,
    merge_legends = FALSE,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE,
    padding = grid::unit(c(4, 0, 4, left_padding_mm), "mm")
  )
)

dev.off()

#Safety fallbacks
if (!is.finite(m_lg$w) || m_lg$w < 1) m_lg$w <- 120
if (!is.finite(m_lg$h) || m_lg$h < 1) m_lg$h <- 240
if (!is.finite(m_ht$h) || m_ht$h < 1) m_ht$h <- 240
if (!is.finite(label_w_mm) || label_w_mm < 1) label_w_mm <- 60

#Heatmap block width MUST include full label width (do NOT subtract “tighten” here)
hm_block_mm <- left_padding_mm + label_w_mm + label_pad_mm + hm_body_mm

#Final device size (mm -> inches)
need_w_mm <- margins_mm[4] + hm_block_mm + gap_mm_pdf + m_lg$w + margins_mm[2]
need_h_mm <- margins_mm[1] + max(m_ht$h, m_lg$h) + margins_mm[3]

pdf(out_pdf, width = need_w_mm/25.4, height = need_h_mm/25.4, useDingbats = FALSE)
grid::grid.newpage()

#Outer margin viewport (top-left anchored)
grid::pushViewport(grid::viewport(
  x = grid::unit(0, "npc"), y = grid::unit(1, "npc"),
  just = c("left", "top")
))

#Inner content viewport (inside margins, fixed mm size)
grid::pushViewport(grid::viewport(
  x = grid::unit(margins_mm[4], "mm"),
  y = grid::unit(need_h_mm - margins_mm[1], "mm"),
  just = c("left", "top"),
  width  = grid::unit(need_w_mm - margins_mm[4] - margins_mm[2], "mm"),
  height = grid::unit(need_h_mm - margins_mm[1] - margins_mm[3], "mm")
))

#3-column layout: heatmap block | gap | legend (all fixed mm, all >= 0)
lay <- grid::grid.layout(
  nrow = 1, ncol = 3,
  widths = grid::unit.c(
    grid::unit(hm_block_mm, "mm"),
    grid::unit(gap_mm_pdf, "mm"),
    grid::unit(m_lg$w, "mm")
  )
)
grid::pushViewport(grid::viewport(layout = lay))

#Column 1: Heatmap (vector)
grid::pushViewport(grid::viewport(layout.pos.col = 1))
ComplexHeatmap::draw(
  ht_all,
  newpage = FALSE,
  merge_legends = FALSE,
  show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE,
  padding = grid::unit(c(4, 0, 4, left_padding_mm), "mm")
)
grid::popViewport()

#Column 3: Legends (vector) — shift LEFT inside its own column to remove whitespace
grid::pushViewport(grid::viewport(layout.pos.col = 3))
lg_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(legend_2col))

#Clamp shift so we never over-shift past the column
legend_shift_mm <- max(0, min(legend_shift_mm, m_lg$w - 5))

grid::pushViewport(grid::viewport(
  x = grid::unit(0, "npc") - grid::unit(legend_shift_mm, "mm"),
  y = grid::unit(0.5, "npc"),
  just = c("left", "center")
))
grid::grid.draw(lg_grob)
grid::popViewport(2)

#Cleanup
grid::popViewport(3)
dev.off()

message("Saved: ", normalizePath(out_pdf))



##==============================================================================
###-----22. PGR & 11q22.1 METHYLATION — RESTRICTED BAYLOR REFERENCE COHORTS-----###
##==============================================================================
#
# Reviewer 3 requested PGR / 11q22.1 methylation analysis restricted to 
# NF2-wildtype and TRAF7-mutant reference subgroups. This section generates:
#
#   Panel A: PGR butterfly        — DMPA vs. Baylor WHO I NF2-WT   (n~69)
#   Panel B: PGR butterfly        — DMPA vs. Baylor WHO I TRAF7-mut (n~20)
#   Panel C: 11q22.1 butterfly    — DMPA vs. Baylor WHO I NF2-WT   (n~69)
#   Panel D: 11q22.1 butterfly    — DMPA vs. Baylor WHO I TRAF7-mut (n~20)
#
# Per-CpG statistics (PGR) and per-gene statistics (11q22.1) saved to Excel.
# Butterfly plots saved to PNG and PDF.
#
# Assumes beta_combat and pd_df are in memory from Section 17 (or loaded below).
##==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(openxlsx)
})


###-----Verify/load beta_combat and pd_df-----###
#If running standalone after Section 17, these should already be in memory.
#Uncomment the block below if they are not:
#
# beta_path <- file.path(base_dir, "Combined Analysis_WHO I_Final/beta_combat_Combined Analysis_All Grades.csv")
# pd_path   <- file.path(base_dir, "Combined Analysis_WHO I_Final/pd_Combined Analysis_All Grades.csv")
# beta_combat <- load_beta_table(beta_path) |> as.matrix()
# storage.mode(beta_combat) <- "double"
# pd_df <- read.csv(pd_path, row.names = 1, check.names = FALSE)

stopifnot(exists("beta_combat"), exists("pd_df"))
beta_combat <- as.matrix(beta_combat); storage.mode(beta_combat) <- "double"


###-----Output directory-----###
out_root <- file.path(ResultsDir, "PGR and 11q22 DMR", "Restricted_Reference_Cohorts")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

###-----Load Baylor clinical Excel and extract mutation status-----###
meta_path  <- file.path(base_dir, "Baylor/GSE189521_Clinical_data (Bayley et al).xlsx")
meta_sheet <- "2. Clinical and genomic dataSH"
meta_bay   <- readxl::read_excel(meta_path, sheet = meta_sheet, skip = 2)

#Extract matching key from Idat column
#pd_df rownames for Baylor = full IDAT basename (e.g. "GSM5702884_203293440013_R01C01")
meta_bay <- meta_bay %>%
  mutate(
    pd_match     = `Idat file`,
    NF2_status   = ifelse(is.na(NF2) | trimws(NF2) == "", "WT", "Mutant"),
    TRAF7_status = ifelse(is.na(TRAF7) | trimws(TRAF7) == "", "WT", "Mutant"),
    WHO_norm     = trimws(`WHO grade`)
  )

#Match to pd_df rownames for Baylor samples
baylor_rows <- rownames(pd_df)[pd_df$Cohort == "Baylor"]
matched_n   <- sum(meta_bay$pd_match %in% baylor_rows)
message(sprintf("Baylor Excel → pd_df match: %d of %d rows", matched_n, nrow(meta_bay)))

#Filter to WHO I and build subgroup ID vectors
meta_whoi    <- meta_bay %>% filter(WHO_norm == "WHO I")
nf2wt_gsm    <- meta_whoi %>% filter(NF2_status == "WT")      %>% pull(pd_match)
traf7mut_gsm <- meta_whoi %>% filter(TRAF7_status == "Mutant") %>% pull(pd_match)

#Intersect with samples actually present in pd_df
nf2wt_ids    <- intersect(nf2wt_gsm, baylor_rows)
traf7mut_ids <- intersect(traf7mut_gsm, baylor_rows)
dmpa_ids     <- rownames(pd_df)[pd_df$Cohort == "DMPA"]

message(sprintf(
  "Sample counts — DMPA: n=%d | Baylor WHO I NF2-WT: n=%d | Baylor WHO I TRAF7-mut: n=%d",
  length(dmpa_ids), length(nf2wt_ids), length(traf7mut_ids)
))

#Sanity checks
stopifnot(length(dmpa_ids) > 0, length(nf2wt_ids) > 0, length(traf7mut_ids) > 0)
stopifnot(all(dmpa_ids %in% colnames(beta_combat)))
stopifnot(all(nf2wt_ids %in% colnames(beta_combat)))
stopifnot(all(traf7mut_ids %in% colnames(beta_combat)))


###-----Load 450k annotation-----###
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_df  <- as.data.frame(ann450k) %>%
  transmute(
    CpG   = Name,
    chr   = as.character(chr),
    pos   = as.numeric(pos),
    Gene  = as.character(UCSC_RefGene_Name),
    Group = as.character(UCSC_RefGene_Group)
  )


###-----Identify PGR CpGs-----###
has_PGR  <- function(x) grepl("(^|;)PGR(;|$)", ifelse(is.na(x), "", x), perl = TRUE)
pgr_cpgs <- ann_df$CpG[has_PGR(ann_df$Gene)]
pgr_cpgs <- intersect(pgr_cpgs, rownames(beta_combat))
stopifnot(length(pgr_cpgs) > 0)
message(sprintf("PGR CpGs in beta_combat: %d", length(pgr_cpgs)))

pgr_ann <- ann_df %>%
  filter(CpG %in% pgr_cpgs) %>%
  arrange(chr, pos)


###-----Identify 11q22.1 probes and build gene-level mapping-----###
band_rng <- c(start = 100000000, end = 102000000)   # hg19 window around PGR
cpg_band <- ann_df %>%
  filter(chr == "chr11", pos >= band_rng["start"], pos <= band_rng["end"]) %>%
  pull(CpG) %>%
  intersect(rownames(beta_combat))
stopifnot(length(cpg_band) > 0)

probe_gene_tbl <- ann_df %>%
  filter(CpG %in% cpg_band, !is.na(Gene), nzchar(Gene)) %>%
  select(CpG, Gene, pos) %>%
  separate_rows(Gene, sep = ";") %>%
  mutate(Gene = trimws(Gene)) %>%
  filter(nzchar(Gene))

gene_pos <- probe_gene_tbl %>%
  group_by(Gene) %>%
  summarise(gene_pos = median(pos, na.rm = TRUE), .groups = "drop")

idx_by_gene <- split(probe_gene_tbl$CpG, probe_gene_tbl$Gene)


##==============================================================================
###-----REUSABLE ANALYSIS + PLOTTING FUNCTIONS-----###
##==============================================================================

###-----Per-CpG PGR differential methylation statistics-----###
run_pgr_stats <- function(beta, dmpa_ids, ref_ids) {
  
  test_one <- function(cg) {
    x  <- beta[cg, dmpa_ids]
    y  <- beta[cg, ref_ids]
    tt <- t.test(x, y, var.equal = FALSE)
    ww <- wilcox.test(x, y, exact = FALSE)
    tibble(
      CpG        = cg,
      mean_DMPA  = mean(x, na.rm = TRUE),
      mean_Ref   = mean(y, na.rm = TRUE),
      delta_beta = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
      t_stat     = unname(tt$statistic),
      pval_t     = tt$p.value,
      W_stat     = unname(ww$statistic),
      pval_wil   = ww$p.value
    )
  }
  
  res <- bind_rows(lapply(pgr_cpgs, test_one)) %>%
    mutate(FDR_t   = p.adjust(pval_t,   method = "BH"),
           FDR_wil = p.adjust(pval_wil, method = "BH")) %>%
    left_join(select(pgr_ann, CpG, chr, pos, Gene, Group), by = "CpG") %>%
    arrange(FDR_t, desc(abs(delta_beta)))
  
  #Gene-level summary (mean β across all PGR CpGs per sample)
  all_ids  <- c(dmpa_ids, ref_ids)
  pgr_mean <- colMeans(beta[pgr_cpgs, all_ids, drop = FALSE], na.rm = TRUE)
  grp_vec  <- factor(ifelse(names(pgr_mean) %in% dmpa_ids, "DMPA", "Reference"))
  
  gene_t <- t.test(pgr_mean[grp_vec == "DMPA"],
                   pgr_mean[grp_vec == "Reference"], var.equal = FALSE)
  gene_w <- wilcox.test(pgr_mean[grp_vec == "DMPA"],
                        pgr_mean[grp_vec == "Reference"], exact = FALSE)
  
  gene_summary <- tibble(
    group     = c("Reference", "DMPA"),
    mean_beta = c(mean(pgr_mean[grp_vec == "Reference"]),
                  mean(pgr_mean[grp_vec == "DMPA"])),
    sd_beta   = c(sd(pgr_mean[grp_vec == "Reference"]),
                  sd(pgr_mean[grp_vec == "DMPA"])),
    n         = c(sum(grp_vec == "Reference"),
                  sum(grp_vec == "DMPA"))
  )
  gene_tests <- tibble(
    metric     = c("Welch t-test", "Wilcoxon rank-sum"),
    statistic  = c(unname(gene_t$statistic), unname(gene_w$statistic)),
    p_value    = c(gene_t$p.value,           gene_w$p.value),
    delta_beta = mean(pgr_mean[grp_vec == "DMPA"]) -
      mean(pgr_mean[grp_vec == "Reference"])
  )
  
  list(per_cpg = res, gene_summary = gene_summary, gene_tests = gene_tests)
}


###-----Per-gene 11q22.1 differential methylation statistics-----###
run_11q22_stats <- function(beta, dmpa_ids, ref_ids) {
  
  all_ids <- c(dmpa_ids, ref_ids)
  mat_sub <- beta[unique(probe_gene_tbl$CpG), all_ids, drop = FALSE]
  
  gene_mat <- vapply(idx_by_gene, function(cpgs) {
    colMeans(mat_sub[intersect(cpgs, rownames(mat_sub)), , drop = FALSE], na.rm = TRUE)
  }, FUN.VALUE = numeric(length(all_ids)))
  gene_mat <- t(gene_mat)
  colnames(gene_mat) <- all_ids
  
  grp_vec <- factor(ifelse(all_ids %in% dmpa_ids, "DMPA", "Reference"))
  
  res <- bind_rows(lapply(rownames(gene_mat), function(g) {
    x  <- gene_mat[g, grp_vec == "DMPA"]
    y  <- gene_mat[g, grp_vec == "Reference"]
    tt <- t.test(x, y, var.equal = FALSE)
    ww <- wilcox.test(x, y, exact = FALSE)
    tibble(
      Gene       = g,
      mean_DMPA  = mean(x, na.rm = TRUE),
      mean_Ref   = mean(y, na.rm = TRUE),
      delta_beta = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
      t_stat     = unname(tt$statistic),
      pval_t     = tt$p.value,
      W_stat     = unname(ww$statistic),
      pval_wil   = ww$p.value
    )
  })) %>%
    mutate(FDR_t   = p.adjust(pval_t,   method = "BH"),
           FDR_wil = p.adjust(pval_wil, method = "BH")) %>%
    left_join(gene_pos, by = "Gene") %>%
    arrange(gene_pos)
  
  list(per_gene = res)
}


###-----Butterfly plot — PGR (CpG-level)-----###
butterfly_pgr <- function(beta, dmpa_ids, ref_ids, ref_label) {
  
  all_ids <- c(dmpa_ids, ref_ids)
  grp     <- factor(ifelse(all_ids %in% dmpa_ids, "DMPA", "Reference"),
                    levels = c("DMPA", "Reference"))
  names(grp) <- all_ids
  
  #CpG order: genomic position (pgr_ann already sorted by chr, pos)
  cpg_order <- pgr_ann$CpG
  
  per_cpg_long <- beta[pgr_cpgs, all_ids, drop = FALSE] %>%
    as.data.frame(check.names = FALSE) %>%
    rownames_to_column("CpG") %>%
    pivot_longer(-CpG, names_to = "Sample", values_to = "Beta") %>%
    mutate(Group = grp[Sample]) %>%
    filter(!is.na(Group))
  per_cpg_long$CpG <- factor(per_cpg_long$CpG, levels = cpg_order)
  
  summ <- per_cpg_long %>%
    group_by(CpG, Group) %>%
    summarise(mean_beta = mean(Beta, na.rm = TRUE), .groups = "drop")
  
  y_key   <- tibble(CpG = cpg_order, y = seq_along(cpg_order))
  summ    <- left_join(summ, y_key, by = "CpG")
  leftDF  <- filter(summ, Group == "DMPA")
  rightDF <- filter(summ, Group == "Reference")
  
  COL_DMPA <- "#582C83"
  COL_REF  <- "#555555"
  
  ggplot() +
    geom_segment(aes(x = 0, xend = 0, y = 0.5, yend = nrow(y_key) + 0.5),
                 color = "grey65", linewidth = 0.6) +
    geom_segment(data = y_key,
                 aes(x = -0.02, xend = 0.02, y = y, yend = y),
                 color = "grey70", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed",
               linewidth = 0.6, color = "grey45") +
    #DMPA wing (left)
    geom_segment(data = leftDF,
                 aes(x = 0, xend = -mean_beta, y = y, yend = y),
                 color = COL_DMPA, linewidth = 2.2, lineend = "round") +
    geom_point(data = leftDF, aes(x = -mean_beta, y = y),
               color = COL_DMPA, size = 2.2) +
    #Reference wing (right)
    geom_segment(data = rightDF,
                 aes(x = 0, xend = mean_beta, y = y, yend = y),
                 color = COL_REF, linewidth = 2.2, lineend = "round") +
    geom_point(data = rightDF, aes(x = mean_beta, y = y),
               color = COL_REF, size = 2.2) +
    #Top cohort labels
    annotate("text", x = -0.65, y = nrow(y_key) + 1.3, label = "DMPA",
             color = COL_DMPA, fontface = "bold", size = 5.2, hjust = 0.5) +
    annotate("text", x =  0.65, y = nrow(y_key) + 1.3, label = ref_label,
             color = COL_REF, fontface = "bold", size = 5.2, hjust = 0.5) +
    scale_y_continuous(breaks = y_key$y, labels = y_key$CpG,
                       expand = expansion(mult = c(0.02, 0.10))) +
    scale_x_continuous(limits = c(-1, 1),
                       breaks = seq(-1, 1, by = 0.25),
                       labels = c("1.00","0.75","0.50","0.25","0",
                                  "0.25","0.50","0.75","1.00")) +
    labs(
      title = "Progesterone receptor methylation",
      x     = expression("Mean "*beta*" magnitude"),
      y     = "CpG sites along PGR gene"
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title        = element_text(hjust = 0.5, face = "bold"),
      axis.line.x       = element_line(linewidth = 1.1, color = "black"),
      axis.line.y       = element_line(linewidth = 1.1, color = "black"),
      axis.ticks        = element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length = unit(5, "pt"),
      axis.text.x       = element_text(size = 13),
      axis.text.y       = element_text(size = 12),
      axis.title.x      = element_text(size = 15, margin = ggplot2::margin(t = 10)),
      axis.title.y      = element_text(size = 15, margin = ggplot2::margin(r = 16)),
      panel.grid        = element_blank()
    )
}


###-----Butterfly plot — 11q22.1 cytoband (gene-level)-----###
butterfly_11q22 <- function(beta, dmpa_ids, ref_ids, ref_label) {
  
  all_ids <- c(dmpa_ids, ref_ids)
  grp     <- factor(ifelse(all_ids %in% dmpa_ids, "DMPA", "Reference"),
                    levels = c("DMPA", "Reference"))
  names(grp) <- all_ids
  
  #Build gene-level mean β matrix
  mat_sub  <- beta[unique(probe_gene_tbl$CpG), all_ids, drop = FALSE]
  gene_mat <- vapply(idx_by_gene, function(cpgs) {
    colMeans(mat_sub[intersect(cpgs, rownames(mat_sub)), , drop = FALSE], na.rm = TRUE)
  }, FUN.VALUE = numeric(length(all_ids)))
  gene_mat <- t(gene_mat)
  colnames(gene_mat) <- all_ids
  
  #Gene order: genomic position
  gene_order <- gene_pos %>% arrange(gene_pos) %>% pull(Gene)
  
  gene_long <- gene_mat %>%
    as.data.frame(check.names = FALSE) %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Beta") %>%
    mutate(Group = grp[Sample]) %>%
    filter(!is.na(Group))
  gene_long$Gene <- factor(gene_long$Gene, levels = gene_order)
  
  summ <- gene_long %>%
    group_by(Gene, Group) %>%
    summarise(mean_beta = mean(Beta, na.rm = TRUE), .groups = "drop")
  
  y_key   <- tibble(Gene = gene_order, y = seq_along(gene_order))
  summ    <- left_join(summ, y_key, by = "Gene")
  leftDF  <- filter(summ, Group == "DMPA")
  rightDF <- filter(summ, Group == "Reference")
  
  COL_DMPA <- "#582C83"
  COL_REF  <- "#555555"
  
  ggplot() +
    geom_segment(aes(x = 0, xend = 0, y = 0.5, yend = nrow(y_key) + 0.5),
                 color = "grey65", linewidth = 0.6) +
    geom_segment(data = y_key,
                 aes(x = -0.02, xend = 0.02, y = y, yend = y),
                 color = "grey70", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed",
               linewidth = 0.6, color = "grey45") +
    #DMPA wing (left)
    geom_segment(data = leftDF,
                 aes(x = 0, xend = -mean_beta, y = y, yend = y),
                 color = COL_DMPA, linewidth = 2.2, lineend = "round") +
    geom_point(data = leftDF, aes(x = -mean_beta, y = y),
               color = COL_DMPA, size = 2.2) +
    #Reference wing (right)
    geom_segment(data = rightDF,
                 aes(x = 0, xend = mean_beta, y = y, yend = y),
                 color = COL_REF, linewidth = 2.2, lineend = "round") +
    geom_point(data = rightDF, aes(x = mean_beta, y = y),
               color = COL_REF, size = 2.2) +
    #Top cohort labels
    annotate("text", x = -0.65, y = nrow(y_key) + 1.3, label = "DMPA",
             color = COL_DMPA, fontface = "bold", size = 5.2, hjust = 0.5) +
    annotate("text", x =  0.65, y = nrow(y_key) + 1.3, label = ref_label,
             color = COL_REF, fontface = "bold", size = 5.2, hjust = 0.5) +
    scale_y_continuous(breaks = y_key$y, labels = y_key$Gene,
                       expand = expansion(mult = c(0.02, 0.10))) +
    scale_x_continuous(limits = c(-1, 1),
                       breaks = seq(-1, 1, by = 0.25),
                       labels = c("1.00","0.75","0.50","0.25","0",
                                  "0.25","0.50","0.75","1.00")) +
    labs(
      title = "11q22.1 cytoband methylation",
      x     = expression("Mean "*beta*" magnitude"),
      y     = "Genes within 11q22.1 cytoband"
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title        = element_text(hjust = 0.5, face = "bold"),
      axis.line.x       = element_line(linewidth = 1.1, color = "black"),
      axis.line.y       = element_line(linewidth = 1.1, color = "black"),
      axis.ticks        = element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length = unit(5, "pt"),
      axis.text.x       = element_text(size = 13),
      axis.text.y       = element_text(size = 12),
      axis.title.x      = element_text(size = 15, margin = ggplot2::margin(t = 10)),
      axis.title.y      = element_text(size = 15, margin = ggplot2::margin(r = 16)),
      panel.grid        = element_blank()
    )
}


##==============================================================================
###-----RUN ANALYSES FOR EACH RESTRICTED REFERENCE COHORT-----###
##==============================================================================

#Define reference subgroups
ref_subgroups <- list(
  list(label = "Baylor NF2-WT",   tag = "NF2_WT",    ids = nf2wt_ids),
  list(label = "Baylor TRAF7-mut", tag = "TRAF7_mut", ids = traf7mut_ids)
)

#Master stats workbook (all comparisons combined into one Excel file)
wb_stats <- createWorkbook()

for (sg in ref_subgroups) {
  
  ref_label <- sg$label
  ref_tag   <- sg$tag
  ref_ids   <- sg$ids
  n_ref     <- length(ref_ids)
  
  message(sprintf("\n=== %s (n=%d) vs DMPA (n=%d) ===", ref_label, n_ref, length(dmpa_ids)))
  
  
  ###-----PGR analysis-----###
  pgr_res <- run_pgr_stats(beta_combat, dmpa_ids, ref_ids)
  p_pgr   <- butterfly_pgr(beta_combat, dmpa_ids, ref_ids, ref_label)
  
  #Add PGR per-CpG stats sheet
  sheet_cpg <- paste0("PGR_CpG_", ref_tag)
  addWorksheet(wb_stats, sheet_cpg)
  writeData(wb_stats, sheet_cpg, pgr_res$per_cpg)
  freezePane(wb_stats, sheet_cpg, firstRow = TRUE)
  
  #Add PGR gene-level summary sheet
  sheet_gene <- paste0("PGR_Gene_", ref_tag)
  addWorksheet(wb_stats, sheet_gene)
  writeData(wb_stats, sheet_gene, pgr_res$gene_summary, startRow = 1)
  writeData(wb_stats, sheet_gene, pgr_res$gene_tests,
            startRow = nrow(pgr_res$gene_summary) + 3)
  
  #Save PGR butterfly plots
  ggsave(file.path(out_root, sprintf("Fig_PGR_butterfly_%s.png", ref_tag)),
         p_pgr, width = 9, height = 7.5, dpi = 300)
  ggsave(file.path(out_root, sprintf("Fig_PGR_butterfly_%s.pdf", ref_tag)),
         p_pgr, width = 9, height = 7.5)
  
  print(p_pgr)
  message(sprintf("  PGR: saved plots + stats for %s", ref_label))
  
  
  ###-----11q22.1 analysis-----###
  cyto_res <- run_11q22_stats(beta_combat, dmpa_ids, ref_ids)
  p_cyto   <- butterfly_11q22(beta_combat, dmpa_ids, ref_ids, ref_label)
  
  #Add 11q22.1 per-gene stats sheet
  sheet_11q <- paste0("11q22_Gene_", ref_tag)
  addWorksheet(wb_stats, sheet_11q)
  writeData(wb_stats, sheet_11q, cyto_res$per_gene)
  freezePane(wb_stats, sheet_11q, firstRow = TRUE)
  
  #Save 11q22.1 butterfly plots
  ggsave(file.path(out_root, sprintf("Fig_11q22_butterfly_%s.png", ref_tag)),
         p_cyto, width = 9, height = 7.5, dpi = 300)
  ggsave(file.path(out_root, sprintf("Fig_11q22_butterfly_%s.pdf", ref_tag)),
         p_cyto, width = 9, height = 7.5)
  
  print(p_cyto)
  message(sprintf("  11q22.1: saved plots + stats for %s", ref_label))
}


###-----Save combined stats workbook-----###
for (sh in names(wb_stats)) setColWidths(wb_stats, sh, cols = 1:50, widths = "auto")
saveWorkbook(wb_stats,
             file.path(out_root, "Restricted_Reference_DiffMeth_Stats.xlsx"),
             overwrite = TRUE)

message("\nAll restricted-reference analyses complete.")
message("Output directory: ", normalizePath(out_root))
message(sprintf("  4 butterfly plots (PNG + PDF) and 1 stats workbook with %d sheets saved.",
                length(names(wb_stats))))

##==============================================================================
###-----23. PROGESTERONE SIGNALING PATHWAY DIFFERENTIAL METHYLATION-----###
##==============================================================================
# Approach:
#   1. Retrieve curated gene set from GO:0032570 ("response to progesterone")
#   2. Map pathway genes to 450k CpG probes present in beta_combat
#   3. Compute per-gene mean beta values across mapped probes
#   4. Per-gene differential methylation (Welch t + Wilcoxon, FDR-corrected)
#   5. Composite pathway methylation score (median-centered) with box plot
#   6. Per-gene butterfly plot matching Figure 3G/H styling
#
# Reference cohort: full Baylor + Heidelberg WHO Grade I (matching Fig 3G/H)
#
# Assumes beta_combat and pd_df are in memory from Section 17.
##==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(openxlsx)
})

BiocManager::install("org.Hs.eg.db")

###-----Verify/load beta_combat and pd_df-----###
#If running standalone, uncomment the load block from Section 17 header.
stopifnot(exists("beta_combat"), exists("pd_df"))
beta_combat <- as.matrix(beta_combat); storage.mode(beta_combat) <- "double"


###-----Output directory-----###
out_root <- file.path(getwd(), "Progesterone_Pathway_DiffMeth")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)


###-----Build group factor (same as Section 17)-----###
grp <- factor(ifelse(pd_df$Cohort == "DMPA", "DMPA", "Reference"),
              levels = c("DMPA", "Reference"))
names(grp) <- rownames(pd_df)
stopifnot(identical(colnames(beta_combat), names(grp)))

dmpa_ids <- names(grp)[grp == "DMPA"]
ref_ids  <- names(grp)[grp == "Reference"]
message(sprintf("Cohort sizes — DMPA: n=%d | Reference: n=%d", length(dmpa_ids), length(ref_ids)))


###-----Retrieve GO:0032570 ("response to progesterone") gene list-----###
#GOALL includes genes annotated to GO:0032570 and all child terms
go_result <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = "GO:0032570",
  keytype = "GOALL",
  columns = c("SYMBOL", "GENENAME")
)
pathway_genes <- sort(unique(go_result$SYMBOL[!is.na(go_result$SYMBOL)]))
message(sprintf("GO:0032570 'response to progesterone': %d unique gene symbols retrieved",
                length(pathway_genes)))


###-----Map pathway genes to 450k CpG probes-----###
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_df  <- as.data.frame(ann450k) %>%
  transmute(
    CpG  = Name,
    chr  = as.character(chr),
    pos  = as.numeric(pos),
    Gene = as.character(UCSC_RefGene_Name)
  )

#Expand semicolon-separated gene names and intersect with pathway genes
probe_gene <- ann_df %>%
  filter(!is.na(Gene), nzchar(Gene)) %>%
  separate_rows(Gene, sep = ";") %>%
  mutate(Gene = trimws(Gene)) %>%
  filter(Gene %in% pathway_genes) %>%
  filter(CpG %in% rownames(beta_combat)) %>%
  distinct(CpG, Gene)

mapped_genes <- sort(unique(probe_gene$Gene))
unmapped     <- setdiff(pathway_genes, mapped_genes)

message(sprintf("Mapped to 450k probes in beta_combat: %d of %d pathway genes (%d CpG probes total)",
                length(mapped_genes), length(pathway_genes), n_distinct(probe_gene$CpG)))
if (length(unmapped) > 0) {
  message(sprintf("Unmapped genes (no 450k probes): %s",
                  paste(head(unmapped, 15), collapse = ", ")))
}

stopifnot(length(mapped_genes) >= 5)  # need at least a handful of genes
idx_by_gene <- split(probe_gene$CpG, probe_gene$Gene)


###-----Compute per-gene mean beta per sample-----###
gene_beta <- vapply(idx_by_gene, function(cpgs) {
  colMeans(beta_combat[cpgs, , drop = FALSE], na.rm = TRUE)
}, FUN.VALUE = numeric(ncol(beta_combat)))
gene_beta <- t(gene_beta)   # genes × samples
colnames(gene_beta) <- colnames(beta_combat)


###-----Per-gene differential methylation-----###
per_gene_stats <- bind_rows(lapply(rownames(gene_beta), function(g) {
  x  <- gene_beta[g, dmpa_ids]
  y  <- gene_beta[g, ref_ids]
  tt <- t.test(x, y, var.equal = FALSE)
  ww <- wilcox.test(x, y, exact = FALSE)
  tibble(
    Gene       = g,
    n_CpGs     = length(idx_by_gene[[g]]),
    mean_DMPA  = mean(x, na.rm = TRUE),
    mean_Ref   = mean(y, na.rm = TRUE),
    delta_beta = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
    t_stat     = unname(tt$statistic),
    pval_t     = tt$p.value,
    W_stat     = unname(ww$statistic),
    pval_wil   = ww$p.value
  )
})) %>%
  mutate(FDR_t   = p.adjust(pval_t,   method = "BH"),
         FDR_wil = p.adjust(pval_wil, method = "BH")) %>%
  arrange(FDR_t, desc(abs(delta_beta)))

#Report top hits
message("\nTop 10 genes by FDR-corrected p-value:")
print(head(select(per_gene_stats, Gene, n_CpGs, delta_beta, pval_t, FDR_t), 10))

sig_genes <- per_gene_stats %>% filter(FDR_t < 0.05)
message(sprintf("\nGenes significant at FDR < 0.05: %d", nrow(sig_genes)))


###-----Composite pathway methylation score (median-centered)-----###
#For each gene, subtract the median across all samples
#Then average centered values per sample → one pathway score per sample
gene_medians  <- apply(gene_beta, 1, median, na.rm = TRUE)
gene_centered <- sweep(gene_beta, 1, gene_medians, "-")

pathway_score <- colMeans(gene_centered, na.rm = TRUE)
pathway_df <- tibble(
  Sample = names(pathway_score),
  Score  = pathway_score,
  Group  = grp[names(pathway_score)]
)

#Pathway-level statistical tests
pw_t <- t.test(Score ~ Group, data = pathway_df, var.equal = FALSE)
pw_w <- wilcox.test(Score ~ Group, data = pathway_df, exact = FALSE)

message(sprintf(
  "\nPathway score summary:\n  DMPA mean = %.5f (SD %.5f)\n  Ref  mean = %.5f (SD %.5f)\n  Delta = %.5f\n  Welch t p = %.4f | Wilcoxon p = %.4f",
  mean(pathway_df$Score[pathway_df$Group == "DMPA"]),
  sd(pathway_df$Score[pathway_df$Group == "DMPA"]),
  mean(pathway_df$Score[pathway_df$Group == "Reference"]),
  sd(pathway_df$Score[pathway_df$Group == "Reference"]),
  mean(pathway_df$Score[pathway_df$Group == "DMPA"]) - mean(pathway_df$Score[pathway_df$Group == "Reference"]),
  pw_t$p.value, pw_w$p.value
))


##==============================================================================
###-----FIGURES-----###
##==============================================================================

COL_DMPA <- "#582C83"
COL_REF  <- "#555555"


###-----Panel A: Pathway score box plot-----###


p_box <- ggplot(pathway_df, aes(x = Group, y = Score, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.8) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA, color = "black", linewidth = 0.6) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.8, aes(color = Group)) +
  scale_fill_manual(values  = c(DMPA = COL_DMPA, Reference = COL_REF)) +
  scale_color_manual(values = c(DMPA = COL_DMPA, Reference = COL_REF)) +
  
  labs(
    title = "Composite progesterone response\npathway methylation score",
    x     = NULL,
    y     = expression("Relative pathway methylation ("*Delta*beta*" from median)"),
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    axis.line.x       = element_line(linewidth = 1.1, color = "black"),
    axis.line.y       = element_line(linewidth = 1.1, color = "black"),
    axis.ticks        = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(5, "pt"),
    axis.text.x       = element_text(size = 14, face = "bold"),
    axis.text.y       = element_text(size = 13),
    axis.title.y      = element_text(size = 15, margin = ggplot2::margin(r = 10)),
    legend.position   = "none",
    panel.grid        = element_blank()
  )

print(p_box)

ggsave(file.path(out_root, "Fig_Progesterone_Pathway_Score_Boxplot.png"),
       p_box, width = 5, height = 7.5, dpi = 300)
ggsave(file.path(out_root, "Fig_Progesterone_Pathway_Score_Boxplot.pdf"),
       p_box, width = 9, height = 7.5)
message("Saved pathway score box plot.")


###-----Panel B: Per-gene butterfly plot-----###

#Sort genes by chromosome and genomic position
gene_positions <- probe_gene %>%
  left_join(select(ann_df, CpG, chr, pos), by = "CpG") %>%
  group_by(Gene) %>%
  summarise(
    chr_mode = names(sort(table(chr), decreasing = TRUE))[1],
    med_pos  = median(pos, na.rm = TRUE),
    .groups  = "drop"
  )

chr_levels <- paste0("chr", c(1:22, "X", "Y"))
gene_order <- gene_positions %>%
  mutate(chr_num = factor(chr_mode, levels = chr_levels)) %>%
  arrange(chr_num, med_pos) %>%
  pull(Gene)

#Build per-group means
summ_long <- gene_beta %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Beta") %>%
  mutate(Group = grp[Sample]) %>%
  filter(!is.na(Group)) %>%
  group_by(Gene, Group) %>%
  summarise(mean_beta = mean(Beta, na.rm = TRUE), .groups = "drop")

summ_long$Gene <- factor(summ_long$Gene, levels = gene_order)

y_key   <- tibble(Gene = gene_order, y = seq_along(gene_order))
summ_long <- left_join(summ_long, y_key, by = "Gene")
leftDF  <- filter(summ_long, Group == "DMPA")
rightDF <- filter(summ_long, Group == "Reference")

#Adaptive sizing based on gene count
n_genes    <- length(gene_order)
y_text_sz  <- 7.5
bar_lw     <- 1.5
pt_sz      <- 1.8
fig_height <- max(7.5, n_genes * 0.18 + 2)
label_gap  <- 1.5

HILITE_LABEL <- "DMPA"
REF_LABEL    <- "Baylor/Heidelberg"

p_butterfly <- ggplot() +
  #Center spine and y ticks
  geom_segment(aes(x = 0, xend = 0, y = 0.5, yend = nrow(y_key) + 0.5),
               color = "grey65", linewidth = 0.6) +
  geom_segment(data = y_key,
               aes(x = -0.02, xend = 0.02, y = y, yend = y),
               color = "grey70", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed",
             linewidth = 0.6, color = "grey45") +
  #DMPA wing (left)
  geom_segment(data = leftDF,
               aes(x = 0, xend = -mean_beta, y = y, yend = y),
               color = COL_DMPA, linewidth = bar_lw, lineend = "round") +
  geom_point(data = leftDF, aes(x = -mean_beta, y = y),
             color = COL_DMPA, size = pt_sz) +
  #Reference wing (right)
  geom_segment(data = rightDF,
               aes(x = 0, xend = mean_beta, y = y, yend = y),
               color = COL_REF, linewidth = bar_lw, lineend = "round") +
  geom_point(data = rightDF, aes(x = mean_beta, y = y),
             color = COL_REF, size = pt_sz) +
  #Top cohort labels
  annotate("text", x = -0.65, y = nrow(y_key) + label_gap, label = HILITE_LABEL,
           color = COL_DMPA, fontface = "bold", size = 5.2, hjust = 0.5) +
  annotate("text", x =  0.65, y = nrow(y_key) + label_gap, label = REF_LABEL,
           color = COL_REF, fontface = "bold", size = 5.2, hjust = 0.5) +
  scale_y_continuous(breaks = y_key$y, labels = y_key$Gene,
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_x_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, by = 0.25),
                     labels = c("1.00","0.75","0.50","0.25","0",
                                "0.25","0.50","0.75","1.00")) +
  labs(
    title = "Per-gene progesterone response\npathway methylation",
    x     = expression("Mean "*beta*" magnitude"),
    y     = "Genes within progesterone response pathway"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    axis.line.x       = element_line(linewidth = 1.1, color = "black"),
    axis.line.y       = element_line(linewidth = 1.1, color = "black"),
    axis.ticks        = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(5, "pt"),
    axis.text.x       = element_text(size = 13),
    axis.text.y       = element_text(size = y_text_sz),
    axis.title.x      = element_text(size = 15, margin = ggplot2::margin(t = 10)),
    axis.title.y      = element_text(size = 15, margin = ggplot2::margin(r = 16)),
    panel.grid        = element_blank()
  )

print(p_butterfly)

ggsave(file.path(out_root, "Fig_Progesterone_Pathway_Butterfly.png"),
       p_butterfly, width = 5, height = fig_height, dpi = 300)
ggsave(file.path(out_root, "Fig_Progesterone_Pathway_Butterfly.pdf"),
       p_butterfly, width = 9, height = fig_height)
message(sprintf("Saved butterfly plot (%d genes, height=%.1f in).", n_genes, fig_height))

p_box_small <- p_box + theme(text = element_text(size = 12))

ggsave(file.path(out_root, "Fig_Supp5A_Response_Boxplot.pdf"),
       p_box_small, width = 4, height = 5)

##==============================================================================
###-----SAVE STATISTICS TO EXCEL-----###
##==============================================================================

wb <- createWorkbook()

#Sheet 1: Per-gene stats
addWorksheet(wb, "Per_Gene_Stats")
writeData(wb, "Per_Gene_Stats", per_gene_stats)
freezePane(wb, "Per_Gene_Stats", firstRow = TRUE)

#Sheet 2: Pathway score summary
pathway_summary <- tibble(
  group      = c("DMPA", "Reference"),
  mean_score = c(mean(pathway_df$Score[pathway_df$Group == "DMPA"]),
                 mean(pathway_df$Score[pathway_df$Group == "Reference"])),
  sd_score   = c(sd(pathway_df$Score[pathway_df$Group == "DMPA"]),
                 sd(pathway_df$Score[pathway_df$Group == "Reference"])),
  n          = c(sum(pathway_df$Group == "DMPA"),
                 sum(pathway_df$Group == "Reference"))
)
pathway_tests <- tibble(
  metric    = c("Welch t-test", "Wilcoxon rank-sum"),
  statistic = c(unname(pw_t$statistic), unname(pw_w$statistic)),
  p_value   = c(pw_t$p.value, pw_w$p.value)
)

addWorksheet(wb, "Pathway_Score_Summary")
writeData(wb, "Pathway_Score_Summary", pathway_summary, startRow = 1)
writeData(wb, "Pathway_Score_Summary", pathway_tests,
          startRow = nrow(pathway_summary) + 3)

#Sheet 3: Gene list documentation (for Methods transparency)
gene_doc <- tibble(
  Gene    = mapped_genes,
  GO_Term = "GO:0032570",
  GO_Name = "response to progesterone",
  n_CpGs  = sapply(idx_by_gene[mapped_genes], length)
)
addWorksheet(wb, "Pathway_Gene_List")
writeData(wb, "Pathway_Gene_List", gene_doc)
freezePane(wb, "Pathway_Gene_List", firstRow = TRUE)

#Sheet 4: Per-sample pathway scores (for reproducibility)
addWorksheet(wb, "Per_Sample_Scores")
writeData(wb, "Per_Sample_Scores", pathway_df)

for (sh in names(wb)) setColWidths(wb, sh, cols = 1:50, widths = "auto")
saveWorkbook(wb, file.path(out_root, "Progesterone_Pathway_DiffMeth_Stats.xlsx"),
             overwrite = TRUE)

message("\nSection 23 complete — progesterone signaling pathway analysis.")
message("Output directory: ", normalizePath(out_root))
message(sprintf("  %d pathway genes analyzed | Box plot + butterfly plot saved (PNG + PDF)",
                n_genes))
message(sprintf("  Stats workbook: %d sheets", length(names(wb))))


###-----GO:0050847 — PROGESTERONE RECEPTOR SIGNALING PATHWAY-----###
# Separate analysis of the PGR signaling machinery (GO:0050847), which includes
# PGR itself, coactivators, corepressors, and receptor turnover regulators.
# Assumes all Section 23 objects are still in memory.

###-----Retrieve GO:0050847 gene list-----###
go_sig_result <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = "GO:0050847",
  keytype = "GOALL",
  columns = c("SYMBOL", "GENENAME")
)
sig_pathway_genes <- sort(unique(go_sig_result$SYMBOL[!is.na(go_sig_result$SYMBOL)]))
message(sprintf("GO:0050847 'progesterone receptor signaling pathway': %d unique gene symbols",
                length(sig_pathway_genes)))


###-----Map to 450k CpG probes-----###
probe_gene_sig <- ann_df %>%
  filter(!is.na(Gene), nzchar(Gene)) %>%
  separate_rows(Gene, sep = ";") %>%
  mutate(Gene = trimws(Gene)) %>%
  filter(Gene %in% sig_pathway_genes) %>%
  filter(CpG %in% rownames(beta_combat)) %>%
  distinct(CpG, Gene)

mapped_genes_sig <- sort(unique(probe_gene_sig$Gene))
unmapped_sig     <- setdiff(sig_pathway_genes, mapped_genes_sig)

message(sprintf("Mapped to 450k probes: %d of %d genes (%d CpG probes)",
                length(mapped_genes_sig), length(sig_pathway_genes), n_distinct(probe_gene_sig$CpG)))
if (length(unmapped_sig) > 0) {
  message(sprintf("Unmapped: %s", paste(unmapped_sig, collapse = ", ")))
}

stopifnot(length(mapped_genes_sig) >= 3)
idx_by_gene_sig <- split(probe_gene_sig$CpG, probe_gene_sig$Gene)


###-----Compute per-gene mean beta per sample-----###
gene_beta_sig <- vapply(idx_by_gene_sig, function(cpgs) {
  colMeans(beta_combat[cpgs, , drop = FALSE], na.rm = TRUE)
}, FUN.VALUE = numeric(ncol(beta_combat)))
gene_beta_sig <- t(gene_beta_sig)
colnames(gene_beta_sig) <- colnames(beta_combat)


###-----Per-gene differential methylation-----###
per_gene_stats_sig <- bind_rows(lapply(rownames(gene_beta_sig), function(g) {
  x  <- gene_beta_sig[g, dmpa_ids]
  y  <- gene_beta_sig[g, ref_ids]
  tt <- t.test(x, y, var.equal = FALSE)
  ww <- wilcox.test(x, y, exact = FALSE)
  tibble(
    Gene       = g,
    n_CpGs     = length(idx_by_gene_sig[[g]]),
    mean_DMPA  = mean(x, na.rm = TRUE),
    mean_Ref   = mean(y, na.rm = TRUE),
    delta_beta = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
    t_stat     = unname(tt$statistic),
    pval_t     = tt$p.value,
    W_stat     = unname(ww$statistic),
    pval_wil   = ww$p.value
  )
})) %>%
  mutate(FDR_t   = p.adjust(pval_t,   method = "BH"),
         FDR_wil = p.adjust(pval_wil, method = "BH")) %>%
  arrange(FDR_t, desc(abs(delta_beta)))

message("\nPer-gene results (GO:0050847):")
print(select(per_gene_stats_sig, Gene, n_CpGs, delta_beta, pval_t, FDR_t))

sig_genes_sig <- per_gene_stats_sig %>% filter(FDR_t < 0.05)
message(sprintf("Genes significant at FDR < 0.05: %d", nrow(sig_genes_sig)))


###-----Composite pathway score (median-centered)-----###
gene_medians_sig  <- apply(gene_beta_sig, 1, median, na.rm = TRUE)
gene_centered_sig <- sweep(gene_beta_sig, 1, gene_medians_sig, "-")

pathway_score_sig <- colMeans(gene_centered_sig, na.rm = TRUE)
pathway_df_sig <- tibble(
  Sample = names(pathway_score_sig),
  Score  = pathway_score_sig,
  Group  = grp[names(pathway_score_sig)]
)

pw_t_sig <- t.test(Score ~ Group, data = pathway_df_sig, var.equal = FALSE)
pw_w_sig <- wilcox.test(Score ~ Group, data = pathway_df_sig, exact = FALSE)

message(sprintf(
  "\nPathway score (GO:0050847):\n  DMPA mean = %.5f (SD %.5f)\n  Ref  mean = %.5f (SD %.5f)\n  Delta = %.5f\n  Welch t p = %.4f | Wilcoxon p = %.4f",
  mean(pathway_df_sig$Score[pathway_df_sig$Group == "DMPA"]),
  sd(pathway_df_sig$Score[pathway_df_sig$Group == "DMPA"]),
  mean(pathway_df_sig$Score[pathway_df_sig$Group == "Reference"]),
  sd(pathway_df_sig$Score[pathway_df_sig$Group == "Reference"]),
  mean(pathway_df_sig$Score[pathway_df_sig$Group == "DMPA"]) - mean(pathway_df_sig$Score[pathway_df_sig$Group == "Reference"]),
  pw_t_sig$p.value, pw_w_sig$p.value
))


###-----Panel C: Pathway score box plot (GO:0050847)-----###
p_box_sig <- ggplot(pathway_df_sig, aes(x = Group, y = Score, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.8) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA, color = "black", linewidth = 0.6) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.8, aes(color = Group)) +
  scale_fill_manual(values  = c(DMPA = COL_DMPA, Reference = COL_REF)) +
  scale_color_manual(values = c(DMPA = COL_DMPA, Reference = COL_REF)) +
  labs(
    title = "Composite PGR signaling\npathway methylation score",
    x     = NULL,
    y     = expression("Relative pathway methylation ("*Delta*beta*" from median)")
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    axis.line.x       = element_line(linewidth = 1.1, color = "black"),
    axis.line.y       = element_line(linewidth = 1.1, color = "black"),
    axis.ticks        = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(5, "pt"),
    axis.text.x       = element_text(size = 14, face = "bold"),
    axis.text.y       = element_text(size = 13),
    axis.title.y      = element_text(size = 15, margin = ggplot2::margin(r = 10)),
    legend.position   = "none",
    panel.grid        = element_blank()
  )

print(p_box_sig)

ggsave(file.path(out_root, "Fig_PGR_Signaling_Pathway_Score_Boxplot.png"),
       p_box_sig, width = 5, height = 7.5, dpi = 300)
ggsave(file.path(out_root, "Fig_PGR_Signaling_Pathway_Score_Boxplot.pdf"),
       p_box_sig, width = 9, height = 7.5)
message("Saved PGR signaling pathway box plot.")


###-----Panel D: Per-gene butterfly plot (GO:0050847)-----###

#Sort genes by genomic position
gene_positions_sig <- probe_gene_sig %>%
  left_join(select(ann_df, CpG, chr, pos), by = "CpG") %>%
  group_by(Gene) %>%
  summarise(
    chr_mode = names(sort(table(chr), decreasing = TRUE))[1],
    med_pos  = median(pos, na.rm = TRUE),
    .groups  = "drop"
  )

chr_levels <- paste0("chr", c(1:22, "X", "Y"))
gene_order_sig <- gene_positions_sig %>%
  mutate(chr_num = factor(chr_mode, levels = chr_levels)) %>%
  arrange(chr_num, med_pos) %>%
  pull(Gene)

#Build per-group means
summ_long_sig <- gene_beta_sig %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Beta") %>%
  mutate(Group = grp[Sample]) %>%
  filter(!is.na(Group)) %>%
  group_by(Gene, Group) %>%
  summarise(mean_beta = mean(Beta, na.rm = TRUE), .groups = "drop")

summ_long_sig$Gene <- factor(summ_long_sig$Gene, levels = gene_order_sig)

y_key_sig   <- tibble(Gene = gene_order_sig, y = seq_along(gene_order_sig))
summ_long_sig <- left_join(summ_long_sig, y_key_sig, by = "Gene")
leftDF_sig  <- filter(summ_long_sig, Group == "DMPA")
rightDF_sig <- filter(summ_long_sig, Group == "Reference")

#Standard sizing — 12 genes is comparable to the main figure panels
HILITE_LABEL <- "DMPA"
REF_LABEL    <- "Baylor/Heidelberg"

p_butterfly_sig <- ggplot() +
  geom_segment(aes(x = 0, xend = 0, y = 0.5, yend = nrow(y_key_sig) + 0.5),
               color = "grey65", linewidth = 0.6) +
  geom_segment(data = y_key_sig,
               aes(x = -0.02, xend = 0.02, y = y, yend = y),
               color = "grey70", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed",
             linewidth = 0.6, color = "grey45") +
  #DMPA wing (left)
  geom_segment(data = leftDF_sig,
               aes(x = 0, xend = -mean_beta, y = y, yend = y),
               color = COL_DMPA, linewidth = 2.2, lineend = "round") +
  geom_point(data = leftDF_sig, aes(x = -mean_beta, y = y),
             color = COL_DMPA, size = 2.2) +
  #Reference wing (right)
  geom_segment(data = rightDF_sig,
               aes(x = 0, xend = mean_beta, y = y, yend = y),
               color = COL_REF, linewidth = 2.2, lineend = "round") +
  geom_point(data = rightDF_sig, aes(x = mean_beta, y = y),
             color = COL_REF, size = 2.2) +
  #Top cohort labels
  annotate("text", x = -0.65, y = nrow(y_key_sig) + 1.3, label = HILITE_LABEL,
           color = COL_DMPA, fontface = "bold", size = 5.2, hjust = 0.5) +
  annotate("text", x =  0.65, y = nrow(y_key_sig) + 1.3, label = REF_LABEL,
           color = COL_REF, fontface = "bold", size = 5.2, hjust = 0.5) +
  scale_y_continuous(breaks = y_key_sig$y, labels = y_key_sig$Gene,
                     expand = expansion(mult = c(0.02, 0.10))) +
  scale_x_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, by = 0.25),
                     labels = c("1.00","0.75","0.50","0.25","0",
                                "0.25","0.50","0.75","1.00")) +
  labs(
    title = "Per-gene PGR signaling\npathway methylation",
    x     = expression("Mean "*beta*" magnitude"),
    y     = "Genes within PGR signaling pathway"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    axis.line.x       = element_line(linewidth = 1.1, color = "black"),
    axis.line.y       = element_line(linewidth = 1.1, color = "black"),
    axis.ticks        = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(5, "pt"),
    axis.text.x       = element_text(size = 13),
    axis.text.y       = element_text(size = 12),
    axis.title.x      = element_text(size = 15, margin = ggplot2::margin(t = 10)),
    axis.title.y      = element_text(size = 15, margin = ggplot2::margin(r = 16)),
    panel.grid        = element_blank()
  )

print(p_butterfly_sig)

ggsave(file.path(out_root, "Fig_PGR_Signaling_Pathway_Butterfly.png"),
       p_butterfly_sig, width = 9, height = 7.5, dpi = 300)
ggsave(file.path(out_root, "Fig_PGR_Signaling_Pathway_Butterfly.pdf"),
       p_butterfly_sig, width = 9, height = 7.5)
message("Saved PGR signaling pathway butterfly plot.")


###-----Save GO:0050847 stats to Excel-----###
wb_sig <- createWorkbook()

addWorksheet(wb_sig, "Per_Gene_Stats")
writeData(wb_sig, "Per_Gene_Stats", per_gene_stats_sig)
freezePane(wb_sig, "Per_Gene_Stats", firstRow = TRUE)

pathway_summary_sig <- tibble(
  group      = c("DMPA", "Reference"),
  mean_score = c(mean(pathway_df_sig$Score[pathway_df_sig$Group == "DMPA"]),
                 mean(pathway_df_sig$Score[pathway_df_sig$Group == "Reference"])),
  sd_score   = c(sd(pathway_df_sig$Score[pathway_df_sig$Group == "DMPA"]),
                 sd(pathway_df_sig$Score[pathway_df_sig$Group == "Reference"])),
  n          = c(sum(pathway_df_sig$Group == "DMPA"),
                 sum(pathway_df_sig$Group == "Reference"))
)
pathway_tests_sig <- tibble(
  metric    = c("Welch t-test", "Wilcoxon rank-sum"),
  statistic = c(unname(pw_t_sig$statistic), unname(pw_w_sig$statistic)),
  p_value   = c(pw_t_sig$p.value, pw_w_sig$p.value)
)

addWorksheet(wb_sig, "Pathway_Score_Summary")
writeData(wb_sig, "Pathway_Score_Summary", pathway_summary_sig, startRow = 1)
writeData(wb_sig, "Pathway_Score_Summary", pathway_tests_sig,
          startRow = nrow(pathway_summary_sig) + 3)

gene_doc_sig <- tibble(
  Gene    = mapped_genes_sig,
  GO_Term = "GO:0050847",
  GO_Name = "progesterone receptor signaling pathway",
  n_CpGs  = sapply(idx_by_gene_sig[mapped_genes_sig], length)
)
addWorksheet(wb_sig, "Pathway_Gene_List")
writeData(wb_sig, "Pathway_Gene_List", gene_doc_sig)
freezePane(wb_sig, "Pathway_Gene_List", firstRow = TRUE)

addWorksheet(wb_sig, "Per_Sample_Scores")
writeData(wb_sig, "Per_Sample_Scores", pathway_df_sig)

for (sh in names(wb_sig)) setColWidths(wb_sig, sh, cols = 1:50, widths = "auto")
saveWorkbook(wb_sig, file.path(out_root, "PGR_Signaling_Pathway_DiffMeth_Stats.xlsx"),
             overwrite = TRUE)

message("\nGO:0050847 analysis complete.")
message(sprintf("  %d genes mapped | Box plot + butterfly saved | Stats workbook saved",
                length(mapped_genes_sig)))



#PATCHWORK
library(patchwork)

pad_left  <- theme(plot.margin = margin(5, 55, 5, 5))
pad_right <- theme(plot.margin = margin(5, 5, 5, 55))

###-----Supplementary Figure 6 (PGR signaling — simple side-by-side)-----###
combined_sig <- (p_box_sig + pad_left) + (p_butterfly_sig + pad_right) + 
  plot_layout(widths = c(1, 2.5))

ggsave(file.path(out_root, "Fig_Supp5_PGR_Signaling_Combined.pdf"),
       combined_sig, width = 15, height = 7.5)

###-----Supplementary Figure 7 (response pathway — taller butterfly)-----###
fig_height_resp <- fig_height * 1.5

left_col_6 <- (p_box + pad_left) / 
  plot_spacer() + 
  plot_layout(heights = c(7.5, fig_height_resp - 7.5))

combined_resp <- left_col_6 | (p_butterfly + pad_right)
combined_resp <- combined_resp + plot_layout(widths = c(1, 2.5))

ggsave(file.path(out_root, "Fig_Supp6_Progesterone_Response_Combined.pdf"),
       combined_resp, width = 15, height = fig_height_resp)






##==============================================================================
###-----SESSION INFO (reproducibility)-----###
##==============================================================================

sessionInfo()
