env_get_chr <- function(name, default = "") {
  x <- Sys.getenv(name, "")
  if (!nzchar(x)) default else x
}

env_get_int <- function(name, default) {
  x <- Sys.getenv(name, "")
  if (!nzchar(x)) return(as.integer(default))
  as.integer(x)
}

env_get_num <- function(name, default) {
  x <- Sys.getenv(name, "")
  if (!nzchar(x)) return(as.numeric(default))
  as.numeric(x)
}

env_get_vec <- function(name, default = character(), sep = ",") {
  x <- Sys.getenv(name, "")
  if (!nzchar(x)) return(default)
  trimws(strsplit(x, sep, fixed = TRUE)[[1]])
}

read_env_kv_file <- function(path) {
  if (!file.exists(path)) return(character())
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines <- lines[!startsWith(lines, "#")]
  if (length(lines) == 0) return(character())

  out <- character()
  for (line in lines) {
    pos <- regexpr("=", line, fixed = TRUE)[1]
    if (pos <= 1) next
    key <- trimws(substr(line, 1, pos - 1))
    val <- trimws(substr(line, pos + 1, nchar(line)))
    if ((startsWith(val, "\"") && endsWith(val, "\"")) ||
        (startsWith(val, "'") && endsWith(val, "'"))) {
      val <- substr(val, 2, nchar(val) - 1)
    }
    if (nzchar(key)) out[[key]] <- val
  }
  out
}

load_env_file_if_present <- function(path, overwrite = FALSE) {
  kv <- read_env_kv_file(path)
  if (length(kv) == 0) return(invisible(FALSE))
  for (nm in names(kv)) {
    existing <- Sys.getenv(nm, "")
    if (overwrite || !nzchar(existing)) {
      do.call(Sys.setenv, stats::setNames(list(kv[[nm]]), nm))
    }
  }
  invisible(TRUE)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

clean_name_for_path <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[/\\\\]", "_", x)
  x
}

resolve_sbc_model_path <- function(sbc_type, sbc_model) {
  candidates <- c(
    file.path("stan", sbc_type, sbc_model),
    file.path("models", sbc_model)
  )

  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop(
      "Could not find Stan model file. Tried: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[[1]]
}

derive_tau2_prior_name <- function(tau2_prior) {
  if (length(tau2_prior) < 3) return("custom")
  if (tau2_prior[3] == 0) return("no")
  if (isTRUE(all.equal(tau2_prior[2], 0.1272))) return("0.1")
  if (isTRUE(all.equal(tau2_prior[2], 1.1013))) return("0.3")
  if (isTRUE(all.equal(tau2_prior[2], 2.8516))) return("0.5")
  "custom"
}

make_sbc_result_stem <- function(SBCmodel, adaptDelta, tau2Prior_name, n_sims) {
  paste0(
    "sbc_",
    clean_name_for_path(SBCmodel), "_",
    clean_name_for_path(as.character(adaptDelta)), "_",
    clean_name_for_path(tau2Prior_name), "_",
    clean_name_for_path(as.character(n_sims))
  )
}

list_sbc_result_files <- function(results_dir = "results") {
  if (!dir.exists(results_dir)) return(character())
  files <- list.files(results_dir, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  keep <- grepl("^sbc_.*\\.(rds|RData)$", basename(files))
  files[keep]
}

read_sbc_result_file <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext == "rds") {
    obj <- readRDS(path)
    if (is.list(obj) && !is.null(obj$results) && !is.null(obj$metadata)) {
      meta <- obj$metadata
      return(list(
        file = normalizePath(path, winslash = "/", mustWork = FALSE),
        name = sub("\\.rds$", "", basename(path), ignore.case = TRUE),
        results = obj$results,
        metadata = list(
          SBCtype = if (!is.null(meta$SBCtype)) meta$SBCtype else NA_character_,
          SBCmodel = if (!is.null(meta$SBCmodel)) meta$SBCmodel else NA_character_,
          nsample = if (!is.null(meta$nsample)) meta$nsample else NA_integer_,
          nchains = if (!is.null(meta$nchains)) meta$nchains else NA_integer_,
          nwarmup = if (!is.null(meta$nwarmup)) meta$nwarmup else NA_integer_,
          adaptDelta = if (!is.null(meta$adaptDelta)) meta$adaptDelta else NA_real_,
          tau2Prior_name = if (!is.null(meta$tau2Prior_name)) meta$tau2Prior_name else NA_character_
        )
      ))
    }
    stop("Unsupported SBC result .rds format: ", path)
  }

  if (ext != "rdata") {
    stop("Unsupported SBC result file extension: ", path)
  }

  e <- new.env(parent = emptyenv())
  load(path, envir = e)

  list(
    file = normalizePath(path, winslash = "/", mustWork = FALSE),
    name = sub("\\.RData$", "", basename(path)),
    results = if (exists("results", envir = e, inherits = FALSE)) get("results", envir = e) else NULL,
    metadata = list(
      SBCtype = if (exists("SBCtype", envir = e, inherits = FALSE)) get("SBCtype", envir = e) else NA_character_,
      SBCmodel = if (exists("SBCmodel", envir = e, inherits = FALSE)) get("SBCmodel", envir = e) else NA_character_,
      nsample = if (exists("nsample", envir = e, inherits = FALSE)) get("nsample", envir = e) else NA_integer_,
      nchains = if (exists("nchains", envir = e, inherits = FALSE)) get("nchains", envir = e) else NA_integer_,
      nwarmup = if (exists("nwarmup", envir = e, inherits = FALSE)) get("nwarmup", envir = e) else NA_integer_,
      adaptDelta = if (exists("adaptDelta", envir = e, inherits = FALSE)) get("adaptDelta", envir = e) else NA_real_,
      tau2Prior_name = if (exists("tau2Prior_name", envir = e, inherits = FALSE)) get("tau2Prior_name", envir = e) else NA_character_
    )
  )
}

build_results_manifest <- function(results_dir = "results") {
  files <- list_sbc_result_files(results_dir)
  entries <- lapply(files, read_sbc_result_file)
  names(entries) <- vapply(entries, `[[`, "", "name")
  entries
}

collect_fitcsv_manifest <- function(root = "fitcsv") {
  if (!dir.exists(root)) return(data.frame())

  files <- list.files(root, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(data.frame())

  rel <- sub(paste0("^", normalizePath(root, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(files, winslash = "/", mustWork = FALSE))

  data.frame(
    file = normalizePath(files, winslash = "/", mustWork = FALSE),
    rel_path = rel,
    size_bytes = file.info(files)$size,
    mtime = file.info(files)$mtime,
    stringsAsFactors = FALSE
  )
}
