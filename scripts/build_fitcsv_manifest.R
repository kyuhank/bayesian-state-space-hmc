source("utils/pipeline_utils.R")

fit_root <- env_get_chr("FITCSV_DIR", "fitcsv")
out_file <- env_get_chr("FITCSV_MANIFEST_FILE", file.path("results", "fitcsv_manifest.rds"))

invisible(ensure_dir(dirname(out_file)))

manifest <- collect_fitcsv_manifest(root = fit_root)
saveRDS(manifest, out_file)

cat("Saved fitcsv manifest:", out_file, "\n")
cat("Rows:", nrow(manifest), "\n")
