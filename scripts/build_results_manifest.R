source("utils/pipeline_utils.R")

results_dir <- env_get_chr("RESULTS_DIR", "results")
out_file <- env_get_chr("RESULTS_MANIFEST_FILE", file.path(results_dir, "results_manifest.rds"))

invisible(ensure_dir(dirname(out_file)))

manifest <- build_results_manifest(results_dir = results_dir)
saveRDS(manifest, out_file)

cat("Saved results manifest:", out_file, "\n")
cat("Entries:", length(manifest), "\n")
