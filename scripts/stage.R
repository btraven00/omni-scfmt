# Shared stage() helper for R modules.
#
# Usage:
#   source(file.path(Sys.getenv("SCFMT_ROOT", "."), "scripts", "stage.R"))
#   stage("load", { ad <- read_h5ad(path) })
#   stage("access", { colsum <- colSums(ad$X) })
#
# Between stages we force a full GC so the RSS baseline for the next stage
# isn't contaminated by the previous stage's residuals. This matters when
# inner denet is attributing per-stage memory peaks. Wall time is reported
# excluding the gc() call itself so the stage timing reflects work only.

if (!exists(".scfmt_stage_timings")) .scfmt_stage_timings <- list()

stage <- function(name, block) {
  env <- parent.frame()
  emit(name, "start")
  t0 <- proc.time()[["elapsed"]]
  eval(substitute(block), envir = env)
  elapsed <- proc.time()[["elapsed"]] - t0
  emit(name, "end", attrs = list(elapsed_s = elapsed))
  .scfmt_stage_timings[[name]] <<- elapsed
  # Full GC *after* emitting end so the stage boundary in the timeline
  # lines up with the actual work boundary, not the cleanup.
  invisible(gc(full = TRUE, verbose = FALSE))
}

stage_timings_df <- function(dataset = NA_character_) {
  data.frame(
    stage     = names(.scfmt_stage_timings),
    elapsed_s = unlist(.scfmt_stage_timings),
    dataset   = dataset,
    row.names = NULL
  )
}
