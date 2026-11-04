# Combine the results of the parallel run of steps 02_01a/b (workflows/02_01_tas_parallel.sh)
# into the tables the serial run of the scripts writes, under the same names:
# {temp_var}_outcome_temp.csv and {temp_var}_outcome_siteyear_temp.csv (02_01a),
# {temp_var}_outcome_temp_water_gpp.csv and {temp_var}_outcome_siteyear_temp_water_gpp.csv
# (02_01b). As in the serial run, the outcome table has row i for site i of data/site_info.csv
# (all NA for a site that was skipped), up to the last site with a result, and the site-year
# table has the sites' rows in site_info order.
#
# A task without results (its script stopped with an error, Slurm stopped it, or it has not run)
# stops the combining of its step and temperature variable, with the reason from the task's
# log.txt. With TAS_ALLOW_FAILED=1, the tasks with results are combined anyway, the failed
# sites have no outcome row (all NA, like a skipped site); the failed tasks and their
# reasons are listed, and written to failed_tasks.csv in the tasks' output directory.
#
# Usage (from the repository root; normally through `bash workflows/02_01_tas_parallel.sh combine`):
#   Rscript workflows/02_01_tas_parallel_combine.R <tasks.tsv> <tasks' output directory> <output directory>
stopifnot(requireNamespace("readr"))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
tasks_file <- args[1]
tasks_dir <- args[2]
combine_dir <- args[3]

tasks <- readr::read_tsv(tasks_file, col_names = c("step", "temp_var", "site"),
                         col_types = readr::cols(.default = readr::col_character()))
site_info <- readr::read_csv(file.path("data", "site_info.csv"), col_types = readr::cols(.default = readr::col_character()))
suffix <- c(a = "_temp", b = "_temp_water_gpp")
allow_failed <- Sys.getenv("TAS_ALLOW_FAILED") == "1"

# a task's table; NULL where the script skipped the site (the site-year table is then written
# from an empty data.frame(), a file with no header)
read_task_table <- function(path, col_types) {
  if (length(readLines(path, n = 2)) < 2) return(NULL)
  as.data.frame(readr::read_csv(path, col_types = col_types, na = "NA"))
}

# why a task has no results: the error that stopped its script (the last error in its log.txt;
# errors caught by try() in the scripts come before it), or what its log says otherwise
task_failure <- function(task_dir) {
  log <- file.path(task_dir, "log.txt")
  if (!file.exists(log)) return("not run (no log.txt)")
  lines <- readLines(log, warn = FALSE)
  i <- utils::tail(grep("^Error", lines), 1)
  if (length(i) == 0) {
    return("no error in log.txt: still running, or stopped by Slurm (time or memory limit; see its Slurm log)")
  }
  # the error message runs up to the call stack or the end of the run
  end <- i
  while (end < length(lines) && !grepl("^(Calls:|In addition:|Execution halted)", lines[end + 1])) end <- end + 1
  paste(trimws(lines[i:end]), collapse = " ")
}

unlink(file.path(tasks_dir, "failed_tasks.csv"))
failed <- data.frame(step = character(), temp_var = character(), site = character(), reason = character())
for (group in split(tasks, paste(tasks$step, tasks$temp_var))) {
  step <- group$step[1]
  temp_var <- group$temp_var[1]
  outcome_name <- paste0(temp_var, "_outcome", suffix[[step]], ".csv")
  siteyear_name <- paste0(temp_var, "_outcome_siteyear", suffix[[step]], ".csv")
  task_dirs <- file.path(tasks_dir, paste0(temp_var, "_", step), group$site)

  done <- file.exists(file.path(task_dirs, outcome_name)) & file.exists(file.path(task_dirs, siteyear_name))
  if (!all(done)) {
    failed <- rbind(failed, data.frame(step = step, temp_var = temp_var, site = group$site[!done],
                                       reason = vapply(task_dirs[!done], task_failure, character(1), USE.NAMES = FALSE)))
    if (!allow_failed) next
    group <- group[done, ]
    task_dirs <- task_dirs[done]
    if (nrow(group) == 0) next
  }

  outcome <- do.call(rbind, lapply(file.path(task_dirs, outcome_name), read_task_table,
                                   col_types = readr::cols(site_ID = readr::col_character(), .default = readr::col_double())))
  # no site with results: the serial run's empty outcome table
  if (is.null(outcome)) {
    outcome <- data.frame(site_ID = character(), RMSE = double(), R2 = double(), control_year = double(),
                          window_size = double(), nwindow = double(), TAS = double(), TASp = double())
  }
  # a one-site run leaves an all-NA row for every site before it in site_info
  outcome <- outcome[!is.na(outcome$site_ID), , drop = FALSE]
  stopifnot(!anyDuplicated(outcome$site_ID))
  idx <- match(site_info$site_ID, outcome$site_ID)
  n_rows <- if (any(!is.na(idx))) max(which(!is.na(idx))) else 0
  outcome <- outcome[idx[seq_len(n_rows)], , drop = FALSE]
  rownames(outcome) <- NULL

  siteyear <- do.call(rbind, lapply(file.path(task_dirs, siteyear_name), read_task_table,
                                    col_types = readr::cols(site_ID = readr::col_character(), window = readr::col_character(),
                                                            .default = readr::col_double())))
  if (is.null(siteyear)) siteyear <- data.frame()
  if (nrow(siteyear) > 0) siteyear <- siteyear[order(match(siteyear$site_ID, site_info$site_ID)), , drop = FALSE]

  utils::write.csv(outcome, file.path(combine_dir, outcome_name), row.names = FALSE)
  utils::write.csv(siteyear, file.path(combine_dir, siteyear_name), row.names = FALSE)
  message(sprintf("02_01%s %s: %d sites with TAS of %d tasks -> %s, %s", step, temp_var, sum(!is.na(outcome$site_ID)),
                  nrow(group), file.path(combine_dir, outcome_name), file.path(combine_dir, siteyear_name)))
}

if (nrow(failed) > 0) {
  failed_file <- file.path(tasks_dir, "failed_tasks.csv")
  utils::write.csv(failed, failed_file, row.names = FALSE)
  listing <- paste0("02_01", failed$step, " ", failed$temp_var, " ", failed$site, ": ", failed$reason, collapse = "\n")
  if (allow_failed) {
    message(nrow(failed), " tasks have no results; the other tasks were combined (TAS_ALLOW_FAILED=1). ",
            "The failed tasks (also in ", failed_file, "):\n", listing)
  } else {
    stop(nrow(failed), " tasks have no results, so their step and temperature variable were not combined ",
         "(TAS_ALLOW_FAILED=1 combines the others; the list is also in ", failed_file, "):\n", listing, call. = FALSE)
  }
}
