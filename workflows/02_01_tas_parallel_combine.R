# Combine the results of the parallel run of steps 02_01a/b (workflows/02_01_tas_parallel.sh)
# into the tables the serial run of the scripts writes, under the same names:
# {temp_var}_outcome_temp.csv and {temp_var}_outcome_siteyear_temp.csv (02_01a),
# {temp_var}_outcome_temp_water_gpp.csv and {temp_var}_outcome_siteyear_temp_water_gpp.csv
# (02_01b). As in the serial run, the outcome table has row i for site i of data/site_info.csv
# (all NA for a site that was skipped), up to the last site with a result, and the site-year
# table has the sites' rows in site_info order.
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

# a task's table; NULL where the script skipped the site (the site-year table is then written
# from an empty data.frame(), a file with no header)
read_task_table <- function(path, col_types) {
  if (length(readLines(path, n = 2)) < 2) return(NULL)
  as.data.frame(readr::read_csv(path, col_types = col_types, na = "NA"))
}

missing <- character()
for (group in split(tasks, paste(tasks$step, tasks$temp_var))) {
  step <- group$step[1]
  temp_var <- group$temp_var[1]
  outcome_name <- paste0(temp_var, "_outcome", suffix[[step]], ".csv")
  siteyear_name <- paste0(temp_var, "_outcome_siteyear", suffix[[step]], ".csv")
  task_dirs <- file.path(tasks_dir, paste0(temp_var, "_", step), group$site)

  done <- file.exists(file.path(task_dirs, outcome_name)) & file.exists(file.path(task_dirs, siteyear_name))
  if (!all(done)) {
    missing <- c(missing, paste0("02_01", step, " ", temp_var, " ", group$site[!done]))
    next
  }

  outcome <- do.call(rbind, lapply(file.path(task_dirs, outcome_name), read_task_table,
                                   col_types = readr::cols(site_ID = readr::col_character(), .default = readr::col_double())))
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

if (length(missing) > 0) {
  stop(length(missing), " tasks have no results (failed or not finished), so their step and temperature variable ",
       "were not combined:\n", paste(missing, collapse = "\n"))
}
