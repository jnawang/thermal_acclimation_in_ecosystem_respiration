# This script demonstrates identifiability of direct and apparent TAS components
# with simulation experiments where the true contributions are known.
# It also includes a confounded case showing when the apparent residual is not
# interpretable as pure mediation via SWC and GPP.
# 
# Authors: Sparkle Malone and Junna Wang in May, 2026.
# It takes about 1-2 minutes to run on a laptop.

# Content of this script: 
# for a virtue site, we need to have T~ER relationship as ER = exp(0.0588 * T)
# this site has 4 windows, 8 years
# for each window, we have reference temperature, SWC, and GPP of each window, we calculate reference ER
# we designed four scenarios: direct effect only, measured indirect effect only, 
# direct + measured indirect, direct + measured indirect + hidden indirect. 
# each scenario has 40 replicates. 

rm(list = ls())

set.seed(20260423)

dir.create("data", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

simulate_scenario <- function(
  scenario,
  beta_direct,
  lambda_swc,
  lambda_gpp,
  gamma_swc,
  gamma_gpp,
  beta_hidden = 0,
  rho_hidden = 0,
  n_sites = 1,
  n_windows = 4,
  n_years = 8,
  n_obs = 50,
  n_rep = 40
) {
  year_signal <- seq(-1.5, 1.5, length.out = n_years)
  out <- vector("list", n_rep)

  true_apparent <- gamma_swc * lambda_swc + gamma_gpp * lambda_gpp
  true_hidden <- beta_hidden * rho_hidden * (lambda_swc + lambda_gpp)
  true_residual <- true_apparent + true_hidden
  true_total <- beta_direct + true_residual

  for (irep in seq_len(n_rep)) {
    summary_list <- vector("list", n_sites * n_windows * n_years)
    idx <- 1

    for (isite in seq_len(n_sites)) {
      site_effect <- rnorm(1, mean = 0, sd = 0.15)

      for (iwindow in seq_len(n_windows)) {
        site_window <- paste0("S", isite, "_W", iwindow)
        ref_swc <- runif(1, min = 0.6, max = 1.8)
        ref_gpp <- runif(1, min = 0.8, max = 2.2)
        window_effect <- rnorm(1, mean = 0, sd = 0.2)

        for (iyear in seq_len(n_years)) {
          t_anom <- year_signal[iyear] + rnorm(1, mean = 0, sd = 0.08)
          swc_mean <- ref_swc + lambda_swc * t_anom + rnorm(1, mean = 0, sd = 0.05)
          gpp_mean <- ref_gpp + lambda_gpp * t_anom + rnorm(1, mean = 0, sd = 0.05)

          ts_center <- rnorm(n_obs, mean = 0, sd = 0.9)
          swc_dev <- rnorm(n_obs, mean = swc_mean - ref_swc, sd = 1.25)  # sd = 0.12
          gpp_dev <- rnorm(n_obs, mean = gpp_mean - ref_gpp, sd = 1.5)  # sd = 0.16
          hidden_dev <- rho_hidden * (swc_dev + gpp_dev) + rnorm(n_obs, mean = 0, sd = 0.10)

          direct_shift <- beta_direct * t_anom
          apparent_shift <- gamma_swc * (swc_mean - ref_swc) + gamma_gpp * (gpp_mean - ref_gpp)
          hidden_shift <- beta_hidden * rho_hidden * ((swc_mean - ref_swc) + (gpp_mean - ref_gpp))

          log_er <- site_effect + window_effect +
            direct_shift +
            0.0588 * ts_center +      # this is that realistic temperature ER relationship. 
            gamma_swc * swc_dev +
            gamma_gpp * gpp_dev +
            beta_hidden * hidden_dev +
            rnorm(n_obs, mean = 0, sd = 0.18)

          dat <- data.frame(
            log_er = log_er,
            ts_center = ts_center,
            swc_dev = swc_dev,
            gpp_dev = gpp_dev
          )

          fit_tot <- lm(log_er ~ ts_center, data = dat)
          fit_dir <- lm(log_er ~ ts_center + swc_dev + gpp_dev, data = dat)

          pred_tot <- unname(predict(fit_tot, newdata = data.frame(ts_center = 0)))
          pred_dir <- unname(predict(fit_dir, newdata = data.frame(ts_center = 0, swc_dev = 0, gpp_dev = 0)))

          summary_list[[idx]] <- data.frame(
            site_window = site_window,
            t_anom = t_anom,
            direct_true = direct_shift,
            apparent_true = apparent_shift,
            hidden_true = hidden_shift,
            residual_true = apparent_shift + hidden_shift,
            total_true = direct_shift + apparent_shift + hidden_shift,
            total_hat = pred_tot,
            direct_hat = pred_dir
          )
          idx <- idx + 1
        }
      }
    }

    df <- do.call(rbind, summary_list)
    df$apparent_hat <- df$total_hat - df$direct_hat

    fit_total <- summary(lm(total_hat ~ t_anom + site_window, data = df))
    fit_direct <- summary(lm(direct_hat ~ t_anom + site_window, data = df))
    fit_app <- summary(lm(apparent_hat ~ t_anom + site_window, data = df))

    out[[irep]] <- data.frame(
      scenario = scenario,
      replicate = irep,
      total_est = fit_total$coefficients["t_anom", "Estimate"],
      total_se = fit_total$coefficients["t_anom", "Std. Error"],
      direct_est = fit_direct$coefficients["t_anom", "Estimate"],
      direct_se = fit_direct$coefficients["t_anom", "Std. Error"],
      apparent_est = fit_app$coefficients["t_anom", "Estimate"],
      apparent_se = fit_app$coefficients["t_anom", "Std. Error"],
      total_true = true_total,
      direct_true = beta_direct,
      apparent_true = true_apparent,
      residual_true = true_residual,
      hidden_true = true_hidden
    )
  }

  do.call(rbind, out)
}

summarise_component <- function(df, component, truth_col) {
  est_col <- paste0(component, "_est")
  se_col <- paste0(component, "_se")
  truth <- df[[truth_col]]
  est <- df[[est_col]]
  se <- df[[se_col]]

  data.frame(
    component = component,
    target = truth_col,
    truth = unique(truth),
    mean_est = mean(est),
    sd_est = sd(est),
    bias = mean(est - truth),
    rmse = sqrt(mean((est - truth)^2)),
    coverage = mean((est - 1.96 * se) <= truth & (est + 1.96 * se) >= truth)
  )
}

scenario_defs <- list(
  list(
    scenario = "Direct only",
    beta_direct = 0.040,
    lambda_swc = 0.00,
    lambda_gpp = 0.00,
    gamma_swc = 0.06,
    gamma_gpp = 0.07,
    beta_hidden = 0.00,
    rho_hidden = 0.00
  ),
  list(
    scenario = "Measured indirect only",
    beta_direct = 0.000,
    lambda_swc = 0.18,
    lambda_gpp = 0.13,
    gamma_swc = 0.06,
    gamma_gpp = 0.07,
    beta_hidden = 0.00,
    rho_hidden = 0.00
  ),
  list(
    scenario = "Direct + measured indirect",
    beta_direct = 0.030,
    lambda_swc = 0.12,
    lambda_gpp = 0.11,
    gamma_swc = 0.06,
    gamma_gpp = 0.07,
    beta_hidden = 0.00,
    rho_hidden = 0.00
  ),
  list(
    scenario = "Direct + measured indirect + hidden indirect",
    beta_direct = 0.030,
    lambda_swc = 0.12,
    lambda_gpp = 0.11,
    gamma_swc = 0.06,
    gamma_gpp = 0.07,
    beta_hidden = 0.10,
    rho_hidden = 0.80
  )
)

replicate_results <- do.call(
  rbind,
  lapply(scenario_defs, function(x) {
    cat("Running scenario:", x$scenario, "\n")
    do.call(simulate_scenario, x)
  })
)

summary_results <- do.call(
  rbind,
  lapply(split(replicate_results, replicate_results$scenario), function(df) {
    rbind(
      summarise_component(df, "total", "total_true"),
      summarise_component(df, "direct", "direct_true"),
      summarise_component(df, "apparent", "apparent_true"),
      summarise_component(df, "apparent", "residual_true")
    )
  })
)
summary_results$scenario <- rep(names(split(replicate_results, replicate_results$scenario)), each = 4)
summary_results <- summary_results[, c("scenario", "component", "target", "truth", "mean_est", "sd_est", "bias", "rmse", "coverage")]

write.csv(replicate_results, file = "data/simulation_identifiability_replications.csv", row.names = FALSE)
write.csv(summary_results, file = "data/simulation_identifiability_summary.csv", row.names = FALSE)

plot_panel <- function(df, title_text, apparent_note = FALSE) {
  vals <- list(df$total_est, df$direct_est, df$apparent_est)
  cols <- c("#4C78A8", "#54A24B", "#F58518")
  boxplot(
    vals,
    names = c("Total", "Direct", "Indirect"),
    col = cols,
    border = NA,
    ylab = "Estimated TRS slope",
    ylim = range(c(unlist(vals), df$total_true, df$direct_true, df$apparent_true, df$residual_true)) + c(-0.02, 0.02)
  )
  abline(h = unique(df$total_true), col = cols[1], lwd = 2, lty = 2)
  abline(h = unique(df$direct_true), col = cols[2], lwd = 2, lty = 2)
  abline(h = unique(df$apparent_true), col = cols[3], lwd = 2, lty = 2)
  if (apparent_note) {
    abline(h = unique(df$residual_true), col = "#B22222", lwd = 2, lty = 3)
    legend(
      "topleft",
      legend = c("true total", "true direct", "true measured indirect", "true residual = measured + hidden indirect"),
      col = c(cols[1], cols[2], cols[3], "#B22222"),
      lty = c(2, 2, 2, 3),
      lwd = 2,
      bty = "n",
      cex = 0.82
    )
  } else {
    legend(
      "topleft",
      legend = c("true total", "true direct", "true measured indirect"),
      col = cols,
      lty = 2,
      lwd = 2,
      bty = "n",
      cex = 0.82
    )
  }
  mtext(title_text, side = 3, line = 0.8, font = 2, cex = 0.95)
  text(
    x = c(1, 2, 3),
    y = sapply(vals, function(v) quantile(v, 0.95)) + 0.001,
    labels = sprintf("%.3f", c(mean(df$total_est), mean(df$direct_est), mean(df$apparent_est))),
    cex = 0.78
  )
}

png(
  filename = "figures/tas_identifiability_simulation.png",
  width = 1800,
  height = 1600,
  res = 180
)
op <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), mar = c(4.8, 4.8, 3.5, 1.2), oma = c(0, 0, 2.4, 0))

plot_panel(subset(replicate_results, scenario == "Direct only"), "A. Direct only")
plot_panel(subset(replicate_results, scenario == "Measured indirect only"), "B. Measured indirect only")
plot_panel(subset(replicate_results, scenario == "Direct + measured indirect"), "C. Direct + measured indirect")
plot_panel(subset(replicate_results, scenario == "Direct + measured indirect + hidden indirect"), "D. Direct + measured indirect + hidden indirect", apparent_note = TRUE)

mtext(
  "Simulation-based identifiability of direct and indirect TRS components",
  outer = TRUE,
  cex = 1.25,
  font = 2
)
par(op)
dev.off()

cat("\nSaved replicate-level results to data/simulation_identifiability_replications.csv\n")
cat("Saved summary results to data/simulation_identifiability_summary.csv\n")
cat("Saved figure to figures/tas_identifiability_simulation.png\n\n")

print(summary_results)
