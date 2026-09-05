test_that("the three estimands are the ones the three functions produce on their own", {
  d    <- two_arm("Stomach", "Small intestine")
  arms <- c("Stomach", "Small intestine")
  tau  <- c(12, 36, 60)
  ce   <- compare_estimands(d, group = "site", arms = arms, tau = tau)

  # this is the test that matters for a function whose job is composition:
  # if the wiring were wrong, every column would still look plausible
  d1 <- d
  d1$site <- droplevels(d1$site)   # the third site is not in this subset
  cx <- fit_cox(d1, covariates = "site")
  row <- cx[cx$level == "Small intestine", ]
  expect_equal(unique(ce$hr),     row$hr,      tolerance = 1e-9)
  expect_equal(unique(ce$hr_lcl), row$ci_low,  tolerance = 1e-9)
  expect_equal(unique(ce$hr_ucl), row$ci_high, tolerance = 1e-9)
  expect_equal(unique(ce$hr_p),   row$p_value, tolerance = 1e-9)

  rm2 <- calc_rmst(d, group = "site", arms = arms, tau = tau)
  expect_equal(ce$rmst_diff,     rm2$rmst_diff,     tolerance = 1e-9)
  expect_equal(ce$rmst_diff_lcl, rm2$rmst_diff_lcl, tolerance = 1e-9)
  expect_equal(ce$rmst_diff_ucl, rm2$rmst_diff_ucl, tolerance = 1e-9)
  expect_equal(ce$rmst_diff_p,   rm2$rmst_diff_p,   tolerance = 1e-9)

  km <- fit_km(d, group = "site", times = tau)
  kt <- attr(km, "km_times")
  for (i in seq_along(tau)) {
    expect_equal(ce$km_surv_arm0[i],
                 kt$surv[kt$level == arms[1] & kt$time == tau[i]], tolerance = 1e-9)
    expect_equal(ce$km_surv_arm1[i],
                 kt$surv[kt$level == arms[2] & kt$time == tau[i]], tolerance = 1e-9)
  }
  expect_equal(ce$km_surv_diff, ce$km_surv_arm1 - ce$km_surv_arm0, tolerance = 1e-12)
  expect_equal(ce$km_surv_diff_pp, 100 * ce$km_surv_diff, tolerance = 1e-12)
})

test_that("the concordance verdict follows from the three significance flags", {
  d <- prep_sim()

  # a contrast where the three do not agree, and one where they do. Both must
  # exist in sim_gist, or the function's whole reason for being is untested.
  disc <- compare_estimands(d, group = "sex", arms = c("Male", "Female"),
                            tau = c(12, 36, 60))
  conc <- compare_estimands(d, group = "site", arms = c("Stomach", "Colorectal"),
                            tau = c(12, 36, 60))

  at60 <- function(x) x[x$tau == 60, ]
  a <- at60(disc); b <- at60(conc)

  # the flags must be the ones the intervals imply
  expect_equal(a$sig_hr,   a$hr_lcl > 1 | a$hr_ucl < 1)
  expect_equal(a$sig_km,   a$km_surv_diff_lcl > 0 | a$km_surv_diff_ucl < 0)
  expect_equal(a$sig_rmst, a$rmst_diff_lcl > 0 | a$rmst_diff_ucl < 0)

  # and the verdict must be the one those flags describe
  expect_equal(a$n_signif, sum(a$sig_hr, a$sig_km, a$sig_rmst))
  expect_equal(b$n_signif, sum(b$sig_hr, b$sig_km, b$sig_rmst))
  expect_lt(a$n_signif, 3L)
  expect_gt(a$n_signif, 0L)
  expect_equal(b$n_signif, 3L)
  expect_match(a$verdict, "DISCORDANT")
  expect_match(b$verdict, "concordant")
  expect_false(a$concordant)
  expect_true(b$concordant)
})

test_that("contrast-level columns are repeated so any single row can be quoted", {
  d   <- two_arm("Stomach", "Small intestine")
  tau <- c(60, 12, 36)                       # deliberately out of order
  ce  <- compare_estimands(d, group = "site",
                           arms = c("Stomach", "Small intestine"), tau = tau)

  expect_equal(ce$tau, sort(tau))            # returned sorted regardless

  # the hazard ratio does not depend on tau, so it must be identical on every
  # row; a per-row value would mean the model was refitted per tau
  for (col in c("hr", "hr_lcl", "hr_ucl", "hr_p", "ph_p", "logrank_p",
                "n_arm0", "n_arm1", "events_arm0", "events_arm1", "tau_upper"))
    expect_equal(length(unique(ce[[col]])), 1L, info = col)

  # the tau-specific ones must not be
  for (col in c("km_surv_arm0", "rmst_arm0", "rmst_diff"))
    expect_equal(length(unique(ce[[col]])), length(tau), info = col)

  # tau past the follow-up is refused here too, by the same guard
  edge <- unique(ce$tau_upper)
  expect_error(compare_estimands(d, group = "site",
                                 arms = c("Stomach", "Small intestine"),
                                 tau = edge + 1), "follow-up")
})
