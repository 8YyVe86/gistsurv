test_that("the medians and counts agree with survfit() called directly", {
  d <- prep_sim()
  k <- fit_km(d, group = "site", times = c(12, 36, 60))

  ref <- survival::survfit(survival::Surv(time_mo, event_os) ~ site, data = d,
                           conf.type = "log")
  s <- summary(ref)$table

  # the strata come back in the order survfit() names them
  expect_equal(k$stratum, rownames(s))
  expect_equal(k$n, unname(s[, "records"]))
  expect_equal(k$events, unname(s[, "events"]))
  expect_equal(k$median, unname(s[, "median"]))
  expect_equal(k$median_lcl, unname(s[, "0.95LCL"]))
  expect_equal(k$median_ucl, unname(s[, "0.95UCL"]))

  # survival probabilities at the requested times, from the same fit
  times <- attr(k, "km_times")
  for (i in seq_len(nrow(times))) {
    lev <- times$level[i]
    sub <- d[d$site == lev, ]
    f1 <- survival::survfit(survival::Surv(time_mo, event_os) ~ 1, data = sub,
                            conf.type = "log")
    expect_equal(times$surv[i], summary(f1, times = times$time[i])$surv,
                 tolerance = 1e-10)
  }

  # the log-rank test, repeated on every row, is the one survdiff computes
  sd <- survival::survdiff(survival::Surv(time_mo, event_os) ~ site, data = d)
  expect_equal(unique(k$logrank_chisq), unname(sd$chisq), tolerance = 1e-10)
  expect_equal(unique(k$logrank_df), length(sd$n) - 1L)
})

test_that("the reverse-Kaplan-Meier follow-up is the reverse fit, not the forward one", {
  d <- prep_sim()
  k <- fit_km(d, group = "site")

  for (lev in levels(d$site)) {
    sub <- d[d$site == lev, ]
    rev <- survival::survfit(survival::Surv(time_mo, 1L - event_os) ~ 1, data = sub)
    want <- unname(summary(rev)$table["median"])
    expect_equal(k$median_fu_revkm[k$level == lev], want, tolerance = 1e-10)
  }

  # person-time is the plain sum of follow-up, and the strata partition it
  expect_equal(sum(k$total_fu), sum(d$time_mo))
  expect_equal(k$max_fu,
               as.vector(tapply(d$time_mo, d$site, max)[levels(d$site)]))
})

test_that("a time point past a stratum's follow-up is refused, or returned as NA", {
  d <- prep_sim()
  worst <- min(tapply(d$time_mo, d$site, max))

  # the closed interval: the largest admissible time point is allowed
  expect_no_error(fit_km(d, group = "site", times = worst))
  expect_error(fit_km(d, group = "site", times = worst + 1), "follow-up")

  na_ok <- fit_km(d, group = "site", times = worst + 1, times_beyond = "na")
  tm <- attr(na_ok, "km_times")
  short <- names(which.min(tapply(d$time_mo, d$site, max)))
  expect_false(tm$estimable[tm$level == short])
  expect_true(is.na(tm$surv[tm$level == short]))
})
