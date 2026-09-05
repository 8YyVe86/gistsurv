test_that("the RMST estimates are the ones survRM2 produces", {
  d   <- two_arm("Stomach", "Small intestine")
  tau <- c(12, 36, 60)
  r   <- calc_rmst(d, group = "site", arms = c("Stomach", "Small intestine"),
                   tau = tau)

  expect_equal(r$tau, sort(tau))          # rows are ordered by tau

  # survRM2 takes arm as 0/1; arm1 is the second element of `arms`
  arm <- as.integer(d$site == "Small intestine")
  for (i in seq_along(r$tau)) {
    ref <- survRM2::rmst2(time = d$time_mo, status = d$event_os, arm = arm,
                          tau = r$tau[i])
    expect_equal(r$rmst_arm0[i], unname(ref$RMST.arm0$rmst["Est."]),
                 tolerance = 1e-9)
    expect_equal(r$rmst_arm1[i], unname(ref$RMST.arm1$rmst["Est."]),
                 tolerance = 1e-9)
    diff_row <- ref$unadjusted.result["RMST (arm=1)-(arm=0)", ]
    expect_equal(r$rmst_diff[i],     unname(diff_row[1]), tolerance = 1e-9)
    expect_equal(r$rmst_diff_lcl[i], unname(diff_row[2]), tolerance = 1e-9)
    expect_equal(r$rmst_diff_ucl[i], unname(diff_row[3]), tolerance = 1e-9)
    expect_equal(r$rmst_diff_p[i],   unname(diff_row[4]), tolerance = 1e-9)
  }

  # RMST is a mean survival time over [0, tau], so it cannot exceed tau, and
  # a longer horizon cannot make it smaller
  expect_true(all(r$rmst_arm0 <= r$tau, r$rmst_arm1 <= r$tau))
  expect_true(all(diff(r$rmst_arm0) > 0, diff(r$rmst_arm1) > 0))

  # the direction columns must agree with the sign of the difference, rather
  # than needing to be re-derived by eye
  expect_equal(r$arm1_worse, r$rmst_diff < 0)
  expect_equal(r$diff_significant, r$rmst_diff_lcl > 0 | r$rmst_diff_ucl < 0)
})

test_that("tau past the follow-up of an arm is refused, and tau at the edge is not", {
  d    <- two_arm("Stomach", "Colorectal")
  arms <- c("Stomach", "Colorectal")
  edge <- min(tapply(d$time_mo, droplevels(d$site), max))

  # the guard is a closed interval: tau == tau_upper is admissible
  ok <- calc_rmst(d, group = "site", arms = arms, tau = edge)
  expect_equal(ok$tau_upper[1], edge)
  expect_error(calc_rmst(d, group = "site", arms = arms, tau = edge + 1),
               "follow-up")

  # one bad value in a vector stops the whole call rather than being dropped
  expect_error(calc_rmst(d, group = "site", arms = arms,
                         tau = c(12, 36, edge + 1)), "follow-up")

  # this is what the guard is for: survRM2 on its own extrapolates instead of
  # complaining, so the wrapper has to be the one that refuses
  expect_no_error(survRM2::rmst2(time = d$time_mo, status = d$event_os,
                                 arm = as.integer(d$site == "Colorectal"),
                                 tau = edge))
})

test_that("an arm too small to estimate warns but still returns its numbers", {
  d <- two_arm("Stomach", "Colorectal")
  arms <- c("Stomach", "Colorectal")

  # thin the second arm below min_n. This is a warning, not an error: the
  # large-sample interval is unreliable, not undefined, and refusing to show
  # the estimate would not make the user's data any bigger.
  keep <- c(which(d$site == "Stomach"), utils::head(which(d$site == "Colorectal"), 5))
  thin <- d[keep, ]
  expect_warning(r <- calc_rmst(thin, group = "site", arms = arms, tau = 12,
                                min_n = 10L, min_events = 0L),
                 "subjects")
  expect_equal(nrow(r), 1L)
  expect_false(is.na(r$rmst_diff))
  expect_equal(r$n_arm1, 5L)

  # and the documented escape hatch really silences it
  expect_silent(calc_rmst(thin, group = "site", arms = arms, tau = 12,
                          min_n = 0L, min_events = 0L))

  # a level that is not in the data at all is a different matter: there is
  # nothing to estimate, so that one stops
  expect_error(calc_rmst(d, group = "site", arms = c("Stomach", "Pancreas"),
                         tau = 12), "Pancreas")

  # so is a grouping column with missing values, unless told what to do
  na_site <- d
  na_site$site[1:5] <- NA
  expect_error(calc_rmst(na_site, group = "site", arms = arms, tau = 12),
               "missing")
  expect_no_error(calc_rmst(na_site, group = "site", arms = arms, tau = 12,
                            group_missing = "drop"))
})
