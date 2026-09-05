test_that("the cumulative incidences are the ones cmprsk::cuminc produces", {
  d  <- prep_sim()
  fc <- fit_competing(d, group = "site", times = c(12, 36, 60),
                      cause_specific = "skip")

  ref <- cmprsk::cuminc(ftime = d$time_mo, fstatus = d$event_cr, group = d$site,
                        cencode = 0)
  at <- function(nm, t) {
    e <- ref[[nm]]
    e$est[max(which(e$time <= t))]
  }
  full <- fc[fc$analysis_set == "full", ]
  for (i in seq_len(nrow(full))) {
    nm <- paste(full$group[i], full$cause[i])
    expect_equal(full$cif[i], at(nm, full$time[i]), tolerance = 1e-9,
                 info = paste(nm, full$time[i]))
  }

  # Gray's test, repeated on every row of its cause
  expect_equal(unique(full$gray_stat[full$cause == 1]),
               unname(ref$Tests[1, "stat"]), tolerance = 1e-8)
  expect_equal(unique(full$gray_p[full$cause == 1]),
               unname(ref$Tests[1, "pv"]), tolerance = 1e-8)

  # a cumulative incidence is non-decreasing, and the two causes together
  # cannot exceed 1
  for (g in levels(d$site)) {
    for (cs in 1:2)
      expect_false(is.unsorted(full$cif[full$group == g & full$cause == cs]),
                   info = paste(g, cs))
    tot <- tapply(full$cif[full$group == g], full$time[full$group == g], sum)
    expect_true(all(tot <= 1))
  }
  expect_true(all(full$cif_lcl <= full$cif, full$cif >= 0,
                  full$cif_ucl >= full$cif, full$cif_ucl <= 1))
})

test_that("the Fine-Gray coefficients are the ones cmprsk::crr produces", {
  d    <- prep_sim()
  vars <- c("site", "age_grp", "surgery")     # surgery has missing values
  fc   <- fit_competing(d, covariates = vars, times = 60,
                        cause_specific = "skip")
  fg   <- attr(fc, "finegray")

  cc <- d[stats::complete.cases(d[, c("time_mo", "event_cr", vars)]), ]
  X  <- stats::model.matrix(~ site + age_grp + surgery, data = cc)[, -1, drop = FALSE]
  for (cs in 1:2) {
    ref <- cmprsk::crr(ftime = cc$time_mo, fstatus = cc$event_cr, cov1 = X,
                       failcode = cs, cencode = 0)
    got <- fg[fg$failcode == cs & !fg$is_reference, ]
    expect_equal(got$term, colnames(X))
    expect_equal(got$coef, unname(ref$coef), tolerance = 1e-6)
    expect_equal(got$shr,  unname(exp(ref$coef)), tolerance = 1e-6)
    expect_equal(got$se,   unname(sqrt(diag(ref$var))), tolerance = 1e-6)
    expect_true(all(got$converged))
    expect_equal(unique(got$n_model), nrow(cc))
  }

  # reference rows carry no estimate, like everywhere else in the package
  refs <- fg[fg$is_reference, ]
  expect_true(all(refs$shr == 1))
  expect_true(all(is.na(refs$ci_low), is.na(refs$ci_high), is.na(refs$p_value)))

  # the listwise deletion is reported, not quietly absorbed
  drop <- attr(fc, "dropna")
  tot  <- drop[drop$is_total, ]
  expect_equal(nrow(tot), 1L)
  expect_equal(tot$n_dropped, nrow(d) - nrow(cc))
  expect_gt(tot$n_dropped, 0L)          # or this proves nothing
  expect_equal(tot$n, nrow(d))
})

test_that("an event_os that contradicts event_cr stops the call", {
  d <- prep_sim()

  # flip one subject's overall-survival indicator only. Both columns remain
  # individually valid; only the pair is wrong, which is exactly what a check
  # on either column alone would miss.
  bad <- d
  bad$event_os[which(bad$event_os == 1L & bad$event_cr == 1L)[1]] <- 0L
  expect_error(fit_competing(bad, group = "site", times = 12), "event")

  # and the other direction: a censored subject given a cause of death
  bad2 <- d
  bad2$event_cr[which(bad2$event_os == 0L)[1]] <- 2L
  expect_error(fit_competing(bad2, group = "site", times = 12), "event")

  # on good data the check reports what it verified rather than only failing
  ok  <- fit_competing(d, group = "site", times = 12, cause_specific = "skip")
  chk <- attr(ok, "outcome_check")
  expect_true(all(chk$checks))
  expect_equal(length(chk$checks), 4L)
  expect_equal(as.vector(chk$table), as.vector(table(d$event_os, d$event_cr)))
  expect_equal(sum(chk$counts$n), nrow(d))
  # the off-diagonal cells are the contradictions, and must be empty
  expect_equal(chk$counts$n[chk$counts$event_os == 0L & chk$counts$event_cr != 0L],
               c(0L, 0L))
  expect_equal(chk$counts$n[chk$counts$event_os == 1L & chk$counts$event_cr == 0L], 0L)
})

test_that("the outcome codes are arguments, not assumptions", {
  d <- prep_sim()

  # Re-code the outcome away from the 0/1/2 default entirely: censoring 9,
  # the disease 4, other causes 7. Nothing about the estimates should change.
  # Without this the codes could be hard-wired anywhere downstream and every
  # other test would still pass, because every other test uses the default.
  e <- d
  e$cr9 <- c(9L, 4L, 7L)[d$event_cr + 1L]

  a <- fit_competing(d, group = "site", times = c(12, 60),
                     cause_specific = "skip")
  b <- fit_competing(e, group = "site", event_cr = "cr9",
                     causes = c(4, 7), censor_code = 9,
                     times = c(12, 60), cause_specific = "skip")

  expect_equal(b$cif,     a$cif,     tolerance = 1e-12)
  expect_equal(b$cif_lcl, a$cif_lcl, tolerance = 1e-12)
  expect_equal(b$cif_ucl, a$cif_ucl, tolerance = 1e-12)
  expect_equal(b$gray_p,  a$gray_p,  tolerance = 1e-12)
  expect_equal(b$n_event_cause, a$n_event_cause)
  expect_equal(b$cause, rep(c(4, 7), each = nrow(b) / 2)[order(order(b$cause))])

  # the labels travel with the codes
  lab <- fit_competing(e, group = "site", event_cr = "cr9",
                       causes = c(4, 7), censor_code = 9,
                       cause_labels = c("GIST death", "Other death"),
                       times = 12, cause_specific = "skip")
  expect_setequal(unique(lab$cause_label), c("GIST death", "Other death"))

  # and a code that is declared but absent, or present but undeclared, is an
  # error rather than a silently empty cause
  expect_error(fit_competing(e, group = "site", event_cr = "cr9",
                             causes = c(4, 7), censor_code = 0, times = 12),
               "0|censor")
})
