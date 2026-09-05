test_that("the hazard ratios are the ones coxph() produces", {
  d    <- prep_sim()
  vars <- c("site", "age_grp", "sex", "size_grp", "surgery")
  cx   <- fit_cox(d, covariates = vars)

  ref <- survival::coxph(
    survival::Surv(time_mo, event_os) ~ site + age_grp + sex + size_grp + surgery,
    data = d, ties = "efron")
  s  <- summary(ref)
  ci <- stats::confint(ref)

  non_ref <- cx[!cx$is_reference, ]
  expect_equal(non_ref$term, rownames(s$coefficients))
  expect_equal(non_ref$coef,    unname(s$coefficients[, "coef"]),  tolerance = 1e-10)
  expect_equal(non_ref$se_coef, unname(s$coefficients[, "se(coef)"]), tolerance = 1e-10)
  expect_equal(non_ref$hr,      unname(s$coefficients[, "exp(coef)"]), tolerance = 1e-10)
  expect_equal(non_ref$p_value, unname(s$coefficients[, "Pr(>|z|)"]), tolerance = 1e-10)
  expect_equal(non_ref$ci_low,  unname(exp(ci[, 1])), tolerance = 1e-10)
  expect_equal(non_ref$ci_high, unname(exp(ci[, 2])), tolerance = 1e-10)

  expect_equal(unique(cx$n_model), ref$n)
  expect_equal(unique(cx$n_event_model), ref$nevent)
})

test_that("cox.zph is reported, and the transform argument actually changes it", {
  d    <- prep_sim()
  vars <- c("site", "age_grp", "sex")
  fml  <- survival::Surv(time_mo, event_os) ~ site + age_grp + sex
  ref  <- survival::coxph(fml, data = d, ties = "efron")

  usable <- c("km", "rank", "identity")
  for (tr in usable) {
    cx <- fit_cox(d, covariates = vars, zph_transform = tr)
    z  <- attr(cx, "cox_zph")
    want <- survival::cox.zph(ref, transform = tr)$table
    expect_equal(z$term, rownames(want))
    expect_equal(z$chisq, unname(want[, "chisq"]), tolerance = 1e-8)
    expect_equal(z$df,    unname(want[, "df"]))
    expect_equal(z$p,     unname(want[, "p"]), tolerance = 1e-8)
    expect_equal(sum(z$is_global), 1L)
  }

  # the transforms are not interchangeable; a test that passed whichever one
  # was used would not be testing anything
  p_by <- vapply(usable, function(tr) {
    z <- attr(fit_cox(d, covariates = vars, zph_transform = tr), "cox_zph")
    z$p[z$is_global]
  }, numeric(1))
  expect_gt(diff(range(p_by)), 0)

  # transform = "log" is not usable here, and must say so rather than return a
  # missing test. sim_gist has events at month 0, as a registry extract does,
  # and log(0) is not finite -- cox.zph() would hand back NA and a caller who
  # did not look would read "no PH problem" off a test that never ran.
  expect_error(fit_cox(d, covariates = vars, zph_transform = "log"),
               "transform")
  expect_gt(sum(d$time_mo == 0 & d$event_os == 1L), 0L)

  # on data without time-0 events it works, and agrees with cox.zph()
  d2  <- d[d$time_mo > 0, ]
  cx2 <- fit_cox(d2, covariates = vars, zph_transform = "log")
  ref2 <- survival::coxph(fml, data = d2, ties = "efron")
  want <- survival::cox.zph(ref2, transform = "log")$table
  expect_equal(attr(cx2, "cox_zph")$p, unname(want[, "p"]), tolerance = 1e-8)
})

test_that("reference rows are marked, and the deletion counts add up", {
  d  <- prep_sim()
  cx <- fit_cox(d, covariates = c("site", "age_grp", "size_grp", "surgery"))

  refs <- cx[cx$is_reference, ]
  expect_equal(nrow(refs), 4L)                     # one per factor
  expect_true(all(refs$hr == 1))
  expect_true(all(is.na(refs$ci_low), is.na(refs$ci_high)))
  expect_true(all(refs$hr_txt == "1.00 (reference)"))
  expect_true(all(is.na(refs$p_value)))

  # the first level of each factor is the reference
  for (v in c("site", "age_grp", "size_grp", "surgery"))
    expect_equal(cx$level[cx$variable == v][1], levels(d[[v]])[1])

  # listwise deletion is accounted for, not silently absorbed
  cnt <- attr(cx, "fit_cox")$counts
  g <- function(k) as.integer(cnt[[k]])
  expect_equal(g("rows_in"), nrow(d))
  expect_equal(g("rows_in"), g("rows_used") + g("rows_dropped"))
  expect_equal(g("events_in"), g("events_used") + g("events_dropped"))
  expect_equal(g("rows_dropped"),
               sum(is.na(d$size_grp) | is.na(d$surgery)))
  expect_equal(g("rows_used"), unique(cx$n_model))

  # a declared level that no row has: dropped with a warning and counted,
  # rather than left as an all-zero design column
  sub <- d[d$site != "Colorectal", ]           # level kept, rows gone
  expect_warning(cx2 <- fit_cox(sub, covariates = c("site", "sex")), "level")
  expect_false("Colorectal" %in% cx2$level)
  st <- attr(cx2, "fit_cox")$settings
  expect_equal(unname(st$levels_dropped[["site"]]), 1L)
  expect_equal(unname(st$levels_dropped[["sex"]]), 0L)
  expect_false("Colorectal" %in% st$levels[["site"]])

  # and the reference level is untouched by that
  expect_equal(cx2$level[cx2$variable == "site"][1], "Stomach")
  expect_true(cx2$is_reference[cx2$variable == "site"][1])
})
