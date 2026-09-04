# ---------------------------------------------------------------------------
# fit_competing(): competing-risk analysis -- cumulative incidence functions,
# Gray's test and Fine-Gray subdistribution hazard models.
#
# Design notes (deliberate, and different from the study script this was
# distilled from):
#   * no column name and no cause code is hard-coded: the time, the two
#     outcome columns, the grouping variable, the covariates, the censoring
#     code and the cause codes are all arguments;
#   * nothing is printed, nothing is plotted and nothing is written to disk.
#     The study script prints table(event_os, event_cr) to satisfy a project
#     rule that the two outcome columns be cross-checked. A package function
#     cannot discharge that rule by printing: printed output is not available
#     to the caller and disappears from a non-interactive run. The check is
#     therefore *enforced* (a disagreement is an error, not a message) and the
#     cross-tabulation is *returned*, in attr(x, "outcome_check");
#   * two analysis sets, and they are labelled. The cumulative incidence
#     functions need only time, cause and group, so they are estimated on
#     every usable row. A Fine-Gray model needs the covariates, so it is
#     fitted on the complete cases. The two n differ by construction; every
#     row carries an `analysis_set` label and the rows lost to casewise
#     deletion are counted, per group level, in attr(x, "dropna");
#   * confidence intervals for the cumulative incidence are built on the
#     complementary log-log scale, not by a Wald interval truncated into
#     [0, 1]. A truncated Wald interval reads 0.000 (0.000-0.000) in a cell
#     with no events, and [0, 0] is not a 95 per cent confidence interval.
#     Cells where the transform is undefined are returned as NA and labelled,
#     never as a degenerate interval;
#   * cmprsk::crr() has no formula interface, so the design matrix is built
#     here. That is where dummy coding goes wrong quietly -- a shifted column,
#     a dropped reference level, a factor level that lost all its rows to
#     casewise deletion -- so the matrix is checked against the levels it is
#     supposed to encode, column by column, before it is passed to crr();
#   * a time point past the largest follow-up time of a group is an error, not
#     a silently extrapolated number, which is the rule fit_km() and the RMST
#     tau guard already use;
#   * no fitted object is returned. A crr object carries its own call and the
#     unique failure times of the fitted data; a survfit object carries the
#     caller's data inside the environment of its formula. Both would put
#     individual records into anything the result is saved into.
# The file is kept ASCII-only so that it behaves the same under a UTF-8 and
# under a non-UTF-8 locale.
#
# The small formatting helpers .ps_values() and .ps_counts() are shared with
# R/prep_surv.R; they are internal to the package.
# ---------------------------------------------------------------------------

#' Competing risks: cumulative incidence, Gray's test and Fine-Gray models
#'
#' @description
#' `fit_competing()` estimates the cumulative incidence function (CIF) of each
#' cause of failure, optionally within the levels of one grouping variable,
#' compares those curves with Gray's test, and fits one Fine-Gray
#' subdistribution hazard model per cause. It returns the numbers a
#' competing-risk paper reports: the CIF and its confidence interval at a set
#' of time points, Gray's test, the subdistribution hazard ratios with their
#' intervals and p-values, the cause-specific hazard ratios to read beside
#' them, the rows lost to casewise deletion, and the cross-tabulation of the
#' overall and the competing-risk outcome column.
#'
#' The estimators are [cmprsk::cuminc()] and [cmprsk::crr()]. They are not
#' interchangeable with `riskRegression::FGR()` or with
#' [survival::finegray()] followed by [survival::coxph()]: those are different
#' programmes, and a number that came out of one cannot be traced back to the
#' other.
#'
#' The function prints nothing, plots nothing and writes nothing. The CIF
#' table comes back as a plain data frame and everything else is attached to
#' it as attributes.
#'
#' @details
#' # Two analysis sets, and why their n differ
#'
#' A cumulative incidence function needs three columns: the follow-up time,
#' the cause of failure and the grouping variable. A Fine-Gray model needs the
#' covariates as well, and covariates go missing. The two are therefore
#' estimated on different sets of rows, deliberately:
#'
#' * the CIF, Gray's test and the counts in the returned data frame use every
#'   usable row (`analysis_set = "full"`);
#' * the Fine-Gray models, and the cause-specific Cox models beside them, use
#'   the complete cases (`analysis_set = "complete_case"`).
#'
#' The two n are not meant to match, and a table that quotes both must say
#' which is which -- hence the `analysis_set` column on every row of every
#' table this function returns. The rows lost are counted in
#' `attr(x, "dropna")`, in total and per level of `group`, split by cause,
#' because a loss that falls unevenly across groups is not missing completely
#' at random and belongs in the limitations of the paper rather than in a
#' footnote. With `cif_complete_case = TRUE` the CIF is estimated a second
#' time on the complete cases, so that the question "did casewise deletion
#' change the cause-of-death structure?" has an answer instead of an
#' assumption.
#'
#' # Confidence intervals for the cumulative incidence
#'
#' [cmprsk::timepoints()] returns the CIF and its variance. A Wald interval
#' built from them, `F +/- z * SE(F)`, has to be truncated into `[0, 1]`, and
#' in a cell with few events the lower limit is pinned at exactly 0; in a cell
#' with no events the interval collapses to `[0, 0]`. Declaring such an
#' interval truncated does not make it a confidence interval.
#'
#' The interval used here is built on the complementary log-log scale:
#'
#' \deqn{\nu = \log(-\log(1 - F)), \qquad
#'       SE(\nu) = \frac{SE(F)}{|(1 - F)\log(1 - F)|}, \qquad
#'       CI = 1 - \exp(-\exp(\nu \mp z \cdot SE(\nu))).}
#'
#' The transform maps `(0, 1)` onto the whole real line, so the
#' back-transformed limits fall strictly inside `(0, 1)` without any
#' truncation, and the interval is asymmetric where the estimate is small.
#' That the limits lie strictly inside `(0, 1)` and bracket the point estimate
#' is checked on every call.
#'
#' # Cells where no interval exists
#'
#' The transform is undefined at `F = 0` and at `F = 1`, and there is nothing
#' to transform when the estimate itself is missing. Such a cell gets
#' `cif_lcl = cif_ucl = NA`, `estimable = FALSE`, and a `ci_method` that says
#' which of the three it is:
#'
#' * `"not estimable (no events at this horizon)"` -- the group has had no
#'   event of that cause by that time, so `cif` is 0. This is the cell to
#'   report as "0 events, CI not estimable", with the number still at risk
#'   from `n_risk` beside it. It is never reported as `0.00 (0.00-0.00)`;
#' * `"not estimable (cumulative incidence is 1)"` -- everyone in the group
#'   has failed from that cause by that time;
#' * `"not estimable (the interval degenerates at 0 or 1)"` -- the transform is
#'   defined, but in floating point the back-transformed limits land on the
#'   boundary: an estimate the product-limit arithmetic returns as `1 - 2e-16`,
#'   or one whose standard error dwarfs it, drives `exp(-exp(.))` to 0 or to 1.
#'   The resulting `0.000-1.000` carries no information and is the degenerate
#'   interval this transform was chosen to avoid, so it is withdrawn rather
#'   than printed;
#' * `"not estimable (beyond the largest follow-up time in this group)"` --
#'   only reachable with `times_beyond = "na"`; see below.
#'
#' `cif == 0` and `n_event_cause_by_time == 0` are checked against each other
#' on every call, in both directions, so a cell labelled "no events" is one.
#' That an estimable interval lies strictly inside `(0, 1)` and brackets its
#' estimate is then true by construction, and is checked as well.
#'
#' # Time points, and the follow-up guard
#'
#' A cumulative incidence function is not estimable past the largest follow-up
#' time in a group. `fit_competing()` refuses to report one: a time point past
#' the largest follow-up of any group stops the call, and `times_beyond =
#' "na"` keeps the row, marks it `estimable = FALSE` and leaves the estimate
#' and its interval missing. The number still at risk is reported for such a
#' row, because an empty risk set is a fact rather than an extrapolation. This
#' is the same rule as the tau guard of a restricted-mean analysis, and the
#' largest usable time point comes back as `follow_up["tau_upper"]` in the
#' diagnostics.
#'
#' # The design matrix that goes into crr()
#'
#' [cmprsk::crr()] takes a numeric matrix, not a formula, so the dummy coding
#' is the caller's responsibility and its failure modes are silent: a column
#' in the wrong position attaches a coefficient to the wrong level, a
#' reference level dropped by [droplevels()] moves the reference without
#' changing any code, and a level emptied by casewise deletion leaves an
#' all-zero column. The matrix built here is checked before it is used:
#'
#' 1. it is built with `contr.treatment` contrasts, forced for the duration of
#'    the call, so the reference is the first level whatever
#'    `options("contrasts")` is set to session-wide;
#' 2. the intercept is dropped and the remaining column names must equal, in
#'    order, the non-reference levels of each covariate;
#' 3. no reference level may appear as a column, and no column name may be
#'    duplicated;
#' 4. every dummy column must be 0/1, and the mean of every column must equal
#'    the proportion of the complete cases at that level -- a continuous
#'    covariate must reproduce its own mean -- to a tolerance of 1e-12. A
#'    shifted column fails this immediately;
#' 5. after the fit, the names of the coefficients returned by `crr()` must
#'    equal the column names of the matrix, and `crr()` must report the same
#'    number of rows, no rows dropped for missing values, and convergence.
#'
#' # Subdistribution and cause-specific hazards
#'
#' The Fine-Gray model estimates a subdistribution hazard ratio. It keeps
#' subjects who died of the competing cause in the risk set, which is what
#' makes it a model for the cumulative incidence -- and what makes its
#' coefficients easy to misread as risks. When the competing cause is strong,
#' a covariate that raises death from cause 1 can show a subdistribution
#' hazard ratio below 1 for cause 2 simply because those subjects are being
#' used up by cause 1.
#'
#' The remedy is not to hide such coefficients but to put the cause-specific
#' hazard ratio beside them, so `fit_competing()` fits both by default. The
#' cause-specific model is an ordinary Cox model with the competing causes
#' treated as censoring; it is fitted by [fit_cox()] on exactly the same
#' complete cases and with the same reference levels, which is why it is
#' fitted here rather than left to the caller -- rebuilding the same
#' complete-case set by hand is where the two stop being comparable. Going
#' through [fit_cox()] also brings the [survival::cox.zph()]
#' proportional-hazards test with it, which the project rules require of every
#' Cox model. Pass `cause_specific = "skip"` to leave it out; any message
#' raised while fitting it names [fit_cox()], because that is the function
#' that raised it.
#'
#' The two models answer different questions and do not have to agree.
#' `attr(x, "fg_vs_csh")` puts them side by side and flags the terms where
#' both are significant and the directions differ; that flag marks a term to
#' explain, not a bug to fix.
#'
#' # The outcome cross-check
#'
#' `event_os` (0 alive or censored, 1 dead of any cause) and `event_cr` (the
#' censoring code, or a cause of death) are two encodings of the same fact and
#' must agree, in both directions:
#'
#' * every censored row of `event_os` must carry the censoring code in
#'   `event_cr`, and every row carrying the censoring code must be censored in
#'   `event_os`;
#' * every death in `event_os` must carry a cause in `event_cr`, and every row
#'   carrying a cause must be a death in `event_os`.
#'
#' A disagreement stops the call and the offending counts are listed. The
#' cross-tabulation itself, a tidy version of it and the four checks come back
#' in `attr(x, "outcome_check")`. Pass `event_os = NULL` to skip the
#' cross-check; that is a deliberate act and is recorded in the result.
#'
#' # Reproducibility
#'
#' Neither [cmprsk::cuminc()] nor [cmprsk::crr()] uses random numbers:
#' `cuminc()` is a closed-form estimator and `crr()` is a Newton-Raphson
#' iteration started from zero. Nothing in this function draws from the random
#' number generator, so the result does not depend on the seed and the state
#' of the generator is left untouched.
#'
#' # Input checks
#'
#' The call stops, with the offending columns, values or counts listed, when
#'
#' 1. a required column is not in `data`, or appears in `data` more than once;
#' 2. `time` is not numeric, or holds missing, negative or infinite values;
#' 3. `event_cr` holds a value that is neither `censor_code` nor one of
#'    `causes`, missing values included. The offending values are listed with
#'    their counts;
#' 4. a cause in `causes` never occurs in `data`. `cuminc()` drops such a
#'    cause from its output altogether and a Fine-Gray model for it has
#'    nothing to fit;
#' 5. `event_os` and `event_cr` disagree in either direction;
#' 6. `group` has missing values and `group_missing = "error"` (the default),
#'    or has fewer than two levels with rows in them;
#' 7. a covariate is of a class that cannot enter the model, is constant, is
#'    an ordered factor, or is missing in every row;
#' 8. casewise deletion empties a level of a covariate, or leaves fewer events
#'    of some cause than the model has coefficients;
#' 9. the design matrix does not encode the levels it is supposed to encode,
#'    or `crr()` does not converge;
#' 10. an element of `times` is past the largest follow-up time of some group
#'     and `times_beyond = "error"` (the default), or `times` holds
#'     duplicated, missing, negative or infinite values.
#'
#' @param data A data frame (or tibble) with one row per subject.
#' @param time Name of the follow-up-time column, as a single string. Must be
#'   numeric. Factors and character digits are refused: use [prep_surv()],
#'   which converts them and returns a numeric `time_mo`.
#' @param event_cr Name of the competing-risk outcome column, as a single
#'   string. Must hold `censor_code` for a censored subject and one of
#'   `causes` for a death, and nothing else.
#' @param event_os Name of the overall-survival event indicator (`0`
#'   censored, `1` dead of any cause), as a single string, used only for the
#'   cross-check against `event_cr`. `NULL` skips the cross-check; the column
#'   is not used for anything else. Defaults to `"event_os"`, so the
#'   cross-check happens unless it is switched off on purpose.
#' @param group Name of the grouping variable for the cumulative incidence
#'   functions and Gray's test, as a single string, or `NULL` (the default)
#'   for one set of curves over the whole cohort. To cross two variables,
#'   build the interaction first, e.g.
#'   `data$site_age <- interaction(data$site, data$age_grp, sep = " / ")`.
#'   `group` does not enter the Fine-Gray model unless it is also named in
#'   `covariates`.
#' @param covariates Names of the Fine-Gray covariates, as a character vector,
#'   in the order they should appear in the result, or `NULL` (the default) to
#'   estimate cumulative incidence only and fit no model. A factor enters as
#'   its levels against the first level; a numeric column enters as one term
#'   per unit.
#' @param causes The cause codes used in `event_cr`, as a numeric vector of at
#'   least two whole numbers, none of them equal to `censor_code`. Defaults to
#'   `c(1, 2)`. Every cause listed must occur in the data, and every death in
#'   the data must carry one of them.
#' @param cause_labels Labels for `causes`, filling the `cause_label` column:
#'   a character vector as long as `causes` and in the same order, or one
#'   named by the cause codes. Defaults to `"Cause 1"`, `"Cause 2"`, and so
#'   on.
#' @param censor_code The value of `event_cr` that means censored. Defaults to
#'   `0`.
#' @param times Time points, in the unit of `time`, at which the cumulative
#'   incidence, its confidence interval and the number at risk are reported.
#'   Defaults to 12, 36 and 60.
#' @param times_beyond What to do with a time point past the largest follow-up
#'   time of a group: `"error"` (the default) to stop, or `"na"` to keep the
#'   row, mark it `estimable = FALSE` and leave the estimate and its interval
#'   missing. There is no setting that extrapolates.
#' @param group_missing What to do with rows whose `group` value is missing:
#'   `"error"` (the default) to stop, `"level"` to give them their own group
#'   labelled `"(Missing)"`, or `"drop"` to leave them out. The number of such
#'   rows is reported in the diagnostics in all three cases.
#' @param missing_text Label used for the group of rows with a missing `group`
#'   value when `group_missing = "level"`.
#' @param labels Labels for the covariates, filling the `var_label` column: a
#'   named character vector or named list, `covariate = "label"`. Covariates
#'   not named keep their column name.
#' @param na_action What to do with rows that have a missing covariate value:
#'   `"omit"` (the default) to fit the models on the complete cases, or
#'   `"fail"` to stop. The counts are reported either way, and the cumulative
#'   incidence is estimated on every usable row in both cases.
#' @param empty_levels What to do with a factor level of a covariate that no
#'   row has: `"drop"` (the default) to drop it with a warning and a count, or
#'   `"error"` to stop. A level emptied by casewise deletion is an error
#'   whatever this is set to.
#' @param cif_complete_case Whether to estimate the cumulative incidence a
#'   second time on the complete cases, so that casewise deletion can be shown
#'   to have changed the cause structure or not to have changed it. Ignored
#'   when `covariates` is `NULL`. Defaults to `TRUE`.
#' @param cause_specific Whether to fit the cause-specific Cox model beside
#'   each Fine-Gray model: `"fit"` (the default) or `"skip"`. Ignored when
#'   `covariates` is `NULL`.
#' @param conf_level Confidence level of every interval returned, and, through
#'   `1 - conf_level`, the level at which `attr(x, "fg_vs_csh")` calls a term
#'   significant. Defaults to `0.95`.
#' @param epv_warn Smallest number of events of a cause per model coefficient
#'   that does not draw a warning. Fewer events than coefficients is an error
#'   whatever this is set to. Set to `0` to silence the warning. Defaults to
#'   `10`.
#' @param digits_hr Decimal places used when a hazard ratio and its interval
#'   are pasted into `shr_txt` and `hr_txt`. The numeric columns are returned
#'   unrounded, so rounding stays a presentation choice and never a stored
#'   one.
#' @param digits_p Decimal places used for the formatted p-value columns; at
#'   the default a p-value below 0.001 is written `"<0.001"`.
#'
#' @return
#' A plain data frame of cumulative incidence estimates, one row per analysis
#' set, cause, group and time point, with the columns
#'
#' * `grouping` -- the name of the grouping variable, or `"overall"`;
#'   `analysis_set` -- `"full"` or `"complete_case"`; `group` -- the level, or
#'   `"ALL"`;
#' * `cause`, `cause_label` -- the cause code and its label;
#' * `time` -- the time point;
#' * `n_group`, `n_risk`, `n_event_cause`, `n_event_cause_by_time` -- subjects
#'   in the group, subjects still at risk at `time`, events of this cause in
#'   the group at any time, and events of this cause by `time`;
#' * `cif`, `cif_se`, `cif_lcl`, `cif_ucl` -- the cumulative incidence, its
#'   standard error and its confidence interval;
#' * `ci_method`, `estimable` -- how the interval was built, or why there is
#'   none;
#' * `gray_stat`, `gray_df`, `gray_p`, `gray_p_fmt` -- Gray's test for that
#'   cause across all groups, repeated on every row of that cause, `NA` when
#'   there is a single group.
#'
#' Six attributes are attached:
#'
#' * `finegray` -- the Fine-Gray results: one row per cause and covariate
#'   level, reference levels included, with `shr`, `ci_low`, `ci_high`,
#'   `shr_txt`, `coef`, `se`, `z`, `p_value`, `p_fmt`, the subjects and events
#'   at that level, and the model-level `converged`, `loglik`, `n_model`,
#'   `n_full` and `n_dropped`. `NULL` when `covariates` is `NULL`;
#' * `cause_specific` -- the same rows for the cause-specific Cox models, with
#'   `hr` and the [survival::cox.zph()] columns that [fit_cox()] returns.
#'   `NULL` when not fitted;
#' * `fg_vs_csh` -- the two side by side for every non-reference term, with
#'   `both_sig`, `same_side` and `direction_flip`. `NULL` when not fitted;
#' * `dropna` -- casewise deletion by level of `group` plus a total row
#'   (`is_total`), with `n`, `n_dropped`, `pct_dropped`, `n_dropped_events`
#'   and one `n_dropped_cause<k>` column per cause. `NULL` when `covariates`
#'   is `NULL`;
#' * `outcome_check` -- `table`, the `table(event_os, event_cr)`
#'   cross-tabulation; `counts`, the same as a data frame; `checks`, the four
#'   agreement invariants; and `note`;
#' * `fit_competing` -- a list with `counts`, `follow_up`, `gray`, `missing`,
#'   `checks`, `cause_specific_zph`, `settings` and `call`.
#'
#' Attributes are dropped by most data-frame verbs, so read them off the
#' object returned by `fit_competing()` before piping it further.
#'
#' @seealso [fit_cox()], which fits the cause-specific models here and is the
#'   right function for an overall-survival Cox model; [fit_km()] and
#'   [calc_rmst()] for the all-cause view of the same cohort; [prep_surv()],
#'   which builds `time_mo`, `event_os` and `event_cr` from vital-status
#'   columns.
#'
#' @examples
#' ## A small simulated cohort. No real patient records are used anywhere in
#' ## this package.
#' set.seed(20260901)
#' n    <- 600
#' site <- factor(rep(c("Stomach", "Small intestine", "Colorectal"),
#'                    length.out = n),
#'                levels = c("Stomach", "Small intestine", "Colorectal"))
#' age  <- factor(sample(c("<65", "65+"), n, TRUE), levels = c("<65", "65+"))
#' sex  <- factor(sample(c("Male", "Female"), n, TRUE),
#'                levels = c("Male", "Female"))
#' r1   <- c(Stomach = 0.006, "Small intestine" = 0.011,
#'           Colorectal = 0.016)[as.character(site)]
#' r2   <- ifelse(age == "65+", 0.012, 0.003)
#' t1   <- stats::rexp(n, r1)
#' t2   <- stats::rexp(n, r2)
#' tc   <- stats::runif(n, 6, 130)
#' tm   <- pmin(t1, t2, tc)
#' toy  <- data.frame(
#'   site     = site,
#'   age      = age,
#'   sex      = sex,
#'   time_mo  = round(tm, 1),
#'   event_os = as.integer(tm < tc),
#'   event_cr = ifelse(tm == tc, 0L, ifelse(t1 <= t2, 1L, 2L))
#' )
#' ## a covariate with missing values, to exercise casewise deletion
#' toy$stage <- factor(sample(c("Localised", "Advanced"), n, TRUE),
#'                     levels = c("Localised", "Advanced"))
#' toy$stage[sample(n, 60)] <- NA
#'
#' ## Cumulative incidence only, one set of curves for the whole cohort.
#' ci <- fit_competing(toy, cause_labels = c("Death from GIST",
#'                                           "Death from other causes"))
#' ci[, c("cause_label", "time", "cif", "cif_lcl", "cif_ucl", "n_risk")]
#'
#' ## The outcome cross-check is returned, not printed.
#' attr(ci, "outcome_check")$table
#' attr(ci, "outcome_check")$checks
#'
#' ## By primary site, with Gray's test repeated on every row of a cause.
#' cif <- fit_competing(toy, group = "site",
#'                      cause_labels = c("Death from GIST",
#'                                       "Death from other causes"))
#' cif[cif$time == 60, c("group", "cause_label", "cif", "cif_lcl", "cif_ucl",
#'                       "gray_p_fmt")]
#'
#' ## With covariates: Fine-Gray models, the cause-specific models beside
#' ## them, and the rows casewise deletion costs.
#' fg <- fit_competing(toy, group = "site",
#'                     covariates = c("site", "age", "sex", "stage"),
#'                     cause_labels = c("Death from GIST",
#'                                      "Death from other causes"),
#'                     labels = c(site = "Primary site", age = "Age group"))
#' attr(fg, "finegray")[, c("cause_label", "variable", "level", "n",
#'                          "n_event_cause", "shr_txt", "p_fmt")]
#' attr(fg, "dropna")
#'
#' ## The two hazards side by side; a flagged term is one to explain.
#' attr(fg, "fg_vs_csh")[, c("cause_label", "term", "shr_txt", "hr_txt",
#'                           "both_sig", "same_side", "direction_flip")]
#'
#' ## Diagnostics travel with the result instead of being printed.
#' attr(fg, "fit_competing")$counts
#' attr(fg, "fit_competing")$follow_up
#' attr(fg, "fit_competing")$gray
#' attr(fg, "fit_competing")$checks
#'
#' ## A cell with no events of a cause has no interval, and says so rather
#' ## than reporting 0.00 (0.00-0.00).
#' recoded <- toy
#' recoded$event_cr[recoded$site == "Colorectal" & recoded$event_cr == 2L] <- 1L
#' z <- fit_competing(recoded, group = "site", times = c(12, 60))
#' z[z$cause == 2 & z$group == "Colorectal",
#'   c("group", "time", "n_risk", "n_event_cause_by_time", "cif", "cif_lcl",
#'     "ci_method")]
#'
#' ## A time point past the end of follow-up is refused, not extrapolated.
#' try(fit_competing(toy, group = "site", times = c(12, 60, 200)))
#'
#' ## An impossible cause code is listed, not silently read as censored.
#' bad <- toy
#' bad$event_cr[c(1, 5, 9)] <- 3L
#' try(fit_competing(bad, group = "site"))
#'
#' ## So is a disagreement between the two outcome columns.
#' bad2 <- toy
#' bad2$event_os[bad2$event_cr == 2L][1:4] <- 0L
#' try(fit_competing(bad2, group = "site"))
#'
#' @export
fit_competing <- function(data,
                          time              = "time_mo",
                          event_cr          = "event_cr",
                          event_os          = "event_os",
                          group             = NULL,
                          covariates        = NULL,
                          causes            = c(1, 2),
                          cause_labels      = NULL,
                          censor_code       = 0,
                          times             = c(12, 36, 60),
                          times_beyond      = c("error", "na"),
                          group_missing     = c("error", "level", "drop"),
                          missing_text      = "Missing",
                          labels            = NULL,
                          na_action         = c("omit", "fail"),
                          empty_levels      = c("drop", "error"),
                          cif_complete_case = TRUE,
                          cause_specific    = c("fit", "skip"),
                          conf_level        = 0.95,
                          epv_warn          = 10,
                          digits_hr         = 2L,
                          digits_p          = 3L) {

  cl <- match.call()

  ## -- 0. data ---------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("fit_competing(): `data` must be a data frame, got an object of ",
         "class ", .ps_values(class(data)), ".", call. = FALSE)
  }
  n_in <- nrow(data)
  if (n_in == 0L) {
    stop("fit_competing(): `data` has 0 rows, there is nothing to estimate.",
         call. = FALSE)
  }

  ## -- 1. the argument values themselves --------------------------------------
  .cr_name_arg(time,     "time")
  .cr_name_arg(event_cr, "event_cr")
  if (!is.null(event_os)) .cr_name_arg(event_os, "event_os")
  if (!is.null(group))    .cr_name_arg(group,    "group")
  .cr_name_arg(missing_text, "missing_text")
  times_beyond   <- match.arg(times_beyond)
  group_missing  <- match.arg(group_missing)
  na_action      <- match.arg(na_action)
  empty_levels   <- match.arg(empty_levels)
  cause_specific <- match.arg(cause_specific)

  if (!is.logical(cif_complete_case) || length(cif_complete_case) != 1L ||
      is.na(cif_complete_case)) {
    stop("fit_competing(): `cif_complete_case` must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop(sprintf(
      paste0("fit_competing(): `conf_level` must be a single number strictly ",
             "between 0 and 1, got %s."),
      .ps_values(conf_level)),
      call. = FALSE)
  }
  if (!is.numeric(epv_warn) || length(epv_warn) != 1L || is.na(epv_warn) ||
      epv_warn < 0) {
    stop(sprintf(
      paste0("fit_competing(): `epv_warn` must be a single number >= 0 (0 ",
             "silences the warning), got %s."),
      .ps_values(epv_warn)),
      call. = FALSE)
  }
  .cr_count_arg(digits_hr, "digits_hr", min = 0L)
  .cr_count_arg(digits_p,  "digits_p",  min = 1L)

  z_q   <- stats::qnorm((1 + conf_level) / 2)
  alpha <- 1 - conf_level
  times <- .cr_times_arg(times)

  ## censoring code and causes
  if (!is.numeric(censor_code) || length(censor_code) != 1L ||
      is.na(censor_code) || !is.finite(censor_code) ||
      censor_code != as.integer(censor_code)) {
    stop(sprintf(
      paste0("fit_competing(): `censor_code` must be a single whole number, ",
             "got %s."),
      .ps_values(censor_code)),
      call. = FALSE)
  }
  censor_code <- as.integer(censor_code)
  if (!is.numeric(causes) || length(causes) < 2L || anyNA(causes) ||
      any(!is.finite(causes)) || any(causes != as.integer(causes))) {
    stop(sprintf(
      paste0("fit_competing(): `causes` must be a numeric vector of at least ",
             "two whole numbers, got %s of length %d.\n",
             "  With a single cause there is no competing risk: use fit_km() ",
             "and fit_cox() on that one event instead."),
      .ps_values(class(causes)), length(causes)),
      call. = FALSE)
  }
  causes <- as.integer(causes)
  if (anyDuplicated(causes) != 0L) {
    stop(sprintf(
      "fit_competing(): `causes` lists %s more than once.",
      .ps_values(unique(causes[duplicated(causes)]))),
      call. = FALSE)
  }
  if (censor_code %in% causes) {
    stop(sprintf(
      paste0("fit_competing(): `censor_code` (%d) is also listed in `causes` ",
             "(%s). A code cannot mean both censored and a cause of failure."),
      censor_code, .ps_values(causes)),
      call. = FALSE)
  }
  n_cause <- length(causes)

  cause_labels <- .cr_cause_labels(cause_labels, causes)

  ## -- 2. required columns present, and present once --------------------------
  needed <- c(time = time, event_cr = event_cr)
  if (!is.null(event_os)) needed <- c(needed, event_os = event_os)
  if (!is.null(group))    needed <- c(needed, group = group)
  if (!is.null(covariates)) {
    if (!is.character(covariates) || length(covariates) == 0L ||
        anyNA(covariates) || any(!nzchar(covariates))) {
      stop(sprintf(
        paste0("fit_competing(): `covariates` must be a character vector of ",
               "at least one non-missing column name, or NULL to fit no ",
               "model, got %s of length %d."),
        .ps_values(class(covariates)), length(covariates)),
        call. = FALSE)
    }
    if (anyDuplicated(covariates) != 0L) {
      stop(sprintf(
        paste0("fit_competing(): `covariates` names %s more than once. The ",
               "same column twice makes the design matrix singular."),
        .ps_values(unique(covariates[duplicated(covariates)]))),
        call. = FALSE)
    }
    needed <- c(needed, stats::setNames(covariates,
                                        paste0("covariate", seq_along(covariates))))
  }

  miss <- needed[!(needed %in% names(data))]
  if (length(miss) > 0L) {
    hint <- if (!is.null(event_os) && event_os %in% unname(miss)) {
      paste0("\n  `event_os` is used only to cross-check `event_cr`. Pass ",
             "event_os = NULL to skip that check on purpose; it is not ",
             "skipped by accident.")
    } else {
      ""
    }
    stop(sprintf(
      paste0("fit_competing(): %d required column%s not found in `data`:\n",
             "  %s\n  `data` has %d columns. Check the spelling, or pass the ",
             "column names explicitly.%s"),
      length(miss), if (length(miss) > 1L) "s" else "",
      paste0(sprintf('%s = "%s"', names(miss), unname(miss)), collapse = "\n  "),
      ncol(data), hint),
      call. = FALSE)
  }
  dup <- needed[needed %in% names(data)[duplicated(names(data))]]
  if (length(dup) > 0L) {
    stop(sprintf(
      paste0("fit_competing(): %s appear%s more than once among the columns ",
             "of `data`; fit_competing() would silently use the first one. ",
             "Make the column names unique first."),
      paste0(sprintf('%s = "%s"', names(dup), unname(dup)), collapse = ", "),
      if (length(dup) > 1L) "" else "s"),
      call. = FALSE)
  }
  ## internal column names must not collide with the caller's
  reserved <- c(".cr_time", ".cr_cause", ".cr_event", ".cr_group")
  clash    <- intersect(reserved, unname(needed))
  if (length(clash) > 0L) {
    stop(sprintf(
      paste0("fit_competing(): %s %s reserved for the internal columns of ",
             "fit_competing(). Rename the column%s before calling."),
      .ps_values(clash), if (length(clash) > 1L) "are" else "is",
      if (length(clash) > 1L) "s" else ""),
      call. = FALSE)
  }

  ## -- 3. follow-up time ------------------------------------------------------
  tv <- data[[time]]
  if (is.factor(tv) || is.character(tv)) {
    stop(sprintf(
      paste0('fit_competing(): the `time` column "%s" is %s. Factor levels ',
             "and character digits are not numbers, so fit_competing() will ",
             "not convert them for you: convert the column explicitly, or ",
             "build it with prep_surv(), which returns a numeric `time_mo`."),
      time, if (is.factor(tv)) "a factor" else "character"),
      call. = FALSE)
  }
  if (!is.numeric(tv)) {
    stop(sprintf(
      paste0('fit_competing(): the `time` column "%s" must be numeric, got an ',
             "object of class %s."),
      time, .ps_values(class(tv))),
      call. = FALSE)
  }
  tnum   <- as.numeric(tv)
  n_t_na <- sum(is.na(tnum))
  if (n_t_na > 0L) {
    stop(sprintf(
      paste0('fit_competing(): the `time` column "%s" has %d missing value%s ',
             "out of %d rows. fit_competing() neither drops nor imputes them; ",
             "decide what to do with these rows before calling it."),
      time, n_t_na, if (n_t_na > 1L) "s" else "", n_in),
      call. = FALSE)
  }
  n_t_inf <- sum(!is.finite(tnum))
  if (n_t_inf > 0L) {
    stop(sprintf(
      paste0('fit_competing(): the `time` column "%s" has %d infinite ',
             "value%s. A cumulative incidence cannot be estimated at an ",
             "infinite follow-up time."),
      time, n_t_inf, if (n_t_inf > 1L) "s" else ""),
      call. = FALSE)
  }
  if (any(tnum < 0)) {
    n_neg <- sum(tnum < 0)
    stop(sprintf(
      paste0('fit_competing(): the `time` column "%s" has %d negative ',
             "value%s (minimum %s). A negative follow-up time is not ",
             "analysable, so fit_competing() stops instead of computing on ",
             "it."),
      time, n_neg, if (n_neg > 1L) "s" else "", format(min(tnum))),
      call. = FALSE)
  }

  ## -- 4. the competing-risk outcome column -----------------------------------
  crr_raw <- data[[event_cr]]
  if (is.factor(crr_raw) || is.character(crr_raw)) {
    stop(sprintf(
      paste0('fit_competing(): the `event_cr` column "%s" is %s. ',
             "fit_competing() will not guess which level is the censoring ",
             "code and which levels are causes: recode it to the numbers in ",
             "`censor_code` and `causes` first.\n  values present : %s"),
      event_cr, if (is.factor(crr_raw)) "a factor" else "character",
      .ps_counts(crr_raw)),
      call. = FALSE)
  }
  if (is.logical(crr_raw)) {
    stop(sprintf(
      paste0('fit_competing(): the `event_cr` column "%s" is logical. A ',
             "competing-risk outcome has at least three states (censored and ",
             "two causes) and cannot be TRUE/FALSE."),
      event_cr),
      call. = FALSE)
  }
  if (!is.numeric(crr_raw)) {
    stop(sprintf(
      paste0('fit_competing(): the `event_cr` column "%s" must be numeric, ',
             "got an object of class %s."),
      event_cr, .ps_values(class(crr_raw))),
      call. = FALSE)
  }
  allowed <- c(censor_code, causes)
  bad_cr  <- is.na(crr_raw) | !(crr_raw %in% allowed)
  if (any(bad_cr)) {
    stop(sprintf(
      paste0('fit_competing(): the `event_cr` column "%s" has %d of %d rows ',
             "whose value is neither the censoring code (%d) nor one of the ",
             "causes (%s).\n  offending values : %s\n",
             "  A missing value counts as an offending value: ",
             "fit_competing() never reads an unknown code as censored, ",
             "because that would move those subjects out of the numerator ",
             "and leave them in the denominator. Recode them, or list the ",
             "code in `causes`."),
      event_cr, sum(bad_cr), n_in, censor_code, .ps_values(causes),
      .ps_counts(crr_raw[bad_cr])),
      call. = FALSE)
  }
  cr <- as.integer(crr_raw)

  n_by_cause <- vapply(causes, function(k) sum(cr == k), integer(1))
  if (any(n_by_cause == 0L)) {
    stop(sprintf(
      paste0("fit_competing(): %d of the %d cause%s in `causes` never occur%s ",
             "in `data` (%s).\n",
             "  cmprsk::cuminc() leaves a cause with no events out of its ",
             "output altogether, so a table built on `causes` would silently ",
             "lose a row, and cmprsk::crr() would have nothing to fit for ",
             "it.\n",
             "  Drop the cause from `causes`, or check that `data` was not ",
             "already filtered."),
      sum(n_by_cause == 0L), n_cause, if (n_cause > 1L) "s" else "",
      if (sum(n_by_cause == 0L) > 1L) "" else "s",
      .ps_values(causes[n_by_cause == 0L])),
      call. = FALSE)
  }

  ## -- 5. the outcome cross-check ---------------------------------------------
  ## The project rule is that table(event_os, event_cr) be cross-checked. A
  ## package function cannot discharge that by printing: it enforces it and
  ## returns the table.
  if (is.null(event_os)) {
    os_tab    <- NULL
    os_counts <- NULL
    os_checks <- c(`event_os and event_cr agree in both directions` = NA)
    os_note   <- paste0("not performed: event_os = NULL was passed, so the ",
                        "competing-risk outcome was not cross-checked ",
                        "against an overall-survival indicator")
  } else {
    ev_raw <- data[[event_os]]
    if (is.factor(ev_raw) || is.character(ev_raw)) {
      stop(sprintf(
        paste0('fit_competing(): the `event_os` column "%s" is %s. ',
               "fit_competing() will not guess which level means a death: ",
               "convert it to 0 (censored) / 1 (dead) explicitly, or build it ",
               "with prep_surv().\n  values present : %s"),
        event_os, if (is.factor(ev_raw)) "a factor" else "character",
        .ps_counts(ev_raw)),
        call. = FALSE)
    }
    if (is.logical(ev_raw)) {
      ev_os <- as.integer(ev_raw)
    } else if (is.numeric(ev_raw)) {
      ev_os <- as.numeric(ev_raw)
    } else {
      stop(sprintf(
        paste0('fit_competing(): the `event_os` column "%s" must be numeric ',
               "0/1 or logical, got an object of class %s."),
        event_os, .ps_values(class(ev_raw))),
        call. = FALSE)
    }
    bad_os <- is.na(ev_os) | !(ev_os %in% c(0, 1))
    if (any(bad_os)) {
      stop(sprintf(
        paste0('fit_competing(): the `event_os` column "%s" has %d of %d rows ',
               "whose value is neither 0 (censored) nor 1 (dead).\n",
               "  offending values : %s"),
        event_os, sum(bad_os), n_in, .ps_counts(ev_raw[bad_os])),
        call. = FALSE)
    }
    ev_os <- as.integer(ev_os)

    n_a <- sum(ev_os == 0L & cr != censor_code)
    n_b <- sum(ev_os == 1L & cr == censor_code)
    if (n_a > 0L || n_b > 0L) {
      lines <- character(0)
      if (n_a > 0L) {
        lines <- c(lines, sprintf(
          paste0("%d row%s censored in \"%s\" (0) but carrying a cause of ",
                 "death in \"%s\" (%s)"),
          n_a, if (n_a > 1L) "s are" else " is", event_os, event_cr,
          .ps_counts(cr[ev_os == 0L & cr != censor_code])))
      }
      if (n_b > 0L) {
        lines <- c(lines, sprintf(
          paste0("%d row%s a death in \"%s\" (1) but carrying the censoring ",
                 "code in \"%s\" (%d)"),
          n_b, if (n_b > 1L) "s are" else " is", event_os, event_cr,
          censor_code))
      }
      stop(sprintf(
        paste0("fit_competing(): \"%s\" and \"%s\" do not agree:\n  %s\n",
               "  They are two encodings of the same fact and a competing-",
               "risk analysis built on one of them contradicts every ",
               "all-cause number built on the other. Fix the recoding that ",
               "produced them; do not pass event_os = NULL to get past this ",
               "message."),
        event_os, event_cr, paste(lines, collapse = "\n  ")),
        call. = FALSE)
    }

    os_tab <- table(ev_os, cr, dnn = c(event_os, event_cr))
    os_counts <- data.frame(
      event_os         = as.integer(rep(rownames(os_tab), times = ncol(os_tab))),
      event_cr         = as.integer(rep(colnames(os_tab), each  = nrow(os_tab))),
      n                = as.integer(os_tab),
      stringsAsFactors = FALSE,
      row.names        = NULL)
    names(os_counts)[1:2] <- c(event_os, event_cr)
    os_checks <- c(
      `censored in event_os implies the censoring code in event_cr` = n_a == 0L,
      `a death in event_os implies a cause in event_cr`             = n_b == 0L,
      `the censoring code in event_cr implies censored in event_os` = n_b == 0L,
      `a cause in event_cr implies a death in event_os`             = n_a == 0L)
    os_note <- sprintf(
      paste0("table(%s, %s) on all %d rows: %d censored, %d death%s (%s). ",
             "Off-diagonal cells are 0 in both directions."),
      event_os, event_cr, n_in, sum(cr == censor_code), sum(ev_os == 1L),
      if (sum(ev_os == 1L) == 1L) "" else "s",
      paste(sprintf("cause %d: %d", causes, n_by_cause), collapse = ", "))
  }

  ## -- 6. the grouping variable ------------------------------------------------
  n_g_na    <- 0L
  n_g_drop  <- 0L
  n_lv_drop <- 0L
  g         <- NULL

  if (!is.null(group)) {
    g_raw <- data[[group]]
    if (is.factor(g_raw)) {
      g <- g_raw
    } else if (is.character(g_raw)) {
      g <- factor(g_raw)
    } else if (is.logical(g_raw)) {
      g <- factor(g_raw, levels = c(FALSE, TRUE))
    } else if (is.numeric(g_raw)) {
      g <- factor(g_raw, levels = sort(unique(g_raw[!is.na(g_raw)])))
    } else {
      stop(sprintf(
        paste0('fit_competing(): the `group` column "%s" is of class %s, ',
               "which cannot be used as a grouping variable. Convert it to a ",
               "factor first."),
        group, .ps_values(class(g_raw))),
        call. = FALSE)
    }

    n_g_na   <- sum(is.na(g))
    n_lv_all <- nlevels(g)
    g        <- droplevels(g)
    n_lv_drop <- n_lv_all - nlevels(g)

    if (n_g_na > 0L) {
      if (group_missing == "error") {
        stop(sprintf(
          paste0('fit_competing(): the `group` column "%s" has %d missing ',
                 "value%s out of %d rows.\n",
                 "  fit_competing() will not decide on its own what those ",
                 "rows are: they are neither a group nor nothing, and ",
                 "dropping them would change the denominator of every ",
                 "cumulative incidence it returns.\n",
                 '  Pass group_missing = "level" to give them their own group ',
                 'labelled "(%s)", or group_missing = "drop" to leave them ',
                 "out (either way the count comes back in the diagnostics), ",
                 "or filter them out before calling."),
          group, n_g_na, if (n_g_na > 1L) "s" else "", n_in, missing_text),
          call. = FALSE)
      } else if (group_missing == "level") {
        na_lab <- paste0("(", missing_text, ")")
        if (na_lab %in% levels(g)) {
          stop(sprintf(
            paste0('fit_competing(): group_missing = "level" would add a ',
                   'group labelled "%s", but the `group` column "%s" already ',
                   "has a level with that name. Rename the level, or change ",
                   "`missing_text`."),
            na_lab, group),
            call. = FALSE)
        }
        g <- factor(ifelse(is.na(g), na_lab, as.character(g)),
                    levels = c(levels(g), na_lab))
      } else {
        n_g_drop <- n_g_na
      }
    }
  }

  keep <- if (is.null(g)) rep(TRUE, n_in) else !is.na(g)
  if (!any(keep)) {
    stop(sprintf(
      paste0('fit_competing(): every row has a missing `group` value ("%s"), ',
             "so there is nothing left to estimate."),
      group),
      call. = FALSE)
  }
  row_id <- which(keep)
  tnum   <- tnum[keep]
  cr     <- cr[keep]
  g      <- if (is.null(g)) factor(rep("ALL", sum(keep))) else droplevels(g[keep])
  n_used <- length(tnum)
  lv     <- levels(g)

  if (!is.null(group) && length(lv) < 2L) {
    stop(sprintf(
      paste0('fit_competing(): the `group` column "%s" has %d level%s with ',
             "rows in it (%s) in the %d rows used; at least 2 are needed.\n",
             "  With one group there is nothing to compare: Gray's test is ",
             "undefined. Call fit_competing() with group = NULL if a single ",
             "set of curves is what you want.\n",
             "  Check that the data were not already filtered down to one ",
             "group, and that an empty group is not hiding behind an unused ",
             "factor level (%d unused level%s dropped here)."),
      group, length(lv), if (length(lv) == 1L) "" else "s",
      if (length(lv) == 0L) "none" else .ps_values(lv), n_used,
      n_lv_drop, if (n_lv_drop == 1L) " was" else "s were"),
      call. = FALSE)
  }

  ## every cause must still occur after the group filter
  n_by_cause_used <- vapply(causes, function(k) sum(cr == k), integer(1))
  if (any(n_by_cause_used == 0L)) {
    stop(sprintf(
      paste0("fit_competing(): after dropping %d row%s with a missing `group` ",
             "value, cause%s %s no longer occur%s in the data. cmprsk::cuminc",
             "() would leave %s out of its output. Use group_missing = ",
             '"level" instead of "drop", or drop the cause from `causes`.'),
      n_g_drop, if (n_g_drop == 1L) "" else "s",
      if (sum(n_by_cause_used == 0L) > 1L) "s" else "",
      .ps_values(causes[n_by_cause_used == 0L]),
      if (sum(n_by_cause_used == 0L) > 1L) "" else "s",
      if (sum(n_by_cause_used == 0L) > 1L) "them" else "it"),
      call. = FALSE)
  }

  ## -- 7. the follow-up guard --------------------------------------------------
  max_fu    <- vapply(split(tnum, g), max, numeric(1))
  tau_upper <- min(max_fu)
  beyond    <- outer(times, unname(max_fu), ">")   # times x groups
  if (any(beyond) && times_beyond == "error") {
    bad_g <- lv[apply(beyond, 2L, any)]
    lines <- vapply(bad_g, function(s) {
      tt <- times[beyond[, match(s, lv)]]
      sprintf("%s : time point%s %s past a largest follow-up of %s",
              s, if (length(tt) > 1L) "s" else "",
              paste(format(tt, trim = TRUE), collapse = ", "),
              format(unname(max_fu[match(s, lv)])))
    }, character(1))
    stop(sprintf(
      paste0("fit_competing(): %d of %d requested time point%s cannot be ",
             "estimated in %d of %d group%s:\n  %s\n",
             "  A cumulative incidence function is not estimable past the ",
             "largest follow-up time of a group. This is the rule an RMST tau ",
             "guard uses as well.\n",
             "  The largest time point every group can be evaluated at is %s. ",
             'Either keep `times` within it, or pass times_beyond = "na" to ',
             "keep those rows and mark them estimable = FALSE."),
      sum(apply(beyond, 1L, any)), length(times),
      if (length(times) > 1L) "s" else "",
      length(bad_g), length(lv), if (length(lv) > 1L) "s" else "",
      paste(lines, collapse = "\n  "), format(tau_upper)),
      call. = FALSE)
  }

  ## -- 8. cumulative incidence on the full analysis set -------------------------
  cif_full <- .cr_cif_set(tnum = tnum, cr = cr, g = g, causes = causes,
                          cause_labels = cause_labels,
                          censor_code = censor_code, times = times, z_q = z_q,
                          grouping = if (is.null(group)) "overall" else group,
                          analysis_set = "full", single_group = is.null(group))
  gray_full <- attr(cif_full, "gray")

  ## -- 9. the covariates --------------------------------------------------------
  ## Every covariate is copied into an internally named column, so that no
  ## other column of `data` reaches model.matrix(), and so that a column name
  ## that is not a syntactic R name cannot reach the model formula. The
  ## caller's names are put back on the way out.
  has_cov <- !is.null(covariates)
  k       <- if (has_cov) length(covariates) else 0L

  var_label <- covariates
  if (has_cov) {
    var_label <- .cr_labels_arg(labels, covariates)
  }

  xs        <- vector("list", k)
  lev_list  <- vector("list", k)
  is_cat    <- logical(k)
  n_lv_cov  <- integer(k)
  ref_lev   <- rep(NA_character_, k)
  if (k > 0L) names(n_lv_cov) <- names(ref_lev) <- covariates

  for (i in seq_len(k)) {
    v  <- covariates[i]
    xr <- data[[v]][keep]

    if (!is.null(dim(xr))) {
      stop(sprintf(
        paste0('fit_competing(): the covariate "%s" is not a plain column (it ',
               "has dimensions %s). A matrix column cannot be reported level ",
               "by level; split it into columns first."),
        v, paste(dim(xr), collapse = " x ")),
        call. = FALSE)
    }
    if (all(is.na(xr))) {
      stop(sprintf(
        paste0('fit_competing(): the covariate "%s" is missing in all %d rows ',
               "used. Casewise deletion would empty the cohort, and crr() ",
               "would then fail for a reason that has nothing to do with the ",
               "real cause. Drop it from `covariates`, or clean the column ",
               "first."),
        v, n_used),
        call. = FALSE)
    }

    if (is.factor(xr)) {
      if (is.ordered(xr)) {
        stop(sprintf(
          paste0('fit_competing(): the covariate "%s" is an ordered factor. R ',
                 "fits ordered factors with polynomial contrasts, which give ",
                 ".L and .Q terms instead of one subdistribution hazard ratio ",
                 "per level against a reference, so fit_competing() refuses ",
                 "it rather than reporting terms that do not mean what the ",
                 "row labels would say.\n",
                 '  Use data$%s <- factor(data$%s, ordered = FALSE) to fit it ',
                 "as a categorical covariate; the level order, and therefore ",
                 "the reference level, is kept."),
          v, v, v),
          call. = FALSE)
      }
      if (!is.null(attr(xr, "contrasts"))) {
        stop(sprintf(
          paste0('fit_competing(): the covariate "%s" carries a `contrasts` ',
                 "attribute, so its coefficients would not be the levels ",
                 'against the first level. Remove it with attr(data$%s, ',
                 '"contrasts") <- NULL, or recode the column, before calling ',
                 "fit_competing()."),
          v, v),
          call. = FALSE)
      }
      x <- xr
      is_cat[i] <- TRUE
    } else if (is.character(xr)) {
      x <- factor(xr)
      is_cat[i] <- TRUE
    } else if (is.logical(xr)) {
      x <- factor(xr, levels = c(FALSE, TRUE))
      is_cat[i] <- TRUE
    } else if (is.numeric(xr)) {
      x <- as.numeric(xr)
      n_x_inf <- sum(!is.finite(x) & !is.na(x))
      if (n_x_inf > 0L) {
        stop(sprintf(
          paste0('fit_competing(): the covariate "%s" has %d infinite ',
                 "value%s. An infinite covariate makes the linear predictor ",
                 "infinite for that subject; recode those values before ",
                 "calling."),
          v, n_x_inf, if (n_x_inf > 1L) "s" else ""),
          call. = FALSE)
      }
      is_cat[i] <- FALSE
    } else {
      stop(sprintf(
        paste0('fit_competing(): the covariate "%s" is of class %s, which ',
               "cannot enter a Fine-Gray model. Convert it to a factor ",
               "(categorical) or to a number (continuous) first, so that it ",
               "is on the record which of the two it is."),
        v, .ps_values(class(xr))),
        call. = FALSE)
    }

    if (is_cat[i]) {
      cnt   <- table(x, useNA = "no")
      empty <- names(cnt)[cnt == 0L]
      n_ok  <- sum(cnt > 0L)
      if (n_ok < 2L) {
        stop(sprintf(
          paste0('fit_competing(): the covariate "%s" has %d level%s with ',
                 "rows in it (%s); at least 2 are needed.\n",
                 "  A covariate that is the same for everybody explains ",
                 "nothing and has no subdistribution hazard ratio to report: ",
                 "its design column would be constant and its coefficient ",
                 "missing.\n",
                 "  Check that the data were not already filtered down to one ",
                 "group, and that empty groups are not hiding behind unused ",
                 "factor levels (%d declared level%s, %d with no rows: %s)."),
          v, n_ok, if (n_ok == 1L) "" else "s",
          if (n_ok == 0L) "none" else .ps_values(names(cnt)[cnt > 0L]),
          length(cnt), if (length(cnt) == 1L) "" else "s",
          length(empty),
          if (length(empty) == 0L) "none" else .ps_values(empty)),
          call. = FALSE)
      }
      if (length(empty) > 0L) {
        if (names(cnt)[1] %in% empty) {
          stop(sprintf(
            paste0('fit_competing(): the reference level of the covariate ',
                   '"%s" is "%s", and no row has it.\n',
                   '  Dropping it would move the reference to "%s" without ',
                   "saying so, and every subdistribution hazard ratio of this ",
                   "covariate would change meaning while the code stayed the ",
                   "same. fit_competing() will not do that.\n",
                   "  Relevel the factor explicitly, e.g. data$%s <- ",
                   'stats::relevel(droplevels(data$%s), ref = "%s").'),
            v, names(cnt)[1], names(cnt)[cnt > 0L][1], v, v,
            names(cnt)[cnt > 0L][1]),
            call. = FALSE)
        }
        msg <- sprintf(
          paste0('the covariate "%s" has %d declared level%s that no row has ',
                 "(%s). An empty level has an all-zero design column and no ",
                 "subdistribution hazard ratio to report."),
          v, length(empty), if (length(empty) > 1L) "s" else "",
          .ps_values(empty))
        if (empty_levels == "error") {
          stop(sprintf(
            paste0("fit_competing(): %s\n  Drop them with droplevels(), or ",
                   'pass empty_levels = "drop" to have fit_competing() drop ',
                   "them and report how many."),
            msg),
            call. = FALSE)
        }
        warning(sprintf(
          paste0("fit_competing(): %s\n  Dropped, and counted in ",
                 'attr(x, "fit_competing")$settings$levels_dropped. The ',
                 'reference level ("%s") is not affected.'),
          msg, names(cnt)[1]),
          call. = FALSE)
        n_lv_cov[i] <- length(empty)
        x <- droplevels(x)
      }
      lev_list[[i]] <- levels(x)
      ref_lev[i]    <- levels(x)[1]
    } else {
      if (length(unique(x[!is.na(x)])) < 2L) {
        stop(sprintf(
          paste0('fit_competing(): the continuous covariate "%s" takes the ',
                 "same value (%s) in every row that has one. A constant ",
                 "covariate explains nothing and its coefficient would be ",
                 "missing."),
          v, .ps_values(unique(x[!is.na(x)]))),
          call. = FALSE)
      }
    }
    xs[[i]] <- x
  }

  ## -- 10. casewise deletion ----------------------------------------------------
  vn <- if (k > 0L) paste0(".cr_v", seq_len(k), "_") else character(0)
  md <- data.frame(.cr_time = tnum, .cr_cause = cr, stringsAsFactors = FALSE)
  for (i in seq_len(k)) md[[vn[i]]] <- xs[[i]]

  n_miss <- vapply(md, function(z) sum(is.na(z)), integer(1))
  missing_tab <- data.frame(
    variable    = c(time, event_cr, covariates),
    role        = c("time", "event_cr", rep("covariate", k)),
    n_missing   = unname(n_miss),
    pct_missing = 100 * unname(n_miss) / n_used,
    stringsAsFactors = FALSE,
    row.names        = NULL)
  miss_txt <- if (any(missing_tab$n_missing > 0L)) {
    paste(sprintf("%s (n = %d)",
                  missing_tab$variable[missing_tab$n_missing > 0L],
                  missing_tab$n_missing[missing_tab$n_missing > 0L]),
          collapse = ", ")
  } else {
    "none"
  }

  cc         <- stats::complete.cases(md)
  n_cc       <- sum(cc)
  n_dropped  <- n_used - n_cc
  ev_by_cause_cc <- vapply(causes, function(x) sum(cr[cc] == x), integer(1))

  dropna_tab <- NULL
  if (has_cov) {
    dropna_tab <- .cr_dropna_table(g = g, lv = lv, cc = cc, cr = cr,
                                   causes = causes)

    if (n_dropped > 0L && na_action == "fail") {
      stop(sprintf(
        paste0('fit_competing(): na_action = "fail" and %d of %d row%s ',
               "(%.1f%%) have a missing covariate value; %d of the %d death%s ",
               "would go with them.\n  missing by column : %s\n",
               '  Pass na_action = "omit" to fit the models on the %d ',
               "complete case%s (the counts come back in ",
               'attr(x, "dropna") either way), or handle the missing values ',
               "before calling."),
        n_dropped, n_used, if (n_dropped > 1L) "s" else "",
        100 * n_dropped / n_used,
        sum(cr != censor_code) - sum(ev_by_cause_cc), sum(cr != censor_code),
        if (sum(cr != censor_code) == 1L) "" else "s", miss_txt,
        n_cc, if (n_cc == 1L) "" else "s"),
        call. = FALSE)
    }
    if (n_cc == 0L) {
      stop(sprintf(
        paste0("fit_competing(): casewise deletion leaves 0 of %d rows, so ",
               "there is nothing to fit.\n  missing by column : %s"),
        n_used, miss_txt),
        call. = FALSE)
    }
    if (any(ev_by_cause_cc == 0L)) {
      stop(sprintf(
        paste0("fit_competing(): after casewise deletion cause%s %s ha%s no ",
               "events left in the %d complete case%s (there were %s in the ",
               "%d rows used).\n",
               "  The %d row%s deleted for missing covariate values took all ",
               "of them; cmprsk::crr() cannot fit a model for a cause with no ",
               "events.\n  missing by column : %s"),
        if (sum(ev_by_cause_cc == 0L) > 1L) "s" else "",
        .ps_values(causes[ev_by_cause_cc == 0L]),
        if (sum(ev_by_cause_cc == 0L) > 1L) "ve" else "s",
        n_cc, if (n_cc == 1L) "" else "s",
        paste(n_by_cause_used[ev_by_cause_cc == 0L], collapse = ", "), n_used,
        n_dropped, if (n_dropped == 1L) "" else "s", miss_txt),
        call. = FALSE)
    }
    md <- md[cc, , drop = FALSE]
    row.names(md) <- NULL
  }

  ## -- 11. the complete cases must still support the models ---------------------
  n_coef      <- 0L
  no_event_lv <- character(0)
  if (has_cov) {
    n_coef <- sum(vapply(seq_len(k), function(i) {
      if (is_cat[i]) length(lev_list[[i]]) - 1L else 1L
    }, integer(1)))

    for (i in seq_len(k)) {
      if (!is_cat[i]) {
        if (length(unique(md[[vn[i]]])) < 2L) {
          stop(sprintf(
            paste0("fit_competing(): after casewise deletion the continuous ",
                   'covariate "%s" is constant in the %d remaining row%s. It ',
                   "was not constant in the %d rows used, so the deletion, ",
                   "not the column, is what emptied it."),
            covariates[i], n_cc, if (n_cc > 1L) "s" else "", n_used),
            call. = FALSE)
        }
        next
      }
      cnt  <- table(md[[vn[i]]])
      gone <- names(cnt)[cnt == 0L]
      if (length(gone) > 0L) {
        stop(sprintf(
          paste0('fit_competing(): after casewise deletion the covariate "%s" ',
                 "has no rows left in %d of its %d levels (%s).\n",
                 "  Those levels had rows in the %d rows used; the %d row%s ",
                 "deleted for missing covariate values took all of them. The ",
                 "design matrix would carry an all-zero column and crr() ",
                 "would report the rest of the model as if nothing had ",
                 "happened, so fit_competing() stops.\n",
                 "  Either drop the covariate that is causing the deletion, ",
                 "or collapse the emptied levels, and say in the methods ",
                 "which of the two you did."),
          covariates[i], length(gone), length(cnt), .ps_values(gone),
          n_used, n_dropped, if (n_dropped == 1L) "" else "s"),
          call. = FALSE)
      }
      for (j in seq_along(causes)) {
        ev_lv <- vapply(split(as.integer(md$.cr_cause == causes[j]),
                              md[[vn[i]]]), sum, integer(1))
        if (any(ev_lv == 0L)) {
          no_event_lv <- c(no_event_lv, sprintf(
            'cause %d, "%s" level%s %s', causes[j], covariates[i],
            if (sum(ev_lv == 0L) > 1L) "s" else "",
            .ps_values(names(ev_lv)[ev_lv == 0L])))
        }
      }
    }
    if (length(no_event_lv) > 0L) {
      warning(sprintf(
        paste0("fit_competing(): %s ha%s rows but no events.\n",
               "  The subdistribution hazard ratio of such a level is driven ",
               "towards 0 with a confidence interval that carries no ",
               'information; read it as "no events observed", not as a strong ',
               "protective effect. cmprsk::crr() often fails to converge on ",
               "such a model, in which case this call will stop."),
        paste(no_event_lv, collapse = "; "),
        if (length(no_event_lv) > 1L) "ve" else "s"),
        call. = FALSE)
    }

    short <- ev_by_cause_cc < n_coef
    if (any(short)) {
      stop(sprintf(
        paste0("fit_competing(): the model has %d coefficient%s and cause%s ",
               "%s ha%s only %s event%s in the %d complete case%s.\n",
               "  A Fine-Gray model cannot identify more coefficients than ",
               "the cause has events; the estimates would be arbitrary. Fit ",
               "fewer covariates, or collapse levels."),
        n_coef, if (n_coef > 1L) "s" else "",
        if (sum(short) > 1L) "s" else "", .ps_values(causes[short]),
        if (sum(short) > 1L) "ve" else "s",
        paste(ev_by_cause_cc[short], collapse = ", "),
        if (all(ev_by_cause_cc[short] == 1L)) "" else "s",
        n_cc, if (n_cc > 1L) "s" else ""),
        call. = FALSE)
    }
    epv <- ev_by_cause_cc / n_coef
    if (epv_warn > 0 && any(epv < epv_warn)) {
      warning(sprintf(
        paste0("fit_competing(): cause%s %s ha%s %s event%s for %d ",
               "coefficient%s, that is %s events per coefficient, below the ",
               "%g this call asked about. The models are estimable but their ",
               "hazard ratios and intervals are unstable; report the events ",
               "per coefficient with them, or fit fewer covariates. Pass ",
               "epv_warn = 0 to silence this."),
        if (sum(epv < epv_warn) > 1L) "s" else "",
        .ps_values(causes[epv < epv_warn]),
        if (sum(epv < epv_warn) > 1L) "ve" else "s",
        paste(ev_by_cause_cc[epv < epv_warn], collapse = ", "),
        if (all(ev_by_cause_cc[epv < epv_warn] == 1L)) "" else "s",
        n_coef, if (n_coef > 1L) "s" else "",
        paste(sprintf("%.1f", epv[epv < epv_warn]), collapse = ", "),
        epv_warn),
        call. = FALSE)
    }
  }

  ## -- 12. cumulative incidence on the complete cases ---------------------------
  cif <- cif_full
  if (has_cov && cif_complete_case) {
    g_cc <- droplevels(g[cc])
    if (nlevels(g_cc) == nlevels(g) &&
        all(vapply(causes, function(x) sum(cr[cc] == x) > 0L, logical(1)))) {
      cif_cc <- .cr_cif_set(tnum = tnum[cc], cr = cr[cc], g = g_cc,
                            causes = causes, cause_labels = cause_labels,
                            censor_code = censor_code, times = times,
                            z_q = z_q,
                            grouping = if (is.null(group)) "overall" else group,
                            analysis_set = "complete_case",
                            single_group = is.null(group))
      cif <- rbind(cif_full, .cr_strip(cif_cc))
      row.names(cif) <- NULL
    } else {
      warning(sprintf(
        paste0("fit_competing(): casewise deletion emptied a level of the ",
               "grouping variable or a cause, so the cumulative incidence was ",
               "not re-estimated on the complete cases. Only the rows with ",
               'analysis_set = "full" are returned.')),
        call. = FALSE)
    }
  }

  ## -- 13. the design matrix, and the Fine-Gray models --------------------------
  fg_tab   <- NULL
  csh_tab  <- NULL
  side_tab <- NULL
  csh_zph  <- NULL
  design   <- NULL

  if (has_cov) {
    op <- options(contrasts = c(unordered = "contr.treatment",
                                ordered   = "contr.poly"))
    on.exit(options(op), add = TRUE)

    frm <- stats::as.formula(paste0("~ ", paste(vn, collapse = " + ")))
    X_full <- stats::model.matrix(frm, data = md)
    if (!identical(colnames(X_full)[1], "(Intercept)")) {
      stop(paste0("fit_competing(): the first column of the design matrix is ",
                  "not the intercept, which should be impossible by ",
                  "construction. Do not use this result; report it as a bug ",
                  "in fit_competing()."),
           call. = FALSE)
    }
    X <- X_full[, -1L, drop = FALSE]

    expected_cols <- unlist(lapply(seq_len(k), function(i) {
      if (is_cat[i]) paste0(vn[i], lev_list[[i]][-1]) else vn[i]
    }), use.names = FALSE)
    ref_cols <- unlist(lapply(seq_len(k), function(i) {
      if (is_cat[i]) paste0(vn[i], lev_list[[i]][1]) else character(0)
    }), use.names = FALSE)

    if (anyDuplicated(expected_cols) != 0L) {
      stop(sprintf(
        paste0("fit_competing(): two covariate levels would produce the same ",
               "design-matrix column (%s), so their coefficients could not be ",
               "told apart. Rename the levels concerned."),
        .ps_values(unique(expected_cols[duplicated(expected_cols)]))),
        call. = FALSE)
    }
    ## The columns are named internally, so the check above can only fail in
    ## contrived cases. The names the caller sees, variable pasted onto level,
    ## can collide much more easily -- covariate "a" at level "b1" and
    ## covariate "ab" at level "1" both read "ab1" -- and those are the keys
    ## the Fine-Gray and cause-specific tables are joined on.
    term_cols <- unlist(lapply(seq_len(k), function(i) {
      if (is_cat[i]) paste0(covariates[i], lev_list[[i]][-1]) else covariates[i]
    }), use.names = FALSE)
    if (anyDuplicated(term_cols) != 0L) {
      stop(sprintf(
        paste0("fit_competing(): two covariates and levels give the same term ",
               "name (%s), so their rows could not be told apart in the ",
               "returned tables. Rename the level or the column concerned."),
        .ps_values(unique(term_cols[duplicated(term_cols)]))),
        call. = FALSE)
    }
    if (!identical(colnames(X), expected_cols)) {
      stop(sprintf(
        paste0("fit_competing(): the design matrix does not encode the levels ",
               "it should, which should be impossible by construction.\n",
               "  model.matrix() : %s\n  expected       : %s\n",
               "  Do not use this result; report it as a bug in ",
               "fit_competing()."),
        .ps_values(colnames(X)), .ps_values(expected_cols)),
        call. = FALSE)
    }
    if (any(ref_cols %in% colnames(X))) {
      stop(sprintf(
        paste0("fit_competing(): a reference level (%s) has a column of its ",
               "own in the design matrix, which would make the matrix ",
               "singular. Do not use this result; report it as a bug in ",
               "fit_competing()."),
        .ps_values(intersect(ref_cols, colnames(X)))),
        call. = FALSE)
    }
    if (nrow(X) != n_cc || anyNA(X)) {
      stop(sprintf(
        paste0("fit_competing(): the design matrix has %d rows and %d missing ",
               "values where the complete cases are %d and 0. Do not use this ",
               "result; report it as a bug in fit_competing()."),
        nrow(X), sum(is.na(X)), n_cc),
        call. = FALSE)
    }

    ## Column by column: the mean of a dummy column must be the proportion of
    ## the complete cases at that level, and the mean of a continuous column
    ## must be the mean of the covariate. A shifted column fails here.
    design <- data.frame(
      column           = colnames(X),
      variable         = character(ncol(X)),
      level            = NA_character_,
      col_mean         = unname(colMeans(X)),
      expected_mean    = NA_real_,
      is_dummy         = FALSE,
      stringsAsFactors = FALSE,
      row.names        = NULL)
    pos <- 0L
    for (i in seq_len(k)) {
      if (is_cat[i]) {
        for (lvl in lev_list[[i]][-1]) {
          pos <- pos + 1L
          cn  <- paste0(vn[i], lvl)
          got <- mean(X[, cn])
          exp_ <- mean(as.character(md[[vn[i]]]) == lvl)
          if (!all(X[, cn] %in% c(0, 1))) {
            stop(sprintf(
              paste0('fit_competing(): the design-matrix column for "%s" = ',
                     '"%s" is not a 0/1 dummy. Do not use this result; report ',
                     "it as a bug in fit_competing()."),
              covariates[i], lvl),
              call. = FALSE)
          }
          if (abs(got - exp_) > 1e-12) {
            stop(sprintf(
              paste0("fit_competing(): the design-matrix column in position ",
                     '%d has mean %.12f, but "%s" == "%s" holds in %.12f of ',
                     "the %d complete cases. The columns are out of step with ",
                     "the levels they are supposed to encode, so every ",
                     "subdistribution hazard ratio would be attached to the ",
                     "wrong level. Do not use this result; report it as a bug ",
                     "in fit_competing()."),
              pos, got, covariates[i], lvl, exp_, n_cc),
              call. = FALSE)
          }
          design$variable[pos]      <- covariates[i]
          design$level[pos]         <- lvl
          design$expected_mean[pos] <- exp_
          design$is_dummy[pos]      <- TRUE
        }
      } else {
        pos  <- pos + 1L
        got  <- mean(X[, vn[i]])
        exp_ <- mean(md[[vn[i]]])
        if (abs(got - exp_) > 1e-12 * max(1, abs(exp_))) {
          stop(sprintf(
            paste0("fit_competing(): the design-matrix column in position %d ",
                   'has mean %.12f, but the covariate "%s" has mean %.12f in ',
                   "the %d complete cases. Do not use this result; report it ",
                   "as a bug in fit_competing()."),
            pos, got, covariates[i], exp_, n_cc),
            call. = FALSE)
        }
        design$variable[pos]      <- covariates[i]
        design$expected_mean[pos] <- exp_
      }
    }
    design$column <- unlist(lapply(seq_len(k), function(i) {
      if (is_cat[i]) paste0(covariates[i], lev_list[[i]][-1]) else covariates[i]
    }), use.names = FALSE)

    ## the fits, one per cause
    fg_list <- lapply(causes, function(fc) {
      .cr_fit_crr(X = X, md = md, failcode = fc, censor_code = censor_code,
                  n_cc = n_cc, z_q = z_q, no_event_lv = no_event_lv)
    })

    ## a level table with the reference rows, in the order of `covariates`
    lvl_tab <- do.call(rbind, lapply(seq_len(k), function(i) {
      if (is_cat[i]) {
        data.frame(variable = covariates[i], var_label = var_label[i],
                   level = lev_list[[i]],
                   int_term = paste0(vn[i], lev_list[[i]]),
                   is_reference = seq_along(lev_list[[i]]) == 1L,
                   stringsAsFactors = FALSE, row.names = NULL)
      } else {
        data.frame(variable = covariates[i], var_label = var_label[i],
                   level = NA_character_, int_term = vn[i],
                   is_reference = FALSE,
                   stringsAsFactors = FALSE, row.names = NULL)
      }
    }))
    lvl_tab$term <- ifelse(is.na(lvl_tab$level), lvl_tab$variable,
                           paste0(lvl_tab$variable, lvl_tab$level))

    fg_tab <- do.call(rbind, lapply(seq_along(causes), function(j) {
      .cr_fg_rows(fit = fg_list[[j]], lvl_tab = lvl_tab, md = md, vn = vn,
                  is_cat = is_cat, covariates = covariates, causes = causes,
                  j = j, cause_labels = cause_labels, censor_code = censor_code,
                  n_cc = n_cc, n_used = n_used, n_dropped = n_dropped,
                  digits_hr = digits_hr, digits_p = digits_p)
    }))
    row.names(fg_tab) <- NULL

    ## -- 14. cause-specific hazards, fitted by fit_cox() -----------------------
    if (cause_specific == "fit") {
      csh <- .cr_fit_csh(md = md, vn = vn, covariates = covariates,
                         var_label = var_label, causes = causes,
                         cause_labels = cause_labels, conf_level = conf_level,
                         epv_warn = epv_warn, digits_hr = digits_hr,
                         digits_p = digits_p)
      csh_tab <- csh$table
      csh_zph <- csh$zph
      side_tab <- .cr_side_by_side(fg_tab, csh_tab, alpha)
    }
  }

  ## -- 15. invariants ------------------------------------------------------------
  inv <- c(
    "every cumulative incidence is between 0 and 1" =
      all(cif$cif >= -1e-10 & cif$cif <= 1 + 1e-10, na.rm = TRUE),
    "the causes never add up to more than 1 at a time point" =
      .cr_sum_ok(cif),
    "an estimable interval lies strictly inside (0, 1)" =
      all(cif$cif_lcl[cif$estimable] > 0) &&
      all(cif$cif_ucl[cif$estimable] < 1),
    "an estimable interval brackets its point estimate" =
      all(cif$cif_lcl[cif$estimable] <= cif$cif[cif$estimable] + 1e-12) &&
      all(cif$cif[cif$estimable] <= cif$cif_ucl[cif$estimable] + 1e-12),
    "a cell with no interval has both limits missing and says why" =
      all(is.na(cif$cif_lcl[!cif$estimable])) &&
      all(is.na(cif$cif_ucl[!cif$estimable])) &&
      all(grepl("^not estimable \\(", cif$ci_method[!cif$estimable])),
    "a cumulative incidence of 0 means no events of that cause by then" =
      .cr_zero_ok(cif),
    "the cumulative incidences of all causes equal 1 minus all-cause survival" =
      .cr_aalen_johansen_ok(cif, tnum, cr, g, censor_code, cc, has_cov),
    "the groups add up to the rows used" =
      sum(cif$n_group[cif$analysis_set == "full" &
                        cif$cause == causes[1] &
                        cif$time == times[1]]) == n_used)
  if (has_cov) {
    inv <- c(inv,
      "every non-reference Fine-Gray row has a hazard ratio" =
        !anyNA(fg_tab$shr[!fg_tab$is_reference]),
      "every reference Fine-Gray row is 1 with no p-value" =
        all(fg_tab$shr[fg_tab$is_reference] == 1) &&
        all(is.na(fg_tab$p_value[fg_tab$is_reference])),
      "the levels of each covariate add up to the complete cases" =
        .cr_levels_add_up(fg_tab, n_cc),
      "the design matrix reproduces the level proportions" = TRUE,
      "every Fine-Gray model converged" = all(fg_tab$converged))
  }
  if (!all(inv)) {
    stop(sprintf(
      paste0("fit_competing(): the returned tables do not add up, which ",
             "should be impossible by construction.\n  failed check%s : %s\n",
             "  Do not use this result; report it as a bug in ",
             "fit_competing()."),
      if (sum(!inv) > 1L) "s" else "",
      paste(names(inv)[!inv], collapse = "; ")),
      call. = FALSE)
  }

  ## -- 16. assemble ---------------------------------------------------------------
  counts <- c(rows_in            = n_in,
              rows_used          = n_used,
              group_missing      = n_g_na,
              group_missing_dropped = n_g_drop,
              rows_model         = if (has_cov) n_cc else NA_integer_,
              rows_model_dropped = if (has_cov) n_dropped else NA_integer_,
              censored           = sum(cr == censor_code),
              deaths             = sum(cr != censor_code),
              coefficients       = if (has_cov) n_coef else NA_integer_)
  counts <- c(counts,
              stats::setNames(as.numeric(n_by_cause_used),
                              paste0("events_cause", causes)))
  if (has_cov) {
    counts <- c(counts,
                stats::setNames(as.numeric(ev_by_cause_cc),
                                paste0("events_cause", causes, "_model")),
                stats::setNames(as.numeric(ev_by_cause_cc / n_coef),
                                paste0("epv_cause", causes)))
  }

  ## strip the working attributes so that the result carries exactly the ones
  ## the documentation lists, and nothing else
  out <- .cr_strip(cif)
  attr(out, "finegray")       <- fg_tab
  attr(out, "cause_specific") <- csh_tab
  attr(out, "fg_vs_csh")      <- side_tab
  attr(out, "dropna")         <- dropna_tab
  attr(out, "outcome_check")  <- list(table  = os_tab,
                                      counts = os_counts,
                                      checks = os_checks,
                                      note   = os_note)
  attr(out, "fit_competing")  <- list(
    counts    = counts,
    follow_up = c(min       = min(tnum),
                  median    = stats::median(tnum),
                  max       = max(tnum),
                  total     = sum(tnum),
                  tau_upper = unname(tau_upper)),
    gray      = gray_full,
    missing   = missing_tab,
    design    = design,
    checks    = inv,
    cause_specific_zph = csh_zph,
    settings  = list(
      time              = time,
      event_cr          = event_cr,
      event_os          = event_os,
      group             = group,
      group_levels      = lv,
      group_levels_dropped = n_lv_drop,
      group_missing     = if (is.null(group)) NULL else group_missing,
      covariates        = covariates,
      labels            = if (has_cov) stats::setNames(var_label, covariates)
                          else NULL,
      reference         = if (has_cov) stats::setNames(ref_lev, covariates)
                          else NULL,
      levels            = if (has_cov) stats::setNames(lev_list, covariates)
                          else NULL,
      levels_dropped    = if (has_cov) n_lv_cov else NULL,
      causes            = causes,
      cause_labels      = stats::setNames(cause_labels, causes),
      censor_code       = censor_code,
      times             = times,
      times_beyond      = times_beyond,
      na_action         = na_action,
      empty_levels      = empty_levels,
      cif_complete_case = cif_complete_case,
      cause_specific    = cause_specific,
      conf_level        = conf_level,
      alpha             = alpha,
      epv_warn          = epv_warn,
      digits_hr         = digits_hr,
      digits_p          = digits_p,
      estimator_cif     = "cmprsk::cuminc (Aalen-Johansen)",
      estimator_test    = "cmprsk::cuminc, Gray's test, rho = 0",
      estimator_model   = "cmprsk::crr (Fine-Gray subdistribution hazard)",
      estimator_csh     = if (has_cov && cause_specific == "fit")
        "fit_cox() -> survival::coxph, competing causes censored" else NULL,
      ci_cif            = "complementary log-log transform",
      ci_model          = "Wald on the log scale, exp(coef +/- z * se)",
      contrasts         = "contr.treatment, reference = first level",
      random_numbers    = "none: cuminc() and crr() are deterministic"),
    call = cl)
  out
}


# -- internal helpers, not exported ------------------------------------------
# .ps_values() and .ps_counts() are defined in R/prep_surv.R and reused here.

#' @noRd
.cr_name_arg <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf(
      "fit_competing(): `%s` must be a single non-missing string.", arg),
      call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.cr_count_arg <- function(x, arg, min = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != as.integer(x) ||
      x < min) {
    stop(sprintf(
      "fit_competing(): `%s` must be a single whole number >= %d, got %s.",
      arg, min, .ps_values(x)),
      call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.cr_times_arg <- function(times) {
  if (!is.numeric(times) || length(times) == 0L) {
    stop(sprintf(
      paste0("fit_competing(): `times` must be a numeric vector of at least ",
             "one time point, got %s of length %d."),
      .ps_values(class(times)), length(times)),
      call. = FALSE)
  }
  bad <- is.na(times) | !is.finite(times)
  if (any(bad)) {
    stop(sprintf(
      paste0("fit_competing(): `times` has %d missing or infinite element%s ",
             "(position%s %s)."),
      sum(bad), if (sum(bad) > 1L) "s" else "",
      if (sum(bad) > 1L) "s" else "", paste(which(bad), collapse = ", ")),
      call. = FALSE)
  }
  if (any(times < 0)) {
    stop(sprintf(
      paste0("fit_competing(): `times` has %d negative element%s (%s). A ",
             "cumulative incidence before the origin does not exist."),
      sum(times < 0), if (sum(times < 0) > 1L) "s" else "",
      .ps_values(times[times < 0])),
      call. = FALSE)
  }
  d <- unique(times[duplicated(times)])
  if (length(d) > 0L) {
    stop(sprintf(
      paste0("fit_competing(): `times` lists %s more than once, which would ",
             "put the same row in the result twice. Keep one copy."),
      .ps_values(d)),
      call. = FALSE)
  }
  sort(as.numeric(times))
}

#' @noRd
.cr_cause_labels <- function(cause_labels, causes) {
  if (is.null(cause_labels)) return(paste("Cause", causes))
  if (!is.character(cause_labels) || anyNA(cause_labels)) {
    stop(sprintf(
      paste0("fit_competing(): `cause_labels` must be a character vector ",
             "without missing values, got %s."),
      .ps_values(class(cause_labels))),
      call. = FALSE)
  }
  if (!is.null(names(cause_labels))) {
    j <- match(as.character(causes), names(cause_labels))
    if (anyNA(j)) {
      stop(sprintf(
        paste0("fit_competing(): `cause_labels` is named, so it must have one ",
               "name per cause. No label for cause%s %s.\n",
               "  names given : %s"),
        if (sum(is.na(j)) > 1L) "s" else "", .ps_values(causes[is.na(j)]),
        .ps_values(names(cause_labels))),
        call. = FALSE)
    }
    return(unname(cause_labels[j]))
  }
  if (length(cause_labels) != length(causes)) {
    stop(sprintf(
      paste0("fit_competing(): `cause_labels` has %d element%s and `causes` ",
             "has %d. Give one label per cause, in the same order, or name ",
             "the labels by the cause codes."),
      length(cause_labels), if (length(cause_labels) > 1L) "s" else "",
      length(causes)),
      call. = FALSE)
  }
  cause_labels
}

#' @noRd
.cr_labels_arg <- function(labels, covariates) {
  out <- covariates
  if (is.null(labels)) return(out)
  if (is.list(labels)) labels <- unlist(labels, use.names = TRUE)
  if (!is.character(labels) || is.null(names(labels)) ||
      any(!nzchar(names(labels)))) {
    stop(paste0("fit_competing(): `labels` must be a named character vector ",
                "or named list, covariate = \"label\"."),
         call. = FALSE)
  }
  unknown <- setdiff(names(labels), covariates)
  if (length(unknown) > 0L) {
    stop(sprintf(
      paste0("fit_competing(): `labels` names %s, which %s not in ",
             "`covariates` (%s). A label for a covariate that is not in the ",
             "model is usually a typo in one of the two."),
      .ps_values(unknown), if (length(unknown) > 1L) "are" else "is",
      .ps_values(covariates)),
      call. = FALSE)
  }
  j <- match(covariates, names(labels))
  out[!is.na(j)] <- unname(labels[j[!is.na(j)]])
  out
}

#' @noRd
.cr_fmt_p <- function(p, digits) {
  thr <- 10^(-digits)
  ifelse(is.na(p), NA_character_,
         ifelse(p < thr,
                paste0("<", formatC(thr, format = "f", digits = digits)),
                formatC(p, format = "f", digits = digits)))
}

#' @noRd
.cr_strip <- function(x) {
  for (a in setdiff(names(attributes(x)), c("names", "row.names", "class"))) {
    attr(x, a) <- NULL
  }
  x
}

# Estimate the cumulative incidence of every cause in every group at `times`,
# with a complementary log-log confidence interval, and Gray's test.
#
# The element names cmprsk::cuminc() gives its output are "<group> <cause>".
# They are matched here by building the expected names and looking them up,
# rather than by chopping characters off the end, so that a two-digit cause
# code or a group label ending in a digit cannot silently mis-assign a curve.
#' @noRd
.cr_cif_set <- function(tnum, cr, g, causes, cause_labels, censor_code, times,
                        z_q, grouping, analysis_set, single_group) {
  lv <- levels(g)
  ci <- tryCatch(
    cmprsk::cuminc(ftime = tnum, fstatus = cr, group = g,
                   cencode = censor_code, na.action = stats::na.fail),
    error = function(e) {
      stop(sprintf(
        paste0("fit_competing(): cmprsk::cuminc() failed on the %s analysis ",
               "set.\n  cuminc() said : %s"),
        analysis_set, conditionMessage(e)),
        call. = FALSE)
    })

  nm_have   <- setdiff(names(ci), "Tests")
  nm_expect <- paste(rep(lv, times = length(causes)),
                     rep(causes, each = length(lv)))
  if (anyDuplicated(nm_expect) != 0L) {
    stop(sprintf(
      paste0("fit_competing(): two (group, cause) pairs would carry the same ",
             "cmprsk::cuminc() element name (%s), so their curves could not ",
             "be told apart. Rename the group levels concerned."),
      .ps_values(unique(nm_expect[duplicated(nm_expect)]))),
      call. = FALSE)
  }
  if (!setequal(nm_have, nm_expect)) {
    stop(sprintf(
      paste0("fit_competing(): cmprsk::cuminc() returned curves for a ",
             "different set of (group, cause) pairs than the data imply, on ",
             "the %s analysis set.\n  cuminc()  : %s\n  expected  : %s\n",
             "  Do not use this result; report it as a bug in ",
             "fit_competing()."),
      analysis_set, .ps_values(nm_have), .ps_values(nm_expect)),
      call. = FALSE)
  }

  tp <- cmprsk::timepoints(ci, times = times)
  if (!isTRUE(all.equal(as.numeric(colnames(tp$est)), times))) {
    stop(paste0("fit_competing(): cmprsk::timepoints() did not return the ",
                "requested time points, which should be impossible by ",
                "construction. Do not use this result; report it as a bug in ",
                "fit_competing()."),
         call. = FALSE)
  }
  i_row <- match(nm_expect, rownames(tp$est))
  if (anyNA(i_row)) {
    stop(paste0("fit_competing(): cmprsk::timepoints() dropped a curve that ",
                "cmprsk::cuminc() returned. Do not use this result; report it ",
                "as a bug in fit_competing()."),
         call. = FALSE)
  }

  nt <- length(times)
  ng <- length(lv)
  nc <- length(causes)
  ## one row per (cause, group, time), cause-major, exactly the order in which
  ## cuminc() lists its elements
  idx_gc  <- rep(seq_len(ng * nc), each = nt)
  grp_i   <- rep(rep(seq_len(ng), times = nc), each = nt)
  cause_i <- rep(rep(seq_len(nc), each = ng), each = nt)
  time_i  <- rep(seq_len(nt), times = ng * nc)

  est <- tp$est[i_row, , drop = FALSE]
  vr  <- tp$var[i_row, , drop = FALSE]
  cifv <- est[cbind(idx_gc, time_i)]
  varv <- vr[cbind(idx_gc, time_i)]
  sev  <- sqrt(varv)

  ## the complementary log-log interval; see the details section
  ok  <- !is.na(cifv) & !is.na(sev) & cifv > 0 & cifv < 1 & sev > 0
  lcl <- rep(NA_real_, length(cifv))
  ucl <- rep(NA_real_, length(cifv))
  j   <- which(ok)
  if (length(j) > 0L) {
    nu    <- log(-log(1 - cifv[j]))
    se_nu <- sev[j] / abs((1 - cifv[j]) * log(1 - cifv[j]))
    lcl[j] <- 1 - exp(-exp(nu - z_q * se_nu))
    ucl[j] <- 1 - exp(-exp(nu + z_q * se_nu))
  }
  ## In exact arithmetic the back-transformed limits are strictly inside
  ## (0, 1). In floating point they are not always: an estimate that the
  ## product-limit arithmetic returns as 1 - 2e-16, or one whose standard
  ## error dwarfs it, drives exp(-exp(.)) to 0 or to 1 and the interval
  ## collapses onto the boundary. Such an interval carries no information and
  ## is exactly the degenerate thing this transform was chosen to avoid, so it
  ## is withdrawn and labelled rather than printed as 0.000-1.000.
  deg <- ok & (!is.finite(lcl) | !is.finite(ucl) | lcl <= 0 | ucl >= 1)
  ci_method <- ifelse(
    ok & !deg, "complementary log-log transform",
    ifelse(is.na(cifv),
           "not estimable (beyond the largest follow-up time in this group)",
           ifelse(cifv >= 1 - 1e-12,
                  "not estimable (cumulative incidence is 1)",
                  ifelse(deg,
                         "not estimable (the interval degenerates at 0 or 1)",
                         "not estimable (no events at this horizon)"))))
  lcl[deg] <- NA_real_
  ucl[deg] <- NA_real_
  ok       <- ok & !deg

  ## counts that belong beside the estimate
  n_group  <- as.integer(table(g)[lv])[grp_i]
  n_risk   <- integer(length(cifv))
  n_ev_all <- integer(length(cifv))
  n_ev_by  <- integer(length(cifv))
  for (r in seq_along(cifv)) {
    in_g <- g == lv[grp_i[r]]
    n_risk[r]   <- sum(in_g & tnum >= times[time_i[r]])
    n_ev_all[r] <- sum(in_g & cr == causes[cause_i[r]])
    n_ev_by[r]  <- sum(in_g & cr == causes[cause_i[r]] &
                         tnum <= times[time_i[r]])
  }

  ## Gray's test: a property of a cause across all groups
  gray <- .cr_gray(ci, causes, cause_labels, grouping, analysis_set,
                   single_group, ng)

  out <- data.frame(
    grouping     = grouping,
    analysis_set = analysis_set,
    group        = if (single_group) "ALL" else lv[grp_i],
    cause        = causes[cause_i],
    cause_label  = cause_labels[cause_i],
    time         = times[time_i],
    n_group      = n_group,
    n_risk       = n_risk,
    n_event_cause         = n_ev_all,
    n_event_cause_by_time = n_ev_by,
    cif          = unname(cifv),
    cif_se       = unname(sev),
    cif_lcl      = lcl,
    cif_ucl      = ucl,
    ci_method    = ci_method,
    estimable    = ok,
    gray_stat    = gray$gray_stat[match(causes[cause_i], gray$cause)],
    gray_df      = gray$gray_df[match(causes[cause_i], gray$cause)],
    gray_p       = gray$gray_p[match(causes[cause_i], gray$cause)],
    gray_p_fmt   = gray$gray_p_fmt[match(causes[cause_i], gray$cause)],
    stringsAsFactors = FALSE,
    row.names        = NULL)
  attr(out, "gray") <- gray
  out
}

#' @noRd
.cr_gray <- function(ci, causes, cause_labels, grouping, analysis_set,
                     single_group, ng) {
  if (single_group || ng < 2L) {
    return(data.frame(
      grouping = grouping, analysis_set = analysis_set,
      cause = causes, cause_label = cause_labels,
      gray_stat = NA_real_, gray_df = NA_integer_, gray_p = NA_real_,
      gray_p_fmt = NA_character_,
      note = paste0("not computed: there is a single group, so there is ",
                    "nothing to compare"),
      stringsAsFactors = FALSE, row.names = NULL))
  }
  tt <- ci$Tests
  if (is.null(tt) || !all(c("stat", "pv", "df") %in% colnames(tt))) {
    stop(paste0("fit_competing(): cmprsk::cuminc() did not return the Gray ",
                "test table this version of fit_competing() knows how to ",
                "read. Do not use this result; report it as a bug in ",
                "fit_competing()."),
         call. = FALSE)
  }
  j <- match(as.character(causes), rownames(tt))
  if (anyNA(j)) {
    stop(sprintf(
      paste0("fit_competing(): cmprsk::cuminc() returned no Gray test for ",
             "cause%s %s. Do not use this result; report it as a bug in ",
             "fit_competing()."),
      if (sum(is.na(j)) > 1L) "s" else "", .ps_values(causes[is.na(j)])),
      call. = FALSE)
  }
  data.frame(
    grouping = grouping, analysis_set = analysis_set,
    cause = causes, cause_label = cause_labels,
    gray_stat = unname(tt[j, "stat"]),
    gray_df   = as.integer(unname(tt[j, "df"])),
    gray_p    = unname(tt[j, "pv"]),
    gray_p_fmt = .cr_fmt_p(unname(tt[j, "pv"]), 3L),
    note = sprintf("Gray's test across the %d groups of \"%s\"", ng, grouping),
    stringsAsFactors = FALSE, row.names = NULL)
}

#' @noRd
.cr_fit_crr <- function(X, md, failcode, censor_code, n_cc, z_q, no_event_lv) {
  fit <- tryCatch(
    cmprsk::crr(ftime = md$.cr_time, fstatus = md$.cr_cause, cov1 = X,
                failcode = failcode, cencode = censor_code,
                na.action = stats::na.fail),
    error = function(e) {
      stop(sprintf(
        paste0("fit_competing(): cmprsk::crr() failed for cause %d.\n",
               "  crr() said : %s"),
        failcode, conditionMessage(e)),
        call. = FALSE)
    })

  if (!isTRUE(fit$converged)) {
    stop(sprintf(
      paste0("fit_competing(): the Fine-Gray model for cause %d did not ",
             "converge, so its coefficients are wherever the iteration ",
             "stopped and must not be reported.\n%s",
             "  Collapse the levels that carry few events of this cause, or ",
             "fit fewer covariates."),
      failcode,
      if (length(no_event_lv) > 0L)
        sprintf("  Levels with rows but no events: %s\n",
                paste(no_event_lv, collapse = "; ")) else ""),
      call. = FALSE)
  }
  if (fit$n != n_cc || fit$n.missing != 0L) {
    stop(sprintf(
      paste0("fit_competing(): cmprsk::crr() used %d rows and dropped %d for ",
             "missing values, where the complete cases are %d and 0. Do not ",
             "use this result; report it as a bug in fit_competing()."),
      fit$n, fit$n.missing, n_cc),
      call. = FALSE)
  }
  if (!identical(names(fit$coef), colnames(X))) {
    stop(sprintf(
      paste0("fit_competing(): the coefficients cmprsk::crr() returned are ",
             "not the columns of the design matrix, so every subdistribution ",
             "hazard ratio would be attached to the wrong level.\n",
             "  crr()    : %s\n  expected : %s\n",
             "  Do not use this result; report it as a bug in ",
             "fit_competing()."),
      .ps_values(names(fit$coef)), .ps_values(colnames(X))),
      call. = FALSE)
  }
  se <- sqrt(diag(fit$var))
  if (anyNA(fit$coef) || anyNA(se) || any(!is.finite(se)) || any(se <= 0)) {
    stop(sprintf(
      paste0("fit_competing(): the Fine-Gray model for cause %d returned a ",
             "missing or non-positive standard error, so its confidence ",
             "intervals do not exist. This is what an unidentified ",
             "coefficient looks like: two covariates carrying the same ",
             "information, or a level with no events of this cause."),
      failcode),
      call. = FALSE)
  }
  list(coef = unname(fit$coef), se = unname(se), names = names(fit$coef),
       converged = isTRUE(fit$converged), loglik = unname(fit$loglik),
       loglik_null = unname(fit$loglik.null), z_q = z_q)
}

#' @noRd
.cr_fg_rows <- function(fit, lvl_tab, md, vn, is_cat, covariates, causes, j,
                        cause_labels, censor_code, n_cc, n_used, n_dropped,
                        digits_hr, digits_p) {
  fc  <- causes[j]
  pos <- match(lvl_tab$int_term, fit$names)
  pos[lvl_tab$is_reference] <- NA_integer_

  coef_v <- ifelse(is.na(pos), NA_real_, fit$coef[pos])
  se_v   <- ifelse(is.na(pos), NA_real_, fit$se[pos])
  z_v    <- coef_v / se_v
  p_v    <- 2 * stats::pnorm(-abs(z_v))
  shr_v  <- ifelse(lvl_tab$is_reference, 1, exp(coef_v))
  lo_v   <- exp(coef_v - fit$z_q * se_v)
  hi_v   <- exp(coef_v + fit$z_q * se_v)
  fmt    <- function(x) sprintf("%.*f", digits_hr, x)

  n_lv    <- integer(nrow(lvl_tab))
  ev_lv   <- integer(nrow(lvl_tab))
  evo_lv  <- integer(nrow(lvl_tab))
  cens_lv <- integer(nrow(lvl_tab))
  for (r in seq_len(nrow(lvl_tab))) {
    i <- match(lvl_tab$variable[r], covariates)
    if (is_cat[i]) {
      sel <- as.character(md[[vn[i]]]) == lvl_tab$level[r]
    } else {
      sel <- rep(TRUE, nrow(md))
    }
    n_lv[r]    <- sum(sel)
    ev_lv[r]   <- sum(sel & md$.cr_cause == fc)
    evo_lv[r]  <- sum(sel & md$.cr_cause != fc & md$.cr_cause != censor_code)
    cens_lv[r] <- sum(sel & md$.cr_cause == censor_code)
  }

  data.frame(
    failcode      = fc,
    cause_label   = cause_labels[j],
    analysis_set  = "complete_case",
    variable      = lvl_tab$variable,
    var_label     = lvl_tab$var_label,
    level         = lvl_tab$level,
    term          = lvl_tab$term,
    is_reference  = lvl_tab$is_reference,
    n             = n_lv,
    n_event_cause = ev_lv,
    n_event_competing = evo_lv,
    n_censored    = cens_lv,
    shr           = shr_v,
    ci_low        = lo_v,
    ci_high       = hi_v,
    shr_txt       = ifelse(lvl_tab$is_reference,
                           paste0(fmt(1), " (reference)"),
                           sprintf("%s (%s-%s)", fmt(shr_v), fmt(lo_v),
                                   fmt(hi_v))),
    coef          = coef_v,
    se            = se_v,
    z             = z_v,
    p_value       = p_v,
    p_fmt         = .cr_fmt_p(p_v, digits_p),
    converged     = fit$converged,
    loglik        = fit$loglik,
    loglik_null   = fit$loglik_null,
    n_model       = n_cc,
    n_full        = n_used,
    n_dropped     = n_dropped,
    stringsAsFactors = FALSE,
    row.names        = NULL)
}

# The cause-specific hazard models. They go through fit_cox() so that they are
# fitted on exactly the same complete cases, with the same reference levels,
# and so that cox.zph() comes with them.
#' @noRd
.cr_fit_csh <- function(md, vn, covariates, var_label, causes, cause_labels,
                        conf_level, epv_warn, digits_hr, digits_p) {
  tab <- vector("list", length(causes))
  zph <- vector("list", length(causes))
  for (j in seq_along(causes)) {
    d <- md[, vn, drop = FALSE]
    names(d) <- covariates
    d$.cr_time  <- md$.cr_time
    d$.cr_event <- as.integer(md$.cr_cause == causes[j])
    fit <- tryCatch(
      fit_cox(d, covariates = covariates, time = ".cr_time",
              event = ".cr_event",
              labels = stats::setNames(as.list(var_label), covariates),
              conf_level = conf_level, na_action = "fail",
              empty_levels = "error", epv_warn = epv_warn,
              digits_hr = digits_hr, digits_p = digits_p),
      error = function(e) {
        stop(sprintf(
          paste0("fit_competing(): the cause-specific Cox model for cause %d ",
                 "could not be fitted.\n  fit_cox() said : %s\n",
                 '  Pass cause_specific = "skip" to return the Fine-Gray ',
                 "models without it, but do not quote a subdistribution ",
                 "hazard ratio without saying that it is one."),
          causes[j], conditionMessage(e)),
          call. = FALSE)
      })
    z <- attr(fit, "cox_zph")
    fit <- .cr_strip(fit)
    tab[[j]] <- cbind(
      data.frame(failcode = causes[j], cause_label = cause_labels[j],
                 analysis_set = "complete_case",
                 stringsAsFactors = FALSE, row.names = NULL),
      fit)
    z$failcode <- causes[j]
    zph[[j]]   <- z
  }
  out <- do.call(rbind, tab)
  row.names(out) <- NULL
  list(table = out, zph = stats::setNames(zph, paste0("cause", causes)))
}

#' @noRd
.cr_side_by_side <- function(fg_tab, csh_tab, alpha) {
  a <- fg_tab[!fg_tab$is_reference,
              c("failcode", "cause_label", "variable", "level", "term", "shr",
                "shr_txt", "ci_low", "ci_high", "p_value")]
  names(a)[names(a) == "p_value"] <- "shr_p"
  names(a)[names(a) == "ci_low"]  <- "shr_ci_low"
  names(a)[names(a) == "ci_high"] <- "shr_ci_high"
  b <- csh_tab[!csh_tab$is_reference,
               c("failcode", "term", "hr", "hr_txt", "ci_low", "ci_high",
                 "p_value")]
  names(b)[names(b) == "p_value"] <- "csh_p"
  names(b)[names(b) == "ci_low"]  <- "csh_ci_low"
  names(b)[names(b) == "ci_high"] <- "csh_ci_high"

  key_a <- paste(a$failcode, a$term, sep = "\r")
  key_b <- paste(b$failcode, b$term, sep = "\r")
  j <- match(key_a, key_b)
  if (anyNA(j) || anyDuplicated(key_b) > 0L) {
    stop(paste0("fit_competing(): the Fine-Gray and the cause-specific tables ",
                "do not have the same terms, which should be impossible by ",
                "construction. Do not use this result; report it as a bug in ",
                "fit_competing()."),
         call. = FALSE)
  }
  out <- cbind(a, b[j, setdiff(names(b), c("failcode", "term")), drop = FALSE])
  row.names(out) <- NULL
  out$both_sig       <- out$shr_p < alpha & out$csh_p < alpha
  out$same_side      <- (out$shr > 1) == (out$hr > 1)
  out$direction_flip <- out$both_sig & !out$same_side
  out$alpha          <- alpha
  out
}

#' @noRd
.cr_dropna_table <- function(g, lv, cc, cr, causes) {
  per <- do.call(rbind, lapply(lv, function(s) {
    sel <- g == s
    row <- data.frame(
      group            = s,
      is_total         = FALSE,
      n                = sum(sel),
      n_dropped        = sum(sel & !cc),
      n_dropped_events = sum(sel & !cc & cr %in% causes),
      stringsAsFactors = FALSE, row.names = NULL)
    for (x in causes) {
      row[[paste0("n_dropped_cause", x)]] <- sum(sel & !cc & cr == x)
    }
    row
  }))
  tot <- data.frame(
    group            = "ALL",
    is_total         = TRUE,
    n                = length(cc),
    n_dropped        = sum(!cc),
    n_dropped_events = sum(!cc & cr %in% causes),
    stringsAsFactors = FALSE, row.names = NULL)
  for (x in causes) {
    tot[[paste0("n_dropped_cause", x)]] <- sum(!cc & cr == x)
  }
  out <- rbind(per, tot)
  out$pct_dropped <- 100 * out$n_dropped / out$n
  out <- out[, c("group", "is_total", "n", "n_dropped", "pct_dropped",
                 "n_dropped_events", paste0("n_dropped_cause", causes))]
  row.names(out) <- NULL
  out
}

#' @noRd
.cr_sum_ok <- function(cif) {
  d   <- cif[!is.na(cif$cif), , drop = FALSE]
  if (nrow(d) == 0L) return(TRUE)
  key <- paste(d$analysis_set, d$group, d$time, sep = "\r")
  all(vapply(split(d$cif, key), sum, numeric(1)) <= 1 + 1e-10)
}

#' @noRd
.cr_zero_ok <- function(cif) {
  d <- cif[!is.na(cif$cif), , drop = FALSE]
  if (nrow(d) == 0L) return(TRUE)
  all((d$cif == 0) == (d$n_event_cause_by_time == 0))
}

# The Aalen-Johansen identity: the cumulative incidences of all causes add up
# to 1 minus the all-cause Kaplan-Meier survival, at every time point of every
# group. The two sides use the same rows and the same censoring, so any gap
# means the causes and the all-cause indicator have come apart.
#' @noRd
.cr_aalen_johansen_ok <- function(cif, tnum, cr, g, censor_code, cc, has_cov) {
  lv  <- levels(g)
  ev  <- as.integer(cr != censor_code)
  gap <- 0
  for (aset in unique(cif$analysis_set)) {
    sel_rows <- if (aset == "full") rep(TRUE, length(tnum)) else cc
    d <- cif[cif$analysis_set == aset, , drop = FALSE]
    for (s in unique(d$group)) {
      in_g <- if (identical(s, "ALL") && nlevels(g) == 1L) sel_rows else
        sel_rows & g == s
      if (!any(in_g)) return(FALSE)
      tt <- sort(unique(d$time[d$group == s]))
      sf <- survival::survfit(
        survival::Surv(tnum[in_g], ev[in_g]) ~ 1)
      km <- summary(sf, times = tt, extend = TRUE)$surv
      for (i in seq_along(tt)) {
        v <- d$cif[d$group == s & d$time == tt[i]]
        if (anyNA(v)) next
        gap <- max(gap, abs(sum(v) - (1 - km[i])))
      }
    }
  }
  gap < 1e-8
}

#' @noRd
.cr_levels_add_up <- function(fg_tab, n_cc) {
  key <- paste(fg_tab$failcode, fg_tab$variable, sep = "\r")
  s   <- vapply(split(fg_tab$n[!is.na(fg_tab$level)],
                      key[!is.na(fg_tab$level)]), sum, numeric(1))
  length(s) == 0L || all(s == n_cc)
}
