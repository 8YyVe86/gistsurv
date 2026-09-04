# ---------------------------------------------------------------------------
# fit_cox(): a Cox proportional hazards model, returned as a data frame,
# always together with the proportional-hazards test.
#
# Design notes (deliberate, and different from the study script this was
# distilled from):
#   * no column name is hard-coded: the time, the event indicator and the
#     covariates are all arguments;
#   * nothing is printed, nothing is plotted and nothing is written to disk:
#     one row per covariate level comes back as a data frame, and the
#     proportional-hazards test, the goodness-of-fit statistics and the
#     bookkeeping as columns and attributes;
#   * the cox.zph() test is not optional. It is computed on every call, it is
#     part of the returned data frame (as columns, which survive the
#     data-frame verbs that drop attributes) and there is no argument that
#     switches it off. If cox.zph() fails, fit_cox() fails: this function does
#     not hand back a hazard ratio without the assumption it rests on;
#   * a failed proportional-hazards test is reported, never acted on.
#     fit_cox() will not silently stratify, add a time-dependent coefficient
#     or truncate follow-up to make the test pass; those change the research
#     question and are the caller's decision;
#   * missing covariate values are deleted casewise, and the number of rows
#     and events deleted is a returned number rather than a footnote inside
#     the fitted object;
#   * no fitted object is attached. A coxph object stores the response
#     (coxph(y = TRUE) is the default), the linear predictor and the
#     residuals -- one value per subject -- and drags the caller's data along
#     inside the environment of its formula; a cox.zph object stores the
#     scaled Schoenfeld residuals, one row per event time. Both would put
#     individual records into anything the result is saved into. Everything
#     needed to reproduce the fit is recorded in the settings diagnostic
#     instead.
# The file is kept ASCII-only so that it behaves the same under a UTF-8 and
# under a non-UTF-8 locale.
#
# The small formatting helpers .ps_values() and .ps_counts() are shared with
# R/prep_surv.R; they are internal to the package.
# ---------------------------------------------------------------------------

#' Cox proportional hazards model with the proportional-hazards test
#'
#' @description
#' `fit_cox()` fits a Cox proportional hazards model and returns the numbers a
#' survival paper reports: one row per covariate level with its hazard ratio,
#' confidence interval and p-value, the reference level marked as such, the
#' subjects and events behind each level, and -- on every row, not as an
#' option -- the [survival::cox.zph()] test of the proportional-hazards
#' assumption. The concordance with its standard error, the likelihood ratio
#' test, the number of subjects and events, and the number of rows deleted for
#' missing covariate values come back with the result.
#'
#' The time, event and covariate columns are arguments, so the function is not
#' tied to any particular registry export. The function prints nothing, plots
#' nothing and writes nothing.
#'
#' @details
#' # The proportional-hazards test is not optional
#'
#' A hazard ratio from a Cox model is a single number standing for the whole
#' of follow-up. It only means what it is usually taken to mean if the hazards
#' really are proportional, so the test of that assumption is reported with
#' the estimate rather than alongside it:
#'
#' * `cox.zph()` runs on every call. There is no argument that turns it off;
#' * its per-covariate chi-square, degrees of freedom and p-value are
#'   **columns** of the returned data frame (`zph_chisq_var`, `zph_df_var`,
#'   `zph_p_var`), repeated on every level of that covariate, and the global
#'   test is in `zph_chisq_global`, `zph_df_global` and `zph_p_global`. They
#'   are columns and not only attributes because most data-frame verbs drop
#'   attributes, and a hazard ratio should not be able to outlive its
#'   assumption check by being piped somewhere;
#' * the full test, including the `GLOBAL` row, is also attached as the
#'   `cox_zph` attribute;
#' * if `cox.zph()` fails, or returns a test with a missing chi-square or
#'   p-value, `fit_cox()` fails. A missing test is not a passed test, and the
#'   hazard ratios are not returned with the test left empty.
#'
#' `ph_alpha` (default `0.05`) only sets the `ph_violated_var` and
#' `ph_violated_global` flags. A flagged model is still returned exactly as
#' fitted: `fit_cox()` never stratifies on the offending covariate, adds a
#' time-dependent coefficient or shortens follow-up on its own. Those are
#' different models answering different questions, and choosing between them
#' -- or reporting the unmodified hazard ratios and describing the time
#' dimension with a restricted mean survival time instead -- is the analyst's
#' decision.
#'
#' `zph_transform` is passed to `cox.zph()` as `transform` and decides the
#' function of time the Schoenfeld residuals are regressed on. The default,
#' `"km"`, is the default of `cox.zph()` itself, and therefore the test a
#' script calling `cox.zph(fit)` without naming `transform` has been reporting
#' all along. Changing it changes the chi-square and the p-value of a test
#' whose estimates do not move at all, so the value used is recorded in the
#' `settings` diagnostic and is never chosen on your behalf.
#'
#' The test is computed with `terms = TRUE`, i.e. one row per covariate, on as
#' many degrees of freedom as that covariate has coefficients, so that a
#' factor is tested as one variable rather than as a set of unrelated
#' contrasts.
#'
#' # What is returned, and what is not
#'
#' One row per level of each categorical covariate, in the order of the levels
#' of the factor, and one row per continuous covariate. The first level of a
#' factor is the reference: `is_reference` is `TRUE`, `hr` is `1`, the
#' confidence limits and the p-value are `NA` and `hr_txt` reads
#' `"1.00 (reference)"`. The hazard ratio of a reference row is written as `1`
#' rather than left missing so that the column can be plotted or sorted
#' without a special case, and the missing interval, the missing p-value and
#' the `is_reference` flag are all there to stop it being read as an estimate.
#'
#' No `coxph` object is attached to the result. A fitted Cox model stores the
#' response -- `coxph()` keeps `y` by default -- together with the linear
#' predictor and the residuals, one value per subject, and its formula carries
#' the environment the model was fitted in; a `cox.zph` object stores the
#' scaled Schoenfeld residuals, one row per event time. Attaching either would
#' put individual follow-up times into every file the result is written to.
#' Every argument that defines the fit is recorded in the `settings`
#' diagnostic instead, so the same fit can be reproduced exactly.
#'
#' # Missing values
#'
#' `na_action = "omit"` (the default) fits the model on the complete cases,
#' which is what [survival::coxph()] does through `na.action` anyway; the
#' difference is that the rows and the events lost are counted and returned
#' (`counts["rows_dropped"]`, `counts["events_dropped"]`, and the
#' per-variable `missing` table) instead of being left inside the fitted
#' object. That count is what a methods section has to quote, and what tells
#' you whether the complete-case analysis is defensible. Nothing is ever
#' imputed.
#'
#' `na_action = "fail"` stops instead, with the same counts in the message,
#' for an analysis that has already decided that no row may be dropped.
#'
#' A missing follow-up time or a missing event indicator is always an error,
#' whatever `na_action` says: those are outcome variables, and a subject with
#' no outcome is not a complete case under any definition.
#'
#' # Contrasts and reference levels
#'
#' The reference level of a factor is its first level; relevel the factor
#' before calling if you want another one. For the duration of the fit
#' `fit_cox()` sets `options(contrasts = c("contr.treatment", "contr.poly"))`
#' and restores the previous setting on exit, so that the reference level is
#' the first level whatever the session-wide option happens to be. An ordered
#' factor is refused rather than fitted with polynomial contrasts, which would
#' produce `.L` and `.Q` terms instead of one hazard ratio per level.
#'
#' Only main effects of existing columns are supported. To fit an interaction
#' or a stratified model, build the column first, e.g.
#' `data$site_age <- interaction(data$site, data$age_grp, sep = " / ")`.
#'
#' # Input checks
#'
#' The call stops, with the offending columns, values or counts listed, when
#'
#' 1. `covariates` is empty, lists a column twice, or lists the `time` or the
#'    `event` column;
#' 2. a required column is not in `data`, or appears in `data` more than once;
#' 3. `time` is not numeric, or holds missing, negative or infinite values, or
#'    `event` holds a value other than `0` and `1`, missing values included.
#'    Rows are never dropped or imputed on your behalf for an outcome
#'    variable;
#' 4. a covariate is entirely missing, which would empty the cohort at the
#'    complete-case step and produce an error that looks like something else;
#' 5. a covariate is of a class that cannot enter a model, is an ordered
#'    factor, or carries a non-default `contrasts` attribute;
#' 6. a categorical covariate has fewer than two levels with rows in it, or a
#'    continuous covariate is constant. The message says how many empty factor
#'    levels were found, because a covariate that looks like it has three
#'    groups and has data in only one is usually a filtering mistake;
#' 7. an empty factor level is the reference level. Dropping it would silently
#'    move the reference to another level and change every hazard ratio of
#'    that covariate, so `fit_cox()` stops and asks you to relevel. Other
#'    empty levels are dropped with a warning, or with an error when
#'    `empty_levels = "error"`;
#' 8. `na_action = "fail"` and any row has a missing covariate value;
#' 9. casewise deletion leaves no rows, no events, a level of a covariate with
#'    no rows at all, or fewer events than the model has coefficients. A level
#'    that survives with rows but no events, and a model with fewer than
#'    `epv_warn` events per coefficient, are warnings: the estimates exist but
#'    should not be read as if they were stable;
#' 10. the fitted model has a missing coefficient, which is what a collinear
#'    pair of covariates or an empty design column produces. Such a model
#'    reports hazard ratios for some covariates and silent blanks for others,
#'    so it is refused rather than returned.
#'
#' @param data A data frame (or tibble) with one row per subject.
#' @param covariates Names of the covariates, as a character vector, in the
#'   order they are to appear in the result. Main effects of existing columns
#'   only; must not contain `time` or `event`, and must not repeat a name.
#'   Factors and character and logical columns are fitted as categorical, with
#'   the first level as reference; numeric columns are fitted as continuous,
#'   so that the hazard ratio is per one unit of that column.
#' @param time Name of the follow-up-time column, as a single string. Must be
#'   numeric, non-missing and non-negative. Factors and character digits are
#'   refused: use [prep_surv()], which converts them and returns a numeric
#'   `time_mo`.
#' @param event Name of the event-indicator column, as a single string. Must
#'   hold only `0` (censored) and `1` (event); logical `FALSE`/`TRUE` is
#'   accepted and read as `0`/`1`.
#' @param labels Labels for the covariates: a named character vector or named
#'   list, e.g. `c(site = "Primary site")`. Covariates with no entry are
#'   labelled with their column name. Entries for columns not in `covariates`
#'   are ignored, so one project-wide label list can be passed to every table.
#' @param ties Handling of tied event times, passed to [survival::coxph()] as
#'   `method`: `"efron"` (the default, and the default of `coxph()` itself),
#'   `"breslow"` or `"exact"`.
#' @param conf_level Confidence level of the hazard-ratio interval. Defaults
#'   to `0.95`. The interval is the Wald interval on the log scale,
#'   `exp(coef +/- z * se)`, which is what [survival::coxph()] itself reports.
#' @param na_action What to do with rows that have a missing covariate value:
#'   `"omit"` (the default) to fit on the complete cases and report how many
#'   rows and events that lost, or `"fail"` to stop. There is no setting that
#'   imputes. A missing time or event value is an error either way.
#' @param empty_levels What to do with a factor level that no row has:
#'   `"drop"` (the default) to drop it with a warning naming it, or `"error"`
#'   to stop. An empty *reference* level is an error under both settings,
#'   because dropping it would move the reference silently.
#' @param zph_transform Function of time used by [survival::cox.zph()], passed
#'   as `transform`: `"km"` (the default, and the default of `cox.zph()`
#'   itself), `"rank"`, `"identity"` or `"log"`. It changes the
#'   proportional-hazards test and nothing else. `"log"` is not finite when an
#'   event happens at time 0, and `cox.zph()` then returns `NaN` rather than
#'   an error; `fit_cox()` stops instead of returning hazard ratios with an
#'   empty test.
#' @param ph_alpha Significance level used only to set the `ph_violated_var`
#'   and `ph_violated_global` flags. Defaults to `0.05`. It has no effect on
#'   the model, which is returned exactly as fitted whether the test passes or
#'   not.
#' @param epv_warn Smallest number of events per model coefficient that does
#'   not draw a warning. Defaults to `10`, the usual rule of thumb; `0`
#'   silences the warning. A model with fewer events than coefficients is an
#'   error, and this argument does not affect that.
#' @param digits_hr Decimal places used when the hazard ratio and its interval
#'   are pasted into `hr_txt`. The numeric columns are returned unrounded, so
#'   rounding stays a presentation choice and never a stored one.
#' @param digits_p Decimal places used for `p_fmt` and for the formatted
#'   p-values of the `cox_zph` attribute. A p-value smaller than
#'   `10^-digits_p` is written `"<0.001"` rather than rounded to `"0.000"`.
#'
#' @return
#' A plain data frame with one row per level of each categorical covariate and
#' one row per continuous covariate, and the columns
#'
#' * `variable`, `var_label`, `level`, `term` -- the covariate, its label, the
#'   level (`NA` for a continuous covariate) and the coefficient name
#'   `coxph()` gives that level, `paste0(variable, level)`;
#' * `is_reference` -- `TRUE` on the first level of a factor;
#' * `n`, `n_event` -- subjects and events in that level among the rows the
#'   model was fitted on; for a continuous covariate, the model totals;
#' * `hr`, `ci_low`, `ci_high`, `hr_txt` -- the hazard ratio, its confidence
#'   interval and the formatted `"hr (low-high)"`. On a reference row `hr` is
#'   `1`, the limits are `NA` and `hr_txt` reads `"1.00 (reference)"`;
#' * `coef`, `se_coef`, `z`, `p_value`, `p_fmt` -- the log hazard ratio, its
#'   standard error, the Wald statistic and the p-value, raw and formatted;
#' * `zph_chisq_var`, `zph_df_var`, `zph_p_var`, `ph_violated_var` -- the
#'   proportional-hazards test of that covariate, repeated on each of its
#'   rows;
#' * `zph_chisq_global`, `zph_df_global`, `zph_p_global`,
#'   `ph_violated_global` -- the global proportional-hazards test, repeated on
#'   every row;
#' * `n_model`, `n_event_model` -- subjects and events in the fitted model.
#'
#' Two attributes are attached:
#'
#' * `cox_zph` -- a data frame with one row per covariate plus the `GLOBAL`
#'   row and the columns `variable`, `var_label`, `term`, `is_global`,
#'   `chisq`, `df`, `p`, `p_fmt`, `violates`, `alpha`, `transform`, `n_model`
#'   and `n_event_model`;
#' * `fit_cox` -- a list with `counts` (rows read, rows used, rows and events
#'   dropped, events, censored, coefficients and events per coefficient),
#'   `missing` (a data frame of missing values per model column), `fit`
#'   (concordance and its standard error, the likelihood ratio, Wald and score
#'   tests, the log-likelihoods, and the subjects and events), `ph` (the
#'   global test, which covariates were flagged, and a note stating that a
#'   flagged model was not modified), `checks` (the invariants that were
#'   enforced), `settings` (including `reference`, the reference level of each
#'   covariate) and `call`.
#'
#' Attributes are dropped by most data-frame verbs, so read them off the
#' object returned by `fit_cox()` before piping it further. The
#' proportional-hazards test is deliberately in the columns as well, so that
#' it survives that.
#'
#' @seealso [prep_surv()], which builds `time_mo` and `event_os` from a
#'   vital-status column; [fit_km()], for the same comparison without the
#'   proportional-hazards assumption.
#'
#' @examples
#' ## A small simulated cohort. No real patient records are used anywhere in
#' ## this package.
#' set.seed(20260901)
#' n    <- 300
#' site <- factor(rep(c("Stomach", "Small intestine", "Colorectal"),
#'                    each = 100),
#'                levels = c("Stomach", "Small intestine", "Colorectal"))
#' age_grp <- factor(sample(c("<50", "50-69", ">=70"), n, replace = TRUE),
#'                   levels = c("<50", "50-69", ">=70"))
#' sex <- factor(sample(c("Female", "Male"), n, replace = TRUE),
#'               levels = c("Female", "Male"))
#' lp <- c(Stomach = 0, "Small intestine" = 0.4,
#'         Colorectal = 0.8)[as.character(site)] +
#'       c("<50" = 0, "50-69" = 0.3, ">=70" = 0.9)[as.character(age_grp)]
#' t_event  <- stats::rexp(n, 0.02 * exp(lp))
#' t_censor <- stats::runif(n, 6, 120)
#' toy <- data.frame(
#'   site     = site,
#'   age_grp  = age_grp,
#'   sex      = sex,
#'   time_mo  = round(pmin(t_event, t_censor), 1),
#'   event_os = as.integer(t_event <= t_censor)
#' )
#'
#' ## One covariate: hazard ratios against the first level of the factor.
#' cox1 <- fit_cox(toy, "site")
#' cox1[, c("variable", "level", "n", "n_event", "hr_txt", "p_fmt")]
#'
#' ## Every column the table carries.
#' names(cox1)
#'
#' ## A multivariable model.
#' cox3 <- fit_cox(toy, c("site", "age_grp", "sex"),
#'                 labels = c(site = "Primary site", age_grp = "Age group"))
#' cox3[, c("var_label", "level", "is_reference", "hr", "ci_low", "ci_high",
#'          "p_fmt")]
#'
#' ## The proportional-hazards test is part of the result, in columns that
#' ## survive the data-frame verbs that drop attributes.
#' unique(cox3[, c("variable", "zph_chisq_var", "zph_df_var", "zph_p_var",
#'                 "ph_violated_var")])
#' attr(cox3, "cox_zph")
#' attr(cox3, "fit_cox")$ph
#'
#' ## Goodness of fit, sample-size bookkeeping and the invariants checked.
#' attr(cox3, "fit_cox")$fit
#' attr(cox3, "fit_cox")$counts
#' attr(cox3, "fit_cox")$checks
#' attr(cox3, "fit_cox")$settings$reference
#'
#' ## Missing covariate values are deleted casewise and counted, never
#' ## imputed; the per-variable table is the one a methods section quotes.
#' gappy <- toy
#' gappy$sex[1:20] <- NA
#' attr(fit_cox(gappy, c("site", "sex")), "fit_cox")$counts
#' attr(fit_cox(gappy, c("site", "sex")), "fit_cox")$missing
#'
#' ## Or refuse to lose them.
#' try(fit_cox(gappy, c("site", "sex"), na_action = "fail"))
#'
#' ## A misspelled column name is named in the error message.
#' try(fit_cox(toy, c("site", "grade")))
#'
#' ## A covariate with one observed level cannot be a covariate.
#' try(fit_cox(toy[toy$site == "Stomach", ], "site"))
#'
#' ## The time column cannot also be a covariate.
#' try(fit_cox(toy, c("site", "time_mo")))
#'
#' ## Two copies of the same covariate give a missing coefficient, which is
#' ## refused rather than returned half-filled.
#' collinear <- toy
#' collinear$site_copy <- collinear$site
#' try(fit_cox(collinear, c("site", "site_copy")))
#'
#' ## An ordered factor would be fitted with polynomial contrasts, so it is
#' ## refused instead.
#' ord <- toy
#' ord$age_grp <- factor(ord$age_grp, ordered = TRUE)
#' try(fit_cox(ord, "age_grp"))
#'
#' @importFrom survival Surv
#' @export
fit_cox <- function(data,
                    covariates,
                    time          = "time_mo",
                    event         = "event_os",
                    labels        = NULL,
                    ties          = c("efron", "breslow", "exact"),
                    conf_level    = 0.95,
                    na_action     = c("omit", "fail"),
                    empty_levels  = c("drop", "error"),
                    zph_transform = c("km", "rank", "identity", "log"),
                    ph_alpha      = 0.05,
                    epv_warn      = 10,
                    digits_hr     = 2L,
                    digits_p      = 3L) {

  cl <- match.call()

  ## -- 0. data -------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("fit_cox(): `data` must be a data frame, got an object of class ",
         .ps_values(class(data)), ".", call. = FALSE)
  }
  n_in <- nrow(data)
  if (n_in == 0L) {
    stop("fit_cox(): `data` has 0 rows, there is nothing to fit.",
         call. = FALSE)
  }

  ## -- 1. the argument values themselves -----------------------------------
  if (missing(covariates)) {
    stop(paste0("fit_cox(): `covariates` must be supplied; fit_cox() has no ",
                "default set of covariates. Pass the column names, e.g. ",
                'covariates = c("site", "age_grp").'),
         call. = FALSE)
  }
  .cx_name_arg(time,  "time")
  .cx_name_arg(event, "event")
  if (!is.character(covariates) || length(covariates) == 0L) {
    stop(sprintf(
      paste0("fit_cox(): `covariates` must be a character vector of at least ",
             "one column name, got %s of length %d. A Cox model with no ",
             "covariates has nothing to estimate."),
      .ps_values(class(covariates)), length(covariates)),
      call. = FALSE)
  }
  bad_v <- is.na(covariates) | !nzchar(covariates)
  if (any(bad_v)) {
    stop(sprintf(
      paste0("fit_cox(): `covariates` has %d empty or missing element%s ",
             "(position%s %s)."),
      sum(bad_v), if (sum(bad_v) > 1L) "s" else "",
      if (sum(bad_v) > 1L) "s" else "", paste(which(bad_v), collapse = ", ")),
      call. = FALSE)
  }
  dup_v <- unique(covariates[duplicated(covariates)])
  if (length(dup_v) > 0L) {
    stop(sprintf(
      paste0("fit_cox(): `covariates` lists %s more than once. The same ",
             "column entered twice is perfectly collinear with itself, so the ",
             "model would have a missing coefficient. Keep one copy."),
      .ps_values(dup_v)),
      call. = FALSE)
  }
  clash <- covariates[covariates %in% c(time, event)]
  if (length(clash) > 0L) {
    stop(sprintf(
      paste0('fit_cox(): %s %s also the `time` ("%s") or the `event` ("%s") ',
             "column. An outcome variable cannot be a covariate of the model ",
             "that explains it: it would be a perfect predictor of itself. ",
             "Drop it from `covariates`, or pass the outcome columns ",
             "explicitly if the names were meant to be different."),
      .ps_values(clash), if (length(clash) > 1L) "are" else "is", time, event),
      call. = FALSE)
  }
  ties          <- match.arg(ties)
  na_action     <- match.arg(na_action)
  empty_levels  <- match.arg(empty_levels)
  zph_transform <- match.arg(zph_transform)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop(sprintf(
      paste0("fit_cox(): `conf_level` must be a single number strictly ",
             "between 0 and 1, got %s."),
      .ps_values(conf_level)),
      call. = FALSE)
  }
  if (!is.numeric(ph_alpha) || length(ph_alpha) != 1L || is.na(ph_alpha) ||
      ph_alpha <= 0 || ph_alpha >= 1) {
    stop(sprintf(
      paste0("fit_cox(): `ph_alpha` must be a single number strictly between ",
             "0 and 1, got %s."),
      .ps_values(ph_alpha)),
      call. = FALSE)
  }
  if (!is.numeric(epv_warn) || length(epv_warn) != 1L || is.na(epv_warn) ||
      !is.finite(epv_warn) || epv_warn < 0) {
    stop(sprintf(
      paste0("fit_cox(): `epv_warn` must be a single non-negative number (0 ",
             "silences the events-per-coefficient warning), got %s."),
      .ps_values(epv_warn)),
      call. = FALSE)
  }
  .cx_count_arg(digits_hr, "digits_hr", min = 0L)
  .cx_count_arg(digits_p,  "digits_p",  min = 1L)

  ## -- 2. required columns present, and present once -----------------------
  needed <- c(time, event, covariates)
  role   <- c("time", "event", rep("covariates", length(covariates)))
  hit    <- needed %in% names(data)
  if (!all(hit)) {
    stop(sprintf(
      paste0("fit_cox(): %d required column%s not found in `data`:\n  %s\n",
             "  `data` has %d columns. Check the spelling, or pass the column ",
             "names explicitly."),
      sum(!hit), if (sum(!hit) > 1L) "s" else "",
      paste0(sprintf('%s = "%s"', role[!hit], needed[!hit]), collapse = "\n  "),
      ncol(data)),
      call. = FALSE)
  }
  dup_c <- needed[needed %in% names(data)[duplicated(names(data))]]
  if (length(dup_c) > 0L) {
    stop(sprintf(
      paste0("fit_cox(): %s appear%s more than once among the columns of ",
             "`data`; fit_cox() would silently use the first one. Make the ",
             "column names unique first."),
      .ps_values(dup_c), if (length(dup_c) > 1L) "" else "s"),
      call. = FALSE)
  }

  ## -- 3. labels ------------------------------------------------------------
  if (is.null(labels)) labels <- list()
  if (!(is.character(labels) || is.list(labels))) {
    stop(sprintf(
      paste0("fit_cox(): `labels` must be a named character vector or a named ",
             "list, got %s."),
      .ps_values(class(labels))),
      call. = FALSE)
  }
  if (length(labels) > 0L &&
      (is.null(names(labels)) || any(!nzchar(names(labels))))) {
    stop(paste0("fit_cox(): every element of `labels` must be named with the ",
                "column it labels, e.g. ",
                'labels = c(site = "Primary site").'),
         call. = FALSE)
  }
  var_label <- vapply(covariates, function(v) {
    if (!(v %in% names(labels))) return(v)
    lab <- labels[[v]]
    if (!is.character(lab) || length(lab) != 1L || is.na(lab)) {
      stop(sprintf(
        paste0('fit_cox(): the label for "%s" must be a single non-missing ',
               "string, got %s of length %d."),
        v, .ps_values(class(lab)), length(lab)),
        call. = FALSE)
    }
    lab
  }, character(1), USE.NAMES = FALSE)

  ## -- 4. the outcome columns ------------------------------------------------
  tv <- data[[time]]
  if (is.factor(tv) || is.character(tv)) {
    stop(sprintf(
      paste0('fit_cox(): the `time` column "%s" is %s. Factor levels and ',
             "character digits are not numbers, so fit_cox() will not convert ",
             "them for you: convert the column explicitly, or build it with ",
             "prep_surv(), which returns a numeric `time_mo`."),
      time, if (is.factor(tv)) "a factor" else "character"),
      call. = FALSE)
  }
  if (!is.numeric(tv)) {
    stop(sprintf(
      paste0('fit_cox(): the `time` column "%s" must be numeric, got an ',
             "object of class %s."),
      time, .ps_values(class(tv))),
      call. = FALSE)
  }
  tnum   <- as.numeric(tv)
  n_t_na <- sum(is.na(tnum))
  if (n_t_na > 0L) {
    stop(sprintf(
      paste0('fit_cox(): the `time` column "%s" has %d missing value%s out of ',
             "%d rows. A subject with no follow-up time is not a complete ",
             "case under any definition, so fit_cox() neither drops nor ",
             "imputes those rows whatever `na_action` says; decide what to do ",
             "with them before calling it."),
      time, n_t_na, if (n_t_na > 1L) "s" else "", n_in),
      call. = FALSE)
  }
  n_t_inf <- sum(!is.finite(tnum))
  if (n_t_inf > 0L) {
    stop(sprintf(
      paste0('fit_cox(): the `time` column "%s" has %d infinite value%s. A ',
             "risk set cannot be formed at an infinite follow-up time."),
      time, n_t_inf, if (n_t_inf > 1L) "s" else ""),
      call. = FALSE)
  }
  if (any(tnum < 0)) {
    n_neg <- sum(tnum < 0)
    stop(sprintf(
      paste0('fit_cox(): the `time` column "%s" has %d negative value%s ',
             "(minimum %s). A negative follow-up time is not analysable, so ",
             "fit_cox() stops instead of computing on it."),
      time, n_neg, if (n_neg > 1L) "s" else "", format(min(tnum))),
      call. = FALSE)
  }

  evr <- data[[event]]
  if (is.factor(evr) || is.character(evr)) {
    stop(sprintf(
      paste0('fit_cox(): the `event` column "%s" is %s. fit_cox() will not ',
             "guess which level means an event: convert it to 0 (censored) / ",
             "1 (event) explicitly, or build it with prep_surv(), which ",
             "returns a 0/1 `event_os` from a vital-status column.\n",
             "  values present : %s"),
      event, if (is.factor(evr)) "a factor" else "character",
      .ps_counts(evr)),
      call. = FALSE)
  }
  if (is.logical(evr)) {
    ev <- as.integer(evr)
  } else if (is.numeric(evr)) {
    ev <- as.numeric(evr)
  } else {
    stop(sprintf(
      paste0('fit_cox(): the `event` column "%s" must be numeric 0/1 or ',
             "logical, got an object of class %s."),
      event, .ps_values(class(evr))),
      call. = FALSE)
  }
  bad_ev <- !(ev %in% c(0, 1))
  if (any(bad_ev)) {
    stop(sprintf(
      paste0('fit_cox(): the `event` column "%s" has %d of %d row%s whose ',
             "value is neither 0 (censored) nor 1 (event).\n",
             "  offending values : %s\n",
             "  A missing indicator (<NA>) counts as an offending value: ",
             "fit_cox() never guesses whether a subject had the event, and ",
             "never drops the row for it. Use prep_surv() to build a 0/1 ",
             "indicator from a vital-status column."),
      event, sum(bad_ev), n_in, if (n_in > 1L) "s" else "",
      .ps_counts(evr[bad_ev])),
      call. = FALSE)
  }
  ev <- as.integer(ev)

  ## -- 5. the covariates -----------------------------------------------------
  ## Every covariate is copied into an internally named column, so that no
  ## other column of `data` reaches coxph(), and so that a column name that is
  ## not a syntactic R name cannot reach the model formula. The caller's names
  ## are put back on the way out.
  k         <- length(covariates)
  xs        <- vector("list", k)
  lev_list  <- vector("list", k)
  is_cat    <- logical(k)
  n_lv_drop <- integer(k)
  ref_lev   <- rep(NA_character_, k)
  names(n_lv_drop) <- names(ref_lev) <- covariates

  for (i in seq_len(k)) {
    v  <- covariates[i]
    xr <- data[[v]]

    if (!is.null(dim(xr))) {
      stop(sprintf(
        paste0('fit_cox(): the covariate "%s" is not a plain column (it has ',
               "dimensions %s). A matrix column cannot be reported level by ",
               "level; split it into columns first."),
        v, paste(dim(xr), collapse = " x ")),
        call. = FALSE)
    }
    if (all(is.na(xr))) {
      stop(sprintf(
        paste0('fit_cox(): the covariate "%s" is missing in all %d rows. ',
               "Casewise deletion would empty the cohort, and coxph() would ",
               "then fail for a reason that has nothing to do with the real ",
               "cause. Drop it from `covariates`, or clean the column first."),
        v, n_in),
        call. = FALSE)
    }

    if (is.factor(xr)) {
      if (is.ordered(xr)) {
        stop(sprintf(
          paste0('fit_cox(): the covariate "%s" is an ordered factor. R fits ',
                 "ordered factors with polynomial contrasts, which give .L ",
                 "and .Q terms instead of one hazard ratio per level against ",
                 "a reference, so fit_cox() refuses it rather than reporting ",
                 "terms that do not mean what the row labels would say.\n",
                 '  Use data$%s <- factor(data$%s, ordered = FALSE) to fit it ',
                 "as a categorical covariate; the level order, and therefore ",
                 "the reference level, is kept."),
          v, v, v),
          call. = FALSE)
      }
      if (!is.null(attr(xr, "contrasts"))) {
        stop(sprintf(
          paste0('fit_cox(): the covariate "%s" carries a `contrasts` ',
                 "attribute, so its coefficients would not be the levels ",
                 "against the first level. Remove it with ",
                 'attr(data$%s, "contrasts") <- NULL, or recode the column, ',
                 "before calling fit_cox()."),
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
          paste0('fit_cox(): the covariate "%s" has %d infinite value%s. An ',
                 "infinite covariate makes the linear predictor infinite for ",
                 "that subject; recode those values before calling."),
          v, n_x_inf, if (n_x_inf > 1L) "s" else ""),
          call. = FALSE)
      }
      is_cat[i] <- FALSE
    } else {
      stop(sprintf(
        paste0('fit_cox(): the covariate "%s" is of class %s, which cannot ',
               "enter a Cox model. Convert it to a factor (categorical) or to ",
               "a number (continuous) first, so that it is on the record ",
               "which of the two it is."),
        v, .ps_values(class(xr))),
        call. = FALSE)
    }

    if (is_cat[i]) {
      cnt   <- table(x, useNA = "no")
      empty <- names(cnt)[cnt == 0L]
      n_ok  <- sum(cnt > 0L)
      if (n_ok < 2L) {
        stop(sprintf(
          paste0('fit_cox(): the covariate "%s" has %d level%s with rows in ',
                 "it (%s); at least 2 are needed.\n",
                 "  A covariate that is the same for everybody explains ",
                 "nothing and has no hazard ratio to report: its design ",
                 "column would be constant and its coefficient missing.\n",
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
            paste0('fit_cox(): the reference level of the covariate "%s" is ',
                   '"%s", and no row has it.\n',
                   '  Dropping it would move the reference to "%s" without ',
                   "saying so, and every hazard ratio of this covariate would ",
                   "change meaning while the code stayed the same. fit_cox() ",
                   "will not do that.\n",
                   "  Relevel the factor explicitly, e.g. data$%s <- ",
                   'stats::relevel(droplevels(data$%s), ref = "%s").'),
            v, names(cnt)[1], names(cnt)[cnt > 0L][1], v, v,
            names(cnt)[cnt > 0L][1]),
            call. = FALSE)
        }
        msg <- sprintf(
          paste0('the covariate "%s" has %d declared level%s that no row has ',
                 "(%s). An empty level has an all-zero design column and no ",
                 "hazard ratio to report."),
          v, length(empty), if (length(empty) > 1L) "s" else "",
          .ps_values(empty))
        if (empty_levels == "error") {
          stop(sprintf(
            paste0("fit_cox(): %s\n  Drop them with droplevels(), or pass ",
                   'empty_levels = "drop" to have fit_cox() drop them and ',
                   "report how many."),
            msg),
            call. = FALSE)
        }
        warning(sprintf(
          paste0("fit_cox(): %s\n  Dropped, and counted in ",
                 'attr(x, "fit_cox")$settings$levels_dropped. The reference ',
                 'level ("%s") is not affected.'),
          msg, names(cnt)[1]),
          call. = FALSE)
        n_lv_drop[i] <- length(empty)
        x <- droplevels(x)
      }
      lev_list[[i]] <- levels(x)
      ref_lev[i]    <- levels(x)[1]
    } else {
      if (length(unique(x[!is.na(x)])) < 2L) {
        stop(sprintf(
          paste0('fit_cox(): the continuous covariate "%s" takes the same ',
                 "value (%s) in every row that has one. A constant covariate ",
                 "explains nothing and its coefficient would be missing."),
          v, .ps_values(unique(x[!is.na(x)]))),
          call. = FALSE)
      }
    }
    xs[[i]] <- x
  }

  ## -- 6. casewise deletion --------------------------------------------------
  vn <- paste0(".cox_v", seq_len(k))
  md <- data.frame(.cox_time = tnum, .cox_event = ev, stringsAsFactors = FALSE)
  for (i in seq_len(k)) md[[vn[i]]] <- xs[[i]]

  n_miss <- vapply(md, function(z) sum(is.na(z)), integer(1))
  missing_tab <- data.frame(
    variable    = c(time, event, covariates),
    role        = c("time", "event", rep("covariate", k)),
    n_missing   = unname(n_miss),
    pct_missing = 100 * unname(n_miss) / n_in,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  miss_txt <- if (any(missing_tab$n_missing > 0L)) {
    paste(sprintf("%s (n = %d)",
                  missing_tab$variable[missing_tab$n_missing > 0L],
                  missing_tab$n_missing[missing_tab$n_missing > 0L]),
          collapse = ", ")
  } else {
    "none"
  }

  keep       <- stats::complete.cases(md)
  n_used     <- sum(keep)
  n_dropped  <- n_in - n_used
  ev_in      <- sum(ev)
  ev_used    <- sum(ev[keep])
  ev_dropped <- ev_in - ev_used

  if (n_dropped > 0L && na_action == "fail") {
    stop(sprintf(
      paste0('fit_cox(): na_action = "fail" and %d of %d row%s (%.1f%%) have ',
             "a missing covariate value; %d of the %d event%s would go with ",
             "them.\n",
             "  missing by column : %s\n",
             '  Pass na_action = "omit" to fit on the %d complete case%s (the ',
             "count comes back in the diagnostics either way), or handle the ",
             "missing values before calling."),
      n_dropped, n_in, if (n_dropped > 1L) "s" else "",
      100 * n_dropped / n_in, ev_dropped, ev_in,
      if (ev_in == 1L) "" else "s", miss_txt,
      n_used, if (n_used == 1L) "" else "s"),
      call. = FALSE)
  }
  if (n_used == 0L) {
    stop(sprintf(
      paste0("fit_cox(): casewise deletion leaves 0 of %d rows, so there is ",
             "nothing to fit.\n  missing by column : %s"),
      n_in, miss_txt),
      call. = FALSE)
  }
  md <- md[keep, , drop = FALSE]
  row.names(md) <- NULL

  ## -- 7. the complete cases must still support the model --------------------
  n_coef <- sum(vapply(seq_len(k), function(i) {
    if (is_cat[i]) length(lev_list[[i]]) - 1L else 1L
  }, integer(1)))

  if (ev_used == 0L) {
    stop(sprintf(
      paste0("fit_cox(): the %d row%s the model would be fitted on contain 0 ",
             "events. A Cox model estimates nothing without events."),
      n_used, if (n_used > 1L) "s" else ""),
      call. = FALSE)
  }

  no_event_lv <- character(0)
  for (i in seq_len(k)) {
    if (!is_cat[i]) {
      if (length(unique(md[[vn[i]]])) < 2L) {
        stop(sprintf(
          paste0('fit_cox(): after casewise deletion the continuous covariate ',
                 '"%s" is constant in the %d remaining row%s. It was not ',
                 "constant in the %d rows read, so the deletion, not the ",
                 "column, is what emptied it."),
          covariates[i], n_used, if (n_used > 1L) "s" else "", n_in),
          call. = FALSE)
      }
      next
    }
    cnt  <- table(md[[vn[i]]])
    gone <- names(cnt)[cnt == 0L]
    if (length(gone) > 0L) {
      stop(sprintf(
        paste0('fit_cox(): after casewise deletion the covariate "%s" has no ',
               "rows left in %d of its %d levels (%s).\n",
               "  Those levels had rows in the %d rows read; the %d row%s ",
               "deleted for missing covariate values took all of them. ",
               "coxph() would return a missing coefficient for each and ",
               "report the rest as if nothing had happened, so fit_cox() ",
               "stops.\n",
               "  Either drop the covariate that is causing the deletion, or ",
               "collapse the emptied levels, and say in the methods which of ",
               "the two you did."),
        covariates[i], length(gone), length(cnt), .ps_values(gone),
        n_in, n_dropped, if (n_dropped == 1L) "" else "s"),
        call. = FALSE)
    }
    ev_lv_i <- vapply(split(md$.cox_event, md[[vn[i]]]), sum, integer(1))
    if (any(ev_lv_i == 0L)) {
      no_event_lv <- c(
        no_event_lv,
        sprintf('"%s" level%s %s', covariates[i],
                if (sum(ev_lv_i == 0L) > 1L) "s" else "",
                .ps_values(names(ev_lv_i)[ev_lv_i == 0L])))
    }
  }
  if (length(no_event_lv) > 0L) {
    warning(sprintf(
      paste0("fit_cox(): %s ha%s rows but no events.\n",
             "  The hazard ratio of such a level is driven towards 0 with a ",
             "confidence interval that carries no information; read it as ",
             '"no events observed", not as a strong protective effect.'),
      paste(no_event_lv, collapse = "; "),
      if (length(no_event_lv) > 1L) "ve" else "s"),
      call. = FALSE)
  }

  if (ev_used < n_coef) {
    stop(sprintf(
      paste0("fit_cox(): the model has %d coefficient%s and only %d event%s ",
             "in the %d row%s it would be fitted on. A Cox model cannot ",
             "identify more coefficients than it has events; the estimates ",
             "would be arbitrary. Fit fewer covariates, or collapse levels."),
      n_coef, if (n_coef > 1L) "s" else "", ev_used,
      if (ev_used > 1L) "s" else "", n_used, if (n_used > 1L) "s" else ""),
      call. = FALSE)
  }
  epv <- ev_used / n_coef
  if (epv_warn > 0 && epv < epv_warn) {
    warning(sprintf(
      paste0("fit_cox(): %d event%s for %d coefficient%s is %.1f events per ",
             "coefficient, below the %g this call asked about. The model is ",
             "estimable but its hazard ratios and intervals are unstable; ",
             "report the events per coefficient with them, or fit fewer ",
             "covariates. Pass epv_warn = 0 to silence this."),
      ev_used, if (ev_used > 1L) "s" else "", n_coef,
      if (n_coef > 1L) "s" else "", epv, epv_warn),
      call. = FALSE)
  }

  ## -- 8. the fit ------------------------------------------------------------
  ## Treatment contrasts are forced for the duration of the fit, so that the
  ## reference level is the first level whatever options("contrasts") is set to
  ## session-wide, and restored on exit.
  op <- options(contrasts = c(unordered = "contr.treatment",
                              ordered   = "contr.poly"))
  on.exit(options(op), add = TRUE)

  frm <- stats::as.formula(paste0("Surv(.cox_time, .cox_event) ~ ",
                                  paste(vn, collapse = " + ")))
  ## x = TRUE, y = TRUE keeps the design matrix inside the fit, so that
  ## cox.zph() never has to rebuild the model frame from the call. The fitted
  ## object stays inside this function and is not returned.
  fit <- survival::coxph(frm, data = md, method = ties, x = TRUE, y = TRUE)

  cf_int   <- character(0)
  cf_owner <- integer(0)
  for (i in seq_len(k)) {
    nm       <- if (is_cat[i]) paste0(vn[i], lev_list[[i]][-1]) else vn[i]
    cf_int   <- c(cf_int, nm)
    cf_owner <- c(cf_owner, rep(i, length(nm)))
  }
  if (!identical(names(stats::coef(fit)), cf_int)) {
    stop(sprintf(
      paste0("fit_cox(): the coefficients of the fit are not the ones the ",
             "covariates imply, which should be impossible by ",
             "construction.\n  coxph()  : %s\n  expected : %s\n",
             "  Do not use this result; report it as a bug in fit_cox()."),
      .ps_values(names(stats::coef(fit))), .ps_values(cf_int)),
      call. = FALSE)
  }
  if (fit$n != n_used || fit$nevent != ev_used) {
    stop(sprintf(
      paste0("fit_cox(): coxph() used %d rows and %d events where the ",
             "complete cases are %d and %d, which should be impossible by ",
             "construction. Do not use this result; report it as a bug in ",
             "fit_cox()."),
      fit$n, fit$nevent, n_used, ev_used),
      call. = FALSE)
  }
  na_cf <- is.na(stats::coef(fit))
  if (any(na_cf)) {
    stop(sprintf(
      paste0("fit_cox(): the fitted model has %d missing coefficient%s, in ",
             "%s.\n",
             "  coxph() returns NA for a coefficient it cannot estimate: two ",
             "covariates carrying the same information, a level that no ",
             "longer has rows, or a covariate that is a function of the ",
             "others. The rest of the model is still reported as if nothing ",
             "were wrong, so fit_cox() refuses to return it.\n",
             "  Drop or combine the covariates concerned and fit again."),
      sum(na_cf), if (sum(na_cf) > 1L) "s" else "",
      .ps_values(unique(covariates[cf_owner[na_cf]]))),
      call. = FALSE)
  }

  ## -- 9. cox.zph: computed on every call, and never optional ----------------
  zph <- tryCatch(
    survival::cox.zph(fit, transform = zph_transform, terms = TRUE,
                      global = TRUE),
    error = function(e) {
      stop(sprintf(
        paste0("fit_cox(): cox.zph() failed, so the proportional-hazards ",
               "assumption behind every hazard ratio of this model is ",
               "untested.\n  cox.zph() said : %s\n",
               "  fit_cox() does not return hazard ratios without that test, ",
               "so the call stops here. Try another `zph_transform` (this ",
               'call used "%s"), or check the model.'),
        conditionMessage(e), zph_transform),
        call. = FALSE)
    })

  zt <- zph$table
  if (!identical(rownames(zt), c(vn, "GLOBAL")) ||
      !identical(colnames(zt), c("chisq", "df", "p"))) {
    stop(sprintf(
      paste0("fit_cox(): cox.zph() returned a table this version of fit_cox() ",
             "does not recognise, so the proportional-hazards test cannot be ",
             "matched to the covariates.\n  rows    : %s\n  columns : %s\n",
             "  Do not use this result; report it as a bug in fit_cox()."),
      .ps_values(rownames(zt)), .ps_values(colnames(zt))),
      call. = FALSE)
  }

  zph_var_chisq <- unname(zt[vn, "chisq"])
  zph_var_df    <- unname(zt[vn, "df"])
  zph_var_p     <- unname(zt[vn, "p"])
  zph_g_chisq   <- unname(zt["GLOBAL", "chisq"])
  zph_g_df      <- unname(zt["GLOBAL", "df"])
  zph_g_p       <- unname(zt["GLOBAL", "p"])
  ## A missing test is not a passed test. cox.zph() returns NaN rather than an
  ## error when the transform of the event times is not finite -- transform =
  ## "log" with an event at time 0 is the case that occurs in practice -- and
  ## the hazard ratios would then be returned with an empty assumption check.
  zph_all_p <- c(zph_var_p, zph_g_p)
  if (anyNA(zph_all_p) || anyNA(c(zph_var_chisq, zph_g_chisq))) {
    n_t0 <- sum(md$.cox_time <= 0 & md$.cox_event == 1L)
    stop(sprintf(
      paste0("fit_cox(): cox.zph() returned no proportional-hazards test for ",
             '%s: the chi-square or the p-value is missing with transform = ',
             '"%s".\n',
             "  A missing test is not a passed test, and fit_cox() does not ",
             "return hazard ratios whose proportional-hazards assumption was ",
             "never actually tested.%s\n",
             '  Use the default zph_transform = "km", which does not ',
             "transform the event times themselves, or check the follow-up ",
             "times."),
      .ps_values(c(covariates, "GLOBAL")[is.na(zph_all_p)]), zph_transform,
      if (zph_transform == "log" && n_t0 > 0L)
        sprintf(paste0("\n  transform = \"log\" takes the logarithm of the ",
                       "event times, which is not finite at time 0; %d of the ",
                       "%d event%s in this model happened at time 0."),
                n_t0, ev_used, if (ev_used > 1L) "s" else "")
      else ""),
      call. = FALSE)
  }
  viol_var      <- !is.na(zph_var_p) & zph_var_p < ph_alpha
  viol_global   <- !is.na(zph_g_p)   & zph_g_p   < ph_alpha

  cox_zph <- data.frame(
    variable      = c(covariates, NA_character_),
    var_label     = c(var_label,  NA_character_),
    term          = c(covariates, "GLOBAL"),
    is_global     = c(rep(FALSE, k), TRUE),
    chisq         = c(zph_var_chisq, zph_g_chisq),
    df            = c(zph_var_df,    zph_g_df),
    p             = c(zph_var_p,     zph_g_p),
    p_fmt         = .cx_fmt_p(c(zph_var_p, zph_g_p), digits_p),
    violates      = c(viol_var, viol_global),
    alpha         = ph_alpha,
    transform     = zph$transform,
    n_model       = n_used,
    n_event_model = ev_used,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  ## -- 10. the coefficient table ---------------------------------------------
  s   <- summary(fit)
  cfm <- s$coefficients
  z_q <- stats::qnorm((1 + conf_level) / 2)

  var_i <- unlist(lapply(seq_len(k), function(i)
    rep(i, if (is_cat[i]) length(lev_list[[i]]) else 1L)), use.names = FALSE)
  level <- unlist(lapply(seq_len(k), function(i)
    if (is_cat[i]) lev_list[[i]] else NA_character_), use.names = FALSE)
  is_ref <- unlist(lapply(seq_len(k), function(i)
    if (is_cat[i]) seq_along(lev_list[[i]]) == 1L else FALSE),
    use.names = FALSE)

  key_int <- ifelse(is.na(level), vn[var_i], paste0(vn[var_i], level))
  key_int[is_ref] <- NA_character_
  j <- match(key_int, rownames(cfm))

  n_lv  <- integer(length(var_i))
  ev_lv <- integer(length(var_i))
  for (r in seq_along(var_i)) {
    i <- var_i[r]
    if (is_cat[i]) {
      in_lv    <- md[[vn[i]]] == level[r]
      n_lv[r]  <- sum(in_lv)
      ev_lv[r] <- sum(md$.cox_event[in_lv])
    } else {
      n_lv[r]  <- n_used
      ev_lv[r] <- ev_used
    }
  }

  coef_v <- ifelse(is.na(j), NA_real_, cfm[j, "coef"])
  se_v   <- ifelse(is.na(j), NA_real_, cfm[j, "se(coef)"])
  z_v    <- ifelse(is.na(j), NA_real_, cfm[j, "z"])
  p_v    <- ifelse(is.na(j), NA_real_, cfm[j, "Pr(>|z|)"])
  hr_v   <- ifelse(is_ref, 1, exp(coef_v))
  lo_v   <- exp(coef_v - z_q * se_v)
  hi_v   <- exp(coef_v + z_q * se_v)
  fmt_hr <- function(x) sprintf("%.*f", digits_hr, x)

  out <- data.frame(
    variable           = covariates[var_i],
    var_label          = var_label[var_i],
    level              = level,
    term               = ifelse(is.na(level), covariates[var_i],
                                paste0(covariates[var_i], level)),
    is_reference       = is_ref,
    n                  = n_lv,
    n_event            = ev_lv,
    hr                 = hr_v,
    ci_low             = lo_v,
    ci_high            = hi_v,
    hr_txt             = ifelse(
      is_ref, paste0(fmt_hr(1), " (reference)"),
      sprintf("%s (%s-%s)", fmt_hr(hr_v), fmt_hr(lo_v), fmt_hr(hi_v))),
    coef               = coef_v,
    se_coef            = se_v,
    z                  = z_v,
    p_value            = p_v,
    p_fmt              = .cx_fmt_p(p_v, digits_p),
    zph_chisq_var      = zph_var_chisq[var_i],
    zph_df_var         = zph_var_df[var_i],
    zph_p_var          = zph_var_p[var_i],
    ph_violated_var    = viol_var[var_i],
    zph_chisq_global   = zph_g_chisq,
    zph_df_global      = zph_g_df,
    zph_p_global       = zph_g_p,
    ph_violated_global = viol_global,
    n_model            = n_used,
    n_event_model      = ev_used,
    stringsAsFactors   = FALSE,
    row.names          = NULL
  )

  ## -- 11. diagnostics -------------------------------------------------------
  sum_by_var <- function(x) {
    keep_cat <- !is.na(out$level)
    if (!any(keep_cat)) return(numeric(0))
    vapply(split(as.numeric(x[keep_cat]), out$variable[keep_cat]), sum,
           numeric(1))
  }
  inv <- c(
    "one row per level of every covariate" = nrow(out) == length(var_i),
    "every non-reference row has a hazard ratio" =
      !anyNA(out$hr[!out$is_reference]),
    "every reference row has a hazard ratio of 1 and no p-value" =
      all(out$hr[out$is_reference] == 1) &&
      all(is.na(out$p_value[out$is_reference])),
    "the levels of each covariate add up to the rows fitted" =
      all(sum_by_var(out$n) == n_used),
    "the levels of each covariate add up to the events fitted" =
      all(sum_by_var(out$n_event) == ev_used),
    "subjects and events agree with coxph()" =
      fit$n == n_used && fit$nevent == ev_used,
    "no missing coefficient" = !anyNA(stats::coef(fit)),
    "cox.zph() reported one test per covariate plus the global test" =
      nrow(cox_zph) == k + 1L,
    "the proportional-hazards test is on every row" =
      !anyNA(out$zph_p_var) && !anyNA(out$zph_p_global))
  if (!all(inv)) {
    stop(sprintf(
      paste0("fit_cox(): the returned table does not add up, which should be ",
             "impossible by construction.\n  failed check%s : %s\n",
             "  Do not use this result; report it as a bug in fit_cox()."),
      if (sum(!inv) > 1L) "s" else "",
      paste(names(inv)[!inv], collapse = "; ")),
      call. = FALSE)
  }

  ph_note <- if (viol_global || any(viol_var)) {
    sprintf(
      paste0("the proportional-hazards assumption is rejected at alpha = %g ",
             "(global p = %s; %s). The model is returned exactly as fitted: ",
             "fit_cox() does not stratify, add a time-dependent coefficient ",
             "or shorten follow-up on its own, because each of those answers ",
             "a different question. A hazard ratio from this model is an ",
             "average over follow-up of a ratio that is changing."),
      ph_alpha, .cx_fmt_p(zph_g_p, digits_p),
      if (any(viol_var)) paste0("flagged: ", .ps_values(covariates[viol_var]))
      else "no single covariate flagged")
  } else {
    sprintf(
      paste0("the proportional-hazards assumption is not rejected at alpha = ",
             "%g (global p = %s, every covariate p >= %g). Not rejecting is ",
             "not the same as holding."),
      ph_alpha, .cx_fmt_p(zph_g_p, digits_p), ph_alpha)
  }

  attr(out, "cox_zph") <- cox_zph
  attr(out, "fit_cox") <- list(
    counts = c(rows_in                = n_in,
               rows_used              = n_used,
               rows_dropped           = n_dropped,
               events_in              = ev_in,
               events_used            = ev_used,
               events_dropped         = ev_dropped,
               censored_used          = n_used - ev_used,
               coefficients           = n_coef,
               events_per_coefficient = epv),
    missing = missing_tab,
    fit = c(concordance    = unname(s$concordance[["C"]]),
            concordance_se = unname(s$concordance[["se(C)"]]),
            lrt_chisq      = unname(s$logtest[["test"]]),
            lrt_df         = unname(s$logtest[["df"]]),
            lrt_p          = unname(s$logtest[["pvalue"]]),
            wald_chisq     = unname(s$waldtest[["test"]]),
            wald_df        = unname(s$waldtest[["df"]]),
            wald_p         = unname(s$waldtest[["pvalue"]]),
            score_chisq    = unname(s$sctest[["test"]]),
            score_df       = unname(s$sctest[["df"]]),
            score_p        = unname(s$sctest[["pvalue"]]),
            loglik_null    = unname(fit$loglik[1]),
            loglik_model   = unname(fit$loglik[2]),
            n              = n_used,
            n_event        = ev_used),
    ph = list(test            = "cox.zph (Grambsch-Therneau)",
              transform       = zph$transform,
              global_chisq    = zph_g_chisq,
              global_df       = zph_g_df,
              global_p        = zph_g_p,
              alpha           = ph_alpha,
              violated_global = viol_global,
              violated_vars   = covariates[viol_var],
              note            = ph_note),
    checks   = inv,
    settings = list(
      time           = time,
      event          = event,
      covariates     = covariates,
      labels         = stats::setNames(var_label, covariates),
      reference      = stats::setNames(ref_lev, covariates),
      levels         = stats::setNames(lev_list, covariates),
      levels_dropped = n_lv_drop,
      ties           = ties,
      conf_level     = conf_level,
      na_action      = na_action,
      empty_levels   = empty_levels,
      zph_transform  = zph_transform,
      ph_alpha       = ph_alpha,
      epv_warn       = epv_warn,
      digits_hr      = digits_hr,
      digits_p       = digits_p,
      model          = paste0("Surv(", time, ", ", event, ") ~ ",
                              paste(covariates, collapse = " + ")),
      estimator      = paste0("Cox proportional hazards (survival::coxph, ",
                              "method = \"", ties, "\")"),
      ci             = "Wald on the log scale, exp(coef +/- z * se)",
      ph_test        = paste0("survival::cox.zph, transform = \"",
                              zph_transform, "\", terms = TRUE"),
      contrasts      = "contr.treatment, reference = first level"),
    call = cl)
  out
}


# -- internal helpers, not exported ------------------------------------------
# .ps_values() and .ps_counts() are defined in R/prep_surv.R and reused here.

#' @noRd
.cx_name_arg <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf(
      "fit_cox(): `%s` must be a single non-missing string.", arg),
      call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.cx_count_arg <- function(x, arg, min = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != as.integer(x) ||
      x < min) {
    stop(sprintf(
      "fit_cox(): `%s` must be a single whole number >= %d, got %s.",
      arg, min, .ps_values(x)),
      call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.cx_fmt_p <- function(p, digits) {
  thr <- 10^(-digits)
  ifelse(is.na(p), NA_character_,
         ifelse(p < thr,
                paste0("<", formatC(thr, format = "f", digits = digits)),
                formatC(p, format = "f", digits = digits)))
}
