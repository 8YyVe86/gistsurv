# ---------------------------------------------------------------------------
# fit_km(): Kaplan-Meier estimates for one optional grouping variable.
#
# Design notes (deliberate, and different from the study script this was
# distilled from):
#   * no column name is hard-coded: the time, the event indicator and the
#     grouping variable are all arguments;
#   * nothing is printed, nothing is plotted and nothing is written to disk:
#     the per-stratum estimates come back as a data frame, and the time-point
#     estimates and the diagnostics as attributes;
#   * a time point past the end of a stratum's follow-up is an error, not a
#     silently extrapolated number. summary.survfit(extend = TRUE) carries the
#     last value forward without saying so; the rule enforced here is the one
#     the RMST tau guard uses, tau <= the smallest per-stratum maximum
#     follow-up;
#   * a grouping variable with fewer than two observed levels is an error, and
#     an unused factor level is counted rather than quietly dropped;
#   * only the requested time points are returned, never the whole step
#     function, and no fitted object is attached: a full Kaplan-Meier table
#     has one step per distinct event time and can be turned back into
#     individual follow-up times, and a survfit object drags the caller's data
#     along inside the environment of its formula.
# The file is kept ASCII-only so that it behaves the same under a UTF-8 and
# under a non-UTF-8 locale.
#
# The small formatting helpers .ps_values() and .ps_counts() are shared with
# R/prep_surv.R; they are internal to the package.
# ---------------------------------------------------------------------------

#' Kaplan-Meier estimates by an optional grouping variable
#'
#' @description
#' `fit_km()` fits a Kaplan-Meier curve, optionally within the levels of one
#' grouping variable, and returns the numbers a survival paper reports: the
#' subjects and events per stratum, the median survival time with its
#' confidence interval, the survival probability and the number still at risk
#' at a set of time points, the log-rank test across strata, the
#' reverse-Kaplan-Meier median follow-up, the largest follow-up time in each
#' stratum, and the total person-time.
#'
#' The time, event and grouping columns are arguments, so the function is not
#' tied to any particular registry export. The function prints nothing, plots
#' nothing and writes nothing: the per-stratum estimates come back as a plain
#' data frame, and the time-point estimates, the log-rank test and the
#' remaining diagnostics are attached to it as attributes.
#'
#' @details
#' # What is estimated
#'
#' One row per stratum, in the order of the levels of `group`, or a single row
#' labelled `"ALL"` when `group` is `NULL`. Each row carries the number of
#' subjects and events, the median survival time with its confidence interval,
#' the median follow-up, the largest follow-up time and the person-time of
#' that stratum, and the log-rank test, which is a property of the whole
#' comparison and is therefore repeated on every row.
#'
#' The survival probability, its confidence interval and the number still at
#' risk at each element of `times` are a second, longer table: one row per
#' stratum and time point. It travels back as the `km_times` attribute, so
#' that the two tables are always two halves of the same fit.
#'
#' Only the requested time points are returned, never the whole step function.
#' That is partly tidiness and partly disclosure control: a full Kaplan-Meier
#' table has one step per distinct event time and can be turned back into
#' individual follow-up times. For the same reason no `survfit` object is
#' attached to the result -- a fitted object carries the caller's data along
#' inside the environment of its formula, and would put individual records
#' into anything the result is saved into. Every argument that defines the fit
#' is recorded in the `settings` diagnostic instead, so the same fit can be
#' reproduced exactly.
#'
#' # A median that is not reached
#'
#' The median is the first time the estimated curve falls to 0.5. In a stratum
#' whose curve never gets there it does not exist: `median` is `NA`,
#' `median_reached` is `FALSE` and `median_ci_txt` reads `"not reached"`. It is
#' never replaced by the largest follow-up time, by the end of the x axis of a
#' figure, or by a restricted mean; those are different quantities, and each
#' of them would understate the survival of that stratum.
#'
#' # Confidence intervals
#'
#' `conf_type` is passed to [survival::survfit()] and decides the confidence
#' interval of both the median and the time-point survival probabilities. The
#' default, `"log"`, is the default of `survfit()` itself, and therefore the
#' interval that a script calling `survfit()` without naming `conf.type` has
#' been reporting all along. Changing it changes published intervals without
#' changing a single point estimate, so the value used is recorded in the
#' `settings` diagnostic and is never chosen on your behalf.
#'
#' # Time points, and the follow-up guard
#'
#' A Kaplan-Meier curve is not estimable past the largest follow-up time in a
#' stratum. `summary.survfit(extend = TRUE)` will nevertheless carry the last
#' value forward, silently, and the resulting row looks exactly like an
#' estimate. `fit_km()` refuses to do that: a time point past the largest
#' follow-up of any stratum stops the call, and `times_beyond = "na"` keeps the
#' row but marks it `estimable = FALSE` and leaves the survival probability and
#' its interval missing. The number at risk is still reported for such a row,
#' because an empty risk set is a fact rather than an extrapolation.
#'
#' This is the same rule as the tau guard of a restricted-mean analysis: the
#' largest usable time point is the smallest of the per-stratum maxima, and it
#' comes back as `follow_up["tau_upper"]` in the diagnostics.
#'
#' # Median follow-up
#'
#' `median_fu_revkm` is the reverse-Kaplan-Meier median follow-up: the event
#' indicator is reversed (`1 - event`, so that a censored observation counts as
#' the event and a death is censored) and a Kaplan-Meier median is taken from
#' that fit. Unlike the median of the observed follow-up times it is not pulled
#' down by early deaths, which is why it is the figure to quote as "median
#' follow-up". It is `NA` when fewer than half the subjects are censored, which
#' is the honest answer: the reversed curve never reaches 0.5.
#'
#' # The log-rank test
#'
#' The test is [survival::survdiff()] with `rho = 0`, the log-rank test, on
#' `length(levels) - 1` degrees of freedom, with the p-value taken from the
#' chi-square distribution. With `group = NULL` there is one stratum, nothing
#' to compare and therefore no test: `logrank_chisq`, `logrank_df` and
#' `logrank_p` are `NA`, and the reason is spelled out in
#' `attr(x, "fit_km")$logrank$note` rather than left as an unexplained blank
#' column.
#'
#' The test compares whole curves and answers a different question from the
#' medians and from the time-point probabilities in the same result; a small
#' p-value next to overlapping medians is not a contradiction.
#'
#' # Input checks
#'
#' The call stops, with the offending columns, values or counts listed, when
#'
#' 1. a required column (`time`, `event`, `group`) is not in `data`, or
#'    appears in `data` more than once;
#' 2. `time` is not numeric, or holds missing, negative or infinite values.
#'    Rows are never dropped or imputed on your behalf;
#' 3. `event` holds a value other than `0` and `1`, missing values included.
#'    The offending values are listed with their counts. Use [prep_surv()] to
#'    build a `0`/`1` indicator from a vital-status column;
#' 4. `group` has missing values and `group_missing = "error"` (the default).
#'    The other settings keep those rows as their own stratum or drop them, and
#'    the number of rows involved is reported either way -- rows are never
#'    dropped without a count;
#' 5. `group` has fewer than two observed levels. There is then one curve and
#'    no log-rank test; the message says how many empty factor levels were
#'    dropped, because a grouping variable that looks like it has three groups
#'    and has data in only one is usually a filtering mistake rather than an
#'    intention;
#' 6. an element of `times` is past the largest follow-up time of some stratum
#'    and `times_beyond = "error"` (the default);
#' 7. `times` holds duplicated, missing, negative or infinite values, or
#'    `group` is of a class that cannot be a grouping variable.
#'
#' @param data A data frame (or tibble) with one row per subject.
#' @param time Name of the follow-up-time column, as a single string. Must be
#'   numeric. Factors and character digits are refused: use [prep_surv()],
#'   which converts them and returns a numeric `time_mo`.
#' @param event Name of the event-indicator column, as a single string. Must
#'   hold only `0` (censored) and `1` (event); logical `FALSE`/`TRUE` is
#'   accepted and read as `0`/`1`.
#' @param group Name of the grouping variable, as a single string, or `NULL`
#'   (the default) for a single curve. One stratum is produced per level, in
#'   the order of the levels of the factor. To cross two variables, build the
#'   interaction first, e.g.
#'   `data$site_age <- interaction(data$site, data$age_grp, sep = " / ")`.
#' @param times Time points, in the unit of `time`, at which the survival
#'   probability, its confidence interval and the number at risk are reported.
#'   Defaults to 12, 36 and 60. Include `0` to have the size of the starting
#'   risk set in the same table. `NULL` skips the time-point table altogether.
#' @param times_beyond What to do with a time point past the largest follow-up
#'   time of a stratum: `"error"` (the default) to stop, or `"na"` to keep the
#'   row, mark it `estimable = FALSE` and leave the survival probability and
#'   its interval missing. There is no setting that extrapolates.
#' @param group_missing What to do with rows whose `group` value is missing:
#'   `"error"` (the default) to stop, `"level"` to give them their own stratum
#'   labelled `"(Missing)"`, or `"drop"` to leave them out. The number of such
#'   rows is reported in the diagnostics in all three cases.
#' @param missing_text Label used for the stratum of rows with a missing
#'   `group` value when `group_missing = "level"`.
#' @param conf_type Confidence-interval transformation, passed to
#'   [survival::survfit()]: one of `"log"` (the default, and the default of
#'   `survfit()` itself), `"log-log"`, `"plain"`, `"logit"` or `"arcsin"`.
#'   `"none"` is not offered, because a confidence interval is part of what
#'   this function returns.
#' @param conf_level Confidence level, passed to [survival::survfit()] as
#'   `conf.int`. Defaults to `0.95`.
#' @param time_per_year How many units of `time` make up one year, used only to
#'   turn the total person-time into person-years. Defaults to `12`, i.e.
#'   `time` is in months; pass `365.25` for days or `1` for years.
#' @param digits_median Decimal places used when the median and its interval
#'   are pasted into `median_ci_txt`. The numeric columns are returned
#'   unrounded, so rounding stays a presentation choice and never a stored one.
#'
#' @return
#' A plain data frame with one row per stratum and the columns
#'
#' * `grouping` -- the name of the grouping variable, or `"overall"`;
#' * `stratum` -- the stratum label, `"<group>=<level>"` as `survfit()` writes
#'   it, or `"ALL"`; and `level`, the bare level, `NA` when `group` is `NULL`;
#' * `n`, `events` -- subjects and events in the stratum;
#' * `median`, `median_lcl`, `median_ucl` -- the median survival time and its
#'   confidence interval, `NA` when the curve never reaches 0.5;
#' * `median_reached`, `median_ci_txt` -- whether the median exists, and the
#'   formatted `"median (lcl-ucl)"` or `"not reached"`;
#' * `median_fu_revkm` -- the reverse-Kaplan-Meier median follow-up;
#' * `max_fu`, `total_fu`, `total_fu_yr` -- the largest follow-up time in the
#'   stratum, and its person-time in the unit of `time` and in years;
#' * `logrank_chisq`, `logrank_df`, `logrank_p` -- the log-rank test across all
#'   strata, repeated on every row, `NA` when there is a single stratum.
#'
#' Two attributes are attached:
#'
#' * `km_times` -- a data frame with one row per stratum and time point and the
#'   columns `grouping`, `stratum`, `level`, `time`, `n_risk`,
#'   `n_event_interval`, `n_censor_interval` (events and censorings since the
#'   previous requested time point, not cumulative), `surv`, `lcl`, `ucl` and
#'   `estimable`. `NULL` when `times` is `NULL`;
#' * `fit_km` -- a list with `counts` (rows read, rows used, strata, events,
#'   censored, and rows with a missing group), `follow_up` (minimum, median,
#'   maximum, person-time, person-years and `tau_upper`, the largest time point
#'   every stratum can still be evaluated at), `logrank` (the test, with a
#'   `note` saying why it was not computed when it was not), `checks` (the
#'   invariants that were enforced), `settings` and `call`.
#'
#' Attributes are dropped by most data-frame verbs, so read them off the object
#' returned by `fit_km()` before piping it further.
#'
#' @seealso [prep_surv()], which builds `time_mo` and `event_os` from a
#'   vital-status column.
#'
#' @examples
#' ## A small simulated cohort. No real patient records are used anywhere in
#' ## this package.
#' set.seed(20260901)
#' site <- factor(rep(c("Stomach", "Small intestine", "Colorectal"), each = 60),
#'                levels = c("Stomach", "Small intestine", "Colorectal"))
#' rate <- c(Stomach = 0.015, "Small intestine" = 0.025,
#'           Colorectal = 0.040)[as.character(site)]
#' t_event  <- stats::rexp(180, rate)
#' t_censor <- stats::runif(180, 6, 120)
#' toy <- data.frame(
#'   site     = site,
#'   time_mo  = round(pmin(t_event, t_censor), 1),
#'   event_os = as.integer(t_event <= t_censor)
#' )
#'
#' ## One curve for the whole cohort.
#' fit_km(toy)
#'
#' ## By primary site, with the log-rank test repeated on every row.
#' km <- fit_km(toy, group = "site")
#' km[, c("level", "n", "events", "median_ci_txt", "logrank_p")]
#'
#' ## Survival and number at risk at the requested time points.
#' attr(km, "km_times")
#'
#' ## Diagnostics travel with the result instead of being printed.
#' attr(km, "fit_km")$counts
#' attr(km, "fit_km")$follow_up
#' attr(km, "fit_km")$logrank
#' attr(km, "fit_km")$checks
#'
#' ## With no grouping variable there is no log-rank test, and the reason is
#' ## part of the result.
#' attr(fit_km(toy), "fit_km")$logrank$note
#'
#' ## A curve that never falls to 0.5 has no median, and says so.
#' censored <- toy
#' censored$event_os[censored$site == "Stomach"] <- 0L
#' fit_km(censored, group = "site")[, c("level", "events", "median_ci_txt")]
#'
#' ## A time point past the end of follow-up is refused, not extrapolated.
#' try(fit_km(toy, group = "site", times = c(12, 60, 96)))
#'
#' ## Unless you ask for it to be marked as not estimable instead.
#' beyond <- fit_km(toy, group = "site", times = c(60, 96), times_beyond = "na")
#' attr(beyond, "km_times")[, c("stratum", "time", "n_risk", "surv", "estimable")]
#'
#' ## A grouping variable with one observed level has nothing to compare.
#' try(fit_km(toy[toy$site == "Stomach", ], group = "site"))
#'
#' ## A misspelled column name is named in the error message.
#' try(fit_km(toy, time = "time_months", group = "site"))
#'
#' @importFrom survival Surv
#' @export
fit_km <- function(data,
                   time          = "time_mo",
                   event         = "event_os",
                   group         = NULL,
                   times         = c(12, 36, 60),
                   times_beyond  = c("error", "na"),
                   group_missing = c("error", "level", "drop"),
                   missing_text  = "Missing",
                   conf_type     = c("log", "log-log", "plain", "logit",
                                     "arcsin"),
                   conf_level    = 0.95,
                   time_per_year = 12,
                   digits_median = 1L) {

  cl <- match.call()

  ## -- 0. data -------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("fit_km(): `data` must be a data frame, got an object of class ",
         .ps_values(class(data)), ".", call. = FALSE)
  }
  n_in <- nrow(data)
  if (n_in == 0L) {
    stop("fit_km(): `data` has 0 rows, there is nothing to fit.", call. = FALSE)
  }

  ## -- 1. the argument values themselves -----------------------------------
  .km_name_arg(time,  "time")
  .km_name_arg(event, "event")
  if (!is.null(group)) .km_name_arg(group, "group")
  .km_name_arg(missing_text, "missing_text")
  times_beyond  <- match.arg(times_beyond)
  group_missing <- match.arg(group_missing)
  conf_type     <- match.arg(conf_type)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop(sprintf(
      paste0("fit_km(): `conf_level` must be a single number strictly between ",
             "0 and 1, got %s."),
      .ps_values(conf_level)),
      call. = FALSE)
  }
  if (!is.numeric(time_per_year) || length(time_per_year) != 1L ||
      is.na(time_per_year) || !is.finite(time_per_year) || time_per_year <= 0) {
    stop(sprintf(
      paste0("fit_km(): `time_per_year` must be a single positive number (12 ",
             "for months, 365.25 for days, 1 for years), got %s."),
      .ps_values(time_per_year)),
      call. = FALSE)
  }
  .km_count_arg(digits_median, "digits_median", min = 0L)
  times <- .km_times_arg(times)

  ## -- 2. required columns present, and present once -----------------------
  needed <- c(time = time, event = event)
  if (!is.null(group)) needed <- c(needed, group = group)

  miss <- needed[!(needed %in% names(data))]
  if (length(miss) > 0L) {
    stop(sprintf(
      paste0("fit_km(): %d required column%s not found in `data`:\n  %s\n",
             "  `data` has %d columns. Check the spelling, or pass the column ",
             "names explicitly."),
      length(miss), if (length(miss) > 1L) "s" else "",
      paste0(sprintf('%s = "%s"', names(miss), unname(miss)), collapse = "\n  "),
      ncol(data)),
      call. = FALSE)
  }
  dup <- needed[needed %in% names(data)[duplicated(names(data))]]
  if (length(dup) > 0L) {
    stop(sprintf(
      paste0("fit_km(): %s appear%s more than once among the columns of ",
             "`data`; fit_km() would silently use the first one. Make the ",
             "column names unique first."),
      paste0(sprintf('%s = "%s"', names(dup), unname(dup)), collapse = ", "),
      if (length(dup) > 1L) "" else "s"),
      call. = FALSE)
  }

  ## -- 3. follow-up time ---------------------------------------------------
  tv <- data[[time]]
  if (is.factor(tv) || is.character(tv)) {
    stop(sprintf(
      paste0('fit_km(): the `time` column "%s" is %s. Factor levels and ',
             "character digits are not numbers, so fit_km() will not convert ",
             "them for you: convert the column explicitly, or build it with ",
             "prep_surv(), which returns a numeric `time_mo`."),
      time, if (is.factor(tv)) "a factor" else "character"),
      call. = FALSE)
  }
  if (!is.numeric(tv)) {
    stop(sprintf(
      'fit_km(): the `time` column "%s" must be numeric, got an object of class %s.',
      time, .ps_values(class(tv))),
      call. = FALSE)
  }
  tnum <- as.numeric(tv)
  n_t_na <- sum(is.na(tnum))
  if (n_t_na > 0L) {
    stop(sprintf(
      paste0('fit_km(): the `time` column "%s" has %d missing value%s out of ',
             "%d rows. fit_km() neither drops nor imputes them; decide what to ",
             "do with these rows before calling it."),
      time, n_t_na, if (n_t_na > 1L) "s" else "", n_in),
      call. = FALSE)
  }
  n_t_inf <- sum(!is.finite(tnum))
  if (n_t_inf > 0L) {
    stop(sprintf(
      paste0('fit_km(): the `time` column "%s" has %d infinite value%s. A ',
             "curve cannot be estimated at an infinite follow-up time."),
      time, n_t_inf, if (n_t_inf > 1L) "s" else ""),
      call. = FALSE)
  }
  if (any(tnum < 0)) {
    n_neg <- sum(tnum < 0)
    stop(sprintf(
      paste0('fit_km(): the `time` column "%s" has %d negative value%s ',
             "(minimum %s). A negative follow-up time is not analysable, so ",
             "fit_km() stops instead of computing on it."),
      time, n_neg, if (n_neg > 1L) "s" else "", format(min(tnum))),
      call. = FALSE)
  }

  ## -- 4. event indicator ---------------------------------------------------
  evr <- data[[event]]
  if (is.factor(evr) || is.character(evr)) {
    stop(sprintf(
      paste0('fit_km(): the `event` column "%s" is %s. fit_km() will not guess ',
             "which level means an event: convert it to 0 (censored) / 1 ",
             "(event) explicitly, or build it with prep_surv(), which returns ",
             "a 0/1 `event_os` from a vital-status column.\n",
             "  values present : %s"),
      event, if (is.factor(evr)) "a factor" else "character", .ps_counts(evr)),
      call. = FALSE)
  }
  if (is.logical(evr)) {
    ev <- as.integer(evr)
  } else if (is.numeric(evr)) {
    ev <- as.numeric(evr)
  } else {
    stop(sprintf(
      paste0('fit_km(): the `event` column "%s" must be numeric 0/1 or ',
             "logical, got an object of class %s."),
      event, .ps_values(class(evr))),
      call. = FALSE)
  }
  bad_ev <- !(ev %in% c(0, 1))
  if (any(bad_ev)) {
    stop(sprintf(
      paste0('fit_km(): the `event` column "%s" has %d of %d row%s whose value ',
             "is neither 0 (censored) nor 1 (event).\n",
             "  offending values : %s\n",
             "  A missing indicator (<NA>) counts as an offending value: ",
             "fit_km() never guesses whether a subject had the event. Use ",
             "prep_surv() to build a 0/1 indicator from a vital-status column."),
      event, sum(bad_ev), n_in, if (n_in > 1L) "s" else "",
      .ps_counts(evr[bad_ev])),
      call. = FALSE)
  }
  ev <- as.integer(ev)

  ## -- 5. the grouping variable --------------------------------------------
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
        paste0('fit_km(): the `group` column "%s" is of class %s, which cannot ',
               "be used as a grouping variable. Convert it to a factor first."),
        group, .ps_values(class(g_raw))),
        call. = FALSE)
    }

    n_g_na    <- sum(is.na(g))
    n_lv_all  <- nlevels(g)
    g         <- droplevels(g)
    n_lv_drop <- n_lv_all - nlevels(g)

    if (n_g_na > 0L) {
      if (group_missing == "error") {
        stop(sprintf(
          paste0('fit_km(): the `group` column "%s" has %d missing value%s out ',
                 "of %d rows.\n",
                 "  fit_km() will not decide on its own what those rows are: ",
                 "they are neither a stratum nor nothing, and dropping them ",
                 "would change the denominator of every number it returns.\n",
                 '  Pass group_missing = "level" to give them their own ',
                 'stratum labelled "(%s)", or group_missing = "drop" to leave ',
                 "them out (either way the count comes back in the ",
                 "diagnostics), or filter them out before calling."),
          group, n_g_na, if (n_g_na > 1L) "s" else "", n_in, missing_text),
          call. = FALSE)
      } else if (group_missing == "level") {
        na_lab <- paste0("(", missing_text, ")")
        if (na_lab %in% levels(g)) {
          stop(sprintf(
            paste0('fit_km(): group_missing = "level" would add a stratum ',
                   'labelled "%s", but the `group` column "%s" already has a ',
                   "level with that name. Rename the level, or change ",
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
  tnum <- tnum[keep]
  ev   <- ev[keep]
  if (!is.null(g)) g <- droplevels(g[keep])
  n_used <- length(tnum)

  if (n_used == 0L) {
    stop(sprintf(
      paste0('fit_km(): every row has a missing `group` value ("%s"), so there ',
             "is nothing left to fit."),
      group),
      call. = FALSE)
  }

  lv <- if (is.null(g)) NULL else levels(g)
  if (!is.null(g) && length(lv) < 2L) {
    stop(sprintf(
      paste0('fit_km(): the `group` column "%s" has %d observed level%s (%s) ',
             "in the %d rows used; at least 2 are needed.\n",
             "  With one stratum there is nothing to compare: the log-rank ",
             "test is undefined and the result would be a single curve. Call ",
             "fit_km() with group = NULL if that single curve is what you ",
             "want.\n",
             "  Check that the data were not already filtered down to one ",
             "group, and that an empty group is not hiding behind an unused ",
             "factor level (%d unused level%s dropped here)."),
      group, length(lv), if (length(lv) == 1L) "" else "s",
      if (length(lv) == 0L) "none" else .ps_values(lv), n_used,
      n_lv_drop, if (n_lv_drop == 1L) " was" else "s were"),
      call. = FALSE)
  }

  ## -- 6. the fit -----------------------------------------------------------
  ## The fitted data frame is built here, with fixed internal column names, so
  ## that no column of `data` other than the three named ones ever reaches
  ## survfit(), and so that a column name that is not a syntactic R name cannot
  ## reach the model formula: survdiff() fails on a back-quoted term that
  ## survfit() accepts. The stratum labels are then written the way survfit()
  ## writes them, but with the caller's column name in place of the internal
  ## one.
  df <- data.frame(.km_time = tnum, .km_event = ev, stringsAsFactors = FALSE)
  if (!is.null(g)) df[[".km_group"]] <- g

  frm <- stats::as.formula(
    if (is.null(group)) "Surv(.km_time, .km_event) ~ 1" else
      "Surv(.km_time, .km_event) ~ .km_group")

  fit <- survival::survfit(frm, data = df,
                           conf.type = conf_type, conf.int = conf_level)

  lab_int <- if (is.null(group)) "ALL" else paste0(".km_group=", lv)
  lab     <- if (is.null(group)) "ALL" else .km_squish(paste0(group, "=", lv))
  tb      <- .km_fit_table(fit, lab_int[1])

  if (!identical(.km_squish(rownames(tb)), lab_int)) {
    stop(sprintf(
      paste0("fit_km(): the strata of the fit do not match the levels of ",
             '"%s", which should be impossible by construction.\n',
             "  survfit() : %s\n  expected  : %s\n",
             "  Do not use this result; report it as a bug in fit_km()."),
      group, .ps_values(rownames(tb)), .ps_values(lab_int)),
      call. = FALSE)
  }

  lcl_i <- grep("LCL$", colnames(tb))
  ucl_i <- grep("UCL$", colnames(tb))
  if (length(lcl_i) != 1L || length(ucl_i) != 1L) {
    stop(paste0("fit_km(): summary(survfit)$table did not carry exactly one ",
                "pair of median confidence limits. Do not use this result; ",
                "report it as a bug in fit_km()."),
         call. = FALSE)
  }

  ## -- 7. per-stratum quantities computed from the data ---------------------
  idx    <- if (is.null(g)) list(seq_len(n_used)) else split(seq_len(n_used), g)
  n_by   <- vapply(idx, length, integer(1))
  ev_by  <- vapply(idx, function(i) sum(ev[i]), integer(1))
  max_fu <- vapply(idx, function(i) max(tnum[i]), numeric(1))
  tot_fu <- vapply(idx, function(i) sum(tnum[i]), numeric(1))

  n_fit  <- as.integer(tb[, "records"])
  ev_fit <- as.integer(tb[, "events"])
  if (!identical(unname(n_by), n_fit) || !identical(unname(ev_by), ev_fit)) {
    stop(paste0("fit_km(): the subjects and events counted from `data` do not ",
                "match those counted by survfit(), which should be impossible ",
                "by construction. Do not use this result; report it as a bug ",
                "in fit_km()."),
         call. = FALSE)
  }

  ## -- 8. reverse-Kaplan-Meier median follow-up -----------------------------
  df_rev              <- df
  df_rev[[".km_event"]] <- 1L - df[[".km_event"]]
  fit_rev <- survival::survfit(frm, data = df_rev,
                               conf.type = conf_type, conf.int = conf_level)
  tb_rev  <- .km_fit_table(fit_rev, lab_int[1])
  if (!identical(.km_squish(rownames(tb_rev)), lab_int)) {
    stop(paste0("fit_km(): the strata of the reverse-Kaplan-Meier fit do not ",
                "match those of the survival fit, which should be impossible ",
                "by construction. Do not use this result; report it as a bug ",
                "in fit_km()."),
         call. = FALSE)
  }
  med_fu <- unname(tb_rev[, "median"])

  ## -- 9. log-rank ----------------------------------------------------------
  if (is.null(group)) {
    lr_chisq <- NA_real_
    lr_df    <- NA_integer_
    lr_p     <- NA_real_
    lr_note  <- paste0("not computed: no grouping variable was given, so there ",
                       "is a single stratum and nothing to compare")
    lr_test  <- NA_character_
  } else {
    sd       <- survival::survdiff(frm, data = df, rho = 0)
    lr_chisq <- unname(sd$chisq)
    lr_df    <- length(sd$n) - 1L
    lr_p     <- stats::pchisq(lr_chisq, df = lr_df, lower.tail = FALSE)
    lr_note  <- sprintf(
      "log-rank test (survdiff, rho = 0) across the %d strata of \"%s\"",
      length(lv), group)
    lr_test  <- "log-rank"
  }

  ## -- 10. survival and number at risk at the requested time points ---------
  km_times   <- NULL
  tau_upper  <- min(max_fu)
  all_within <- NA

  if (!is.null(times)) {
    grid <- expand.grid(time = times, stratum = lab,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    grid <- grid[, c("stratum", "time")]
    grid$estimable <- grid$time <= max_fu[match(grid$stratum, lab)]
    all_within     <- all(grid$estimable)

    if (!all_within && times_beyond == "error") {
      bad_s <- unique(grid$stratum[!grid$estimable])
      lines <- vapply(bad_s, function(s) {
        tt <- grid$time[grid$stratum == s & !grid$estimable]
        sprintf("%s : time point%s %s past a largest follow-up of %s",
                s, if (length(tt) > 1L) "s" else "",
                paste(format(tt, trim = TRUE), collapse = ", "),
                format(unname(max_fu[match(s, lab)])))
      }, character(1))
      stop(sprintf(
        paste0("fit_km(): %d of %d requested time point%s cannot be estimated ",
               "in %d of %d strat%s:\n  %s\n",
               "  A Kaplan-Meier curve is not estimable past the largest ",
               "follow-up time of a stratum, and summary.survfit(extend = ",
               "TRUE) would carry the last value forward without saying so. ",
               "This is the rule an RMST tau guard uses as well.\n",
               "  The largest time point every stratum can be evaluated at is ",
               "%s. Either keep `times` within it, or pass times_beyond = ",
               '"na" to keep those rows and mark them estimable = FALSE.'),
        length(unique(grid$time[!grid$estimable])), length(times),
        if (length(times) > 1L) "s" else "",
        length(bad_s), length(lab), if (length(lab) > 1L) "a" else "um",
        paste(lines, collapse = "\n  "), format(tau_upper)),
        call. = FALSE)
    }

    s <- summary(fit, times = times, extend = TRUE)
    s_lab <- if (is.null(group)) rep(lab_int, length(s$time)) else
      .km_squish(as.character(s$strata))
    k_have <- paste(s_lab, .km_key(s$time), sep = "\r")
    k_want <- paste(lab_int[match(grid$stratum, lab)], .km_key(grid$time),
                    sep = "\r")
    j <- match(k_want, k_have)
    if (anyNA(j) || anyDuplicated(k_have) > 0L) {
      stop(paste0("fit_km(): summary(survfit) did not return exactly one row ",
                  "per stratum and requested time point, which should be ",
                  "impossible with extend = TRUE. Do not use this result; ",
                  "report it as a bug in fit_km()."),
           call. = FALSE)
    }

    km_times <- data.frame(
      grouping          = if (is.null(group)) "overall" else group,
      stratum           = grid$stratum,
      level             = if (is.null(group)) NA_character_ else
        lv[match(grid$stratum, lab)],
      time              = grid$time,
      n_risk            = as.integer(s$n.risk[j]),
      n_event_interval  = as.integer(s$n.event[j]),
      n_censor_interval = as.integer(s$n.censor[j]),
      surv              = unname(s$surv[j]),
      lcl               = unname(s$lower[j]),
      ucl               = unname(s$upper[j]),
      estimable         = grid$estimable,
      stringsAsFactors  = FALSE,
      row.names         = NULL
    )
    km_times$surv[!km_times$estimable] <- NA_real_
    km_times$lcl[!km_times$estimable]  <- NA_real_
    km_times$ucl[!km_times$estimable]  <- NA_real_
  }

  ## -- 11. assemble ---------------------------------------------------------
  med <- unname(tb[, "median"])
  lcl <- unname(tb[, lcl_i])
  ucl <- unname(tb[, ucl_i])
  fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.*f", digits_median, x))

  out <- data.frame(
    grouping        = if (is.null(group)) "overall" else group,
    stratum         = lab,
    level           = if (is.null(group)) NA_character_ else lv,
    n               = n_fit,
    events          = ev_fit,
    median          = med,
    median_lcl      = lcl,
    median_ucl      = ucl,
    median_reached  = !is.na(med),
    median_ci_txt   = ifelse(is.na(med), "not reached",
                             sprintf("%s (%s-%s)", fmt(med), fmt(lcl), fmt(ucl))),
    median_fu_revkm = med_fu,
    max_fu          = unname(max_fu),
    total_fu        = unname(tot_fu),
    total_fu_yr     = unname(tot_fu) / time_per_year,
    logrank_chisq   = lr_chisq,
    logrank_df      = lr_df,
    logrank_p       = lr_p,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  ## -- 12. diagnostics ------------------------------------------------------
  inv <- c(
    "stratum sizes add up to the rows used"        = sum(out$n) == n_used,
    "events add up to the events in the rows used" = sum(out$events) == sum(ev),
    "person-time adds up to the follow-up in the rows used" =
      isTRUE(all.equal(sum(out$total_fu), sum(tnum))),
    "subjects and events agree with survfit()"     = TRUE,
    "reverse-Kaplan-Meier strata match the survival strata" = TRUE)
  if (!all(inv)) {
    stop(sprintf(
      paste0("fit_km(): the returned per-stratum totals do not add up, which ",
             "should be impossible by construction.\n",
             "  failed check%s : %s\n",
             "  Do not use this result; report it as a bug in fit_km()."),
      if (sum(!inv) > 1L) "s" else "",
      paste(names(inv)[!inv], collapse = "; ")),
      call. = FALSE)
  }
  checks <- c(inv, "requested time points within each stratum's follow-up" =
                all_within)

  attr(out, "km_times") <- km_times
  attr(out, "fit_km")   <- list(
    counts    = c(rows_in            = n_in,
                  rows_used          = n_used,
                  strata             = nrow(out),
                  events             = sum(ev),
                  censored           = sum(ev == 0L),
                  group_missing      = n_g_na,
                  group_missing_dropped = n_g_drop),
    follow_up = c(min       = min(tnum),
                  median    = stats::median(tnum),
                  max       = max(tnum),
                  total     = sum(tnum),
                  total_yr  = sum(tnum) / time_per_year,
                  tau_upper = unname(tau_upper)),
    logrank   = list(test  = lr_test,
                     chisq = lr_chisq,
                     df    = lr_df,
                     p     = lr_p,
                     note  = lr_note),
    checks    = checks,
    settings  = list(
      time                = time,
      event               = event,
      group               = group,
      group_levels        = lv,
      group_levels_dropped = n_lv_drop,
      group_missing       = if (is.null(group)) NULL else group_missing,
      times               = times,
      times_beyond        = times_beyond,
      conf_type           = conf_type,
      conf_level          = conf_level,
      time_per_year       = time_per_year,
      digits_median       = digits_median,
      estimator           = "Kaplan-Meier (survival::survfit, type = \"kaplan-meier\")",
      median_follow_up    = "reverse Kaplan-Meier on 1 - event",
      logrank             = "survival::survdiff, rho = 0"),
    call      = cl)
  out
}


# -- internal helpers, not exported ------------------------------------------
# .ps_values() and .ps_counts() are defined in R/prep_surv.R and reused here.

#' @noRd
.km_name_arg <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf(
      "fit_km(): `%s` must be a single non-missing string.", arg),
      call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.km_count_arg <- function(x, arg, min = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != as.integer(x) ||
      x < min) {
    stop(sprintf(
      "fit_km(): `%s` must be a single whole number >= %d, got %s.",
      arg, min, .ps_values(x)),
      call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.km_times_arg <- function(times) {
  if (is.null(times)) return(NULL)
  if (!is.numeric(times) || length(times) == 0L) {
    stop(sprintf(
      paste0("fit_km(): `times` must be a numeric vector of at least one time ",
             "point, or NULL to skip the time-point table, got %s of length ",
             "%d."),
      .ps_values(class(times)), length(times)),
      call. = FALSE)
  }
  bad <- is.na(times) | !is.finite(times)
  if (any(bad)) {
    stop(sprintf(
      "fit_km(): `times` has %d missing or infinite element%s (position%s %s).",
      sum(bad), if (sum(bad) > 1L) "s" else "",
      if (sum(bad) > 1L) "s" else "", paste(which(bad), collapse = ", ")),
      call. = FALSE)
  }
  if (any(times < 0)) {
    stop(sprintf(
      paste0("fit_km(): `times` has %d negative element%s (%s). A survival ",
             "probability before the origin does not exist."),
      sum(times < 0), if (sum(times < 0) > 1L) "s" else "",
      .ps_values(times[times < 0])),
      call. = FALSE)
  }
  d <- unique(times[duplicated(times)])
  if (length(d) > 0L) {
    stop(sprintf(
      paste0("fit_km(): `times` lists %s more than once, which would put the ",
             "same row in the time-point table twice. Keep one copy."),
      .ps_values(d)),
      call. = FALSE)
  }
  sort(times)
}

#' @noRd
.km_squish <- function(x) gsub("[[:space:]]+", " ", trimws(x))

#' @noRd
.km_key <- function(x) sprintf("%.12g", x)

#' @noRd
.km_fit_table <- function(fit, label) {
  tb <- summary(fit)$table
  if (is.null(dim(tb))) {
    tb <- matrix(tb, nrow = 1L, dimnames = list(label, names(tb)))
  }
  tb
}
