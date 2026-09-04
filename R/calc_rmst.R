# ---------------------------------------------------------------------------
# calc_rmst(): restricted mean survival time for a two-arm contrast, on one or
# more truncation times, returned as a data frame.
#
# Design notes (deliberate, and different from the study script this was
# distilled from):
#   * no column name is hard-coded: the time, the event indicator, the
#     grouping variable and the two arms are all arguments;
#   * nothing is printed and nothing is written to disk: one row per tau comes
#     back as a plain data frame, and the bookkeeping as attributes;
#   * the tau guard is structural. tau must not exceed the smallest of the two
#     arms' largest follow-up times, and there is no argument, no default and
#     no escape hatch that relaxes it: a caller cannot switch it off, and
#     calc_rmst() does not quietly shorten tau to the largest feasible value
#     either. Silently truncating tau would make the table report a tau other
#     than the one the text says was used, which is worse than an error;
#   * the guard limit is computed by exactly the calculation fit_km() uses for
#     follow_up[["tau_upper"]] -- split the rows by the grouping factor, take
#     the largest follow-up time within each level, take the smallest of those
#     -- so that the two functions can never disagree about where follow-up
#     ends. See .rm_tau_upper() at the bottom of this file;
#   * the arm indicator handed to survRM2 is checked against the grouping
#     column in both directions before anything is estimated. An arm vector
#     that had collapsed to a single value would put both groups in one arm,
#     and the tau guard would then "pass" on data that no longer held two
#     groups;
#   * no fitted object is attached. survRM2::rmst2() keeps a survfit object
#     per arm inside its result, and a survfit object can be turned back into
#     individual follow-up times; only summary statistics are returned.
#
# Why survRM2 rather than an in-house area calculation: the analysis plan
# names survRM2::rmst2() as the estimator, and calling the same function is
# the only way this package and the study scripts are guaranteed to agree to
# the last bit rather than to some tolerance. rmst2() prints nothing when its
# value is assigned (the note it composes is stored in the result and shown
# only by print.rmst2), it uses no random numbers anywhere, and it depends on
# survival alone, which this package already imports. Its own truncation-time
# check is looser than the one required here -- when the largest observation
# in both arms is an event it allows tau up to the larger of the two arms'
# maxima -- so it is never relied on: the guard below runs first and is
# strictly the tighter of the two.
#
# alpha is passed to rmst2() as signif(1 - conf_level, 12) rather than as
# 1 - conf_level: 1 - 0.95 is not the same double as the literal 0.05 that
# rmst2() defaults to and the study scripts use, and rounding to 12 significant
# digits makes the two the same bit pattern.
#
# The file is kept ASCII-only so that it behaves the same under a UTF-8 and
# under a non-UTF-8 locale.
#
# The small formatting helpers .ps_values() and .ps_counts() are shared with
# R/prep_surv.R; they are internal to the package.
# ---------------------------------------------------------------------------

#' Restricted mean survival time for a two-arm contrast
#'
#' @description
#' `calc_rmst()` estimates the restricted mean survival time (RMST) of two
#' arms at one or more truncation times and returns the numbers a survival
#' paper reports: the subjects and events in each arm, each arm's RMST with
#' its confidence interval, the between-arm difference with its confidence
#' interval and p-value, the ratio of the two RMSTs, the ratio of restricted
#' mean time lost, and the largest truncation time the contrast supports.
#'
#' The time, event and grouping columns and the two arms are arguments, so the
#' function is not tied to any particular registry export. The function prints
#' nothing, plots nothing and writes nothing: one row per truncation time comes
#' back as a plain data frame, with the bookkeeping attached as attributes.
#'
#' @details
#' # What is estimated
#'
#' The RMST at truncation time `tau` is the area under the Kaplan-Meier curve
#' between 0 and `tau`: the mean survival time of an arm, computed as if
#' follow-up stopped at `tau`. Unlike a hazard ratio it is a number of months
#' rather than a ratio, it needs no proportional-hazards assumption, and it
#' exists whether or not a median is reached.
#'
#' Three contrasts come back for each `tau`, all of them from
#' [survRM2::rmst2()] with no covariates, i.e. an unadjusted comparison of the
#' two arms:
#'
#' * `rmst_diff`, the difference `RMST(arm 1) - RMST(arm 0)`, in the unit of
#'   `time`. A negative difference means the arm named second lives, on
#'   average and within `tau`, that many months less;
#' * `rmst_ratio`, the ratio `RMST(arm 1) / RMST(arm 0)`, with an interval and
#'   a p-value computed on the log scale;
#' * `rmtl_ratio`, the ratio of restricted mean time lost,
#'   `(tau - RMST(arm 1)) / (tau - RMST(arm 0))`. Time lost is the complement
#'   of the RMST within the window, so this is the contrast that behaves like
#'   a risk ratio: it moves away from 1 much faster than `rmst_ratio` does,
#'   which is why a table that quotes only `rmst_ratio` can make two visibly
#'   different curves look almost identical.
#'
#' All three describe the same two curves over the same window. They cannot
#' disagree about direction, and they are not three independent findings.
#'
#' # The tau guard
#'
#' A Kaplan-Meier curve is not defined past the largest follow-up time in its
#' own arm, so the area under it is not defined past that point either. The
#' largest truncation time a two-arm contrast supports is therefore the
#' smaller of the two arms' largest follow-up times, and `calc_rmst()` stops
#' when any element of `tau` is past it.
#'
#' The guard is inclusive, `tau <= tau_upper`: a `tau` exactly equal to the
#' limit is estimated, because at that point both curves are still defined.
#'
#' There is no argument that relaxes the guard, no `tau_beyond = "na"` and no
#' automatic truncation to the largest feasible value. Shortening `tau` for
#' you would leave the result reporting a truncation time other than the one
#' that was asked for, and a table whose `tau` column disagrees with the
#' sentence citing it is much harder to catch than a stopped call. The limit
#' itself is returned, as the `tau_upper` column and as
#' `attr(x, "calc_rmst")$follow_up[["tau_upper"]]`, so a caller who needs a
#' feasible `tau` can ask for the limit and then choose one.
#'
#' The limit is computed exactly as [fit_km()] computes its
#' `follow_up[["tau_upper"]]`: the rows are split by the grouping factor, the
#' largest follow-up time is taken within each level, and the smallest of
#' those is the limit. On the same two arms the two functions return the same
#' number.
#'
#' The limit belongs to a contrast, not to a data set. A `tau` that is
#' feasible for one pair of sites can be refused for another, and one that is
#' feasible in the whole cohort can be refused inside an age stratum, because
#' each subset has its own shortest arm. Loop over the contrasts and let each
#' call check its own.
#'
#' # Choosing the arms
#'
#' RMST as estimated here is a two-arm contrast. With a grouping variable of
#' exactly two observed levels they are used in the order of the levels. With
#' more than two levels, name the two to compare in `arms`, reference first:
#' `arms = c("Stomach", "Colorectal")` puts Stomach in arm 0 and Colorectal in
#' arm 1, and every contrast is then Colorectal against Stomach. Rows
#' belonging to any other level are left out of that contrast and counted in
#' the diagnostics.
#'
#' For all pairs of a multi-level variable, loop over `combn(levels(x), 2)`
#' and stack the results; each pair then gets its own tau guard, which is the
#' point. For a contrast within a stratum, subset first and then call.
#'
#' # Reproducibility
#'
#' Nothing here is random. [survRM2::rmst2()] contains no call to the random
#' number generator, and neither does `calc_rmst()`: the Kaplan-Meier areas,
#' the Greenwood-type variances and the normal-approximation intervals are all
#' closed form. Repeated calls on the same data return bit-identical numbers,
#' with or without `set.seed()`, and a call does not change `.Random.seed`.
#'
#' # Confidence intervals and small arms
#'
#' The intervals and p-values are large-sample normal approximations: on the
#' RMST and on its difference directly, and on the log scale for the two
#' ratios. They are unreliable in a small arm, so an arm with fewer than
#' `min_n` subjects or fewer than `min_events` events raises a warning that
#' names the arm and the counts, and an arm with no event on or before some
#' `tau` raises a warning too -- its RMST at that `tau` is exactly `tau`, with
#' a standard error of 0 and a zero-width interval, and its restricted mean
#' time lost is 0, so `rmtl_ratio` at that `tau` is 0, infinite or `NaN` with a
#' `NaN` interval. The estimates are returned unchanged in both cases. The
#' small-arm warning is silenced by `min_n = 0, min_events = 0`, the
#' no-event-before-tau warning by `min_events = 0`.
#'
#' An arm with no event anywhere in follow-up is an error rather than a
#' warning, and that one cannot be silenced: no `tau` would be estimable for
#' it on the loss scale, so the whole table would be infinities.
#'
#' # Input checks
#'
#' The call stops, with the offending columns, values or counts listed, when
#'
#' 1. a required column (`time`, `event`, `group`) is not in `data`, appears
#'    in `data` more than once, or two of the three name the same column;
#' 2. `time` is not numeric, or holds missing, negative or infinite values.
#'    Rows are never dropped or imputed on your behalf;
#' 3. `event` holds a value other than `0` and `1`, missing values included.
#'    Use [prep_surv()] to build a `0`/`1` indicator from a vital-status
#'    column;
#' 4. `group` has fewer than two observed levels, or has more than two and
#'    `arms` was not given, or `arms` does not name exactly two levels that
#'    occur in the data;
#' 5. `group` has missing values and `group_missing = "error"` (the default);
#' 6. `tau` is not a positive, finite, non-missing, non-duplicated numeric
#'    vector;
#' 7. any element of `tau` is past the smallest of the two arms' largest
#'    follow-up times;
#' 8. an arm has no rows, or has no event anywhere in follow-up.
#'
#' @param data A data frame (or tibble) with one row per subject.
#' @param group Name of the grouping variable, as a single string. Required:
#'   an RMST contrast needs two arms, and this is the column they come from.
#' @param arms The two levels of `group` to compare, as a vector of length 2,
#'   reference first: `arms[1]` becomes arm 0 and `arms[2]` arm 1, and every
#'   contrast is arm 1 against arm 0. `NULL` (the default) is allowed only
#'   when `group` has exactly two observed levels, which are then used in the
#'   order of the levels.
#' @param time Name of the follow-up-time column, as a single string. Must be
#'   numeric. Factors and character digits are refused: use [prep_surv()],
#'   which converts them and returns a numeric `time_mo`.
#' @param event Name of the event-indicator column, as a single string. Must
#'   hold only `0` (censored) and `1` (event); logical `FALSE`/`TRUE` is
#'   accepted and read as `0`/`1`.
#' @param tau Truncation times, in the unit of `time`, one row of the result
#'   each. Defaults to 12, 36 and 60, i.e. months when `time` is in months.
#'   Every element must be positive and at most `tau_upper`; see the tau guard
#'   above. Returned sorted.
#' @param conf_level Confidence level of every interval, passed to
#'   [survRM2::rmst2()] as `alpha = signif(1 - conf_level, 12)`. Defaults to
#'   `0.95`.
#' @param group_missing What to do with rows whose `group` value is missing:
#'   `"error"` (the default) to stop, or `"drop"` to leave them out. The count
#'   comes back in the diagnostics either way. Rows belonging to a level other
#'   than the two in `arms` are always left out -- that is what naming the arms
#'   means -- and are counted separately.
#' @param min_n Smallest number of subjects in an arm that does not raise a
#'   warning. Defaults to `10`; `0` silences that warning.
#' @param min_events Smallest number of events in an arm that does not raise a
#'   warning. Defaults to `5`; `0` silences that warning and the separate
#'   warning about an arm with no event on or before some `tau`. An arm with no
#'   event anywhere in follow-up is an error whatever this is set to.
#' @param digits_rmst Decimal places used when an RMST or a difference is
#'   pasted into its `_txt` column. Defaults to `1`.
#' @param digits_ratio Decimal places used when a ratio is pasted into its
#'   `_txt` column. Defaults to `3`.
#' @param digits_p Decimal places used for the `_p_fmt` columns, which read
#'   `"<0.001"` below the corresponding threshold. Defaults to `3`. The numeric
#'   columns are returned unrounded, so all three `digits_` arguments stay a
#'   presentation choice and never a stored one.
#'
#' @return
#' A plain data frame with one row per element of `tau`, in increasing order of
#' `tau`, and the columns
#'
#' * `grouping`, `arm0`, `arm1`, `comparison` -- the grouping column, the two
#'   levels compared, and the label `"<arm1> vs <arm0>"`;
#' * `tau` -- the truncation time; and `tau_upper`, the largest truncation
#'   time this contrast supports, repeated on every row;
#' * `n_arm0`, `n_arm1`, `events_arm0`, `events_arm1` -- subjects and events in
#'   each arm over the whole of follow-up; `events_by_tau_arm0`,
#'   `events_by_tau_arm1`, the events on or before `tau`, which are the ones
#'   the variance at that `tau` is built from; `max_fu_arm0`, `max_fu_arm1`,
#'   the largest follow-up time in each arm;
#' * `rmst_arm0`, `rmst_arm0_se`, `rmst_arm0_lcl`, `rmst_arm0_ucl`,
#'   `rmst_arm0_txt`, and the same five for `arm1` -- each arm's RMST, its
#'   standard error, its confidence interval and the formatted
#'   `"est (lcl-ucl)"`;
#' * `rmst_diff`, `rmst_diff_lcl`, `rmst_diff_ucl`, `rmst_diff_p`,
#'   `rmst_diff_p_fmt`, `rmst_diff_txt` -- the difference, arm 1 minus arm 0;
#' * `rmst_ratio`, `rmst_ratio_lcl`, `rmst_ratio_ucl`, `rmst_ratio_p`,
#'   `rmst_ratio_p_fmt`, `rmst_ratio_txt` -- the ratio, arm 1 over arm 0;
#' * `rmtl_ratio`, `rmtl_ratio_lcl`, `rmtl_ratio_ucl`, `rmtl_ratio_p`,
#'   `rmtl_ratio_p_fmt`, `rmtl_ratio_txt` -- the ratio of restricted mean time
#'   lost;
#' * `diff_significant`, whether the confidence interval of the difference
#'   excludes 0, and `arm1_worse`, whether the difference is negative. Both are
#'   columns so that a direction is read off the table rather than re-derived
#'   by eye each time.
#'
#' One attribute is attached:
#'
#' * `calc_rmst` -- a list with `counts` (rows read, rows used, rows left out
#'   for a missing group or for belonging to another level, and the subjects
#'   and events per arm), `follow_up` (the largest follow-up time in each arm
#'   and `tau_upper`), `guard` (the rule enforced, the limit, whether the
#'   boundary is inclusive, and `relaxable`, which is `FALSE`), `checks` (the
#'   invariants enforced on the returned numbers), `settings` and `call`.
#'
#' Attributes are dropped by most data-frame verbs, so read them off the object
#' returned by `calc_rmst()` before piping it further; `tau_upper` is a column
#' as well as a diagnostic for that reason.
#'
#' @seealso [fit_km()], whose `follow_up[["tau_upper"]]` is the same limit
#'   computed the same way; [prep_surv()], which builds `time_mo` and
#'   `event_os` from a vital-status column; [survRM2::rmst2()], the estimator.
#'
#' @examples
#' ## A small simulated cohort. No real patient records are used anywhere in
#' ## this package.
#' set.seed(20260901)
#' site <- factor(rep(c("Stomach", "Small intestine", "Colorectal"), each = 60),
#'                levels = c("Stomach", "Small intestine", "Colorectal"))
#' rate <- c(Stomach = 0.010, "Small intestine" = 0.018,
#'           Colorectal = 0.030)[as.character(site)]
#' t_event  <- stats::rexp(180, rate)
#' t_censor <- stats::runif(180, 24, 140)
#' toy <- data.frame(
#'   site     = site,
#'   time_mo  = round(pmin(t_event, t_censor), 1),
#'   event_os = as.integer(t_event <= t_censor)
#' )
#'
#' ## Two arms picked out of three sites, reference first, at three truncation
#' ## times: one row per tau.
#' r <- calc_rmst(toy, "site", arms = c("Stomach", "Colorectal"))
#' r[, c("comparison", "tau", "rmst_arm0_txt", "rmst_arm1_txt",
#'       "rmst_diff_txt", "rmst_diff_p_fmt")]
#'
#' ## The ratio of restricted mean time lost moves further from 1 than the
#' ## ratio of the RMSTs themselves does.
#' r[, c("tau", "rmst_ratio_txt", "rmtl_ratio_txt")]
#'
#' ## Every column the table carries.
#' names(r)
#'
#' ## A single truncation time, and the bookkeeping that travels with it.
#' r60 <- calc_rmst(toy, "site", arms = c("Stomach", "Colorectal"), tau = 60)
#' attr(r60, "calc_rmst")$counts
#' attr(r60, "calc_rmst")$follow_up
#' attr(r60, "calc_rmst")$guard
#' attr(r60, "calc_rmst")$checks
#'
#' ## The tau limit is the number fit_km() reports as tau_upper, computed the
#' ## same way on the same two arms.
#' pair <- toy[toy$site %in% c("Stomach", "Colorectal"), ]
#' pair$site <- droplevels(pair$site)
#' c(calc_rmst = attr(r60, "calc_rmst")$follow_up[["tau_upper"]],
#'   fit_km    = attr(fit_km(pair, group = "site"),
#'                    "fit_km")$follow_up[["tau_upper"]])
#'
#' ## A cohort followed for at most 30 months: the same call is now refused,
#' ## because the area under a curve that has ended is not an estimate. The
#' ## guard cannot be switched off, and tau is never shortened for you.
#' recent <- toy
#' recent$event_os[recent$time_mo > 30] <- 0L
#' recent$time_mo <- pmin(recent$time_mo, 30)
#' try(calc_rmst(recent, "site", arms = c("Stomach", "Colorectal")))
#'
#' ## The limit itself is estimable: the guard is inclusive.
#' calc_rmst(recent, "site", arms = c("Stomach", "Colorectal"),
#'           tau = c(12, 30))[, c("tau", "tau_upper", "rmst_diff_txt")]
#'
#' ## Three levels and no `arms`: which two were meant is not guessed.
#' try(calc_rmst(toy, "site"))
#'
#' ## A level that does not occur is named in the error message.
#' try(calc_rmst(toy, "site", arms = c("Stomach", "Rectum")))
#'
#' ## A misspelled column name is named in the error message.
#' try(calc_rmst(toy, "site", arms = c("Stomach", "Colorectal"),
#'               time = "time_months"))
#'
#' ## A tau of 0 has no area under it.
#' try(calc_rmst(toy, "site", arms = c("Stomach", "Colorectal"),
#'               tau = c(0, 12)))
#'
#' ## An arm with no event has no restricted mean time lost, so the loss ratio
#' ## would be infinite; that is an error rather than a table of infinities.
#' immortal <- toy
#' immortal$event_os[immortal$site == "Stomach"] <- 0L
#' try(calc_rmst(immortal, "site", arms = c("Stomach", "Colorectal")))
#'
#' @export
calc_rmst <- function(data,
                      group,
                      arms          = NULL,
                      time          = "time_mo",
                      event         = "event_os",
                      tau           = c(12, 36, 60),
                      conf_level    = 0.95,
                      group_missing = c("error", "drop"),
                      min_n         = 10L,
                      min_events    = 5L,
                      digits_rmst   = 1L,
                      digits_ratio  = 3L,
                      digits_p      = 3L) {

  cl <- match.call()

  ## -- 0. data --------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("calc_rmst(): `data` must be a data frame, got an object of class ",
         .ps_values(class(data)), ".", call. = FALSE)
  }
  n_in <- nrow(data)
  if (n_in == 0L) {
    stop("calc_rmst(): `data` has 0 rows, there is nothing to estimate.",
         call. = FALSE)
  }

  ## -- 1. the argument values themselves ------------------------------------
  if (missing(group)) {
    stop(paste0("calc_rmst(): `group` is required. A restricted mean survival ",
                "time is estimated here as a two-arm contrast, and `group` is ",
                "the column the two arms come from; pass its name as a single ",
                'string, e.g. group = "site".'),
         call. = FALSE)
  }
  .rm_name_arg(time,  "time")
  .rm_name_arg(event, "event")
  .rm_name_arg(group, "group")
  group_missing <- match.arg(group_missing)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop(sprintf(
      paste0("calc_rmst(): `conf_level` must be a single number strictly ",
             "between 0 and 1, got %s."),
      .ps_values(conf_level)),
      call. = FALSE)
  }
  .rm_count_arg(min_n,        "min_n")
  .rm_count_arg(min_events,   "min_events")
  .rm_count_arg(digits_rmst,  "digits_rmst")
  .rm_count_arg(digits_ratio, "digits_ratio")
  .rm_count_arg(digits_p,     "digits_p", min = 1L)
  tau <- .rm_tau_arg(tau)

  ## -- 2. required columns present, present once, and distinct --------------
  needed <- c(time = time, event = event, group = group)

  miss <- needed[!(needed %in% names(data))]
  if (length(miss) > 0L) {
    stop(sprintf(
      paste0("calc_rmst(): %d required column%s not found in `data`:\n  %s\n",
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
      paste0("calc_rmst(): %s appear%s more than once among the columns of ",
             "`data`; calc_rmst() would silently use the first one. Make the ",
             "column names unique first."),
      paste0(sprintf('%s = "%s"', names(dup), unname(dup)), collapse = ", "),
      if (length(dup) > 1L) "" else "s"),
      call. = FALSE)
  }
  same <- duplicated(needed) | duplicated(needed, fromLast = TRUE)
  if (any(same)) {
    stop(sprintf(
      paste0("calc_rmst(): %s all name the same column of `data`. The ",
             "follow-up time, the event indicator and the grouping variable ",
             "must be three different columns."),
      paste0(sprintf('%s = "%s"', names(needed)[same], unname(needed)[same]),
             collapse = ", ")),
      call. = FALSE)
  }

  ## -- 3. follow-up time ----------------------------------------------------
  tv <- data[[time]]
  if (is.factor(tv) || is.character(tv)) {
    stop(sprintf(
      paste0('calc_rmst(): the `time` column "%s" is %s. Factor levels and ',
             "character digits are not numbers, so calc_rmst() will not ",
             "convert them for you: convert the column explicitly, or build ",
             "it with prep_surv(), which returns a numeric `time_mo`."),
      time, if (is.factor(tv)) "a factor" else "character"),
      call. = FALSE)
  }
  if (!is.numeric(tv)) {
    stop(sprintf(
      paste0('calc_rmst(): the `time` column "%s" must be numeric, got an ',
             "object of class %s."),
      time, .ps_values(class(tv))),
      call. = FALSE)
  }
  tnum   <- as.numeric(tv)
  n_t_na <- sum(is.na(tnum))
  if (n_t_na > 0L) {
    stop(sprintf(
      paste0('calc_rmst(): the `time` column "%s" has %d missing value%s out ',
             "of %d rows. calc_rmst() neither drops nor imputes them; decide ",
             "what to do with these rows before calling it."),
      time, n_t_na, if (n_t_na > 1L) "s" else "", n_in),
      call. = FALSE)
  }
  n_t_inf <- sum(!is.finite(tnum))
  if (n_t_inf > 0L) {
    stop(sprintf(
      paste0('calc_rmst(): the `time` column "%s" has %d infinite value%s. An ',
             "area under a curve up to tau cannot be built from an infinite ",
             "follow-up time."),
      time, n_t_inf, if (n_t_inf > 1L) "s" else ""),
      call. = FALSE)
  }
  if (any(tnum < 0)) {
    n_neg <- sum(tnum < 0)
    stop(sprintf(
      paste0('calc_rmst(): the `time` column "%s" has %d negative value%s ',
             "(minimum %s). A negative follow-up time is not analysable, so ",
             "calc_rmst() stops instead of computing on it."),
      time, n_neg, if (n_neg > 1L) "s" else "", format(min(tnum))),
      call. = FALSE)
  }

  ## -- 4. event indicator ---------------------------------------------------
  evr <- data[[event]]
  if (is.factor(evr) || is.character(evr)) {
    stop(sprintf(
      paste0('calc_rmst(): the `event` column "%s" is %s. calc_rmst() will ',
             "not guess which level means an event: convert it to 0 ",
             "(censored) / 1 (event) explicitly, or build it with ",
             "prep_surv(), which returns a 0/1 `event_os` from a vital-status ",
             "column.\n  values present : %s"),
      event, if (is.factor(evr)) "a factor" else "character", .ps_counts(evr)),
      call. = FALSE)
  }
  if (is.logical(evr)) {
    ev <- as.integer(evr)
  } else if (is.numeric(evr)) {
    ev <- as.numeric(evr)
  } else {
    stop(sprintf(
      paste0('calc_rmst(): the `event` column "%s" must be numeric 0/1 or ',
             "logical, got an object of class %s."),
      event, .ps_values(class(evr))),
      call. = FALSE)
  }
  bad_ev <- !(ev %in% c(0, 1))
  if (any(bad_ev)) {
    stop(sprintf(
      paste0('calc_rmst(): the `event` column "%s" has %d of %d row%s whose ',
             "value is neither 0 (censored) nor 1 (event).\n",
             "  offending values : %s\n",
             "  A missing indicator (<NA>) counts as an offending value: ",
             "calc_rmst() never guesses whether a subject had the event. Use ",
             "prep_surv() to build a 0/1 indicator from a vital-status ",
             "column."),
      event, sum(bad_ev), n_in, if (n_in > 1L) "s" else "",
      .ps_counts(evr[bad_ev])),
      call. = FALSE)
  }
  ev <- as.integer(ev)

  ## -- 5. the grouping variable and the two arms ----------------------------
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
      paste0('calc_rmst(): the `group` column "%s" is of class %s, which ',
             "cannot be used as a grouping variable. Convert it to a factor ",
             "first."),
      group, .ps_values(class(g_raw))),
      call. = FALSE)
  }

  n_g_na    <- sum(is.na(g))
  n_lv_all  <- nlevels(g)
  g         <- droplevels(g)
  n_lv_drop <- n_lv_all - nlevels(g)

  if (n_g_na > 0L && group_missing == "error") {
    stop(sprintf(
      paste0('calc_rmst(): the `group` column "%s" has %d missing value%s out ',
             "of %d rows.\n",
             "  calc_rmst() will not decide on its own what those rows are: ",
             "they belong to neither arm, and leaving them out silently would ",
             "hide part of the cohort from the counts it reports.\n",
             '  Pass group_missing = "drop" to leave them out (the count comes ',
             "back in the diagnostics), or filter them out before calling."),
      group, n_g_na, if (n_g_na > 1L) "s" else "", n_in),
      call. = FALSE)
  }

  lv <- levels(droplevels(g[!is.na(g)]))
  if (length(lv) < 2L) {
    stop(sprintf(
      paste0('calc_rmst(): the `group` column "%s" has %d observed level%s ',
             "(%s) in the %d rows with a non-missing group; a two-arm RMST ",
             "contrast needs 2.\n",
             "  Check that the data were not already filtered down to one ",
             "group, and that an empty group is not hiding behind an unused ",
             "factor level (%d unused level%s dropped here)."),
      group, length(lv), if (length(lv) == 1L) "" else "s",
      if (length(lv) == 0L) "none" else .ps_values(lv), n_in - n_g_na,
      n_lv_drop, if (n_lv_drop == 1L) " was" else "s were"),
      call. = FALSE)
  }

  if (is.null(arms)) {
    if (length(lv) != 2L) {
      stop(sprintf(
        paste0('calc_rmst(): the `group` column "%s" has %d observed levels ',
               "(%s) and `arms` was not given.\n",
               "  A restricted mean survival time is estimated here as a ",
               "two-arm contrast: which two of the %d levels were meant is ",
               "not guessed, because the answer would change the estimate ",
               "and the tau limit alike.\n",
               "  Name the two levels, reference first:\n",
               '    calc_rmst(data, group = "%s", arms = c("%s", "%s"))\n',
               "  For every pair, loop over combn(levels(data[[%s]]), 2) and ",
               "stack the results; each pair then gets its own tau guard, ",
               "which is the point of checking it per contrast."),
        group, length(lv), .ps_values(lv), length(lv), group, lv[1], lv[2],
        paste0('"', group, '"')),
        call. = FALSE)
    }
    arms <- lv
  } else {
    if (is.factor(arms)) arms <- as.character(arms)
    if (!is.atomic(arms) || length(arms) != 2L) {
      stop(sprintf(
        paste0("calc_rmst(): `arms` must name exactly 2 levels of \"%s\", ",
               "reference first, got %s of length %d.\n",
               "  Two is not a default that can be raised: the difference, ",
               "the ratio and the ratio of time lost are all defined between ",
               "two arms. Loop over combn(levels(x), 2) for all pairs."),
        group, .ps_values(class(arms)), length(arms)),
        call. = FALSE)
    }
    arms <- as.character(arms)
    if (anyNA(arms)) {
      stop(sprintf(
        paste0("calc_rmst(): `arms` has a missing element. Both arms must ",
               "name a level of \"%s\"."),
        group),
        call. = FALSE)
    }
    if (arms[1] == arms[2]) {
      stop(sprintf(
        paste0('calc_rmst(): `arms` names the same level twice ("%s"). An arm ',
               "compared with itself has a difference of 0 by construction ",
               "and is not a contrast."),
        arms[1]),
        call. = FALSE)
    }
    not_found <- arms[!(arms %in% lv)]
    if (length(not_found) > 0L) {
      stop(sprintf(
        paste0("calc_rmst(): %d of the 2 values of `arms` %s not occur in the ",
               "`group` column \"%s\":\n",
               "  not found       : %s\n  observed levels : %s\n",
               "  Level labels are matched exactly, including case and ",
               "spacing."),
        length(not_found), if (length(not_found) > 1L) "do" else "does",
        group, .ps_values(not_found), .ps_values(lv)),
        call. = FALSE)
    }
  }

  ## -- 6. subset to the two arms, and count what that left out --------------
  g_chr    <- as.character(g)
  in_arms  <- !is.na(g_chr) & g_chr %in% arms
  n_other  <- sum(!is.na(g_chr) & !(g_chr %in% arms))
  n_g_drop <- n_g_na

  tnum   <- tnum[in_arms]
  ev     <- ev[in_arms]
  g_used <- factor(g_chr[in_arms], levels = arms)
  n_used <- length(tnum)

  ## -- 7. the arm indicator, checked against the grouping column both ways --
  ## A contrast whose arm vector has collapsed to one value looks like a valid
  ## call all the way through: the tau guard would pass, rmst2() would run, and
  ## the table would report two arms that are in fact one. The check is cheap
  ## and it is not optional.
  arm <- ifelse(g_used == arms[2], 1L, ifelse(g_used == arms[1], 0L, NA_integer_))
  n0  <- sum(g_used == arms[1])
  n1  <- sum(g_used == arms[2])
  arm_ok <- !anyNA(arm) && all(arm %in% c(0L, 1L)) &&
    n0 > 0L && n1 > 0L && n0 + n1 == n_used &&
    all(arm[g_used == arms[1]] == 0L) && all(arm[g_used == arms[2]] == 1L) &&
    all(g_used[arm == 0L] == arms[1]) && all(g_used[arm == 1L] == arms[2])

  if (n0 == 0L || n1 == 0L) {
    stop(sprintf(
      paste0('calc_rmst(): arm %d ("%s") has 0 rows in `data`, so there is ',
             "nothing to compare. Arm 0 has %d rows and arm 1 has %d."),
      if (n0 == 0L) 0L else 1L, if (n0 == 0L) arms[1] else arms[2], n0, n1),
      call. = FALSE)
  }
  if (!arm_ok) {
    stop(paste0("calc_rmst(): the 0/1 arm indicator does not agree with the ",
                "grouping column in both directions, which should be ",
                "impossible by construction. Do not use this result; report ",
                "it as a bug in calc_rmst()."),
         call. = FALSE)
  }

  ev0 <- sum(ev[arm == 0L] == 1L)
  ev1 <- sum(ev[arm == 1L] == 1L)
  if (ev0 == 0L || ev1 == 0L) {
    zero <- which(c(ev0, ev1) == 0L)
    stop(sprintf(
      paste0("calc_rmst(): %s has no event anywhere in follow-up (%s).\n",
             "  Its restricted mean survival time is then exactly tau at ",
             "every tau, with a standard error of 0, and its restricted mean ",
             "time lost is exactly 0, so the ratio of time lost is 0 or ",
             "infinite and its interval and p-value are NaN. calc_rmst() will ",
             "not return a table of infinities, and there is no argument that ",
             "turns this off.\n",
             "  Check that the event indicator is the right column and the ",
             "right coding before concluding that nobody in this arm died."),
      paste(sprintf('arm %d ("%s")', zero - 1L, arms[zero]), collapse = " and "),
      paste(sprintf("%d event%s in %d subject%s",
                    c(ev0, ev1), ifelse(c(ev0, ev1) == 1L, "", "s"),
                    c(n0, n1), ifelse(c(n0, n1) == 1L, "", "s")),
            collapse = "; ")),
      call. = FALSE)
  }

  small <- (c(n0, n1) < min_n) | (c(ev0, ev1) < min_events)
  if (any(small)) {
    warning(sprintf(
      paste0("calc_rmst(): %s. The confidence intervals and p-values here are ",
             "large-sample normal approximations on the area under a ",
             "Kaplan-Meier curve, and below %d subject%s or %d event%s per arm ",
             "they are not reliable. The estimates are returned unchanged; ",
             "pass min_n = 0, min_events = 0 to silence this."),
      paste(sprintf('arm %d ("%s") has %d subject%s and %d event%s',
                    which(small) - 1L, arms[small], c(n0, n1)[small],
                    ifelse(c(n0, n1)[small] == 1L, "", "s"),
                    c(ev0, ev1)[small],
                    ifelse(c(ev0, ev1)[small] == 1L, "", "s")),
            collapse = "; "),
      min_n, if (min_n == 1L) "" else "s",
      min_events, if (min_events == 1L) "" else "s"),
      call. = FALSE)
  }

  ## -- 8. the tau guard -----------------------------------------------------
  ## The limit is min over arms of (largest follow-up in that arm), which is
  ## the calculation fit_km() uses for follow_up[["tau_upper"]]; see
  ## .rm_tau_upper() below. It is deliberately not an argument.
  max_fu    <- .rm_max_fu(tnum, g_used)
  tau_upper <- .rm_tau_upper(tnum, g_used)

  beyond <- tau > tau_upper
  if (any(beyond)) {
    stop(sprintf(
      paste0("calc_rmst(): %d of %d tau value%s past the follow-up of this ",
             "contrast and cannot be estimated.\n",
             "  contrast          : %s (arm 1) vs %s (arm 0), on \"%s\"\n",
             "  tau out of range  : %s\n",
             "  largest follow-up : %s in arm 0, %s in arm 1\n",
             "  largest usable tau: %s, the smaller of the two\n",
             "  The restricted mean survival time at tau is the area under a ",
             "Kaplan-Meier curve up to tau, and neither curve is defined past ",
             "the largest follow-up time of its own arm; past that point the ",
             "area would be an extrapolation rather than an estimate. The ",
             "rule enforced here is tau <= the smallest of the per-arm ",
             "largest follow-up times, and there is no argument that relaxes ",
             "it.\n",
             "  Lower `tau`, or compare arms with longer follow-up. ",
             "calc_rmst() does not shorten `tau` to the largest feasible ",
             "value on your behalf: the table would then report a truncation ",
             "time other than the one that was asked for, which is far harder ",
             "to notice than this error. The limit is available on its own as ",
             'attr(x, "calc_rmst")$follow_up[["tau_upper"]] of a call with a ',
             "feasible tau."),
      sum(beyond), length(tau), if (length(tau) > 1L) "s are" else " is",
      arms[2], arms[1], group,
      paste(.rm_num(tau[beyond]), collapse = ", "),
      .rm_num(unname(max_fu[1])), .rm_num(unname(max_fu[2])),
      .rm_num(tau_upper)),
      call. = FALSE)
  }

  ## -- 9. the estimates, one call to rmst2() per tau ------------------------
  ## alpha is rounded to 12 significant digits so that conf_level = 0.95 hands
  ## rmst2() exactly the double 0.05, the value its own default and the study
  ## scripts use; 1 - 0.95 is a different bit pattern.
  alpha <- signif(1 - conf_level, 12)
  k     <- length(tau)
  num   <- function() rep(NA_real_, k)
  est   <- list(r0 = num(), r0_se = num(), r0_lcl = num(), r0_ucl = num(),
                r1 = num(), r1_se = num(), r1_lcl = num(), r1_ucl = num(),
                d = num(), d_lcl = num(), d_ucl = num(), d_p = num(),
                rr = num(), rr_lcl = num(), rr_ucl = num(), rr_p = num(),
                lr = num(), lr_lcl = num(), lr_ucl = num(), lr_p = num(),
                tau_back = num())
  ev_by_tau0 <- integer(k)
  ev_by_tau1 <- integer(k)

  for (i in seq_len(k)) {
    ev_by_tau0[i] <- sum(ev[arm == 0L] == 1L & tnum[arm == 0L] <= tau[i])
    ev_by_tau1[i] <- sum(ev[arm == 1L] == 1L & tnum[arm == 1L] <= tau[i])

    fit <- survRM2::rmst2(time = tnum, status = ev, arm = arm,
                          tau = tau[i], alpha = alpha)
    .rm_check_rmst2(fit, tau[i])

    a0 <- fit$RMST.arm0$rmst
    a1 <- fit$RMST.arm1$rmst
    ur <- fit$unadjusted.result

    est$r0[i]     <- a0[[1L]]; est$r0_se[i]  <- a0[[2L]]
    est$r0_lcl[i] <- a0[[3L]]; est$r0_ucl[i] <- a0[[4L]]
    est$r1[i]     <- a1[[1L]]; est$r1_se[i]  <- a1[[2L]]
    est$r1_lcl[i] <- a1[[3L]]; est$r1_ucl[i] <- a1[[4L]]
    est$d[i]      <- ur[1L, 1L]; est$d_lcl[i]  <- ur[1L, 2L]
    est$d_ucl[i]  <- ur[1L, 3L]; est$d_p[i]    <- ur[1L, 4L]
    est$rr[i]     <- ur[2L, 1L]; est$rr_lcl[i] <- ur[2L, 2L]
    est$rr_ucl[i] <- ur[2L, 3L]; est$rr_p[i]   <- ur[2L, 4L]
    est$lr[i]     <- ur[3L, 1L]; est$lr_lcl[i] <- ur[3L, 2L]
    est$lr_ucl[i] <- ur[3L, 3L]; est$lr_p[i]   <- ur[3L, 4L]
    est$tau_back[i] <- fit$tau
  }

  no_ev <- (ev_by_tau0 == 0L) | (ev_by_tau1 == 0L)
  if (any(no_ev) && min_events > 0L) {
    warning(sprintf(
      paste0("calc_rmst(): an arm has no event on or before tau at %d of %d ",
             "truncation time%s (tau = %s). Its RMST there is exactly tau ",
             "with a standard error of 0 and a zero-width interval, its ",
             "restricted mean time lost is 0, and the rmtl_ratio columns of ",
             "those rows are 0, infinite or NaN with a NaN interval and ",
             "p-value. The rows are returned as computed; see ",
             "events_by_tau_arm0 and events_by_tau_arm1. Pass min_events = 0 ",
             "to silence this."),
      sum(no_ev), k, if (k > 1L) "s" else "",
      paste(.rm_num(tau[no_ev]), collapse = ", ")),
      call. = FALSE)
  }

  ## -- 10. assemble ---------------------------------------------------------
  f_r <- function(x) sprintf("%.*f", digits_rmst, x)
  f_q <- function(x) sprintf("%.*f", digits_ratio, x)

  out <- data.frame(
    grouping           = group,
    arm0               = arms[1],
    arm1               = arms[2],
    comparison         = paste0(arms[2], " vs ", arms[1]),
    tau                = tau,
    tau_upper          = tau_upper,
    n_arm0             = n0,
    n_arm1             = n1,
    events_arm0        = ev0,
    events_arm1        = ev1,
    events_by_tau_arm0 = ev_by_tau0,
    events_by_tau_arm1 = ev_by_tau1,
    max_fu_arm0        = unname(max_fu[1]),
    max_fu_arm1        = unname(max_fu[2]),
    rmst_arm0          = est$r0,
    rmst_arm0_se       = est$r0_se,
    rmst_arm0_lcl      = est$r0_lcl,
    rmst_arm0_ucl      = est$r0_ucl,
    rmst_arm0_txt      = paste0(f_r(est$r0), " (", f_r(est$r0_lcl), "-",
                                f_r(est$r0_ucl), ")"),
    rmst_arm1          = est$r1,
    rmst_arm1_se       = est$r1_se,
    rmst_arm1_lcl      = est$r1_lcl,
    rmst_arm1_ucl      = est$r1_ucl,
    rmst_arm1_txt      = paste0(f_r(est$r1), " (", f_r(est$r1_lcl), "-",
                                f_r(est$r1_ucl), ")"),
    rmst_diff          = est$d,
    rmst_diff_lcl      = est$d_lcl,
    rmst_diff_ucl      = est$d_ucl,
    rmst_diff_p        = est$d_p,
    rmst_diff_p_fmt    = .rm_fmt_p(est$d_p, digits_p),
    rmst_diff_txt      = paste0(f_r(est$d), " (", f_r(est$d_lcl), " to ",
                                f_r(est$d_ucl), ")"),
    rmst_ratio         = est$rr,
    rmst_ratio_lcl     = est$rr_lcl,
    rmst_ratio_ucl     = est$rr_ucl,
    rmst_ratio_p       = est$rr_p,
    rmst_ratio_p_fmt   = .rm_fmt_p(est$rr_p, digits_p),
    rmst_ratio_txt     = paste0(f_q(est$rr), " (", f_q(est$rr_lcl), "-",
                                f_q(est$rr_ucl), ")"),
    rmtl_ratio         = est$lr,
    rmtl_ratio_lcl     = est$lr_lcl,
    rmtl_ratio_ucl     = est$lr_ucl,
    rmtl_ratio_p       = est$lr_p,
    rmtl_ratio_p_fmt   = .rm_fmt_p(est$lr_p, digits_p),
    rmtl_ratio_txt     = paste0(f_q(est$lr), " (", f_q(est$lr_lcl), "-",
                                f_q(est$lr_ucl), ")"),
    diff_significant   = (est$d_lcl > 0) | (est$d_ucl < 0),
    arm1_worse         = est$d < 0,
    stringsAsFactors   = FALSE,
    row.names          = NULL
  )

  ## -- 11. diagnostics ------------------------------------------------------
  inv <- c(
    "every tau is within the follow-up of both arms" =
      isTRUE(all(out$tau <= tau_upper)),
    "the two arms partition the rows used" =
      n0 + n1 == n_used && n_used + n_g_drop + n_other == n_in,
    "the arm indicator agrees with the grouping column in both directions" =
      isTRUE(arm_ok),
    "survRM2 was given the tau it reported back" =
      .rm_agree(est$tau_back, tau, 1e-12),
    "no RMST is negative or larger than its own tau" =
      isTRUE(all(out$rmst_arm0 >= 0 & out$rmst_arm1 >= 0)) &&
      isTRUE(all(out$rmst_arm0 <= out$tau + 1e-8 &
                   out$rmst_arm1 <= out$tau + 1e-8)),
    "the difference is RMST(arm 1) minus RMST(arm 0)" =
      .rm_agree(out$rmst_diff, out$rmst_arm1 - out$rmst_arm0),
    "the ratio is RMST(arm 1) over RMST(arm 0)" =
      .rm_agree(out$rmst_ratio, out$rmst_arm1 / out$rmst_arm0),
    "the loss ratio is the time lost in arm 1 over the time lost in arm 0" =
      .rm_agree(out$rmtl_ratio,
                (out$tau - out$rmst_arm1) / (out$tau - out$rmst_arm0)),
    "one row per tau, no tau twice" =
      nrow(out) == k && !anyDuplicated(out$tau))
  if (!all(inv)) {
    stop(sprintf(
      paste0("calc_rmst(): the returned table does not add up, which should ",
             "be impossible by construction.\n  failed check%s : %s\n",
             "  Do not use this result; report it as a bug in calc_rmst()."),
      if (sum(!inv) > 1L) "s" else "",
      paste(names(inv)[!inv], collapse = "; ")),
      call. = FALSE)
  }

  attr(out, "calc_rmst") <- list(
    counts = c(rows_in               = n_in,
               rows_used             = n_used,
               rows_group_missing    = n_g_na,
               rows_group_dropped    = if (group_missing == "drop") n_g_drop
               else 0L,
               rows_other_levels     = n_other,
               n_arm0                = n0,
               n_arm1                = n1,
               events_arm0           = ev0,
               events_arm1           = ev1,
               censored_arm0         = n0 - ev0,
               censored_arm1         = n1 - ev1),
    follow_up = c(max_fu_arm0 = unname(max_fu[1]),
                  max_fu_arm1 = unname(max_fu[2]),
                  tau_upper   = tau_upper),
    guard = list(
      rule      = paste0("tau <= min over arms of (largest follow-up time in ",
                         "that arm)"),
      tau_upper = tau_upper,
      boundary  = "inclusive: a tau exactly equal to tau_upper is estimated",
      relaxable = FALSE,
      note      = sprintf(
        paste0("the largest follow-up time is %s in arm 0 (\"%s\") and %s in ",
               "arm 1 (\"%s\"), so this contrast supports tau up to %s. The ",
               "limit is a property of these two arms, not of the data set: ",
               "another pair, or the same pair inside a stratum, has its own. ",
               "This is the same quantity fit_km() returns as ",
               "follow_up[[\"tau_upper\"]], computed the same way."),
        .rm_num(unname(max_fu[1])), arms[1],
        .rm_num(unname(max_fu[2])), arms[2], .rm_num(tau_upper))),
    checks   = inv,
    settings = list(
      time          = time,
      event         = event,
      group         = group,
      arms          = arms,
      group_levels  = lv,
      group_missing = group_missing,
      tau           = tau,
      conf_level    = conf_level,
      alpha         = alpha,
      min_n         = min_n,
      min_events    = min_events,
      digits_rmst   = digits_rmst,
      digits_ratio  = digits_ratio,
      digits_p      = digits_p,
      estimator     = paste0("restricted mean survival time, area under the ",
                             "Kaplan-Meier curve (survRM2::rmst2, ",
                             "unadjusted)"),
      contrasts     = paste0("arm 1 (\"", arms[2], "\") against arm 0 (\"",
                             arms[1], "\"): difference, ratio, and ratio of ",
                             "restricted mean time lost"),
      ci            = paste0("large-sample normal; on the log scale for the ",
                             "two ratios"),
      random        = "none: rmst2() and calc_rmst() use no random numbers"),
    call = cl)
  out
}


# -- internal helpers, not exported ------------------------------------------
# .ps_values() and .ps_counts() are defined in R/prep_surv.R and reused here.

#' @noRd
.rm_name_arg <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("calc_rmst(): `%s` must be a single non-missing string.", arg),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.rm_count_arg <- function(x, arg, min = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != as.integer(x) || x < min) {
    stop(sprintf(
      "calc_rmst(): `%s` must be a single whole number >= %d, got %s.",
      arg, min, .ps_values(x)),
      call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.rm_tau_arg <- function(tau) {
  if (!is.numeric(tau) || length(tau) == 0L) {
    stop(sprintf(
      paste0("calc_rmst(): `tau` must be a numeric vector of at least one ",
             "truncation time, got %s of length %d. There is no default that ",
             "means \"as far as the data go\": a truncation time is part of ",
             "what an RMST reports and has to be chosen."),
      .ps_values(class(tau)), length(tau)),
      call. = FALSE)
  }
  bad <- is.na(tau) | !is.finite(tau)
  if (any(bad)) {
    stop(sprintf(
      paste0("calc_rmst(): `tau` has %d missing or infinite element%s ",
             "(position%s %s)."),
      sum(bad), if (sum(bad) > 1L) "s" else "",
      if (sum(bad) > 1L) "s" else "", paste(which(bad), collapse = ", ")),
      call. = FALSE)
  }
  if (any(tau <= 0)) {
    stop(sprintf(
      paste0("calc_rmst(): `tau` has %d element%s that %s not positive (%s). ",
             "The restricted mean survival time is the area under the curve ",
             "between 0 and tau, and there is no area to either side of ",
             "tau = 0."),
      sum(tau <= 0), if (sum(tau <= 0) > 1L) "s" else "",
      if (sum(tau <= 0) > 1L) "are" else "is", .ps_values(tau[tau <= 0])),
      call. = FALSE)
  }
  d <- unique(tau[duplicated(tau)])
  if (length(d) > 0L) {
    stop(sprintf(
      paste0("calc_rmst(): `tau` lists %s more than once, which would put the ",
             "same estimate in the table twice. Keep one copy."),
      .ps_values(d)),
      call. = FALSE)
  }
  sort(tau)
}

# Numbers inside messages, printed with enough significant digits that a tau
# a hair past the limit does not read as if it were exactly at it.
#' @noRd
.rm_num <- function(x) format(x, trim = TRUE, digits = 15)

#' @noRd
.rm_fmt_p <- function(p, digits) {
  thr <- 10^(-digits)
  ifelse(is.na(p), NA_character_,
         ifelse(p < thr,
                paste0("<", formatC(thr, format = "f", digits = digits)),
                formatC(p, format = "f", digits = digits)))
}

# Do two vectors of estimates agree? Used only on the algebraic identities
# checked at the end of calc_rmst(). An arm with no event on or before tau has
# a restricted mean time lost of exactly 0, so the loss ratio is legitimately
# 0, Inf or NaN there; two identical infinities and two NaNs agree, and
# abs(Inf - Inf) does not.
#' @noRd
.rm_agree <- function(a, b, tol = 1e-8) {
  both_na  <- is.na(a) & is.na(b)
  both_num <- !is.na(a) & !is.na(b)
  both_inf <- both_num & is.infinite(a) & is.infinite(b) & (sign(a) == sign(b))
  close    <- both_num & is.finite(a) & is.finite(b) & abs(a - b) < tol
  isTRUE(all(both_na | both_inf | close))
}

# The largest follow-up time within each level of the grouping factor, and the
# smallest of those. These two lines are the calculation fit_km() runs to fill
# follow_up[["tau_upper"]] -- split the row indices by the grouping factor,
# max() the times within each, min() across levels -- kept here in the same
# form so that the RMST tau guard and the Kaplan-Meier time-point guard can
# never drift apart. Do not "simplify" either one without the other.
#' @noRd
.rm_max_fu <- function(time, g) {
  idx <- split(seq_along(time), g)
  vapply(idx, function(i) max(time[i]), numeric(1))
}

#' @noRd
.rm_tau_upper <- function(time, g) {
  min(.rm_max_fu(time, g))
}

#' @noRd
.rm_check_rmst2 <- function(fit, tau) {
  bug <- function(what) {
    stop(sprintf(
      paste0("calc_rmst(): survRM2::rmst2() returned %s, so its result cannot ",
             "be read the way calc_rmst() reads it. Do not use this result; ",
             "check the installed version of survRM2 and report it as a bug ",
             "in calc_rmst()."),
      what),
      call. = FALSE)
  }
  if (!is.list(fit) ||
      !all(c("tau", "RMST.arm0", "RMST.arm1", "unadjusted.result") %in%
             names(fit))) {
    bug("an object without the elements tau, RMST.arm0, RMST.arm1 and unadjusted.result")
  }
  if (!isTRUE(all.equal(fit$tau, tau))) {
    bug(sprintf("tau = %s when it was given tau = %s",
                format(fit$tau), format(tau)))
  }
  ur <- fit$unadjusted.result
  if (!is.matrix(ur) || nrow(ur) != 3L || ncol(ur) != 4L ||
      !identical(rownames(ur), c("RMST (arm=1)-(arm=0)",
                                 "RMST (arm=1)/(arm=0)",
                                 "RMTL (arm=1)/(arm=0)")) ||
      !identical(colnames(ur)[c(1L, 4L)], c("Est.", "p"))) {
    bug("an unadjusted.result matrix with unexpected dimensions or labels")
  }
  for (a in c("RMST.arm0", "RMST.arm1")) {
    v <- fit[[a]]$rmst
    if (!is.numeric(v) || length(v) != 4L ||
        !identical(names(v)[1:2], c("Est.", "se"))) {
      bug(sprintf("a %s$rmst vector with unexpected length or labels", a))
    }
  }
  invisible(TRUE)
}
