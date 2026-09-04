# ---------------------------------------------------------------------------
# compare_estimands(): the hazard ratio, the Kaplan-Meier survival difference
# at tau and the RMST difference at tau, for one two-arm contrast, side by
# side, with a pre-specified judgement of whether the three agree.
#
# Design notes (deliberate, and different from the study script this was
# distilled from):
#   * no column name is hard-coded: the time, the event indicator, the
#     grouping variable and the two arms are all arguments;
#   * nothing is printed, nothing is plotted and nothing is written to disk:
#     one row per tau comes back as a plain data frame, and the bookkeeping as
#     attributes;
#   * the three estimands are not re-implemented here. fit_cox(), fit_km() and
#     calc_rmst() are called, on one analysis set built once and handed to all
#     three, so that this table and a table assembled from the three functions
#     separately cannot drift apart. The cost is that the same input checks run
#     three times; that cost buys the reconciliation below, which is the point
#     of the function;
#   * the tau guard is calc_rmst()'s and is not duplicated, not re-implemented
#     and not bypassed. calc_rmst() is called first, so a tau past the
#     follow-up of the contrast stops the whole call before anything else is
#     computed, with calc_rmst()'s own message;
#   * "the three estimands describe the same patients" is an enforced
#     invariant, not an assumption. Rows with a missing time, event or group
#     are refused at the door rather than dropped by each estimator on its own
#     terms, the analysis set is handed to the three with their strictest
#     settings (na_action = "fail", group_missing = "error",
#     empty_levels = "error"), and the subjects and events each of them reports
#     back are reconciled against each other before a single estimate is read.
#     A side-by-side table whose three columns rest on different subsets does
#     not compare estimands, it compares analysis sets, and it would do so
#     silently;
#   * the concordance rules are fixed in this file and written out in full in
#     the documentation. They do not depend on the data, there is no argument
#     that tunes them, and the only threshold in them, alpha, is 1 - conf_level
#     so that the interval and the p-value of every estimand can never give
#     different verdicts;
#   * no fitted object is attached. fit_cox(), fit_km() and calc_rmst() attach
#     none either; only summary statistics are returned.
#
# Why the Kaplan-Meier standard error is computed here rather than taken from
# fit_km(): the difference of two survival probabilities needs the standard
# error of each, and fit_km() returns the confidence interval but not the
# standard error it was built from. The interval is on the log scale by
# default, so the standard error cannot be recovered from it exactly. One
# survfit(~ 1) per arm is therefore fitted on the same analysis set and read
# with summary(..., extend = TRUE), exactly as the study script does, and the
# survival probabilities and numbers at risk it produces are checked against
# fit_km()'s to 1e-10 before they are used. The check is what makes the extra
# fit safe: if the two routes ever disagreed, the call stops.
#
# The arm-0 Kaplan-Meier curve is needed a second time, as a step function, to
# integrate the survival curve a constant hazard ratio would imply. That
# integral is self-checked against survRM2's RMST for arm 0 to 1e-6, as in the
# study script: if the same integration cannot reproduce a number survRM2
# computed independently, the implied quantity is wrong and the call stops.
#
# The file is kept ASCII-only so that it behaves the same under a UTF-8 and
# under a non-UTF-8 locale.
#
# The small formatting helpers .ps_values() and .ps_counts() are shared with
# R/prep_surv.R; they are internal to the package.
# ---------------------------------------------------------------------------

#' Hazard ratio, Kaplan-Meier difference and RMST difference side by side
#'
#' @description
#' `compare_estimands()` takes one two-arm contrast and reports, on the same
#' patients and over the same truncation times, the three effect measures a
#' survival paper usually quotes one at a time: the unadjusted hazard ratio
#' with its proportional-hazards test, the difference in Kaplan-Meier survival
#' probability at `tau`, and the difference in restricted mean survival time
#' up to `tau`. Each row also carries how far apart the three have ended up:
#' whether they point at the same arm, whether they agree about significance,
#' what a *constant* hazard ratio would have implied for the RMST difference,
#' and by how many months the observed difference misses that implication.
#'
#' The three estimands are not re-implemented here. [fit_cox()], [fit_km()] and
#' [calc_rmst()] are called on one analysis set, built once and handed to all
#' three, and the subjects and events each of them reports back are reconciled
#' before any estimate is read: a side-by-side table whose columns rest on
#' different subsets compares analysis sets, not estimands.
#'
#' The time, event and grouping columns and the two arms are arguments, so the
#' function is not tied to any particular registry export. The function prints
#' nothing, plots nothing and writes nothing: one row per truncation time comes
#' back as a plain data frame, with the bookkeeping attached as attributes.
#'
#' @details
#' # Why three numbers for one comparison
#'
#' The three estimands answer three different questions about the same two
#' curves, and each is silent about what the others measure:
#'
#' * the **hazard ratio** is one number for the whole of follow-up. It is a
#'   ratio, so it says nothing about how much time the difference is worth,
#'   and it only means what it is usually taken to mean while the hazards stay
#'   proportional. The [survival::cox.zph()] test of that assumption is
#'   returned in the `ph_p` column of every row;
#' * the **Kaplan-Meier difference** is the gap between the two curves at one
#'   instant, `tau`. It is in percentage points and needs no modelling
#'   assumption, but it uses only the two curve heights at `tau` and throws
#'   away everything that happened on the way there;
#' * the **RMST difference** is the area between the two curves from 0 to
#'   `tau`. It is in months, it needs no proportional-hazards assumption, and
#'   it exists whether or not a median is reached, but it depends on `tau` and
#'   averages away the shape of the difference inside the window.
#'
#' Because they measure different things they can, and on real data do,
#' disagree: a hazard ratio far from 1 can sit next to an RMST difference of a
#' few weeks, and two estimands can even point at different arms. The columns
#' below are there to make that disagreement a number rather than an
#' impression.
#'
#' # One analysis set for all three
#'
#' Everything in a row is computed on the same subjects. The analysis set is
#' built once: the rows whose `group` is one of the two `arms` are kept, rows
#' belonging to any other level are counted and left out, and rows with a
#' missing `time`, `event` or `group` are refused (`group_missing = "drop"`
#' leaves the missing-group rows out instead, uniformly, before any estimator
#' sees the data).
#'
#' That set is then handed to [calc_rmst()], [fit_cox()] and [fit_km()] with
#' their strictest settings -- `group_missing = "error"`, `na_action = "fail"`,
#' `empty_levels = "error"` -- so none of them can quietly drop a row of its
#' own accord, and the counts they report back are reconciled against each
#' other and against the analysis set:
#'
#' * rows used, by all three;
#' * subjects in each arm, by all three;
#' * events in each arm, by all three.
#'
#' Any disagreement is an error, not a warning and not a note. The reconciled
#' table travels back as `attr(x, "compare_estimands")$same_patients`, so the
#' invariant can be read off the result rather than taken on trust.
#'
#' A second invariant of the same kind is checked on the estimates themselves:
#' the Kaplan-Meier survival probabilities used for the difference come from
#' one `survfit(~ 1)` per arm (see below), and they must equal the ones
#' [fit_km()] reports at the same `tau` to within 1e-10, numbers at risk
#' included.
#'
#' # The tau guard
#'
#' The truncation times are checked by [calc_rmst()], with the rule it
#' documents: `tau` must not exceed the smaller of the two arms' largest
#' follow-up times, the boundary is inclusive, and there is no argument that
#' relaxes it. The guard is not duplicated, not re-implemented and not
#' bypassed: [calc_rmst()] is called first, before anything else is computed,
#' so a `tau` past the follow-up of the contrast stops the whole call, and the
#' error message that comes back says `calc_rmst()` because that is the
#' function that refused. The limit is returned as the `tau_upper` column.
#'
#' The same `tau` is used for the Kaplan-Meier time points, so the three
#' columns of a row are truncated at the same place. [fit_km()] applies the
#' same rule to its time points, which is why the two can never disagree about
#' where follow-up ends.
#'
#' # What a constant hazard ratio would imply
#'
#' If the hazards really were proportional with ratio `theta`, the comparator
#' arm's survival curve would be `S1(t) = S0(t)^theta`. Taking `S0` to be the
#' Kaplan-Meier curve of arm 0 and `theta` to be the estimated hazard ratio,
#' integrating that implied curve from 0 to `tau` and subtracting the area
#' under `S0` gives `rmst_diff_implied`: the RMST difference one would see at
#' this `tau` if a single hazard ratio told the whole story.
#'
#' `rmst_gap_mo = rmst_diff - rmst_diff_implied` is what is left over -- the
#' months of time-scale information a single hazard ratio cannot carry. It is
#' a difference of two areas, never a ratio, so it can be read whatever the
#' numbers are; `gap_lo_mo` and `gap_hi_mo` repeat it at the smallest and the
#' largest `tau` on every row.
#'
#' The integration is the same step-function sum that a Kaplan-Meier area is:
#' `S` is held at `S(t_i)` on the interval from `t_i` up to but not including
#' `t_i+1`, starting from `S(0) = 1`. It is self-checked at every `tau`: the
#' same code applied to `S0` itself must reproduce the arm-0 RMST
#' [survRM2::rmst2()] computed independently, to within 1e-6. If it does not,
#' the integration is wrong, the implied difference is meaningless, and the
#' call stops rather than returning it.
#'
#' Two further summaries of the same idea are contrast-level columns, computed
#' between the smallest and the largest `tau` and repeated on every row:
#' `ratio_obs` is `d_hi / d_lo`, how much the observed difference grew across
#' the window, and `ratio_implied` is the same ratio for the implied
#' differences; `ratio_dev_pct` is `100 * (ratio_obs / ratio_implied - 1)` and
#' `rel_hi_pct` is `100 * (d_hi / di_hi - 1)`. A ratio whose denominator sits
#' against zero is not informative, so a denominator smaller than `min_denom`
#' in absolute value returns `NA` with the reason in `het_note`, and
#' `gap_hi_mo` is read instead. With a single `tau` the smallest and the
#' largest coincide, `ratio_obs` and `ratio_implied` are 1 by construction, and
#' only the gaps are meaningful.
#'
#' # The concordance rules
#'
#' The rules are fixed here, in this file and in this section. They do not
#' depend on the data, no argument tunes them, and nothing in them is chosen
#' after seeing a result. Read them before reading a `verdict`.
#'
#' **Which arm is worse.** For each estimand, one label, in a column of its
#' own:
#'
#' * `dir_hr` is `"arm1_worse"` when `hr > 1`, else `"arm0_worse"`;
#' * `dir_km` is `"arm1_worse"` when `km_surv_diff < 0`, i.e. arm 1 has the
#'   lower survival probability at `tau`, else `"arm0_worse"`;
#' * `dir_rmst` is `"arm1_worse"` when `rmst_diff < 0`, i.e. arm 1 has the
#'   shorter restricted mean survival time, else `"arm0_worse"`.
#'
#' The comparisons are strict, so an exact tie -- a hazard ratio of exactly 1,
#' a difference of exactly 0 -- reads as `"arm0_worse"`. Continuous estimates
#' do not tie in practice. The column `dir_km_median` is the same label for the
#' difference in median survival; it is reported for reference only, is `NA`
#' when a median is not reached, and never enters the judgement.
#'
#' **Significant or not.** At `alpha = 1 - conf_level`, judged from the
#' confidence interval:
#'
#' * `sig_hr` : the interval of the hazard ratio excludes 1;
#' * `sig_km` : the interval of the Kaplan-Meier difference excludes 0;
#' * `sig_rmst` : the interval of the RMST difference excludes 0.
#'
#' The p-values are returned next to them, and an invariant checks that each
#' interval gives the same verdict as its own `p < alpha`; the interval and the
#' p-value of each estimand come from the same normal approximation, so they
#' cannot legitimately differ.
#'
#' **Do they agree.**
#'
#' * `direction_agree` : all three direction labels are the same;
#' * `signif_agree` : all three significance verdicts are the same;
#' * `n_signif` : how many of the three are significant, 0 to 3;
#' * `concordant` : `direction_agree & signif_agree`;
#' * `inconsistent` : `!concordant`.
#'
#' **What kind of disagreement.** `discordance_type` takes exactly one of four
#' values, tested in this order:
#'
#' 1. `"concordant"` -- `concordant` is `TRUE`;
#' 2. `"direction_conflict_significant"` -- the directions differ **and** some
#'    pair of estimands that point at different arms are **both** significant.
#'    That pair contradicts each other on the data's own terms, and this is the
#'    row a Results section has to address;
#' 3. `"direction_conflict_ci_covers_null"` -- the directions differ, but no
#'    such pair is both significant: every sign flip involves an estimand whose
#'    interval covers the null. A difference sitting against zero has a sign,
#'    and that sign is noise, so this is **not** counted as a contradiction;
#' 4. `"significance_only"` -- the directions all agree and only the
#'    significance verdicts differ.
#'
#' The pair test is the column `dir_conflict_signif`, which is `TRUE` when
#' `(sig_hr & sig_km & dir_hr != dir_km)` or
#' `(sig_hr & sig_rmst & dir_hr != dir_rmst)` or
#' `(sig_km & sig_rmst & dir_km != dir_rmst)`.
#'
#' The boundary case is settled by rule 3 and is worth stating on its own,
#' because it decides how many contradictions a table appears to contain:
#' **two estimands that point at different arms while at least one of the two
#' has an interval covering the null are not a contradiction.** Such a row is
#' still `inconsistent` -- the directions really do differ -- but it is typed
#' apart from a real conflict and its `verdict` begins with `minor:` rather
#' than `DISCORDANT`. Counting sign flips near zero as findings would inflate
#' the disagreement the function exists to measure.
#'
#' **The sentence.** `sig_pattern` is `"HR sig/ns / KM sig/ns / RMST sig/ns"`,
#' and `verdict` is one of
#'
#' * `"concordant: all three significant, same direction"`;
#' * `"concordant: all three non-significant, same direction"`;
#' * `"DISCORDANT (direction): two significant estimands point to different arms (<sig_pattern>)"`;
#' * `"minor: directions differ but the conflicting estimands are not both significant (<sig_pattern>)"`;
#' * `"DISCORDANT (significance): <n> of 3 significant, same direction (<sig_pattern>)"`.
#'
#' A concordant row always has `n_signif` equal to 0 or 3, since `signif_agree`
#' is part of `concordant`; that is checked rather than assumed.
#'
#' # Adjustment
#'
#' All three estimands are unadjusted. That is deliberate and it is what makes
#' the comparison a comparison of estimands: a multivariable hazard ratio put
#' next to an unadjusted RMST difference would differ from it partly because of
#' covariate adjustment, and there would be no way to tell which part is which.
#' If an adjusted hazard ratio is wanted alongside, fit it with [fit_cox()] and
#' add it as a reference column of your own; it does not belong in the
#' concordance judgement.
#'
#' # Input checks
#'
#' The call stops, with the offending columns, values or counts listed, when
#'
#' 1. a required column (`time`, `event`, `group`) is not in `data`, appears in
#'    `data` more than once, or two of the three name the same column;
#' 2. `time` is not numeric, or holds missing, negative or infinite values;
#' 3. `event` holds a value other than `0` and `1`, missing values included;
#' 4. `group` has fewer than two observed levels, or has more than two and
#'    `arms` was not given, or `arms` does not name exactly two levels that
#'    occur in the data;
#' 5. `group` has missing values and `group_missing = "error"` (the default);
#' 6. `tau` is not a positive, finite, non-missing, non-duplicated numeric
#'    vector;
#' 7. any element of `tau` is past the smallest of the two arms' largest
#'    follow-up times. This one is raised by [calc_rmst()] and says so;
#' 8. an arm has no rows, or no event anywhere in follow-up. Also
#'    [calc_rmst()]'s;
#' 9. `conf_level`, `ph_alpha`, `min_denom` or any `digits_` argument is not a
#'    single value of the right kind;
#' 10. the three estimators do not report the same subjects and events, or the
#'    Kaplan-Meier probabilities from the two routes disagree, or the
#'    step-function integral does not reproduce survRM2's arm-0 RMST, or the
#'    assembled table fails one of its own identities. None of these can be
#'    caused by an ordinary input; they are reported as bugs.
#'
#' @param data A data frame (or tibble) with one row per subject.
#' @param group Name of the grouping variable, as a single string. Required:
#'   the comparison needs two arms, and this is the column they come from.
#' @param arms The two levels of `group` to compare, as a vector of length 2,
#'   reference first: `arms[1]` becomes arm 0 and `arms[2]` arm 1, and every
#'   contrast is arm 1 against arm 0. `NULL` (the default) is allowed only when
#'   `group` has exactly two observed levels, which are then used in the order
#'   of the levels.
#' @param time Name of the follow-up-time column, as a single string. Must be
#'   numeric. Factors and character digits are refused: use [prep_surv()],
#'   which converts them and returns a numeric `time_mo`.
#' @param event Name of the event-indicator column, as a single string. Must
#'   hold only `0` (censored) and `1` (event); logical `FALSE`/`TRUE` is
#'   accepted and read as `0`/`1`.
#' @param tau Truncation times, in the unit of `time`, one row of the result
#'   each. Defaults to 12, 36 and 60, i.e. months when `time` is in months. The
#'   same values are used as the Kaplan-Meier time points, so that the three
#'   estimands of a row are truncated at the same place. Every element must be
#'   positive and at most `tau_upper`; the guard is [calc_rmst()]'s. Returned
#'   sorted.
#' @param conf_level Confidence level of every interval, and the source of the
#'   significance threshold `alpha = 1 - conf_level`. Defaults to `0.95`. It is
#'   one argument rather than two on purpose: an interval and a p-value that
#'   are allowed to use different levels can give a row two different verdicts.
#' @param group_missing What to do with rows whose `group` value is missing:
#'   `"error"` (the default) to stop, or `"drop"` to leave them out. Dropping
#'   happens once, before any estimator is called, so all three still see the
#'   same rows. The count comes back in the diagnostics either way. Rows
#'   belonging to a level other than the two in `arms` are always left out and
#'   are counted separately.
#' @param ties Handling of tied event times in the Cox model, passed to
#'   [fit_cox()] and on to [survival::coxph()]. Defaults to `"efron"`.
#' @param zph_transform Function of time used by [survival::cox.zph()], passed
#'   to [fit_cox()]. Defaults to `"km"`, the default of `cox.zph()` itself.
#' @param ph_alpha Level at which the proportional-hazards test is called
#'   violated in the `ph_violated` column. Defaults to `0.05`. It flags a row;
#'   it never changes a model, and it takes no part in the concordance
#'   judgement.
#' @param min_denom Smallest absolute RMST difference, in the unit of `time`,
#'   that may be used as the denominator of `ratio_obs` or `ratio_implied`.
#'   Defaults to `0.25`, about a week when `time` is in months: below that the
#'   ratio is driven by a denominator sitting against zero, so `NA` is returned
#'   with the reason in `het_note` and `gap_hi_mo` is read instead.
#' @param min_n,min_events Smallest number of subjects and of events in an arm
#'   that does not raise a warning, passed to [calc_rmst()]. Default to `10`
#'   and `5`; `0` silences the corresponding warning.
#' @param epv_warn Smallest number of events per model coefficient that does
#'   not raise a warning, passed to [fit_cox()]. Defaults to `10`.
#' @param digits_hr Decimal places for the hazard ratio in `hr_txt`. Defaults
#'   to `2`.
#' @param digits_time Decimal places for a time in months -- RMST, its
#'   difference, the median -- in the `_txt` columns. Defaults to `1`.
#' @param digits_pct Decimal places for the Kaplan-Meier difference in
#'   percentage points in `km_surv_diff_txt`. Defaults to `1`.
#' @param digits_ratio Decimal places for the RMST ratio in `rmst_ratio_txt`.
#'   Defaults to `3`.
#' @param digits_p Decimal places for the `_p_fmt` columns, which read
#'   `"<0.001"` below the corresponding threshold. Defaults to `3`. The numeric
#'   columns are returned unrounded, so every `digits_` argument stays a
#'   presentation choice and never a stored one.
#'
#' @return
#' A plain data frame with one row per element of `tau`, in increasing order of
#' `tau`. Columns that describe the contrast rather than a single `tau` are
#' repeated on every row, so that any single row can be quoted on its own.
#'
#' * the contrast -- `grouping`, `arm0`, `arm1`, `comparison`, `tau`,
#'   `tau_upper`, `n_used`, `n_arm0`, `n_arm1`, `events_arm0`, `events_arm1`;
#' * the hazard ratio -- `hr`, `hr_lcl`, `hr_ucl`, `hr_p`, `hr_p_fmt`,
#'   `hr_txt`, `sig_hr`, `dir_hr`, and the proportional-hazards test `ph_p`,
#'   `ph_p_fmt`, `ph_violated`;
#' * the Kaplan-Meier curves at `tau` -- `km_surv_arm0`, `km_surv_arm0_se`,
#'   `km_surv_arm0_lcl`, `km_surv_arm0_ucl` and the same four for `arm1`;
#'   `km_nrisk_arm0`, `km_nrisk_arm1`;
#' * the Kaplan-Meier difference at `tau` -- `km_surv_diff` (a proportion),
#'   `km_surv_diff_se`, `km_surv_diff_lcl`, `km_surv_diff_ucl`, the same three
#'   in percentage points as `km_surv_diff_pp`, `km_surv_diff_pp_lcl`,
#'   `km_surv_diff_pp_ucl`, then `km_surv_diff_p`, `km_surv_diff_p_fmt`,
#'   `km_surv_diff_txt`, `sig_km`, `dir_km`;
#' * the medians and the log-rank test, for reference -- `km_median_arm0`,
#'   `km_median_arm0_lcl`, `km_median_arm0_ucl`, the same three for `arm1`,
#'   `km_median_diff`, `km_median_txt`, `dir_km_median`, `logrank_chisq`,
#'   `logrank_df`, `logrank_p`, `logrank_p_fmt`;
#' * the RMST -- `rmst_arm0`, `rmst_arm0_lcl`, `rmst_arm0_ucl` and the same
#'   three for `arm1`; `rmst_diff`, `rmst_diff_lcl`, `rmst_diff_ucl`,
#'   `rmst_diff_p`, `rmst_diff_p_fmt`, `rmst_diff_txt`; `rmst_ratio`,
#'   `rmst_ratio_lcl`, `rmst_ratio_ucl`, `rmst_ratio_p`, `rmst_ratio_txt`;
#'   `sig_rmst`, `dir_rmst`;
#' * what a constant hazard ratio would imply -- `rmst_diff_implied` and
#'   `rmst_gap_mo` at this `tau`, then the contrast-level `tau_lo`, `tau_hi`,
#'   `d_lo`, `d_hi`, `di_lo`, `di_hi`, `gap_lo_mo`, `gap_hi_mo`, `ratio_obs`,
#'   `ratio_implied`, `ratio_dev_pct`, `rel_hi_pct` and `het_note`;
#' * the concordance judgement -- `direction_agree`, `signif_agree`,
#'   `n_signif`, `sig_pattern`, `dir_conflict_signif`, `concordant`,
#'   `inconsistent`, `discordance_type`, `verdict`.
#'
#' One attribute is attached:
#'
#' * `compare_estimands` -- a list with `counts` (rows read, rows used, rows
#'   left out for a missing group or for belonging to another level, subjects,
#'   events and censored per arm), `follow_up` (the largest follow-up time in
#'   each arm and `tau_upper`), `same_patients` (the reconciliation: one row
#'   per estimator with the subjects and events it reported, and whether they
#'   all agree), `guard` (whose the tau guard is, the limit and that it is not
#'   relaxable), `rules` (the concordance rules as text, the direction and
#'   significance definitions, the four discordance types and how the boundary
#'   case is settled), `checks` (the invariants enforced on the returned
#'   numbers, including the two cross-route checks), `settings` and `call`.
#'
#' Attributes are dropped by most data-frame verbs, so read them off the object
#' returned by `compare_estimands()` before piping it further; the columns
#' carry everything a table needs, which is why they are columns.
#'
#' @seealso [fit_cox()], [fit_km()] and [calc_rmst()], the three estimators
#'   this function calls and reconciles; [prep_surv()], which builds `time_mo`
#'   and `event_os` from a vital-status column.
#'
#' @examples
#' ## A small simulated cohort. No real patient records are used anywhere in
#' ## this package.
#' set.seed(20260901)
#' site <- factor(rep(c("Stomach", "Small intestine", "Colorectal"), each = 80),
#'                levels = c("Stomach", "Small intestine", "Colorectal"))
#' rate <- c(Stomach = 0.010, "Small intestine" = 0.018,
#'           Colorectal = 0.030)[as.character(site)]
#' t_event  <- stats::rexp(240, rate)
#' t_censor <- stats::runif(240, 24, 140)
#' toy <- data.frame(
#'   site     = site,
#'   time_mo  = round(pmin(t_event, t_censor), 1),
#'   event_os = as.integer(t_event <= t_censor)
#' )
#'
#' ## The three estimands side by side, one row per truncation time.
#' ce <- compare_estimands(toy, "site", arms = c("Stomach", "Colorectal"))
#' ce[, c("comparison", "tau", "hr_txt", "km_surv_diff_txt", "rmst_diff_txt")]
#'
#' ## One hazard ratio stands for the whole of follow-up; the other two do not.
#' ## The proportional-hazards test travels on every row.
#' ce[, c("tau", "hr", "hr_p_fmt", "ph_p_fmt", "ph_violated")]
#'
#' ## The concordance judgement, by the rules in the Details section.
#' ce[, c("tau", "dir_hr", "dir_km", "dir_rmst", "sig_pattern",
#'        "direction_agree", "signif_agree", "discordance_type")]
#' ce$verdict
#'
#' ## What a constant hazard ratio would have implied for the RMST difference,
#' ## and how many months the observed difference misses it by.
#' ce[, c("tau", "rmst_diff", "rmst_diff_implied", "rmst_gap_mo")]
#'
#' ## Every column the table carries.
#' names(ce)
#'
#' ## The diagnostics. same_patients is the invariant that makes the row a
#' ## comparison of estimands rather than of analysis sets.
#' d <- attr(ce, "compare_estimands")
#' d$counts
#' d$follow_up
#' d$same_patients
#' d$guard
#' d$checks
#'
#' ## The concordance rules come back with the result, so a verdict never has
#' ## to be taken on trust.
#' d$rules$direction
#' d$rules$significance
#' d$rules$discordance_type
#' d$rules$boundary_case
#'
#' ## The rules earn their keep when the curves cross. A scenario built so
#' ## that arm B is harmed early and helped later: one hazard ratio stands for
#' ## both halves and fails its own assumption, the Kaplan-Meier difference
#' ## changes sign between 12 and 36 months, and the RMST difference needs a
#' ## long window before it becomes significant.
#' set.seed(20260901)
#' n     <- 150
#' early <- stats::rexp(n, 0.055)
#' late  <- 30 + stats::rexp(n, 0.004)
#' t_b   <- ifelse(stats::runif(n) < 0.35, early, late)
#' t_a   <- stats::rexp(n, 0.020)
#' cens  <- stats::runif(2 * n, 36, 120)
#' cross <- data.frame(trt = factor(rep(c("A", "B"), each = n)),
#'                     raw = c(t_a, t_b))
#' cross$time_mo  <- round(pmin(cross$raw, cens), 1)
#' cross$event_os <- as.integer(cross$raw <= cens)
#' cc <- compare_estimands(cross, "trt")
#' cc[, c("tau", "hr_txt", "km_surv_diff_txt", "rmst_diff_txt", "ph_p_fmt")]
#'
#' ## Three truncation times, three different verdicts, all by the rules
#' ## above: a sign flip whose intervals cover the null is typed apart from a
#' ## real conflict, and never called one.
#' cc[, c("tau", "dir_hr", "dir_km", "dir_rmst", "sig_pattern",
#'        "discordance_type")]
#' cc$verdict
#'
#' ## The same hazard ratio implies a much larger RMST difference than the one
#' ## observed, at every tau. That gap, in months, is what the single number
#' ## cannot carry.
#' cc[, c("tau", "rmst_diff", "rmst_diff_implied", "rmst_gap_mo")]
#'
#' ## A cohort followed for at most 30 months: the tau guard refuses the same
#' ## call. It is calc_rmst()'s guard, not a copy of it, and the message says
#' ## so.
#' recent <- toy
#' recent$event_os[recent$time_mo > 30] <- 0L
#' recent$time_mo <- pmin(recent$time_mo, 30)
#' try(compare_estimands(recent, "site", arms = c("Stomach", "Colorectal")))
#'
#' ## Rows with a missing time, event or group are refused rather than dropped
#' ## by each estimator on its own terms: three columns built on three
#' ## different subsets would not be a comparison of estimands.
#' gappy <- toy
#' gappy$site[c(1, 2, 3)] <- NA
#' try(compare_estimands(gappy, "site", arms = c("Stomach", "Colorectal")))
#'
#' ## Three levels and no `arms`: which two were meant is not guessed.
#' try(compare_estimands(toy, "site"))
#'
#' @export
compare_estimands <- function(data,
                              group,
                              arms          = NULL,
                              time          = "time_mo",
                              event         = "event_os",
                              tau           = c(12, 36, 60),
                              conf_level    = 0.95,
                              group_missing = c("error", "drop"),
                              ties          = c("efron", "breslow", "exact"),
                              zph_transform = c("km", "rank", "identity",
                                                "log"),
                              ph_alpha      = 0.05,
                              min_denom     = 0.25,
                              min_n         = 10L,
                              min_events    = 5L,
                              epv_warn      = 10,
                              digits_hr     = 2L,
                              digits_time   = 1L,
                              digits_pct    = 1L,
                              digits_ratio  = 3L,
                              digits_p      = 3L) {

  cl <- match.call()

  ## -- 0. data --------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("compare_estimands(): `data` must be a data frame, got an object of ",
         "class ", .ps_values(class(data)), ".", call. = FALSE)
  }
  n_in <- nrow(data)
  if (n_in == 0L) {
    stop("compare_estimands(): `data` has 0 rows, there is nothing to compare.",
         call. = FALSE)
  }

  ## -- 1. the argument values themselves ------------------------------------
  if (missing(group)) {
    stop(paste0("compare_estimands(): `group` is required. The three ",
                "estimands are compared on one two-arm contrast, and `group` ",
                "is the column the two arms come from; pass its name as a ",
                'single string, e.g. group = "site".'),
         call. = FALSE)
  }
  .ce_name_arg(time,  "time")
  .ce_name_arg(event, "event")
  .ce_name_arg(group, "group")
  group_missing <- match.arg(group_missing)
  ties          <- match.arg(ties)
  zph_transform <- match.arg(zph_transform)
  .ce_prob_arg(conf_level, "conf_level")
  .ce_prob_arg(ph_alpha,   "ph_alpha")
  if (!is.numeric(min_denom) || length(min_denom) != 1L || is.na(min_denom) ||
      !is.finite(min_denom) || min_denom < 0) {
    stop(sprintf(
      paste0("compare_estimands(): `min_denom` must be a single finite ",
             "non-negative number, got %s."),
      .ps_values(min_denom)),
      call. = FALSE)
  }
  .ce_count_arg(min_n,        "min_n")
  .ce_count_arg(min_events,   "min_events")
  .ce_count_arg(digits_hr,    "digits_hr")
  .ce_count_arg(digits_time,  "digits_time")
  .ce_count_arg(digits_pct,   "digits_pct")
  .ce_count_arg(digits_ratio, "digits_ratio")
  .ce_count_arg(digits_p,     "digits_p", min = 1L)
  if (!is.numeric(epv_warn) || length(epv_warn) != 1L || is.na(epv_warn) ||
      !is.finite(epv_warn) || epv_warn < 0) {
    stop(sprintf(
      paste0("compare_estimands(): `epv_warn` must be a single finite ",
             "non-negative number, got %s."),
      .ps_values(epv_warn)),
      call. = FALSE)
  }
  tau   <- .ce_tau_arg(tau)
  alpha <- 1 - conf_level

  ## -- 2. required columns present, present once, and distinct --------------
  needed <- c(time = time, event = event, group = group)

  miss <- needed[!(needed %in% names(data))]
  if (length(miss) > 0L) {
    stop(sprintf(
      paste0("compare_estimands(): %d required column%s not found in `data`:",
             "\n  %s\n  `data` has %d columns. Check the spelling, or pass ",
             "the column names explicitly."),
      length(miss), if (length(miss) > 1L) "s" else "",
      paste0(sprintf('%s = "%s"', names(miss), unname(miss)), collapse = "\n  "),
      ncol(data)),
      call. = FALSE)
  }
  dup <- needed[needed %in% names(data)[duplicated(names(data))]]
  if (length(dup) > 0L) {
    stop(sprintf(
      paste0("compare_estimands(): %s appear%s more than once among the ",
             "columns of `data`; the three estimators would each silently use ",
             "the first one, and there would be no way to show they used the ",
             "same column. Make the column names unique first."),
      paste0(sprintf('%s = "%s"', names(dup), unname(dup)), collapse = ", "),
      if (length(dup) > 1L) "" else "s"),
      call. = FALSE)
  }
  same <- duplicated(needed) | duplicated(needed, fromLast = TRUE)
  if (any(same)) {
    stop(sprintf(
      paste0("compare_estimands(): %s all name the same column of `data`. The ",
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
      paste0('compare_estimands(): the `time` column "%s" is %s. Factor ',
             "levels and character digits are not numbers, so ",
             "compare_estimands() will not convert them for you: convert the ",
             "column explicitly, or build it with prep_surv(), which returns ",
             "a numeric `time_mo`."),
      time, if (is.factor(tv)) "a factor" else "character"),
      call. = FALSE)
  }
  if (!is.numeric(tv)) {
    stop(sprintf(
      paste0('compare_estimands(): the `time` column "%s" must be numeric, ',
             "got an object of class %s."),
      time, .ps_values(class(tv))),
      call. = FALSE)
  }
  tnum   <- as.numeric(tv)
  n_t_na <- sum(is.na(tnum))
  if (n_t_na > 0L) {
    stop(sprintf(
      paste0('compare_estimands(): the `time` column "%s" has %d missing ',
             "value%s out of %d rows.\n",
             "  compare_estimands() neither drops nor imputes them. Each of ",
             "the three estimators has its own rule for an incomplete row, ",
             "and letting each apply its own would build the hazard ratio, ",
             "the Kaplan-Meier difference and the RMST difference on ",
             "different subsets -- which is a comparison of analysis sets, ",
             "not of estimands. Decide what to do with these rows before ",
             "calling."),
      time, n_t_na, if (n_t_na > 1L) "s" else "", n_in),
      call. = FALSE)
  }
  n_t_inf <- sum(!is.finite(tnum))
  if (n_t_inf > 0L) {
    stop(sprintf(
      paste0('compare_estimands(): the `time` column "%s" has %d infinite ',
             "value%s. Neither an area up to tau nor a survival probability ",
             "at tau can be built from an infinite follow-up time."),
      time, n_t_inf, if (n_t_inf > 1L) "s" else ""),
      call. = FALSE)
  }
  if (any(tnum < 0)) {
    n_neg <- sum(tnum < 0)
    stop(sprintf(
      paste0('compare_estimands(): the `time` column "%s" has %d negative ',
             "value%s (minimum %s). A negative follow-up time is not ",
             "analysable, so compare_estimands() stops instead of computing ",
             "on it."),
      time, n_neg, if (n_neg > 1L) "s" else "", format(min(tnum))),
      call. = FALSE)
  }

  ## -- 4. event indicator ---------------------------------------------------
  evr <- data[[event]]
  if (is.factor(evr) || is.character(evr)) {
    stop(sprintf(
      paste0('compare_estimands(): the `event` column "%s" is %s. ',
             "compare_estimands() will not guess which level means an event: ",
             "convert it to 0 (censored) / 1 (event) explicitly, or build it ",
             "with prep_surv(), which returns a 0/1 `event_os` from a ",
             "vital-status column.\n  values present : %s"),
      event, if (is.factor(evr)) "a factor" else "character", .ps_counts(evr)),
      call. = FALSE)
  }
  if (is.logical(evr)) {
    ev <- as.integer(evr)
  } else if (is.numeric(evr)) {
    ev <- as.numeric(evr)
  } else {
    stop(sprintf(
      paste0('compare_estimands(): the `event` column "%s" must be numeric ',
             "0/1 or logical, got an object of class %s."),
      event, .ps_values(class(evr))),
      call. = FALSE)
  }
  bad_ev <- !(ev %in% c(0, 1))
  if (any(bad_ev)) {
    stop(sprintf(
      paste0('compare_estimands(): the `event` column "%s" has %d of %d ',
             "row%s whose value is neither 0 (censored) nor 1 (event).\n",
             "  offending values : %s\n",
             "  A missing indicator (<NA>) counts as an offending value: ",
             "compare_estimands() never guesses whether a subject had the ",
             "event, and it will not let the three estimands be built on ",
             "three different sets of complete rows. Use prep_surv() to build ",
             "a 0/1 indicator from a vital-status column."),
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
      paste0('compare_estimands(): the `group` column "%s" is of class %s, ',
             "which cannot be used as a grouping variable. Convert it to a ",
             "factor first."),
      group, .ps_values(class(g_raw))),
      call. = FALSE)
  }

  n_g_na    <- sum(is.na(g))
  n_lv_all  <- nlevels(g)
  g         <- droplevels(g)
  n_lv_drop <- n_lv_all - nlevels(g)

  if (n_g_na > 0L && group_missing == "error") {
    stop(sprintf(
      paste0('compare_estimands(): the `group` column "%s" has %d missing ',
             "value%s out of %d rows.\n",
             "  Those rows belong to neither arm. compare_estimands() will ",
             "not leave them out on its own, because the whole point of the ",
             "table is that the hazard ratio, the Kaplan-Meier difference and ",
             "the RMST difference describe the same patients, and a rule ",
             "applied silently is a rule that cannot be checked.\n",
             '  Pass group_missing = "drop" to leave them out -- once, before ',
             "any estimator sees the data, so all three still see the same ",
             "rows, and the count comes back in the diagnostics -- or filter ",
             "them out before calling."),
      group, n_g_na, if (n_g_na > 1L) "s" else "", n_in),
      call. = FALSE)
  }

  lv <- levels(droplevels(g[!is.na(g)]))
  if (length(lv) < 2L) {
    stop(sprintf(
      paste0('compare_estimands(): the `group` column "%s" has %d observed ',
             "level%s (%s) in the %d rows with a non-missing group; a two-arm ",
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
        paste0('compare_estimands(): the `group` column "%s" has %d observed ',
               "levels (%s) and `arms` was not given.\n",
               "  The three estimands are compared on one two-arm contrast: ",
               "which two of the %d levels were meant is not guessed, because ",
               "the answer would change every estimate, the tau limit and the ",
               "concordance verdict alike.\n",
               "  Name the two levels, reference first:\n",
               '    compare_estimands(data, group = "%s", arms = c("%s", "%s"))\n',
               "  For every pair, loop over combn(levels(data[[%s]]), 2) and ",
               "stack the results; each pair then gets its own tau guard and ",
               "its own verdict."),
        group, length(lv), .ps_values(lv), length(lv), group, lv[1], lv[2],
        paste0('"', group, '"')),
        call. = FALSE)
    }
    arms <- lv
  } else {
    if (is.factor(arms)) arms <- as.character(arms)
    if (!is.atomic(arms) || length(arms) != 2L) {
      stop(sprintf(
        paste0('compare_estimands(): `arms` must name exactly 2 levels of ',
               '"%s", reference first, got %s of length %d.\n',
               "  Two is not a default that can be raised: a hazard ratio, a ",
               "difference in survival probability and a difference in ",
               "restricted mean survival time are all defined between two ",
               "arms. Loop over combn(levels(x), 2) for all pairs."),
        group, .ps_values(class(arms)), length(arms)),
        call. = FALSE)
    }
    arms <- as.character(arms)
    if (anyNA(arms)) {
      stop(sprintf(
        paste0("compare_estimands(): `arms` has a missing element. Both arms ",
               'must name a level of "%s".'),
        group),
        call. = FALSE)
    }
    if (arms[1] == arms[2]) {
      stop(sprintf(
        paste0('compare_estimands(): `arms` names the same level twice ',
               '("%s"). An arm compared with itself has a hazard ratio of 1 ',
               "and differences of 0 by construction, and is not a contrast."),
        arms[1]),
        call. = FALSE)
    }
    not_found <- arms[!(arms %in% lv)]
    if (length(not_found) > 0L) {
      stop(sprintf(
        paste0("compare_estimands(): %d of the 2 values of `arms` %s not ",
               'occur in the `group` column "%s":\n',
               "  not found       : %s\n  observed levels : %s\n",
               "  Level labels are matched exactly, including case and ",
               "spacing."),
        length(not_found), if (length(not_found) > 1L) "do" else "does",
        group, .ps_values(not_found), .ps_values(lv)),
        call. = FALSE)
    }
  }

  ## -- 6. the analysis set: built once, handed to all three -----------------
  g_chr   <- as.character(g)
  in_arms <- !is.na(g_chr) & g_chr %in% arms
  n_other <- sum(!is.na(g_chr) & !(g_chr %in% arms))

  ana_t <- tnum[in_arms]
  ana_e <- ev[in_arms]
  ana_g <- factor(g_chr[in_arms], levels = arms)
  n_used <- length(ana_t)

  ana <- data.frame(ana_t, ana_e, ana_g, stringsAsFactors = FALSE)
  names(ana) <- c(time, event, group)

  arm <- ifelse(ana_g == arms[2], 1L, ifelse(ana_g == arms[1], 0L, NA_integer_))
  n0  <- sum(ana_g == arms[1])
  n1  <- sum(ana_g == arms[2])
  arm_ok <- !anyNA(arm) && all(arm %in% c(0L, 1L)) &&
    n0 > 0L && n1 > 0L && n0 + n1 == n_used &&
    all(arm[ana_g == arms[1]] == 0L) && all(arm[ana_g == arms[2]] == 1L) &&
    all(ana_g[arm == 0L] == arms[1]) && all(ana_g[arm == 1L] == arms[2])
  if (!arm_ok) {
    stop(paste0("compare_estimands(): the 0/1 arm indicator does not agree ",
                "with the grouping column in both directions, which should be ",
                "impossible by construction. Do not use this result; report ",
                "it as a bug in compare_estimands()."),
         call. = FALSE)
  }
  ev0 <- sum(ana_e[arm == 0L] == 1L)
  ev1 <- sum(ana_e[arm == 1L] == 1L)

  ## -- 7. RMST first: the tau guard is calc_rmst()'s and is not duplicated ---
  ## calc_rmst() also raises the empty-arm and no-event-in-an-arm errors, so
  ## nothing below has to repeat them either.
  rm_tab <- calc_rmst(ana, group = group, arms = arms, time = time,
                      event = event, tau = tau, conf_level = conf_level,
                      group_missing = "error", min_n = min_n,
                      min_events = min_events, digits_rmst = digits_time,
                      digits_ratio = digits_ratio, digits_p = digits_p)
  rm_diag   <- attr(rm_tab, "calc_rmst")
  tau_upper <- rm_diag$follow_up[["tau_upper"]]
  max_fu0   <- rm_diag$follow_up[["max_fu_arm0"]]
  max_fu1   <- rm_diag$follow_up[["max_fu_arm1"]]

  ## -- 8. the unadjusted Cox model, with its proportional-hazards test ------
  cox_tab <- fit_cox(ana, covariates = group, time = time, event = event,
                     ties = ties, conf_level = conf_level,
                     na_action = "fail", empty_levels = "error",
                     zph_transform = zph_transform, ph_alpha = ph_alpha,
                     epv_warn = epv_warn, digits_hr = digits_hr,
                     digits_p = digits_p)
  cox_diag <- attr(cox_tab, "fit_cox")
  if (nrow(cox_tab) != 2L || sum(cox_tab$is_reference) != 1L ||
      !identical(cox_tab$level, arms)) {
    stop(paste0("compare_estimands(): fit_cox() did not return one reference ",
                "row and one comparator row for the two arms, which should be ",
                "impossible on a two-level factor. Do not use this result; ",
                "report it as a bug in compare_estimands()."),
         call. = FALSE)
  }
  cmp <- which(!cox_tab$is_reference)
  ref <- which(cox_tab$is_reference)

  ## -- 9. the Kaplan-Meier curves, medians and log-rank test ----------------
  km_tab <- fit_km(ana, time = time, event = event, group = group,
                   times = tau, times_beyond = "error",
                   group_missing = "error", conf_type = "log",
                   conf_level = conf_level, digits_median = digits_time)
  km_diag <- attr(km_tab, "fit_km")
  km_pts  <- attr(km_tab, "km_times")
  if (nrow(km_tab) != 2L || !identical(km_tab$level, arms) ||
      is.null(km_pts) || nrow(km_pts) != 2L * length(tau)) {
    stop(paste0("compare_estimands(): fit_km() did not return one row per arm ",
                "and one time-point row per arm and tau, which should be ",
                "impossible on a two-level factor. Do not use this result; ",
                "report it as a bug in compare_estimands()."),
         call. = FALSE)
  }
  k0 <- match(paste(arms[1], tau, sep = "\r"),
              paste(km_pts$level, km_pts$time, sep = "\r"))
  k1 <- match(paste(arms[2], tau, sep = "\r"),
              paste(km_pts$level, km_pts$time, sep = "\r"))
  if (anyNA(k0) || anyNA(k1)) {
    stop(paste0("compare_estimands(): fit_km() did not return a row for every ",
                "arm and tau. Do not use this result; report it as a bug in ",
                "compare_estimands()."),
         call. = FALSE)
  }

  ## -- 10. the invariant: the three estimands describe the same patients ----
  ## Every column below has to hold one value. If it does not, one of the three
  ## estimators used a different subset, and the row would put three numbers
  ## from three cohorts next to each other -- a comparison of analysis sets
  ## wearing the clothes of a comparison of estimands. That is an error here,
  ## and it is checked before any estimate is read.
  same_patients <- data.frame(
    source      = c("compare_estimands", "calc_rmst", "fit_cox", "fit_km"),
    rows_used   = c(n_used,
                    rm_diag$counts[["rows_used"]],
                    cox_diag$counts[["rows_used"]],
                    km_diag$counts[["rows_used"]]),
    n_arm0      = c(n0, rm_tab$n_arm0[1], cox_tab$n[ref], km_tab$n[1]),
    n_arm1      = c(n1, rm_tab$n_arm1[1], cox_tab$n[cmp], km_tab$n[2]),
    events_arm0 = c(ev0, rm_tab$events_arm0[1], cox_tab$n_event[ref],
                    km_tab$events[1]),
    events_arm1 = c(ev1, rm_tab$events_arm1[1], cox_tab$n_event[cmp],
                    km_tab$events[2]),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  sp_cols  <- c("rows_used", "n_arm0", "n_arm1", "events_arm0", "events_arm1")
  sp_one   <- vapply(sp_cols,
                     function(v) length(unique(same_patients[[v]])) == 1L,
                     logical(1))
  same_patients$agrees <- all(sp_one)
  if (!all(sp_one)) {
    lines <- vapply(sp_cols[!sp_one], function(v) sprintf(
      "%-11s : %s", v,
      paste(sprintf("%s = %s", same_patients$source, same_patients[[v]]),
            collapse = ", ")), character(1))
    stop(sprintf(
      paste0("compare_estimands(): the three estimands were not computed on ",
             "the same patients, so they cannot be put side by side.\n  %s\n",
             "  All three are called on one analysis set of %d rows, so this ",
             "should be impossible. Do not use this result; report it as a ",
             "bug in compare_estimands()."),
      paste(lines, collapse = "\n  "), n_used),
      call. = FALSE)
  }

  ## -- 11. one Kaplan-Meier fit per arm, for the standard error and the -----
  ##        step function fit_km() does not return
  s_at <- function(keep) {
    d <- data.frame(.ce_t = ana_t[keep], .ce_e = ana_e[keep])
    f <- survival::survfit(survival::Surv(.ce_t, .ce_e) ~ 1, data = d,
                           conf.type = "log", conf.int = conf_level)
    s <- summary(f, times = tau, extend = TRUE)
    list(fit    = f,
         surv   = unname(s$surv),
         se     = unname(s$std.err),
         lcl    = unname(s$lower),
         ucl    = unname(s$upper),
         n_risk = as.integer(s$n.risk))
  }
  a0 <- s_at(arm == 0L)
  a1 <- s_at(arm == 1L)

  if (anyNA(a0$surv) || anyNA(a1$surv) || anyNA(a0$se) || anyNA(a1$se)) {
    stop(sprintf(
      paste0("compare_estimands(): the Kaplan-Meier survival probability or ",
             "its standard error is missing at one of the requested ",
             "truncation times (tau = %s), so the difference at that tau has ",
             "no interval. Do not use this result; report it as a bug in ",
             "compare_estimands()."),
      paste(.ce_num(tau), collapse = ", ")),
      call. = FALSE)
  }

  ## Cross-route check: the probabilities used here must be the ones fit_km()
  ## reports, or the KM column and the KM difference would come from two
  ## different curves.
  km_route <- .ce_agree(a0$surv, km_pts$surv[k0], 1e-10) &&
    .ce_agree(a1$surv, km_pts$surv[k1], 1e-10) &&
    identical(a0$n_risk, as.integer(km_pts$n_risk[k0])) &&
    identical(a1$n_risk, as.integer(km_pts$n_risk[k1]))
  if (!km_route) {
    stop(paste0("compare_estimands(): the per-arm Kaplan-Meier fit used for ",
                "the standard error does not reproduce the survival ",
                "probabilities fit_km() reports at the same tau. Do not use ",
                "this result; report it as a bug in compare_estimands()."),
         call. = FALSE)
  }

  ## -- 12. the Kaplan-Meier difference at each tau --------------------------
  z_q      <- stats::qnorm(1 - alpha / 2)
  d_surv   <- a1$surv - a0$surv
  se_dsurv <- sqrt(a0$se^2 + a1$se^2)
  if (any(se_dsurv <= 0)) {
    bad <- tau[se_dsurv <= 0]
    stop(sprintf(
      paste0("compare_estimands(): the standard error of the Kaplan-Meier ",
             "difference is 0 at tau = %s, so the difference has no interval ",
             "and no p-value there.\n",
             "  That happens when neither arm has had an event on or before ",
             "that tau: both curves are still at 1 and the difference is ",
             "exactly 0. There is nothing to compare at that truncation time. ",
             "Raise `tau`, or check that the event column is the right one ",
             "(events in arm 0 = %d, in arm 1 = %d over the whole of ",
             "follow-up)."),
      paste(.ce_num(bad), collapse = ", "), ev0, ev1),
      call. = FALSE)
  }
  d_lcl <- d_surv - z_q * se_dsurv
  d_ucl <- d_surv + z_q * se_dsurv
  d_p   <- 2 * stats::pnorm(-abs(d_surv / se_dsurv))

  ## -- 13. the RMST difference a constant hazard ratio would imply ----------
  ## S1(t) = S0(t)^theta, with S0 the arm-0 Kaplan-Meier curve and theta the
  ## estimated hazard ratio. Both curves are integrated with the same step
  ## function, and the integral of S0 itself is checked against the arm-0 RMST
  ## survRM2 computed independently: if the two disagree, the integration is
  ## wrong and everything derived from it is meaningless.
  hr_est    <- cox_tab$hr[cmp]
  k         <- length(tau)
  implied   <- numeric(k)
  rmst0_int <- numeric(k)
  for (i in seq_len(k)) {
    rmst0_int[i] <- .ce_rmst_step(a0$fit$time, a0$fit$surv, tau[i])
    implied[i]   <- .ce_rmst_step(a0$fit$time, a0$fit$surv^hr_est, tau[i]) -
      rmst0_int[i]
  }
  int_ok <- .ce_agree(rmst0_int, rm_tab$rmst_arm0, 1e-6)
  if (!int_ok) {
    stop(sprintf(
      paste0("compare_estimands(): the step-function integral of the arm-0 ",
             "Kaplan-Meier curve does not reproduce the arm-0 restricted mean ",
             "survival time survRM2 computed (largest gap %s at tau = %s, ",
             "tolerance 1e-6).\n",
             "  The same integration is what turns the estimated hazard ratio ",
             "into the RMST difference it would imply, so if it cannot ",
             "reproduce a number computed independently, rmst_diff_implied ",
             "and every gap derived from it are wrong. Do not use this ",
             "result; report it as a bug in compare_estimands()."),
      .ce_num(max(abs(rmst0_int - rm_tab$rmst_arm0))),
      .ce_num(tau[which.max(abs(rmst0_int - rm_tab$rmst_arm0))])),
      call. = FALSE)
  }
  gap_mo <- rm_tab$rmst_diff - implied

  ## -- 14. the concordance judgement, by the rules in the documentation -----
  dir_lab <- function(worse_is_arm1) {
    ifelse(is.na(worse_is_arm1), NA_character_,
           ifelse(worse_is_arm1, "arm1_worse", "arm0_worse"))
  }

  hr_lcl <- cox_tab$ci_low[cmp]
  hr_ucl <- cox_tab$ci_high[cmp]
  hr_p   <- cox_tab$p_value[cmp]
  ph_p   <- cox_tab$zph_p_var[cmp]

  dir_hr   <- rep(dir_lab(hr_est > 1), k)
  dir_km   <- dir_lab(d_surv < 0)
  dir_rmst <- dir_lab(rm_tab$rmst_diff < 0)

  sig_hr   <- rep((hr_lcl > 1) | (hr_ucl < 1), k)
  sig_km   <- (d_lcl > 0) | (d_ucl < 0)
  sig_rmst <- (rm_tab$rmst_diff_lcl > 0) | (rm_tab$rmst_diff_ucl < 0)

  direction_agree <- (dir_hr == dir_km) & (dir_km == dir_rmst)
  signif_agree    <- (sig_hr == sig_km) & (sig_km == sig_rmst)
  n_signif        <- as.integer(sig_hr) + as.integer(sig_km) +
    as.integer(sig_rmst)
  concordant      <- direction_agree & signif_agree
  inconsistent    <- !concordant

  dir_conflict_signif <-
    (sig_hr & sig_km   & dir_hr != dir_km)   |
    (sig_hr & sig_rmst & dir_hr != dir_rmst) |
    (sig_km & sig_rmst & dir_km != dir_rmst)

  discordance_type <- ifelse(
    concordant, "concordant",
    ifelse(!direction_agree & dir_conflict_signif,
           "direction_conflict_significant",
           ifelse(!direction_agree & !dir_conflict_signif,
                  "direction_conflict_ci_covers_null",
                  "significance_only")))

  sig_word    <- function(x) ifelse(x, "sig", "ns")
  sig_pattern <- sprintf("HR %s / KM %s / RMST %s",
                         sig_word(sig_hr), sig_word(sig_km),
                         sig_word(sig_rmst))
  verdict <- ifelse(
    concordant & n_signif == 3L,
    "concordant: all three significant, same direction",
    ifelse(
      concordant & n_signif == 0L,
      "concordant: all three non-significant, same direction",
      ifelse(
        discordance_type == "direction_conflict_significant",
        sprintf(paste0("DISCORDANT (direction): two significant estimands ",
                       "point to different arms (%s)"), sig_pattern),
        ifelse(
          discordance_type == "direction_conflict_ci_covers_null",
          sprintf(paste0("minor: directions differ but the conflicting ",
                         "estimands are not both significant (%s)"),
                  sig_pattern),
          sprintf(paste0("DISCORDANT (significance): %d of 3 significant, ",
                         "same direction (%s)"),
                  n_signif, sig_pattern)))))

  ## -- 15. does one hazard ratio hide the time scale ------------------------
  i_lo   <- which.min(tau)
  i_hi   <- which.max(tau)
  tau_lo <- tau[i_lo]
  tau_hi <- tau[i_hi]
  d_lo   <- rm_tab$rmst_diff[i_lo]
  d_hi   <- rm_tab$rmst_diff[i_hi]
  di_lo  <- implied[i_lo]
  di_hi  <- implied[i_hi]

  ratio_obs     <- if (abs(d_lo)  >= min_denom) d_hi  / d_lo  else NA_real_
  ratio_implied <- if (abs(di_lo) >= min_denom) di_hi / di_lo else NA_real_
  ratio_dev_pct <- 100 * (ratio_obs / ratio_implied - 1)
  rel_hi_pct    <- if (abs(di_hi) >= min_denom) 100 * (d_hi / di_hi - 1) else
    NA_real_
  gap_lo_mo <- d_lo - di_lo
  gap_hi_mo <- d_hi - di_hi

  het_note <- if (is.na(ratio_obs) && is.na(ratio_implied)) {
    sprintf(paste0("the observed and the implied RMST difference at tau = %s ",
                   "are both smaller than %s in absolute value, so both ",
                   "amplification ratios would divide by a number against ",
                   "zero; read gap_hi_mo instead"),
            .ce_num(tau_lo), .ce_num(min_denom))
  } else if (is.na(ratio_obs)) {
    sprintf(paste0("the observed RMST difference at tau = %s is smaller than ",
                   "%s in absolute value, so ratio_obs would divide by a ",
                   "number against zero; read gap_hi_mo instead"),
            .ce_num(tau_lo), .ce_num(min_denom))
  } else if (is.na(ratio_implied)) {
    sprintf(paste0("the RMST difference implied by a constant hazard ratio at ",
                   "tau = %s is smaller than %s in absolute value, so ",
                   "ratio_implied would divide by a number against zero; read ",
                   "gap_hi_mo instead"),
            .ce_num(tau_lo), .ce_num(min_denom))
  } else if (tau_lo == tau_hi) {
    paste0("a single tau was requested, so the smallest and the largest ",
           "coincide and both amplification ratios are 1 by construction; ",
           "read gap_hi_mo instead")
  } else {
    NA_character_
  }

  ## -- 16. assemble ---------------------------------------------------------
  f_t  <- function(x) sprintf("%.*f", digits_time, x)
  f_pp <- function(x) sprintf("%.*f", digits_pct, x)

  med0 <- km_tab$median[1]; med1 <- km_tab$median[2]
  med_diff <- med1 - med0

  out <- data.frame(
    grouping             = group,
    arm0                 = arms[1],
    arm1                 = arms[2],
    comparison           = paste0(arms[2], " vs ", arms[1]),
    tau                  = tau,
    tau_upper            = tau_upper,
    n_used               = n_used,
    n_arm0               = n0,
    n_arm1               = n1,
    events_arm0          = ev0,
    events_arm1          = ev1,

    hr                   = hr_est,
    hr_lcl               = hr_lcl,
    hr_ucl               = hr_ucl,
    hr_p                 = hr_p,
    hr_p_fmt             = cox_tab$p_fmt[cmp],
    hr_txt               = cox_tab$hr_txt[cmp],
    sig_hr               = sig_hr,
    dir_hr               = dir_hr,
    ph_p                 = ph_p,
    ph_p_fmt             = .ce_fmt_p(ph_p, digits_p),
    ph_violated          = ph_p < ph_alpha,

    km_surv_arm0         = a0$surv,
    km_surv_arm0_se      = a0$se,
    km_surv_arm0_lcl     = a0$lcl,
    km_surv_arm0_ucl     = a0$ucl,
    km_surv_arm1         = a1$surv,
    km_surv_arm1_se      = a1$se,
    km_surv_arm1_lcl     = a1$lcl,
    km_surv_arm1_ucl     = a1$ucl,
    km_nrisk_arm0        = a0$n_risk,
    km_nrisk_arm1        = a1$n_risk,
    km_surv_diff         = d_surv,
    km_surv_diff_se      = se_dsurv,
    km_surv_diff_lcl     = d_lcl,
    km_surv_diff_ucl     = d_ucl,
    km_surv_diff_pp      = 100 * d_surv,
    km_surv_diff_pp_lcl  = 100 * d_lcl,
    km_surv_diff_pp_ucl  = 100 * d_ucl,
    km_surv_diff_p       = d_p,
    km_surv_diff_p_fmt   = .ce_fmt_p(d_p, digits_p),
    km_surv_diff_txt     = paste0(f_pp(100 * d_surv), " (", f_pp(100 * d_lcl),
                                  " to ", f_pp(100 * d_ucl), ")"),
    sig_km               = sig_km,
    dir_km               = dir_km,

    km_median_arm0       = med0,
    km_median_arm0_lcl   = km_tab$median_lcl[1],
    km_median_arm0_ucl   = km_tab$median_ucl[1],
    km_median_arm1       = med1,
    km_median_arm1_lcl   = km_tab$median_lcl[2],
    km_median_arm1_ucl   = km_tab$median_ucl[2],
    km_median_diff       = med_diff,
    km_median_txt        = paste0(
      if (is.na(med1)) "not reached" else f_t(med1), " vs ",
      if (is.na(med0)) "not reached" else f_t(med0)),
    dir_km_median        = dir_lab(med_diff < 0),
    logrank_chisq        = km_tab$logrank_chisq[1],
    logrank_df           = km_tab$logrank_df[1],
    logrank_p            = km_tab$logrank_p[1],
    logrank_p_fmt        = .ce_fmt_p(km_tab$logrank_p[1], digits_p),

    rmst_arm0            = rm_tab$rmst_arm0,
    rmst_arm0_lcl        = rm_tab$rmst_arm0_lcl,
    rmst_arm0_ucl        = rm_tab$rmst_arm0_ucl,
    rmst_arm1            = rm_tab$rmst_arm1,
    rmst_arm1_lcl        = rm_tab$rmst_arm1_lcl,
    rmst_arm1_ucl        = rm_tab$rmst_arm1_ucl,
    rmst_diff            = rm_tab$rmst_diff,
    rmst_diff_lcl        = rm_tab$rmst_diff_lcl,
    rmst_diff_ucl        = rm_tab$rmst_diff_ucl,
    rmst_diff_p          = rm_tab$rmst_diff_p,
    rmst_diff_p_fmt      = rm_tab$rmst_diff_p_fmt,
    rmst_diff_txt        = rm_tab$rmst_diff_txt,
    rmst_ratio           = rm_tab$rmst_ratio,
    rmst_ratio_lcl       = rm_tab$rmst_ratio_lcl,
    rmst_ratio_ucl       = rm_tab$rmst_ratio_ucl,
    rmst_ratio_p         = rm_tab$rmst_ratio_p,
    rmst_ratio_txt       = rm_tab$rmst_ratio_txt,
    sig_rmst             = sig_rmst,
    dir_rmst             = dir_rmst,

    rmst_diff_implied    = implied,
    rmst_gap_mo          = gap_mo,
    tau_lo               = tau_lo,
    tau_hi               = tau_hi,
    d_lo                 = d_lo,
    d_hi                 = d_hi,
    di_lo                = di_lo,
    di_hi                = di_hi,
    gap_lo_mo            = gap_lo_mo,
    gap_hi_mo            = gap_hi_mo,
    ratio_obs            = ratio_obs,
    ratio_implied        = ratio_implied,
    ratio_dev_pct        = ratio_dev_pct,
    rel_hi_pct           = rel_hi_pct,
    het_note             = het_note,

    direction_agree      = direction_agree,
    signif_agree         = signif_agree,
    n_signif             = n_signif,
    sig_pattern          = sig_pattern,
    dir_conflict_signif  = dir_conflict_signif,
    concordant           = concordant,
    inconsistent         = inconsistent,
    discordance_type     = discordance_type,
    verdict              = verdict,

    stringsAsFactors     = FALSE,
    row.names            = NULL
  )

  ## -- 17. diagnostics ------------------------------------------------------
  need <- c(out$hr, out$hr_lcl, out$hr_ucl, out$hr_p, out$ph_p,
            out$km_surv_diff, out$km_surv_diff_lcl, out$km_surv_diff_ucl,
            out$km_surv_diff_p, out$rmst_diff, out$rmst_diff_lcl,
            out$rmst_diff_ucl, out$rmst_diff_p, out$rmst_diff_implied)
  inv <- c(
    "every tau is within the follow-up of both arms" =
      isTRUE(all(out$tau <= tau_upper)),
    "the three estimators used the same subjects and the same events" =
      isTRUE(all(sp_one)),
    "the arm indicator agrees with the grouping column in both directions" =
      isTRUE(arm_ok),
    "the two arms partition the rows used" =
      n0 + n1 == n_used && n_used + n_g_na + n_other == n_in,
    "fit_km() and the per-arm fit give the same survival probabilities" =
      isTRUE(km_route),
    "the step-function integral reproduces survRM2's arm-0 RMST" =
      isTRUE(int_ok),
    "no estimate, interval or p-value entering the judgement is missing" =
      !anyNA(need),
    "the Kaplan-Meier difference is arm 1 minus arm 0 at tau" =
      .ce_agree(out$km_surv_diff, out$km_surv_arm1 - out$km_surv_arm0),
    "the RMST difference is arm 1 minus arm 0 up to tau" =
      .ce_agree(out$rmst_diff, out$rmst_arm1 - out$rmst_arm0),
    "the gap is the observed RMST difference minus the implied one" =
      .ce_agree(out$rmst_gap_mo, out$rmst_diff - out$rmst_diff_implied),
    "each confidence interval and its own p-value give the same verdict" =
      isTRUE(all(out$sig_hr   == (out$hr_p           < alpha))) &&
      isTRUE(all(out$sig_km   == (out$km_surv_diff_p < alpha))) &&
      isTRUE(all(out$sig_rmst == (out$rmst_diff_p    < alpha))),
    "a concordant row has 0 or 3 significant estimands" =
      isTRUE(all(!out$concordant | out$n_signif %in% c(0L, 3L))),
    "discordance_type is one of its four values and matches concordant" =
      isTRUE(all(out$discordance_type %in%
                   c("concordant", "direction_conflict_significant",
                     "direction_conflict_ci_covers_null",
                     "significance_only"))) &&
      isTRUE(all((out$discordance_type == "concordant") == out$concordant)),
    "one row per tau, no tau twice" =
      nrow(out) == k && !anyDuplicated(out$tau))
  if (!all(inv)) {
    stop(sprintf(
      paste0("compare_estimands(): the returned table does not add up, which ",
             "should be impossible by construction.\n  failed check%s : %s\n",
             "  Do not use this result; report it as a bug in ",
             "compare_estimands()."),
      if (sum(!inv) > 1L) "s" else "",
      paste(names(inv)[!inv], collapse = "; ")),
      call. = FALSE)
  }

  attr(out, "compare_estimands") <- list(
    counts = c(rows_in            = n_in,
               rows_used          = n_used,
               rows_group_missing = n_g_na,
               rows_group_dropped = if (group_missing == "drop") n_g_na else 0L,
               rows_other_levels  = n_other,
               n_arm0             = n0,
               n_arm1             = n1,
               events_arm0        = ev0,
               events_arm1        = ev1,
               censored_arm0      = n0 - ev0,
               censored_arm1      = n1 - ev1),
    follow_up = c(max_fu_arm0 = max_fu0,
                  max_fu_arm1 = max_fu1,
                  tau_upper   = tau_upper),
    same_patients = same_patients,
    guard = list(
      owner     = "calc_rmst(), called first and not duplicated here",
      rule      = paste0("tau <= min over arms of (largest follow-up time in ",
                         "that arm)"),
      tau_upper = tau_upper,
      boundary  = "inclusive: a tau exactly equal to tau_upper is estimated",
      relaxable = FALSE,
      note      = paste0("the same tau is used for the Kaplan-Meier time ",
                         "points, and fit_km() applies the same rule to them, ",
                         "so the three columns of a row are truncated at the ",
                         "same place")),
    rules = list(
      alpha        = alpha,
      prespecified = paste0("fixed in R/compare_estimands.R and in its ",
                            "documentation; they do not depend on the data ",
                            "and no argument tunes them"),
      direction = c(
        hr             = "dir_hr = arm1_worse if hr > 1, else arm0_worse",
        km             = paste0("dir_km = arm1_worse if km_surv_diff < 0, ",
                                "else arm0_worse"),
        rmst           = paste0("dir_rmst = arm1_worse if rmst_diff < 0, ",
                                "else arm0_worse"),
        ties           = paste0("the comparisons are strict, so an exact tie ",
                                "reads as arm0_worse"),
        reference_only = paste0("dir_km_median is reported for reference and ",
                                "never enters the judgement")),
      significance = c(
        level = sprintf("alpha = 1 - conf_level = %s", .ce_num(alpha)),
        hr    = "sig_hr: the interval of the hazard ratio excludes 1",
        km    = paste0("sig_km: the interval of the Kaplan-Meier difference ",
                       "excludes 0"),
        rmst  = paste0("sig_rmst: the interval of the RMST difference ",
                       "excludes 0"),
        note  = paste0("judged from the interval; an invariant checks that ",
                       "each interval agrees with its own p < alpha")),
      agreement = c(
        direction_agree = "all three direction labels are the same",
        signif_agree    = "all three significance verdicts are the same",
        n_signif        = "how many of the three are significant, 0 to 3",
        concordant      = "direction_agree & signif_agree",
        inconsistent    = "!concordant"),
      discordance_type = c(
        concordant = "concordant is TRUE",
        direction_conflict_significant = paste0(
          "the directions differ and some pair of estimands pointing at ",
          "different arms are both significant"),
        direction_conflict_ci_covers_null = paste0(
          "the directions differ but no such pair is both significant: every ",
          "sign flip involves an estimand whose interval covers the null"),
        significance_only = paste0(
          "the directions all agree and only the significance verdicts ",
          "differ"),
        order = paste0("tested in the order listed; dir_conflict_signif is ",
                       "the pair test")),
      boundary_case = paste0(
        "two estimands that point at different arms while at least one of the ",
        "two has an interval covering the null are NOT counted as a ",
        "contradiction: the row is inconsistent, but its type is ",
        "direction_conflict_ci_covers_null and its verdict begins with ",
        "'minor:'. A difference sitting against zero has a sign, and that ",
        "sign is noise"),
      implied_rmst = paste0(
        "under a constant hazard ratio theta the comparator curve is ",
        "S1(t) = S0(t)^theta, with S0 the arm-0 Kaplan-Meier curve; ",
        "rmst_diff_implied is the area under S0^theta minus the area under ",
        "S0, both integrated from 0 to tau as step functions, and ",
        "rmst_gap_mo is the observed difference minus that")),
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
      ties          = ties,
      zph_transform = zph_transform,
      ph_alpha      = ph_alpha,
      min_denom     = min_denom,
      min_n         = min_n,
      min_events    = min_events,
      epv_warn      = epv_warn,
      digits_hr     = digits_hr,
      digits_time   = digits_time,
      digits_pct    = digits_pct,
      digits_ratio  = digits_ratio,
      digits_p      = digits_p,
      km_conf_type  = "log",
      estimators    = c(
        hr   = paste0("unadjusted Cox model on the two arms, fit_cox(), with ",
                      "survival::cox.zph"),
        km   = paste0("Kaplan-Meier probability at tau, fit_km(); the ",
                      "difference is arm 1 minus arm 0 with a normal ",
                      "interval on the sum of the two Greenwood variances"),
        rmst = paste0("restricted mean survival time to tau, calc_rmst(), ",
                      "i.e. survRM2::rmst2() unadjusted")),
      adjustment    = paste0("none: all three estimands are unadjusted, so a ",
                             "disagreement between them is a property of the ",
                             "estimands and not of a modelling choice"),
      random        = paste0("none: every estimator called here is closed ",
                             "form and uses no random numbers")),
    call = cl)
  out
}


# -- internal helpers, not exported ------------------------------------------
# .ps_values() and .ps_counts() are defined in R/prep_surv.R and reused here.

#' @noRd
.ce_name_arg <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("compare_estimands(): `%s` must be a single non-missing string.",
                 arg),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.ce_prob_arg <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    stop(sprintf(
      paste0("compare_estimands(): `%s` must be a single number strictly ",
             "between 0 and 1, got %s."),
      arg, .ps_values(x)),
      call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.ce_count_arg <- function(x, arg, min = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != as.integer(x) || x < min) {
    stop(sprintf(
      "compare_estimands(): `%s` must be a single whole number >= %d, got %s.",
      arg, min, .ps_values(x)),
      call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.ce_tau_arg <- function(tau) {
  if (!is.numeric(tau) || length(tau) == 0L) {
    stop(sprintf(
      paste0("compare_estimands(): `tau` must be a numeric vector of at least ",
             "one truncation time, got %s of length %d. The same values are ",
             "used for the RMST window and for the Kaplan-Meier time points, ",
             "so that the three estimands of a row are truncated at the same ",
             "place; there is no default that means \"as far as the data ",
             "go\"."),
      .ps_values(class(tau)), length(tau)),
      call. = FALSE)
  }
  bad <- is.na(tau) | !is.finite(tau)
  if (any(bad)) {
    stop(sprintf(
      paste0("compare_estimands(): `tau` has %d missing or infinite element%s ",
             "(position%s %s)."),
      sum(bad), if (sum(bad) > 1L) "s" else "",
      if (sum(bad) > 1L) "s" else "", paste(which(bad), collapse = ", ")),
      call. = FALSE)
  }
  if (any(tau <= 0)) {
    stop(sprintf(
      paste0("compare_estimands(): `tau` has %d element%s that %s not ",
             "positive (%s). There is no area under a curve and no survival ",
             "probability to compare at or before tau = 0."),
      sum(tau <= 0), if (sum(tau <= 0) > 1L) "s" else "",
      if (sum(tau <= 0) > 1L) "are" else "is", .ps_values(tau[tau <= 0])),
      call. = FALSE)
  }
  d <- unique(tau[duplicated(tau)])
  if (length(d) > 0L) {
    stop(sprintf(
      paste0("compare_estimands(): `tau` lists %s more than once, which would ",
             "put the same row in the table twice. Keep one copy."),
      .ps_values(d)),
      call. = FALSE)
  }
  sort(tau)
}

# Numbers inside messages, printed with enough significant digits that a tau a
# hair past the limit does not read as if it were exactly at it.
#' @noRd
.ce_num <- function(x) format(x, trim = TRUE, digits = 15)

#' @noRd
.ce_fmt_p <- function(p, digits) {
  thr <- 10^(-digits)
  ifelse(is.na(p), NA_character_,
         ifelse(p < thr,
                paste0("<", formatC(thr, format = "f", digits = digits)),
                formatC(p, format = "f", digits = digits)))
}

# Do two vectors of estimates agree? Used only on the identities checked at the
# end, and on the two cross-route checks.
#' @noRd
.ce_agree <- function(a, b, tol = 1e-8) {
  both_na  <- is.na(a) & is.na(b)
  both_num <- !is.na(a) & !is.na(b)
  both_inf <- both_num & is.infinite(a) & is.infinite(b) & (sign(a) == sign(b))
  close    <- both_num & is.finite(a) & is.finite(b) & abs(a - b) < tol
  isTRUE(all(both_na | both_inf | close))
}

# The area under a Kaplan-Meier step function between 0 and tau: S is held at
# S(t_i) from t_i up to but not including t_i+1, starting from S(0) = 1. This
# is the calculation survRM2 performs for the RMST, and it is checked against
# survRM2's own arm-0 RMST on every call before it is used on the curve a
# constant hazard ratio implies.
#' @noRd
.ce_rmst_step <- function(times, surv, tau) {
  keep <- times < tau
  t_lo <- c(0, times[keep])
  s_lo <- c(1, surv[keep])
  t_hi <- c(times[keep], tau)
  sum(s_lo * (t_hi - t_lo))
}
