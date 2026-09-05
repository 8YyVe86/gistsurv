# ---------------------------------------------------------------------------
# tab_baseline(): a "Table 1" -- row variables described by group, with the
# standardized mean difference (SMD) between groups.
#
# Design notes (deliberate, and different from the study script this was
# distilled from):
#   * no column name and no label is hard-coded: the row variables, the
#     grouping variable and the labels are all arguments;
#   * nothing is printed and nothing is written to disk: the function returns
#     one data frame, and the diagnostics (missing values per variable per
#     group, group sizes, how each SMD was computed) travel back as
#     attributes;
#   * no table-rendering package is used. median [IQR], n (%) and the SMD are
#     computed with base R and stats, so the result is a plain data frame that
#     the caller can write out, join to, or restyle with whatever package the
#     journal wants;
#   * a grouping variable with fewer than two non-missing levels, a grouping
#     variable with missing values, an all-missing row variable and a
#     categorical variable with an implausible number of levels are errors,
#     not silently handled cases.
# The file is kept ASCII-only so that it behaves the same under a UTF-8 and
# under a non-UTF-8 locale.
#
# The small formatting helpers .ps_values() and .ps_counts() are shared with
# R/prep_surv.R; they are internal to the package.
# ---------------------------------------------------------------------------

#' Baseline characteristics table with standardized mean differences
#'
#' @description
#' `tab_baseline()` builds the table a cohort paper prints first: one block of
#' rows per variable, one column per group, continuous variables summarised as
#' median [IQR] and categorical variables as n (%), plus the standardized mean
#' difference (SMD) between the groups.
#'
#' The row variables, the grouping variable and the row labels are arguments,
#' so the function is not tied to any particular registry export. The function
#' prints nothing and writes nothing: it returns a plain data frame, and the
#' per-variable missing counts, the group sizes and the SMD bookkeeping are
#' attached to it as attributes.
#'
#' @details
#' # What is reported
#'
#' Each row variable contributes a block of rows:
#'
#' * a `"label"` row, which for a continuous variable carries the statistic
#'   `median [p25, p75]` and for a categorical variable is a heading with
#'   empty cells;
#' * one `"level"` row per category, carrying `n (%)`, for categorical
#'   variables only. Factor levels are kept in their declared order, including
#'   levels with no observations; character and logical variables are
#'   tabulated in sorted order;
#' * a `"missing"` row carrying the number of missing values, present when the
#'   variable has at least one missing value anywhere in the data used. That
#'   row deliberately carries a count and no percentage.
#'
#' A variable with two categories gets both rows. It is never collapsed to a
#' single "n (%) of the second level" row, which is what
#' `gtsummary::tbl_summary()` does by default: a reader who can see only one
#' of the two rows cannot see where the denominator went, and cannot tell a
#' variable with no missing values from one with many.
#'
#' # Rounding, and which quartile
#'
#' Two conventions are worth stating, because they decide the last printed
#' digit:
#'
#' * quartiles come from `stats::quantile(..., type = quantile_type)`, and the
#'   default is `type = 2`, the SAS definition, which is what
#'   `gtsummary::tbl_summary()` uses. R's own default is `type = 7`; the two
#'   differ in the quartiles but never in the median. Pass
#'   `quantile_type = 7` for R's default. Whichever is used is recorded in the
#'   `settings` diagnostic, so a table always carries its own definition;
#' * ties are rounded **away from zero**, so `53.25%` prints as `53.3%`. Base
#'   R's `round()` rounds half to even and would print `53.2%`. The rule used
#'   here is `trunc(abs(x) * 10^d + 0.5 + sqrt(.Machine$double.eps)) / 10^d`,
#'   the same rule as `gtsummary`'s `style_number()`, the epsilon being what
#'   makes a decimal literal stored just below a tie, such as `0.145`, round
#'   up rather than down.
#'
#' With these two settings the cells reproduce
#' `gtsummary::tbl_summary(percent = "column", missing = "ifany")` exactly.
#'
#' # Percentages
#'
#' The denominator of a percentage is **the number of non-missing values of
#' that variable in that column**, not the column total. Missing values
#' therefore never dilute the percentages: they are shown on their own row and
#' counted in the `missing` diagnostic. This is the convention of the study
#' script this function was distilled from, and of
#' `gtsummary::tbl_summary(percent = "column", missing = "ifany")`. A column in
#' which every value of a variable is missing has no estimable percentage, so
#' its cells are `NA`.
#'
#' # Standardized mean difference
#'
#' The SMD is an effect size and, unlike a p-value, does not grow with the
#' sample size, which is why it is the statistic to read in an observational
#' Table 1. For two groups:
#'
#' * continuous: `d = (m1 - m2) / sqrt((s1^2 + s2^2) / 2)`;
#' * categorical: the multi-category generalisation of Yang and Dalton (2012),
#'   `d = sqrt(t(p1 - p2) %*% solve(S) %*% (p1 - p2))`, where `p` holds the
#'   proportions of the first K-1 categories, `S = (S1 + S2) / 2` and
#'   `S_k = diag(p_k) - tcrossprod(p_k)`. At K = 2 this reduces to the usual
#'   `|p1 - p2| / sqrt((p1 (1 - p1) + p2 (1 - p2)) / 2)`, and that identity is
#'   re-checked on every call and reported in the `checks` diagnostic.
#'
#' With more than two groups the reported value is the **largest absolute SMD
#' over all pairs of groups**, and `smd_pair` names the pair; taking the
#' maximum keeps one badly imbalanced pair from being diluted by the others.
#' SMDs are computed on non-missing values only and `smd_n` reports how many
#' rows that was. Differential missingness is deliberately not folded into the
#' SMD: a single number mixing "the distributions differ" with "the
#' missingness differs" cannot be read back apart. The missing rows are there
#' for that.
#'
#' # What is deliberately not reported
#'
#' No p-values. In an observational cohort the groups were never randomised,
#' so the null hypothesis of exactly equal distributions is known to be false
#' before the data are seen, and whether it is rejected is driven by the group
#' sizes. If a journal insists, compute the tests separately and join them
#' onto the returned data frame by `variable`; keeping them out of this
#' function means no test is ever chosen on your behalf.
#'
#' # Input checks
#'
#' The call stops, with the offending columns, values or counts listed, when
#'
#' 1. a required column (`by`, or any element of `vars`) is not in `data`, or
#'    appears in `data` more than once;
#' 2. `by` also appears in `vars`, `vars` is empty, or `vars` has duplicates;
#' 3. `by` has missing values and `by_missing = "error"` (the default). The
#'    other settings keep those rows as their own column or drop them, and the
#'    number of rows involved is reported in the returned diagnostics either
#'    way -- rows are never dropped without a count;
#' 4. `by` has fewer than two non-missing levels: there is nothing to compare
#'    and every SMD would be undefined;
#' 5. a row variable is missing in every row, so no statistic is estimable.
#'    Drop it from `vars` and report its missingness on its own;
#' 6. a variable summarised as categorical would produce more than
#'    `max_levels` rows, the usual sign of a continuous or identifier column
#'    reaching the table by mistake;
#' 7. a column has a class the function does not summarise (dates, lists), or
#'    a group level would collide with one of the fixed column names of the
#'    result.
#'
#' @param data A data frame (or tibble) with one row per subject.
#' @param by Name of the grouping variable, as a single string. One column of
#'   the table is produced per level. Factor levels are kept in their declared
#'   order; levels with no rows are dropped, and how many were dropped is
#'   reported in the diagnostics.
#' @param vars Names of the row variables, as a character vector, in the order
#'   they are to appear. Must not contain `by`.
#' @param labels Labels for the row variables: a named character vector or
#'   named list, e.g. `c(age = "Age at diagnosis, years")`. Variables with no
#'   entry are labelled with their column name. Entries for variables not in
#'   `vars` are ignored, so one project-wide label list can be passed to every
#'   table.
#' @param type Optional named character vector overriding how a variable is
#'   summarised, e.g. `c(dx_year = "continuous", stage = "categorical")`.
#'   Allowed values are `"continuous"` and `"categorical"`. Without an
#'   override, numeric variables are continuous and factor, character and
#'   logical variables are categorical. Every name must be in `vars`.
#' @param by_missing What to do with rows whose `by` value is missing:
#'   `"error"` (the default) to stop, `"level"` to give them their own column
#'   labelled `"(Missing)"`, or `"drop"` to leave them out. The number of such
#'   rows is reported in the diagnostics in all three cases.
#' @param missing_text Row label used for the count of missing values.
#' @param overall Whether to prepend a column describing all groups together.
#' @param overall_label Name of that column.
#' @param max_levels Largest number of categories a variable may have before
#'   the call stops. The check is on the number of rows the variable would
#'   produce: the number of factor levels, or the number of distinct
#'   non-missing values for character, logical and numeric columns.
#' @param digits_continuous,digits_pct,digits_smd Decimal places for the
#'   median and quartiles, for percentages, and for the formatted SMD. Ties
#'   are rounded away from zero, not half to even.
#' @param quantile_type Which of the nine algorithms of [stats::quantile()]
#'   computes the quartiles. Defaults to `2`, the SAS definition used by
#'   `gtsummary::tbl_summary()`; pass `7` for R's own default. The median is
#'   the same under both.
#' @param big_mark Thousands separator for counts, passed to
#'   [base::formatC()]. Defaults to `""`, i.e. none; house styles that write
#'   `12 345` or `12,345` pass `" "` or `","` here rather than having the
#'   convention baked into the package.
#' @param stat_prefix Optional prefix for the group column names, e.g. `"n_"`.
#'   Use it when a group level would otherwise collide with one of the fixed
#'   column names.
#'
#' @return
#' A plain data frame with one row per table row and the columns
#'
#' * `variable`, `var_label`, `var_type` -- the row variable, its label, and
#'   whether it was summarised as continuous or categorical;
#' * `row_type` -- `"label"`, `"level"` or `"missing"`;
#' * `label` -- what to print in the first column of the table;
#' * `stat_label` -- what the cells of that row hold (`"Median [IQR]"`,
#'   `"n (%)"`, `"n"`);
#' * one character column per group, named after the group level, preceded by
#'   the overall column when `overall = TRUE`;
#' * `smd`, `smd_fmt`, `smd_pair`, `smd_n` -- the standardized mean
#'   difference, filled on the `"label"` row of each variable and `NA`
#'   elsewhere.
#'
#' Two attributes are attached:
#'
#' * `smd_table` -- one row per variable, with `smd_max`, `smd_pair`, `smd_n`
#'   and the `method` used;
#' * `tab_baseline` -- a list with `counts` (rows read, rows used, rows with a
#'   missing group), `group_n` (the column denominators), `missing` (missing
#'   counts and percentages by variable and group), `smd`, `checks` (the
#'   invariants that were enforced), `settings` and `call`.
#'
#' Attributes are dropped by most data-frame verbs, so read them off the
#' object returned by `tab_baseline()` before piping it further.
#'
#' @references
#' Yang D, Dalton JE. A unified approach to measuring the effect size between
#' two groups using SAS. SAS Global Forum 2012, paper 335-2012.
#'
#' @examples
#' ## A small hand-built cohort. No real patient records are used anywhere
#' ## in this package.
#' set.seed(20260901)
#' toy <- data.frame(
#'   site = factor(rep(c("Stomach", "Small intestine", "Colorectal"),
#'                     times = c(50, 30, 20)),
#'                 levels = c("Stomach", "Small intestine", "Colorectal")),
#'   age  = round(c(rnorm(50, 66, 12), rnorm(30, 62, 13), rnorm(20, 70, 10))),
#'   sex  = factor(sample(c("Male", "Female"), 100, TRUE, prob = c(.55, .45)),
#'                 levels = c("Male", "Female")),
#'   size = round(c(rlnorm(50, 3.6, .6), rlnorm(30, 4, .6), rlnorm(20, 3.8, .6)))
#' )
#' toy$size[c(3, 17, 45, 88)] <- NA          # a few unrecorded tumour sizes
#'
#' tab <- tab_baseline(
#'   toy, by = "site", vars = c("age", "sex", "size"),
#'   labels = c(age = "Age at diagnosis, years", sex = "Sex",
#'              size = "Tumour size, mm")
#' )
#' tab
#'
#' ## Diagnostics travel with the result instead of being printed.
#' attr(tab, "smd_table")
#' attr(tab, "tab_baseline")$group_n
#' subset(attr(tab, "tab_baseline")$missing, variable == "size")
#' attr(tab, "tab_baseline")$checks
#'
#' ## An identifier column is not a categorical variable, and the number of
#' ## rows it would produce is checked before anything is tabulated.
#' toy$record_id <- sprintf("R%03d", seq_len(nrow(toy)))
#' try(tab_baseline(toy, by = "site", vars = c("age", "record_id")))
#'
#' ## A grouping variable with one non-missing level has nothing to compare.
#' try(tab_baseline(toy[toy$site == "Stomach", ], by = "site", vars = "age"))
#'
#' ## A misspelled column name is named in the error message.
#' try(tab_baseline(toy, by = "site", vars = c("age", "sexe")))
#'
#' @export
tab_baseline <- function(data,
                         by,
                         vars,
                         labels            = NULL,
                         type              = NULL,
                         by_missing        = c("error", "level", "drop"),
                         missing_text      = "Missing",
                         overall           = TRUE,
                         overall_label     = "Overall",
                         max_levels        = 20L,
                         digits_continuous = 1L,
                         digits_pct        = 1L,
                         digits_smd        = 3L,
                         quantile_type     = 2L,
                         big_mark          = "",
                         stat_prefix       = "") {

  cl <- match.call()

  ## -- 0. data -------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("tab_baseline(): `data` must be a data frame, got an object of class ",
         .ps_values(class(data)), ".", call. = FALSE)
  }
  n_in <- nrow(data)
  if (n_in == 0L) {
    stop("tab_baseline(): `data` has 0 rows, there is nothing to describe.",
         call. = FALSE)
  }

  ## -- 1. the argument values themselves -----------------------------------
  if (missing(by) || missing(vars)) {
    stop(paste0("tab_baseline(): both `by` (the grouping variable, one column ",
                "per level) and `vars` (the row variables) must be supplied; ",
                "tab_baseline() has no default set of columns."),
         call. = FALSE)
  }
  .tb_name_arg(by, "by")
  if (!is.character(vars) || length(vars) == 0L) {
    stop(sprintf(
      paste0("tab_baseline(): `vars` must be a character vector of at least ",
             "one column name, got %s of length %d. A table with no rows is ",
             "not a table."),
      .ps_values(class(vars)), length(vars)),
      call. = FALSE)
  }
  bad_v <- is.na(vars) | !nzchar(vars)
  if (any(bad_v)) {
    stop(sprintf(
      "tab_baseline(): `vars` has %d empty or missing element%s (position%s %s).",
      sum(bad_v), if (sum(bad_v) > 1L) "s" else "",
      if (sum(bad_v) > 1L) "s" else "", paste(which(bad_v), collapse = ", ")),
      call. = FALSE)
  }
  dup_v <- unique(vars[duplicated(vars)])
  if (length(dup_v) > 0L) {
    stop(sprintf(
      paste0("tab_baseline(): `vars` lists %s more than once, which would put ",
             "the same block of rows in the table twice. Keep one copy."),
      .ps_values(dup_v)),
      call. = FALSE)
  }
  if (by %in% vars) {
    stop(sprintf(
      paste0('tab_baseline(): `by` = "%s" is also in `vars`. The grouping ',
             "variable is the column headings, so it cannot be a row of the ",
             "table as well: every one of its cells would be 100%% or 0%% by ",
             "construction, and its SMD would only measure the grouping. Drop ",
             "it from `vars`, or group by something else."),
      by),
      call. = FALSE)
  }
  by_missing <- match.arg(by_missing)
  .tb_name_arg(missing_text,  "missing_text")
  .tb_name_arg(overall_label, "overall_label")
  .tb_string_arg(stat_prefix, "stat_prefix")
  .tb_string_arg(big_mark,    "big_mark")
  if (!is.logical(overall) || length(overall) != 1L || is.na(overall)) {
    stop("tab_baseline(): `overall` must be TRUE or FALSE.", call. = FALSE)
  }
  .tb_count_arg(max_levels,        "max_levels",        min = 2L)
  .tb_count_arg(digits_continuous, "digits_continuous", min = 0L)
  .tb_count_arg(digits_pct,        "digits_pct",        min = 0L)
  .tb_count_arg(digits_smd,        "digits_smd",        min = 0L)
  .tb_count_arg(quantile_type,     "quantile_type",     min = 1L)
  if (quantile_type > 9L) {
    stop(sprintf(
      paste0("tab_baseline(): `quantile_type` must be one of the 9 algorithms ",
             "of stats::quantile(), got %d."),
      as.integer(quantile_type)),
      call. = FALSE)
  }

  ## -- 2. required columns present, and present once -----------------------
  needed <- c(by, vars)
  role   <- c("by", rep("vars", length(vars)))
  hit    <- needed %in% names(data)
  if (!all(hit)) {
    stop(sprintf(
      paste0("tab_baseline(): %d required column%s not found in `data`:\n  %s\n",
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
      paste0("tab_baseline(): %s appear%s more than once among the columns of ",
             "`data`; tab_baseline() would silently use the first one. Make ",
             "the column names unique first."),
      .ps_values(dup_c), if (length(dup_c) > 1L) "" else "s"),
      call. = FALSE)
  }

  ## -- 3. labels -----------------------------------------------------------
  if (is.null(labels)) labels <- list()
  if (!(is.character(labels) || is.list(labels))) {
    stop(sprintf(
      paste0("tab_baseline(): `labels` must be a named character vector or a ",
             "named list, got %s."),
      .ps_values(class(labels))),
      call. = FALSE)
  }
  if (length(labels) > 0L &&
      (is.null(names(labels)) || any(!nzchar(names(labels))))) {
    stop(paste0("tab_baseline(): every element of `labels` must be named with ",
                "the column it labels, e.g. ",
                "labels = c(age = \"Age at diagnosis, years\")."),
         call. = FALSE)
  }
  var_label <- vapply(vars, function(v) {
    if (!(v %in% names(labels))) return(v)
    lab <- labels[[v]]
    if (!is.character(lab) || length(lab) != 1L || is.na(lab)) {
      stop(sprintf(
        paste0('tab_baseline(): the label for "%s" must be a single ',
               "non-missing string, got %s of length %d."),
        v, .ps_values(class(lab)), length(lab)),
        call. = FALSE)
    }
    lab
  }, character(1), USE.NAMES = FALSE)

  ## -- 4. how each variable is summarised ----------------------------------
  var_type <- vapply(vars, function(v) {
    x <- data[[v]]
    if (is.factor(x) || is.character(x) || is.logical(x)) return("categorical")
    if (is.numeric(x)) return("continuous")
    NA_character_
  }, character(1), USE.NAMES = FALSE)

  bad_cls <- vars[is.na(var_type)]
  if (length(bad_cls) > 0L) {
    stop(sprintf(
      paste0("tab_baseline(): %d row variable%s of a class tab_baseline() ",
             "does not summarise:\n  %s\n",
             "  Only numeric (continuous) and factor, character or logical ",
             "(categorical) columns are handled. Convert %s first, or drop ",
             "%s from `vars`."),
      length(bad_cls), if (length(bad_cls) > 1L) "s are" else " is",
      paste0(sprintf('vars = "%s" (class %s)', bad_cls,
                     vapply(bad_cls,
                            function(v) .ps_values(class(data[[v]])),
                            character(1))),
             collapse = "\n  "),
      if (length(bad_cls) > 1L) "them" else "it",
      if (length(bad_cls) > 1L) "them" else "it"),
      call. = FALSE)
  }

  if (!is.null(type)) {
    if (!is.character(type) || length(type) == 0L || is.null(names(type)) ||
        any(!nzchar(names(type)))) {
      stop(paste0("tab_baseline(): `type` must be a named character vector, ",
                  "e.g. type = c(dx_year = \"continuous\")."),
           call. = FALSE)
    }
    unknown <- setdiff(names(type), vars)
    if (length(unknown) > 0L) {
      stop(sprintf(
        paste0("tab_baseline(): `type` names %s, which %s not in `vars`, so ",
               "the override would be silently ignored. Remove the entr%s, or ",
               "add the variable to `vars`."),
        .ps_values(unknown), if (length(unknown) > 1L) "are" else "is",
        if (length(unknown) > 1L) "ies" else "y"),
        call. = FALSE)
    }
    bad_t <- type[!(type %in% c("continuous", "categorical"))]
    if (length(bad_t) > 0L) {
      stop(sprintf(
        paste0("tab_baseline(): `type` must be \"continuous\" or ",
               "\"categorical\", got %s for %s."),
        .ps_values(unname(bad_t)), .ps_values(names(bad_t))),
        call. = FALSE)
    }
    num_ok <- vapply(names(type), function(v) is.numeric(data[[v]]),
                     logical(1))
    bad_ct <- names(type)[type == "continuous" & !num_ok]
    if (length(bad_ct) > 0L) {
      stop(sprintf(
        paste0("tab_baseline(): `type` asks for %s to be summarised as ",
               "continuous, but %s not numeric. A median of factor levels is ",
               "not a number; convert the column first."),
        .ps_values(bad_ct), if (length(bad_ct) > 1L) "they are" else "it is"),
        call. = FALSE)
    }
    var_type[match(names(type), vars)] <- unname(type)
  }
  names(var_type)  <- vars
  names(var_label) <- vars

  ## -- 5. the grouping variable --------------------------------------------
  g_raw <- data[[by]]
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
      paste0('tab_baseline(): the `by` column "%s" is of class %s, which ',
             "cannot be used as a grouping variable. Convert it to a factor ",
             "first."),
      by, .ps_values(class(g_raw))),
      call. = FALSE)
  }

  n_by_na   <- sum(is.na(g))
  n_lv_all  <- nlevels(g)
  g         <- droplevels(g)
  n_lv_drop <- n_lv_all - nlevels(g)

  if (nlevels(g) > max_levels) {
    stop(sprintf(
      paste0('tab_baseline(): the `by` column "%s" has %d distinct non-missing ',
             "values, more than max_levels = %d, so the table would have %d ",
             "columns.\n",
             "  first values : %s\n",
             "  A continuous variable cannot be a grouping variable: bin it ",
             "first. Raise max_levels only if that many columns really are ",
             "intended."),
      by, nlevels(g), max_levels, nlevels(g) + as.integer(overall),
      .ps_values(levels(g))),
      call. = FALSE)
  }

  n_by_drop <- 0L
  if (n_by_na > 0L) {
    if (by_missing == "error") {
      stop(sprintf(
        paste0('tab_baseline(): the `by` column "%s" has %d missing value%s ',
               "out of %d rows.\n",
               "  tab_baseline() will not decide on its own what those rows ",
               "are: they are neither a group nor nothing.\n",
               "  Pass by_missing = \"level\" to give them their own column ",
               "labelled \"(%s)\", or by_missing = \"drop\" to leave them out ",
               "(either way the count comes back in the diagnostics), or ",
               "filter them out before calling."),
        by, n_by_na, if (n_by_na > 1L) "s" else "", n_in, missing_text),
        call. = FALSE)
    } else if (by_missing == "level") {
      na_lab <- paste0("(", missing_text, ")")
      if (na_lab %in% levels(g)) {
        stop(sprintf(
          paste0('tab_baseline(): by_missing = "level" would add a group ',
                 'labelled "%s", but the `by` column "%s" already has a level ',
                 "with that name. Rename the level, or change `missing_text`."),
          na_lab, by),
          call. = FALSE)
      }
      g <- factor(ifelse(is.na(g), na_lab, as.character(g)),
                  levels = c(levels(g), na_lab))
    } else {
      n_by_drop <- n_by_na
    }
  }

  keep   <- !is.na(g)
  dat    <- data[keep, , drop = FALSE]
  g      <- droplevels(g[keep])
  n_used <- nrow(dat)
  if (n_used == 0L) {
    stop(sprintf(
      paste0('tab_baseline(): every row has a missing `by` value ("%s"), so ',
             "no group is left to describe."),
      by),
      call. = FALSE)
  }

  lv <- levels(g)
  if (length(lv) < 2L) {
    stop(sprintf(
      paste0('tab_baseline(): the `by` column "%s" has %d non-missing level%s ',
             "(%s) in the %d rows used; at least 2 are needed.\n",
             "  With one group there is nothing to compare: every ",
             "standardized mean difference is undefined and the table would ",
             "be a single column.\n",
             "  Check that the data were not already filtered down to one ",
             "group, and that an empty group is not hiding behind an unused ",
             "factor level (%d unused level%s dropped here)."),
      by, length(lv), if (length(lv) == 1L) "" else "s",
      if (length(lv) == 0L) "none" else .ps_values(lv), n_used,
      n_lv_drop, if (n_lv_drop == 1L) " was" else "s were"),
      call. = FALSE)
  }

  ## -- 6. row variables: all-missing, and too many levels ------------------
  n_miss_all <- vapply(vars, function(v) sum(is.na(dat[[v]])), integer(1))
  all_na     <- vars[n_miss_all == n_used]
  if (length(all_na) > 0L) {
    stop(sprintf(
      paste0("tab_baseline(): %d row variable%s missing in all %d rows used:\n",
             "  %s\n",
             "  No median, no percentage and no SMD is estimable from an ",
             "empty column, so tab_baseline() will not print a block of blanks ",
             "for it. Drop %s from `vars` and report the missingness ",
             "separately."),
      length(all_na), if (length(all_na) > 1L) "s are" else " is", n_used,
      .ps_values(all_na), if (length(all_na) > 1L) "them" else "it"),
      call. = FALSE)
  }

  var_levels <- lapply(vars, function(v) {
    if (var_type[[v]] != "categorical") return(character(0))
    x <- dat[[v]]
    if (is.factor(x)) return(levels(x))
    u <- unique(x[!is.na(x)])
    as.character(sort(u))          # numeric codes sort as numbers, not as text
  })
  names(var_levels) <- vars

  n_lv     <- vapply(var_levels, length, integer(1))
  too_many <- vars[var_type == "categorical" & n_lv > max_levels]
  if (length(too_many) > 0L) {
    stop(sprintf(
      paste0("tab_baseline(): %d variable%s summarised as categorical would ",
             "produce more than max_levels = %d rows:\n  %s\n",
             "  That is usually a continuous measurement or an identifier ",
             "column reaching the table by mistake. Either pass ",
             "type = c(%s = \"continuous\") if it is a measurement, or ",
             "collapse the categories first, or raise max_levels if that many ",
             "rows really are intended."),
      length(too_many), if (length(too_many) > 1L) "s" else "", max_levels,
      paste0(sprintf('vars = "%s" (%d levels, class %s)', too_many,
                     n_lv[too_many],
                     vapply(too_many,
                            function(v) .ps_values(class(dat[[v]])),
                            character(1))),
             collapse = "\n  "),
      too_many[1]),
      call. = FALSE)
  }

  ## -- 7. the columns of the table -----------------------------------------
  stat_lv   <- c(if (overall) overall_label, lv)
  stat_cols <- paste0(stat_prefix, stat_lv)
  grp_idx   <- c(if (overall) list(seq_len(n_used)),
                 lapply(lv, function(l) which(g == l)))
  names(grp_idx) <- stat_cols

  fixed <- c("variable", "var_label", "var_type", "row_type", "label",
             "stat_label", "smd", "smd_fmt", "smd_pair", "smd_n")
  clash <- intersect(stat_cols, fixed)
  if (length(clash) > 0L) {
    stop(sprintf(
      paste0("tab_baseline(): the group column%s %s would collide with the ",
             "fixed columns of the returned data frame (%s).\n",
             "  Pass stat_prefix = \"n_\" (or any other prefix), or rename the ",
             "level%s."),
      if (length(clash) > 1L) "s" else "", .ps_values(clash),
      paste(fixed, collapse = ", "),
      if (length(clash) > 1L) "s" else ""),
      call. = FALSE)
  }
  if (anyDuplicated(stat_cols)) {
    stop(sprintf(
      "tab_baseline(): the group column names are not unique: %s.",
      .ps_values(stat_cols[duplicated(stat_cols)])),
      call. = FALSE)
  }

  group_n <- vapply(grp_idx, length, integer(1))
  if (sum(group_n[.tb_group_cols(stat_cols, overall)]) != n_used) {
    stop(paste0("tab_baseline(): the group sizes do not add up to the number ",
                "of rows used, which should be impossible by construction. Do ",
                "not use this result; report it as a bug in tab_baseline()."),
         call. = FALSE)
  }

  ## -- 8. the body of the table --------------------------------------------
  blocks <- lapply(vars, function(v) {
    x     <- dat[[v]]
    ty    <- var_type[[v]]
    has_m <- anyNA(x)

    if (ty == "continuous") {
      rows <- data.frame(
        row_type   = "label",
        label      = var_label[[v]],
        stat_label = "Median [IQR]",
        stringsAsFactors = FALSE)
      cells <- matrix(vapply(grp_idx, function(i)
        .tb_stat_cont(x[i], digits_continuous, big_mark, quantile_type),
        character(1)), nrow = 1L)
    } else {
      levs <- var_levels[[v]]
      rows <- data.frame(
        row_type   = c("label", rep("level", length(levs))),
        label      = c(var_label[[v]], levs),
        stat_label = c(NA_character_, rep("n (%)", length(levs))),
        stringsAsFactors = FALSE)
      cells <- rbind(
        rep(NA_character_, length(grp_idx)),
        vapply(grp_idx, function(i)
          .tb_stat_cat(x[i], levs, digits_pct, big_mark),
          character(length(levs))))
    }

    if (has_m) {
      rows  <- rbind(rows, data.frame(
        row_type = "missing", label = missing_text, stat_label = "n",
        stringsAsFactors = FALSE))
      cells <- rbind(cells, vapply(grp_idx, function(i)
        .tb_int(sum(is.na(x[i])), big_mark), character(1)))
    }

    colnames(cells) <- stat_cols
    cbind(data.frame(variable  = v,
                     var_label = var_label[[v]],
                     var_type  = ty,
                     rows,
                     stringsAsFactors = FALSE),
          as.data.frame(cells, stringsAsFactors = FALSE, check.names = FALSE))
  })
  tab <- do.call(rbind, blocks)
  rownames(tab) <- NULL

  ## -- 9. standardized mean differences ------------------------------------
  chk_smd <- .tb_smd_selfcheck()
  if (!isTRUE(chk_smd)) {
    stop(paste0("tab_baseline(): the categorical SMD does not reproduce the ",
                "closed-form two-group formula, which should be impossible by ",
                "construction. Do not use this result; report it as a bug in ",
                "tab_baseline()."),
         call. = FALSE)
  }

  smd <- do.call(rbind, lapply(vars, function(v) {
    is_cat <- var_type[[v]] == "categorical"
    s <- .tb_smd_summary(dat[[v]], g, is_cat, var_levels[[v]])
    data.frame(variable  = v,
               var_label = var_label[[v]],
               smd_max   = s$smd,
               smd_pair  = s$pair,
               smd_n     = s$n,
               method    = if (is_cat)
                 "categorical: Yang-Dalton (2012) multi-category SMD" else
                   "continuous: (m1 - m2) / sqrt((s1^2 + s2^2) / 2)",
               stringsAsFactors = FALSE)
  }))
  rownames(smd) <- NULL

  i_lab <- tab$row_type == "label"
  m     <- match(tab$variable[i_lab], smd$variable)
  tab$smd      <- NA_real_
  tab$smd_fmt  <- NA_character_
  tab$smd_pair <- NA_character_
  tab$smd_n    <- NA_integer_
  tab$smd[i_lab]      <- smd$smd_max[m]
  tab$smd_fmt[i_lab]  <- .tb_num(smd$smd_max[m], digits_smd, "")
  tab$smd_pair[i_lab] <- smd$smd_pair[m]
  tab$smd_n[i_lab]    <- smd$smd_n[m]

  ## -- 10. diagnostics ------------------------------------------------------
  miss <- do.call(rbind, lapply(vars, function(v) {
    x  <- dat[[v]]
    nm <- vapply(grp_idx, function(i) sum(is.na(x[i])), integer(1))
    nn <- vapply(grp_idx, length, integer(1))
    data.frame(variable  = v,
               var_label = var_label[[v]],
               group     = names(grp_idx),
               n         = unname(nn),
               n_miss    = unname(nm),
               pct_miss  = round(100 * unname(nm) / unname(nn), 2),
               stringsAsFactors = FALSE)
  }))
  rownames(miss) <- NULL

  ok_denom <- all(vapply(vars, function(v) {
    x <- dat[[v]]
    all(vapply(grp_idx,
               function(i) sum(!is.na(x[i])) + sum(is.na(x[i])) == length(i),
               logical(1)))
  }, logical(1)))
  if (!ok_denom) {
    stop(paste0("tab_baseline(): the non-missing and missing counts do not add ",
                "up to the group size, which should be impossible by ",
                "construction. Do not use this result; report it as a bug in ",
                "tab_baseline()."),
         call. = FALSE)
  }

  checks <- c(
    "grouping variable has at least two non-missing levels"     = TRUE,
    "group sizes add up to the rows used"                       = TRUE,
    "no row variable is missing in every row"                   = TRUE,
    "percent denominators are the non-missing count per column" = ok_denom,
    "categorical SMD reproduces the two-group closed form"      = chk_smd)

  attr(tab, "smd_table") <- smd
  attr(tab, "tab_baseline") <- list(
    counts   = c(rows_in            = n_in,
                 rows_used          = n_used,
                 by_missing         = n_by_na,
                 by_missing_dropped = n_by_drop),
    group_n  = group_n,
    missing  = miss,
    smd      = smd,
    checks   = checks,
    settings = list(
      by                   = by,
      vars                 = vars,
      var_label            = var_label,
      var_type             = var_type,
      by_missing           = by_missing,
      by_levels_dropped    = n_lv_drop,
      overall              = overall,
      max_levels           = max_levels,
      percent_denominator  =
        "non-missing values of the variable within the column",
      continuous_statistic = sprintf(
        "median [p25, p75], stats::quantile() type %d, ties rounded away from zero",
        as.integer(quantile_type)),
      smd_rule             = if (length(lv) > 2L)
        "largest absolute SMD over all pairs of groups" else
          "two groups, a single pairwise SMD",
      digits               = c(continuous = digits_continuous,
                               pct        = digits_pct,
                               smd        = digits_smd),
      big_mark             = big_mark),
    call     = cl)
  tab
}


# -- internal helpers, not exported ------------------------------------------
# .ps_values() and .ps_counts() are defined in R/prep_surv.R and reused here.

#' @noRd
.tb_name_arg <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf(
      "tab_baseline(): `%s` must be a single non-missing, non-empty string.",
      arg), call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.tb_string_arg <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("tab_baseline(): `%s` must be a single non-missing string.",
                 arg), call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.tb_count_arg <- function(x, arg, min = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != as.integer(x) ||
      x < min) {
    stop(sprintf(
      "tab_baseline(): `%s` must be a single whole number >= %d, got %s.",
      arg, min, .ps_values(x)), call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
# Round half away from zero (53.25 -> 53.3 at one decimal), with a
# sqrt(.Machine$double.eps) guard so that a decimal literal stored just below
# the tie, such as 0.145, still rounds up. This is the rule clinical tables
# use and the rule behind gtsummary's style_number(); base R's round() would
# give 53.2, because it rounds half to even.
.tb_round5 <- function(x, digits) {
  trunc(abs(x) * 10^digits + 0.5 + sqrt(.Machine$double.eps)) /
    10^digits * sign(x)
}

#' @noRd
.tb_num <- function(x, digits, big_mark) {
  out <- formatC(.tb_round5(x, digits), format = "f", digits = digits,
                 big.mark = big_mark)
  out[is.na(x)] <- NA_character_
  out
}

#' @noRd
.tb_int <- function(x, big_mark) {
  out <- formatC(as.integer(x), format = "d", big.mark = big_mark)
  out[is.na(x)] <- NA_character_
  out
}

#' @noRd
.tb_stat_cont <- function(x, digits, big_mark, q_type) {
  x <- as.numeric(x[!is.na(x)])
  if (length(x) == 0L) return(NA_character_)
  q <- stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE, type = q_type)
  sprintf("%s [%s, %s]",
          .tb_num(q[2], digits, big_mark),
          .tb_num(q[1], digits, big_mark),
          .tb_num(q[3], digits, big_mark))
}

#' @noRd
# n (%) per level. The denominator is the non-missing count in this column.
.tb_stat_cat <- function(x, levs, digits, big_mark) {
  x   <- as.character(x[!is.na(x)])
  den <- length(x)
  if (den == 0L) return(rep(NA_character_, length(levs)))
  cnt <- as.integer(table(factor(x, levels = levs)))
  sprintf("%s (%s%%)",
          .tb_int(cnt, big_mark),
          .tb_num(100 * cnt / den, digits, big_mark))
}

#' @noRd
# Moore-Penrose inverse through the SVD, so that a singular covariance matrix
# does not oblige the package to depend on MASS::ginv().
.tb_ginv <- function(S, tol = sqrt(.Machine$double.eps)) {
  s   <- svd(S)
  pos <- s$d > max(tol * s$d[1L], 0)
  if (!any(pos)) return(matrix(0, nrow = ncol(S), ncol = nrow(S)))
  s$v[, pos, drop = FALSE] %*% ((1 / s$d[pos]) * t(s$u[, pos, drop = FALSE]))
}

#' @noRd
# Yang and Dalton (2012), from two vectors of category proportions.
.tb_smd_cat_props <- function(p1, p2) {
  keep <- (p1 + p2) > 0         # a category empty in both groups makes the
  p1   <- p1[keep]              # covariance matrix singular; drop it first
  p2   <- p2[keep]
  K    <- length(p1)
  if (K <= 1L) return(0)
  idx <- seq_len(K - 1L)        # the last proportion is 1 minus the others
  dv  <- p1[idx] - p2[idx]
  S1  <- diag(p1[idx], nrow = K - 1L) - tcrossprod(p1[idx])
  S2  <- diag(p2[idx], nrow = K - 1L) - tcrossprod(p2[idx])
  S   <- (S1 + S2) / 2
  Sinv <- tryCatch(solve(S), error = function(e) .tb_ginv(S))
  v <- as.numeric(t(dv) %*% Sinv %*% dv)
  if (!is.finite(v) || v < 0) return(NA_real_)
  sqrt(v)
}

#' @noRd
.tb_smd_two_cat <- function(x, g, a, b, levs) {
  x1 <- as.character(x[g == a]); x1 <- x1[!is.na(x1)]
  x2 <- as.character(x[g == b]); x2 <- x2[!is.na(x2)]
  if (length(x1) < 1L || length(x2) < 1L) return(NA_real_)
  p1 <- as.numeric(table(factor(x1, levels = levs))) / length(x1)
  p2 <- as.numeric(table(factor(x2, levels = levs))) / length(x2)
  .tb_smd_cat_props(p1, p2)
}

#' @noRd
.tb_smd_two_cont <- function(x, g, a, b) {
  x1 <- as.numeric(x[g == a]); x1 <- x1[!is.na(x1)]
  x2 <- as.numeric(x[g == b]); x2 <- x2[!is.na(x2)]
  if (length(x1) < 2L || length(x2) < 2L) return(NA_real_)
  den <- sqrt((stats::var(x1) + stats::var(x2)) / 2)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  (mean(x1) - mean(x2)) / den
}

#' @noRd
# The largest absolute SMD over all pairs of groups, and which pair it was.
.tb_smd_summary <- function(x, g, is_cat, levs) {
  lv  <- levels(g)
  K   <- length(lv)
  prs <- vector("list", K * (K - 1L) / 2L)
  k   <- 0L
  for (i in seq_len(K - 1L)) {
    for (j in seq.int(i + 1L, K)) {
      k <- k + 1L
      prs[[k]] <- c(lv[i], lv[j])
    }
  }
  vals <- vapply(prs, function(p) {
    if (is_cat) .tb_smd_two_cat(x, g, p[1], p[2], levs) else
      .tb_smd_two_cont(x, g, p[1], p[2])
  }, numeric(1))

  n_ok <- sum(!is.na(x))
  if (all(is.na(vals))) {
    return(list(smd = NA_real_, pair = NA_character_, n = n_ok))
  }
  w <- which.max(abs(vals))
  list(smd = abs(vals[w]), pair = paste(prs[[w]], collapse = " vs "), n = n_ok)
}

#' @noRd
# At two categories the Yang-Dalton form must equal the textbook binary
# formula. Checked on every call: it is a handful of arithmetic operations,
# and it is the one place where a silent algebra error would leave no trace in
# the printed table.
.tb_smd_selfcheck <- function() {
  p_a    <- 0.6
  p_b    <- 0.4
  closed <- abs(p_a - p_b) / sqrt((p_a * (1 - p_a) + p_b * (1 - p_b)) / 2)
  got    <- .tb_smd_cat_props(c(1 - p_a, p_a), c(1 - p_b, p_b))
  isTRUE(abs(got - closed) < 1e-10)
}

#' @noRd
# The group columns without the overall column; used only by the internal
# "group sizes add up" invariant.
.tb_group_cols <- function(stat_cols, overall) {
  if (overall) stat_cols[-1L] else stat_cols
}
