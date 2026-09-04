# ---------------------------------------------------------------------------
# prep_surv(): build time_mo / event_os / event_cr from one row per subject.
#
# Design notes (deliberate, and different from the study script this was
# distilled from):
#   * no column name is hard-coded: every variable is passed as an argument;
#   * nothing is printed: diagnostics travel back as attributes;
#   * an undeclared value is an error, never a silent recode into
#     "censored" or "death from other causes".
# The file is kept ASCII-only so that it behaves the same under a UTF-8 and
# under a non-UTF-8 locale.
# ---------------------------------------------------------------------------

#' Prepare survival and competing-risk outcome variables
#'
#' @description
#' `prep_surv()` takes a data frame with one row per subject and returns it
#' with the columns a survival analysis needs:
#'
#' * `time_mo` -- follow-up time, numeric, non-missing, non-negative;
#' * `event_os` -- overall-survival indicator, `1` = death, `0` = censored;
#' * `event_cr` -- competing-risk outcome, `0` = censored,
#'   `1` = death from the disease of interest, `2` = death from other causes
#'   (only produced when `cause` is supplied).
#'
#' Column names and coding values are arguments, so the function is not tied
#' to any particular registry export. Every recoding rule is explicit: a
#' value the function has not been told about stops the call instead of being
#' folded into "censored" or "other cause". The function prints nothing; the
#' `event_os` by `event_cr` cross-tabulation and the remaining diagnostics
#' are attached to the returned data frame as attributes.
#'
#' @details
#' # Recoding rules
#'
#' `event_os` is `0` for rows whose `status` value is in `alive_value` and
#' `1` for rows whose `status` value is in `dead_value`. Any other value,
#' including `NA`, is an error.
#'
#' `event_cr` is derived from `status` and `cause` together, so that the two
#' outcome variables cannot contradict each other by construction:
#'
#' * alive -> `0`;
#' * dead and `cause` in `cause_gist_value` -> `1`;
#' * dead and `cause` in `cause_unknown_value` -> `cause_unknown_to`;
#' * dead and `cause` in `cause_other_value` -> `2`;
#' * dead and none of the above -> `2` when `cause_other_value` is `NULL`
#'   (open coding: any undeclared cause is a death from other causes), and an
#'   error when `cause_other_value` is supplied (closed coding: the set of
#'   admissible causes has been declared in full).
#'
#' The four value sets must not overlap. `cause_unknown_value` defaults to
#' `NA`, i.e. a dead subject with a missing cause of death is treated as an
#' unknown cause and sent to `cause_unknown_to`; the number of rows recoded
#' this way is reported in the returned attributes, because it is the figure
#' a limitations paragraph has to quote.
#'
#' # Input checks
#'
#' The call stops, with the offending values listed, when
#'
#' 1. a required column (`time`, `status`, `cause`) is not in `data`;
#' 2. `time` is missing, negative, infinite, or not readable as a number;
#' 3. `status` or `cause` holds a value outside the declared sets;
#' 4. `status` and `cause` contradict each other -- a subject recorded as
#'    alive that carries a declared cause of death, or a subject recorded as
#'    dead that carries a `cause_censor_value`. `cause_unknown_to = 0` is
#'    refused for the same reason: it would censor a dead subject in
#'    `event_cr` while `event_os` records a death.
#'
#' The bidirectional `event_os` / `event_cr` agreement (`event_os == 0` only
#' with `event_cr == 0`; `event_os == 1` only with `event_cr` in `1:2`; and
#' the count of deaths equal to the sum of the two competing-risk events) is
#' verified again on the derived columns before the result is returned.
#'
#' @param data A data frame (or tibble) with one row per subject. The class
#'   of `data` is preserved; the new columns are appended on the right.
#' @param time Name of the follow-up-time column, as a single string. Its
#'   values may be numeric or character digits (a zero-padded registry export
#'   such as `"0012"` is accepted); factors are refused, because factor
#'   levels are not numbers.
#' @param status Name of the vital-status column, as a single string.
#' @param alive_value,dead_value Values of `status` that mean alive
#'   (censored) and dead. Each may be a vector if several codes are used.
#'   The two sets must not overlap, and together they must cover every value
#'   present in the column.
#' @param cause Name of the cause-of-death column, as a single string, or
#'   `NULL` (the default) for overall survival only. When `NULL`, no
#'   `event_cr` column is produced.
#' @param cause_gist_value Values of `cause` that mean death from the disease
#'   of interest (mapped to `event_cr = 1`). Required whenever `cause` is
#'   supplied.
#' @param cause_other_value Values of `cause` that mean death from another
#'   cause (mapped to `event_cr = 2`). Leave `NULL` to accept any undeclared
#'   cause as a death from other causes; supply it to force every admissible
#'   value to be declared.
#' @param cause_unknown_value Values of `cause` that mean the cause of death
#'   is unknown. Defaults to `NA`, i.e. a missing cause of death.
#' @param cause_censor_value Values of `cause` that mean "no event", for a
#'   cause column that already carries its own censoring code. A subject
#'   recorded as dead that carries one of these values is an error.
#' @param cause_unknown_to Category an unknown cause of death is assigned to:
#'   `1L` (death from the disease of interest) or `2L` (death from other
#'   causes). `0L` is refused.
#'
#' @return
#' `data` with `time_mo` and `event_os` appended, plus `event_cr` when
#' `cause` was supplied, and with two attributes:
#'
#' * `os_cr_table` -- the `table(event_os, event_cr)` cross-check, or `NULL`
#'   when no `cause` was supplied;
#' * `prep_surv` -- a list with `counts` (subjects, censored, deaths, and the
#'   `event_cr` categories including how many rows an unknown cause was
#'   recoded into), `follow_up` (minimum, median, maximum and the number of
#'   zero follow-up times), `checks` (the invariants that were enforced),
#'   `settings` (the column names and value sets used) and `call`.
#'
#' Attributes are dropped by most data-frame verbs, so read them off the
#' object returned by `prep_surv()` before piping it further.
#'
#' @examples
#' ## A small hand-built cohort. No real patient records are used anywhere
#' ## in this package.
#' toy <- data.frame(
#'   id    = 1:8,
#'   fu    = c(0, 5, 12, 18, 24, 36, 48, 60),
#'   vital = c("Dead", "Dead", "Alive", "Dead", "Alive", "Dead", "Alive", "Dead"),
#'   cod   = c("GIST", "Heart disease", NA, "Unknown", NA, "GIST", NA, "Stroke"),
#'   stringsAsFactors = FALSE
#' )
#'
#' ## Overall survival only.
#' os <- prep_surv(toy, time = "fu", status = "vital")
#' os[, c("id", "time_mo", "event_os")]
#'
#' ## Competing risks: "GIST" is the event of interest, "Unknown" and a
#' ## missing cause are sent to "death from other causes".
#' cr <- prep_surv(toy, time = "fu", status = "vital",
#'                 cause = "cod", cause_gist_value = "GIST",
#'                 cause_unknown_value = c("Unknown", NA),
#'                 cause_unknown_to = 2L)
#' cr[, c("id", "time_mo", "event_os", "event_cr")]
#'
#' ## The cross-check and the counts travel with the result.
#' attr(cr, "os_cr_table")
#' attr(cr, "prep_surv")$counts
#'
#' ## Closed coding: every admissible cause is declared, so a value that was
#' ## not declared stops the call instead of becoming "other cause".
#' try(prep_surv(toy, time = "fu", status = "vital",
#'               cause = "cod", cause_gist_value = "GIST",
#'               cause_other_value = "Heart disease",
#'               cause_unknown_value = c("Unknown", NA)))
#'
#' ## A misspelled column name is named in the error message.
#' try(prep_surv(toy, time = "fu", status = "vital_status"))
#'
#' @export
prep_surv <- function(data,
                      time                = "time_mo",
                      status              = "vital_status",
                      alive_value         = "Alive",
                      dead_value          = "Dead",
                      cause               = NULL,
                      cause_gist_value    = NULL,
                      cause_other_value   = NULL,
                      cause_unknown_value = NA,
                      cause_censor_value  = NULL,
                      cause_unknown_to    = 2L) {

  cl <- match.call()

  ## -- 0. data ------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("prep_surv(): `data` must be a data frame, got an object of class ",
         .ps_values(class(data)), ".", call. = FALSE)
  }
  n <- nrow(data)
  if (n == 0L) {
    stop("prep_surv(): `data` has 0 rows, there is nothing to prepare.",
         call. = FALSE)
  }

  ## -- 1. the column-name arguments themselves ----------------------------
  .ps_name_arg(time, "time")
  .ps_name_arg(status, "status")
  if (!is.null(cause)) .ps_name_arg(cause, "cause")

  needed <- c(time = time, status = status)
  if (!is.null(cause)) needed <- c(needed, cause = cause)

  ## -- 2. required columns present ----------------------------------------
  miss <- needed[!(needed %in% names(data))]
  if (length(miss) > 0L) {
    stop(sprintf(
      paste0("prep_surv(): %d required column%s not found in `data`:\n  %s\n",
             "  `data` has %d columns. Check the spelling, or pass the ",
             "column names explicitly."),
      length(miss), if (length(miss) > 1L) "s" else "",
      paste0(sprintf('%s = "%s"', names(miss), unname(miss)), collapse = "\n  "),
      ncol(data)),
      call. = FALSE)
  }
  dup <- needed[needed %in% names(data)[duplicated(names(data))]]
  if (length(dup) > 0L) {
    stop(sprintf(
      paste0("prep_surv(): %s appear%s more than once among the columns of ",
             "`data`; prep_surv() would silently use the first one. Make the ",
             "column names unique first."),
      paste0(sprintf('%s = "%s"', names(dup), unname(dup)), collapse = ", "),
      if (length(dup) > 1L) "" else "s"),
      call. = FALSE)
  }

  ## -- 3. no silent overwriting of existing columns ------------------------
  new_cols <- c("time_mo", "event_os", if (!is.null(cause)) "event_cr")
  clash    <- setdiff(intersect(new_cols, names(data)), unname(needed))
  if (length(clash) > 0L) {
    stop(sprintf(
      paste0("prep_surv(): `data` already has the column%s %s, which ",
             "prep_surv() would overwrite. Rename or drop %s first."),
      if (length(clash) > 1L) "s" else "", .ps_values(clash),
      if (length(clash) > 1L) "them" else "it"),
      call. = FALSE)
  }

  ## -- 4. follow-up time ---------------------------------------------------
  tv <- data[[time]]
  if (is.factor(tv)) {
    stop(sprintf(
      paste0('prep_surv(): the `time` column "%s" is a factor. Factor levels ',
             "are not numbers, so prep_surv() will not convert it for you; ",
             "convert it explicitly first."), time),
      call. = FALSE)
  }
  if (is.character(tv)) {
    tnum <- suppressWarnings(as.numeric(tv))
    bad  <- is.na(tnum) & !is.na(tv)
    if (any(bad)) {
      stop(sprintf(
        paste0('prep_surv(): the `time` column "%s" is character and %d of %d ',
               "value%s cannot be read as a number: %s"),
        time, sum(bad), n, if (sum(bad) > 1L) "s" else "",
        .ps_counts(tv[bad])),
        call. = FALSE)
    }
  } else if (is.numeric(tv)) {
    tnum <- as.numeric(tv)
  } else {
    stop(sprintf(
      paste0('prep_surv(): the `time` column "%s" must be numeric or ',
             "character digits, got an object of class %s."),
      time, .ps_values(class(tv))),
      call. = FALSE)
  }

  n_na <- sum(is.na(tnum))
  if (n_na > 0L) {
    stop(sprintf(
      paste0('prep_surv(): the `time` column "%s" has %d missing value%s out ',
             "of %d rows. prep_surv() neither drops nor imputes them; decide ",
             "what to do with these rows before calling it."),
      time, n_na, if (n_na > 1L) "s" else "", n),
      call. = FALSE)
  }
  n_inf <- sum(!is.finite(tnum))
  if (n_inf > 0L) {
    stop(sprintf(
      'prep_surv(): the `time` column "%s" has %d infinite value%s.',
      time, n_inf, if (n_inf > 1L) "s" else ""),
      call. = FALSE)
  }
  neg <- tnum < 0
  if (any(neg)) {
    stop(sprintf(
      paste0('prep_surv(): the `time` column "%s" has %d negative value%s ',
             "(minimum %s). A negative follow-up time is not analysable, so ",
             "prep_surv() stops instead of computing on it."),
      time, sum(neg), if (sum(neg) > 1L) "s" else "",
      format(min(tnum))),
      call. = FALSE)
  }

  ## -- 5. vital status -----------------------------------------------------
  .ps_value_arg(alive_value, "alive_value")
  .ps_value_arg(dead_value, "dead_value")
  both <- intersect(as.character(alive_value), as.character(dead_value))
  if (length(both) > 0L) {
    stop(sprintf(
      paste0("prep_surv(): `alive_value` and `dead_value` must not share ",
             "values, but %s appear%s in both."),
      .ps_values(both), if (length(both) > 1L) "" else "s"),
      call. = FALSE)
  }

  sv <- data[[status]]
  if (is.factor(sv)) sv <- as.character(sv)
  is_alive <- sv %in% alive_value
  is_dead  <- sv %in% dead_value
  unknown  <- !(is_alive | is_dead)
  if (any(unknown)) {
    stop(sprintf(
      paste0('prep_surv(): the `status` column "%s" has %d of %d row%s whose ',
             "value is neither `alive_value` nor `dead_value`.\n",
             "  offending values : %s\n",
             "  alive_value      : %s\n",
             "  dead_value       : %s\n",
             "  A missing status (<NA>) counts as an offending value: ",
             "prep_surv() never guesses whether a subject died."),
      status, sum(unknown), n, if (n > 1L) "s" else "",
      .ps_counts(sv[unknown]), .ps_values(alive_value), .ps_values(dead_value)),
      call. = FALSE)
  }

  event_os <- integer(n)
  event_os[is_dead] <- 1L

  ## -- 6. competing-risk outcome -------------------------------------------
  event_cr  <- NULL
  n_unk_cod <- 0L

  if (!is.null(cause)) {

    if (is.null(cause_gist_value)) {
      stop(paste0("prep_surv(): `cause` was supplied but `cause_gist_value` ",
                  "is NULL. prep_surv() will not decide on its own which ",
                  "values of the cause column count as a death from the ",
                  "disease of interest."),
           call. = FALSE)
    }
    if (length(cause_unknown_to) != 1L || is.na(cause_unknown_to) ||
        !(cause_unknown_to %in% c(1L, 2L))) {
      stop(sprintf(
        paste0("prep_surv(): `cause_unknown_to` must be 1 (death from the ",
               "disease of interest) or 2 (death from other causes), got %s.\n",
               "  0 is refused: it would record a dead subject as censored in ",
               "`event_cr` while `event_os` records a death, so the two ",
               "outcome variables would contradict each other."),
        .ps_values(cause_unknown_to)),
        call. = FALSE)
    }
    cause_unknown_to <- as.integer(cause_unknown_to)

    sets <- list(cause_gist_value    = cause_gist_value,
                 cause_other_value   = cause_other_value,
                 cause_unknown_value = cause_unknown_value,
                 cause_censor_value  = cause_censor_value)
    sets <- sets[!vapply(sets, is.null, logical(1))]
    for (i in seq_along(sets)) {
      .ps_value_arg(sets[[i]], names(sets)[i], allow_na = TRUE)
      for (j in seq_along(sets)) {
        if (j <= i) next
        shared <- intersect(sets[[i]], sets[[j]])
        if (length(shared) > 0L) {
          stop(sprintf(
            paste0("prep_surv(): `%s` and `%s` must not overlap, but %s ",
                   "appear%s in both. One cause of death cannot map to two ",
                   "competing-risk categories."),
            names(sets)[i], names(sets)[j], .ps_values(shared),
            if (length(shared) > 1L) "" else "s"),
            call. = FALSE)
        }
      }
    }

    cv <- data[[cause]]
    if (is.factor(cv)) cv <- as.character(cv)

    hit_gist <- cv %in% cause_gist_value
    hit_unk  <- cv %in% cause_unknown_value
    hit_oth  <- cv %in% cause_other_value
    hit_cens <- cv %in% cause_censor_value

    ## (a) alive, yet carrying a declared cause of death
    death_decl <- c(cause_gist_value, cause_other_value, cause_unknown_value)
    death_decl <- death_decl[!is.na(death_decl)]
    bad_alive  <- is_alive & (cv %in% death_decl)
    if (any(bad_alive)) {
      stop(sprintf(
        paste0("prep_surv(): `status` and `cause` contradict each other. %d ",
               'row%s recorded as alive (event_os = 0) but column "%s" holds ',
               "a declared cause of death.\n",
               "  offending values : %s\n",
               "  A censored subject cannot have died of a specific cause. ",
               "If this value is the placeholder your registry gives living ",
               "subjects, leave it out of the cause_*_value sets."),
        sum(bad_alive), if (sum(bad_alive) > 1L) "s are" else " is", cause,
        .ps_counts(cv[bad_alive])),
        call. = FALSE)
    }

    ## (b) dead, yet the cause column says "no event"
    bad_dead <- is_dead & hit_cens
    if (any(bad_dead)) {
      stop(sprintf(
        paste0("prep_surv(): `status` and `cause` contradict each other. %d ",
               'row%s recorded as dead (event_os = 1) but column "%s" holds a ',
               "`cause_censor_value`, i.e. no event.\n",
               "  offending values : %s"),
        sum(bad_dead), if (sum(bad_dead) > 1L) "s are" else " is", cause,
        .ps_counts(cv[bad_dead])),
        call. = FALSE)
    }

    event_cr <- integer(n)                       # alive -> 0
    event_cr[is_dead & hit_gist] <- 1L
    event_cr[is_dead & hit_unk]  <- cause_unknown_to
    event_cr[is_dead & hit_oth]  <- 2L

    rest <- is_dead & !hit_gist & !hit_unk & !hit_oth
    if (any(rest)) {
      if (is.null(cause_other_value)) {
        event_cr[rest] <- 2L                     # open coding
      } else {
        stop(sprintf(
          paste0('prep_surv(): the `cause` column "%s" has %d dead subject%s ',
                 "whose value is in none of `cause_gist_value`, ",
                 "`cause_other_value`, `cause_unknown_value`.\n",
                 "  offending values : %s\n",
                 "  Add these values to one of the three sets, or set ",
                 "cause_other_value = NULL to let any undeclared cause count ",
                 "as a death from other causes."),
          cause, sum(rest), if (sum(rest) > 1L) "s" else "",
          .ps_counts(cv[rest])),
          call. = FALSE)
      }
    }
    n_unk_cod <- sum(is_dead & hit_unk)

    ## (c) the bidirectional cross-check, on the derived columns
    inv <- c(
      "event_os == 0 only with event_cr == 0" =
        all(event_cr[event_os == 0L] == 0L),
      "event_os == 1 only with event_cr in {1, 2}" =
        all(event_cr[event_os == 1L] %in% c(1L, 2L)),
      "event_cr == 0 only with event_os == 0" =
        all(event_os[event_cr == 0L] == 0L),
      "event_cr in {1, 2} only with event_os == 1" =
        all(event_os[event_cr %in% c(1L, 2L)] == 1L),
      "deaths == event_cr 1 plus event_cr 2" =
        sum(event_os == 1L) == sum(event_cr == 1L) + sum(event_cr == 2L)
    )
    if (!all(inv)) {
      stop(sprintf(
        paste0("prep_surv(): the derived `event_os` and `event_cr` do not ",
               "agree, which should be impossible by construction.\n",
               "  failed check%s : %s\n",
               "  Do not use this result; report it as a bug in prep_surv()."),
        if (sum(!inv) > 1L) "s" else "",
        paste(names(inv)[!inv], collapse = "; ")),
        call. = FALSE)
    }
  }

  ## -- 7. assemble ---------------------------------------------------------
  out <- data
  out[["time_mo"]]  <- tnum
  out[["event_os"]] <- event_os
  if (!is.null(event_cr)) out[["event_cr"]] <- event_cr

  counts <- c(subjects = n,
              censored = sum(event_os == 0L),
              deaths   = sum(event_os == 1L))
  if (!is.null(event_cr)) {
    counts <- c(counts,
                cr_censored           = sum(event_cr == 0L),
                cr_disease            = sum(event_cr == 1L),
                cr_other              = sum(event_cr == 2L),
                cause_unknown_recoded = n_unk_cod)
  }

  checks <- c("time non-missing"       = TRUE,
              "time non-negative"      = TRUE,
              "status values declared" = TRUE)
  if (!is.null(event_cr)) checks <- c(checks, inv)

  attr(out, "os_cr_table") <- if (is.null(event_cr)) NULL else
    table(event_os, event_cr, dnn = c("event_os", "event_cr"))
  attr(out, "prep_surv") <- list(
    counts    = counts,
    follow_up = c(min    = min(tnum),
                  median = stats::median(tnum),
                  max    = max(tnum),
                  n_zero = sum(tnum == 0)),
    checks    = checks,
    settings  = list(time                = time,
                     status              = status,
                     alive_value         = alive_value,
                     dead_value          = dead_value,
                     cause               = cause,
                     cause_gist_value    = cause_gist_value,
                     cause_other_value   = cause_other_value,
                     cause_unknown_value = cause_unknown_value,
                     cause_censor_value  = cause_censor_value,
                     cause_unknown_to    = if (is.null(cause)) NULL else
                       cause_unknown_to),
    call      = cl
  )
  out
}


# -- internal helpers, not exported ------------------------------------------

#' @noRd
.ps_name_arg <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf(
      "prep_surv(): `%s` must be a single non-missing column name (a string).",
      arg), call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.ps_value_arg <- function(x, arg, allow_na = FALSE) {
  if (!is.atomic(x) || length(x) == 0L) {
    stop(sprintf(
      "prep_surv(): `%s` must be an atomic vector of at least one value.",
      arg), call. = FALSE)
  }
  if (!allow_na && anyNA(x)) {
    stop(sprintf("prep_surv(): `%s` must not contain NA.", arg), call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.ps_values <- function(x, max_show = 8L) {
  u       <- unique(x)
  n_extra <- max(0L, length(u) - max_show)
  u       <- u[seq_len(min(length(u), max_show))]
  s       <- as.character(u)
  if (is.character(u) || is.factor(u)) s <- paste0('"', s, '"')
  s[is.na(u)] <- "<NA>"
  s <- paste(s, collapse = ", ")
  if (n_extra > 0L) s <- paste0(s, ", ... (", n_extra, " more)")
  s
}

#' @noRd
.ps_counts <- function(x, max_show = 8L) {
  lab <- as.character(x)
  if (is.character(x) || is.factor(x)) lab <- paste0('"', lab, '"')
  lab[is.na(x)] <- "<NA>"
  tb      <- sort(table(lab), decreasing = TRUE)
  n_extra <- max(0L, length(tb) - max_show)
  tb      <- tb[seq_len(min(length(tb), max_show))]
  s <- paste0(names(tb), " (n = ", as.integer(tb), ")", collapse = ", ")
  if (n_extra > 0L) {
    s <- paste0(s, ", ... (", n_extra, " more distinct values)")
  }
  s
}
