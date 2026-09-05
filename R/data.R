# ---------------------------------------------------------------------------
# Documentation for the data shipped with the package. The data themselves are
# built by data-raw/make-sim-gist.R; this file only describes them.
# ASCII-only, like the rest of the package.
# ---------------------------------------------------------------------------

#' A simulated GIST cohort
#'
#' @description
#' A synthetic cohort of 1200 subjects, used by the examples, the tests and
#' the vignette so that every function in the package can be run end to end
#' without access to registry data.
#'
#' The outcome is supplied in the form a registry export arrives in -- a
#' follow-up time, a vital-status string and a cause-of-death string -- and not
#' as ready-made event indicators. Call [prep_surv()] first to obtain
#' `event_os` and `event_cr`; the other six functions take it from there.
#'
#' @section What this is not:
#' `sim_gist` is **not** a sample, a subset, a perturbation or a summary of
#' SEER data, and it is not fitted to any real cohort. Every parameter is a
#' literal written into `data-raw/make-sim-gist.R`, which opens no file. The
#' marginal distributions and the two cause-specific hazards were chosen to
#' exercise the branches of the seven functions -- both competing causes
#' present in every stratum, follow-up long enough that `tau = 60` is
#' admissible, a covariate whose missingness depends on age, an unknown cause
#' of death -- and not to resemble any population. No epidemiological claim
#' can be read off it.
#'
#' What the dataset does borrow from the study the package was distilled from
#' is its *structure*: the variable names, the factor levels and the reference
#' levels. That is deliberate, so that a user can substitute their own extract
#' for `sim_gist` and have the same calls still fit.
#'
#' @format A data frame with 1200 rows and 14 columns:
#' \describe{
#'   \item{id}{Row identifier, 1 to 1200. Not a patient identifier.}
#'   \item{site}{Primary site, a factor with levels `"Stomach"` (reference),
#'     `"Small intestine"`, `"Colorectal"`.}
#'   \item{age}{Age at diagnosis in whole years, 20 to 90.}
#'   \item{age_grp}{Age group, a factor with levels `"<50"` (reference),
#'     `"50-64"`, `"65-74"`, `"75+"`.}
#'   \item{sex}{A factor with levels `"Male"` (reference), `"Female"`.}
#'   \item{race}{A factor with levels `"NH White"` (reference), `"NH Black"`,
#'     `"NH API"`, `"Hispanic"`, `"Other/Unknown"`.}
#'   \item{size_mm}{Largest tumour dimension in millimetres, 1 to 400.
#'     Missing for 150 subjects, more often in the older groups.}
#'   \item{size_grp}{Size band, a factor with levels `"<=20"` (reference),
#'     `"21-50"`, `"51-100"`, `">100"`. Missing wherever `size_mm` is.}
#'   \item{stage}{A factor with levels `"Localized"` (reference),
#'     `"Regional"`, `"Distant"`, `"Unknown"`.}
#'   \item{surgery}{Surgery of the primary site, a factor with levels `"No"`
#'     (reference) and `"Yes"`. Missing for 25 subjects.}
#'   \item{dx_year}{Year of diagnosis, 2000 to 2019. Determines how much
#'     potential follow-up a subject has, because everyone is censored at the
#'     same notional cutoff.}
#'   \item{time_mo}{Follow-up in whole months, 0 to 299.}
#'   \item{vital_status}{`"Alive"` or `"Dead"`. 564 subjects died.}
#'   \item{cause_of_death}{Cause of death, `NA` for the 636 living subjects and
#'     one of three values for the 564 who died:
#'     `"Dead (attributable to this cancer dx)"` (332),
#'     `"Other cause"` (226), `"Dead (missing/unknown COD)"` (6).
#'     The registry's own label for the middle group is not what you see here;
#'     see the section below.}
#' }
#'
#' @section The one step you must take before prep_surv():
#' The registry codes cause of death with a value -- `"Alive or dead of other
#' cause"` -- that is carried by living subjects **and** by those who died of
#' something else. It cannot be handed to [prep_surv()] as
#' `cause_other_value`: the function refuses a declared cause of death on a
#' living subject, and it is right to, because that same check is what catches
#' a genuine status/cause contradiction everywhere else.
#'
#' Resolve the ambiguity first, in data preparation, where `vital_status` is
#' at hand to resolve it with:
#'
#' ```
#' raw$cause_of_death <- ifelse(
#'   raw$cause_of_death == "Alive or dead of other cause",
#'   ifelse(raw$vital_status == "Dead", "Other cause", NA),
#'   raw$cause_of_death)
#' ```
#'
#' `sim_gist` ships with those two lines already applied -- that is why its
#' cause column reads `"Other cause"` and `NA` rather than the registry
#' string -- so the documented call below runs as it stands. Your own extract
#' will need them.
#'
#' The alternative is to leave `cause_other_value` at its default `NULL` and
#' let open coding sweep the undeclared label into `event_cr = 2`. That gives
#' the same numbers here and costs the whole closed set: a registry code
#' nobody anticipated, or a typo, would also become a death from other causes,
#' silently. Preparing the column instead keeps the closed set intact, and an
#' unexpected value on a dead subject remains an error.
#'
#' Whichever route you take, read `attr(d, "os_cr_table")` and the counts in
#' `attr(d, "prep_surv")` before going further.
#'
#' @source
#' Generated by the script `data-raw/make-sim-gist.R` with
#' `set.seed(20260901)`. That script lives in the package's **source
#' repository**, not in the installed package: `data-raw/` is listed in
#' `.Rbuildignore`, as is conventional, so it is not shipped in the tarball or
#' the installed library. Re-running it rebuilds `data/sim_gist.rda`; the same
#' seed gives the same 1200 rows.
#'
#' @examples
#' str(sim_gist)
#'
#' ## cause of death is NA for the living and declared for every death
#' table(sim_gist$vital_status, sim_gist$cause_of_death, useNA = "ifany")
#'
#' ## prep_surv() turns that into event_os / event_cr. Every admissible value
#' ## is declared, so an unexpected code would stop the call rather than be
#' ## folded into "death from other causes".
#' d <- prep_surv(
#'   sim_gist,
#'   time                = "time_mo",
#'   status              = "vital_status",
#'   cause               = "cause_of_death",
#'   cause_gist_value    = "Dead (attributable to this cancer dx)",
#'   cause_other_value   = "Other cause",
#'   cause_unknown_value = "Dead (missing/unknown COD)",
#'   cause_censor_value  = NA,
#'   cause_unknown_to    = 2L
#' )
#' attr(d, "os_cr_table")
#' attr(d, "prep_surv")$counts
"sim_gist"
