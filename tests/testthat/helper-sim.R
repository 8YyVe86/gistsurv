# Shared fixtures. Sourced by testthat before the test files.
#
# Every test starts from sim_gist, so a change to the example data that breaks
# an invariant breaks the tests too -- which is the point: the generating
# script's stopifnot block and these tests are meant to fail together.
# ASCII-only, like the rest of the package.

GIST_VALUE  <- "Dead (attributable to this cancer dx)"
OTHER_VALUE <- "Other cause"
UNK_VALUE   <- "Dead (missing/unknown COD)"

# sim_gist with event_os / event_cr, in closed coding.
prep_sim <- function(data = gistsurv::sim_gist, unknown_to = 2L) {
  prep_surv(
    data,
    time                = "time_mo",
    status              = "vital_status",
    cause               = "cause_of_death",
    cause_gist_value    = GIST_VALUE,
    cause_other_value   = OTHER_VALUE,
    cause_unknown_value = UNK_VALUE,
    cause_censor_value  = NA,
    cause_unknown_to    = unknown_to
  )
}

# A two-arm subset, for the tests that compare against survRM2 / coxph
# directly. Smaller than the whole cohort so the suite stays quick.
two_arm <- function(a = "Stomach", b = "Small intestine") {
  d <- prep_sim()
  d[d$site %in% c(a, b), , drop = FALSE]
}
