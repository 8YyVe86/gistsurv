# ---------------------------------------------------------------------------
# make-sim-gist.R -- generate data/sim_gist.rda
#
# WHAT THIS IS
#   A synthetic cohort used by the examples, the tests and the vignette, so
#   that the package can be exercised end to end without any registry data.
#
# WHAT THIS IS NOT
#   It is not a sample, a subset, a perturbation or a summary of SEER data,
#   and it is not fitted to any real cohort. Every number below is written
#   into this file by hand. The script opens no file and connects to nothing:
#   its only inputs are the literals you can read here and the seed. The
#   marginal distributions and the hazards were chosen to exercise each
#   branch of the seven functions (both competing causes present in every
#   stratum, follow-up long enough for tau = 60, a covariate with
#   age-dependent missingness, an unknown cause of death), NOT to resemble
#   any particular population. Do not read an epidemiological claim off it.
#
#   The variable names, factor levels and reference levels DO mirror the
#   analysis set of the underlying study. That is the point: a user should be
#   able to swap sim_gist for their own extract and have the calls still fit.
#   Names and level labels are structure, not data.
#
# HOW TO RUN
#   Rscript data-raw/make-sim-gist.R
#   (run from the package root; it rewrites data/sim_gist.rda)
#
# ASCII-only, like the rest of the package.
# ---------------------------------------------------------------------------

set.seed(20260901)

n <- 1200L

## -- 1. covariates ----------------------------------------------------------

SITE_LEVELS  <- c("Stomach", "Small intestine", "Colorectal")
AGE_LEVELS   <- c("<50", "50-64", "65-74", "75+")
RACE_LEVELS  <- c("NH White", "NH Black", "NH API", "Hispanic", "Other/Unknown")
SIZE_LEVELS  <- c("<=20", "21-50", "51-100", ">100")
STAGE_LEVELS <- c("Localized", "Regional", "Distant", "Unknown")

site <- factor(sample(SITE_LEVELS, n, TRUE, prob = c(0.62, 0.30, 0.08)),
               levels = SITE_LEVELS)
age_grp <- factor(sample(AGE_LEVELS, n, TRUE, prob = c(0.18, 0.36, 0.25, 0.21)),
                  levels = AGE_LEVELS)
# age in whole years, drawn inside the band its group names; the study's own
# age variable is a 5-year band lower bound, so integers are the right shape.
age_lo  <- c("<50" = 20L, "50-64" = 50L, "65-74" = 65L, "75+" = 75L)[as.character(age_grp)]
age_hi  <- c("<50" = 49L, "50-64" = 64L, "65-74" = 74L, "75+" = 90L)[as.character(age_grp)]
age     <- as.integer(round(age_lo + (age_hi - age_lo) * stats::runif(n)))

sex  <- factor(sample(c("Male", "Female"), n, TRUE, prob = c(0.51, 0.49)),
               levels = c("Male", "Female"))
race <- factor(sample(RACE_LEVELS, n, TRUE,
                      prob = c(0.55, 0.18, 0.14, 0.12, 0.01)),
               levels = RACE_LEVELS)

# tumour size: lognormal, shifted a little by site, then held to the 400 mm
# ceiling the registry applies.
size_mu <- c(Stomach = 3.5, "Small intestine" = 3.8,
             Colorectal = 3.4)[as.character(site)]
size_mm <- round(pmin(400, exp(stats::rnorm(n, size_mu, 0.75))))
size_mm[size_mm < 1] <- 1

stage <- factor(sample(STAGE_LEVELS, n, TRUE, prob = c(0.58, 0.22, 0.16, 0.04)),
                levels = STAGE_LEVELS)

# surgery: less likely with age and with distant disease
p_surg <- 0.90 -
  c("<50" = 0.00, "50-64" = 0.03, "65-74" = 0.09, "75+" = 0.22)[as.character(age_grp)] -
  c(Localized = 0.00, Regional = 0.04, Distant = 0.30,
    Unknown = 0.12)[as.character(stage)]
surgery <- factor(ifelse(stats::runif(n) < p_surg, "Yes", "No"),
                  levels = c("No", "Yes"))

dx_year <- sample(2000:2019, n, TRUE)

## -- 2. missingness ---------------------------------------------------------
# Size is missing more often in older patients. This is deliberate: it makes
# the complete-case subset differ from the full cohort in a way the functions'
# deletion counters are supposed to reveal. A dataset with MCAR missingness
# would let a broken counter pass unnoticed.
p_size_na <- c("<50" = 0.06, "50-64" = 0.09,
               "65-74" = 0.14, "75+" = 0.20)[as.character(age_grp)]
size_mm[stats::runif(n) < p_size_na] <- NA_integer_

size_grp <- cut(size_mm, breaks = c(0, 20, 50, 100, 400),
                labels = SIZE_LEVELS, right = TRUE)
size_grp <- factor(as.character(size_grp), levels = SIZE_LEVELS)

surgery[stats::runif(n) < 0.02] <- NA   # a handful of unknown surgery codes

## -- 3. two competing causes ------------------------------------------------
# Cause-specific exponential hazards, per month. Chosen so that roughly 45% of
# the cohort dies within the observation window and both causes appear in
# every site and every age group. They are not estimates of anything.

lin_gist <-
  log(0.0014) +
  c(Stomach = 0.00, "Small intestine" = 0.15, Colorectal = 0.25)[as.character(site)] +
  c("<50" = 0.00, "50-64" = 0.10, "65-74" = 0.25, "75+" = 0.45)[as.character(age_grp)] +
  c(Localized = 0.00, Regional = 0.55, Distant = 1.30,
    Unknown = 0.30)[as.character(stage)] +
  c("<=20" = 0.00, "21-50" = 0.20, "51-100" = 0.55,
    ">100" = 0.95)[as.character(size_grp)] +
  ifelse(surgery %in% "Yes", -0.45, 0)
# size and surgery are missing for some rows; an unknown covariate must not
# change the hazard, so the missing contribution is simply zero.
lin_gist[is.na(lin_gist)] <-
  log(0.0014) + c(Stomach = 0.00, "Small intestine" = 0.15,
                  Colorectal = 0.25)[as.character(site)][is.na(lin_gist)]

lin_other <-
  log(0.00040) +
  c("<50" = 0.00, "50-64" = 0.75, "65-74" = 1.60, "75+" = 2.55)[as.character(age_grp)] +
  ifelse(sex == "Female", -0.20, 0) +
  c("NH White" = 0.00, "NH Black" = 0.20, "NH API" = -0.15,
    Hispanic = -0.05, "Other/Unknown" = 0.00)[as.character(race)]

t_gist  <- stats::rexp(n, exp(lin_gist))
t_other <- stats::rexp(n, exp(lin_other))

## -- 4. follow-up -----------------------------------------------------------
# Administrative censoring: everyone is followed to the same calendar cutoff,
# so a 2019 diagnosis has far less potential follow-up than a 2000 one. This
# is what makes tau = 60 admissible in some strata and not others, which is
# exactly the guard calc_rmst() exists to enforce.
STUDY_END <- 2024                     # calendar year of the notional cutoff
t_admin   <- (STUDY_END - dx_year) * 12 + sample(0:11, n, TRUE)
t_ltfu    <- stats::rexp(n, 1 / 900)  # rare loss to follow-up

t_obs   <- pmin(t_gist, t_other, t_admin, t_ltfu)
died    <- t_obs == pmin(t_gist, t_other) & t_obs < pmin(t_admin, t_ltfu)
time_mo <- as.numeric(floor(t_obs))   # whole months, as the registry reports

## -- 5. raw outcome columns, in the registry's own coding -------------------
# The cause-of-death column is deliberately left in its awkward registry form:
# one value covers both the living and those who died of something else, so it
# cannot be read on its own. prep_surv() resolves it against vital_status.
COD_ALIVE_OR_OTHER <- "Alive or dead of other cause"
COD_DEAD_CANCER    <- "Dead (attributable to this cancer dx)"
COD_DEAD_UNKNOWN   <- "Dead (missing/unknown COD)"

vital_status <- ifelse(died, "Dead", "Alive")
of_gist      <- died & (t_gist <= t_other)
cause_of_death <- ifelse(of_gist, COD_DEAD_CANCER, COD_ALIVE_OR_OTHER)

# about 1% of deaths carry no usable cause; the count of these is the figure a
# limitations paragraph has to quote, so the example data must contain some.
unk <- which(died)[stats::runif(sum(died)) < 0.012]
cause_of_death[unk] <- COD_DEAD_UNKNOWN

sim_gist <- data.frame(
  id             = seq_len(n),
  site           = site,
  age            = age,
  age_grp        = age_grp,
  sex            = sex,
  race           = race,
  size_mm        = as.numeric(size_mm),
  size_grp       = size_grp,
  stage          = stage,
  surgery        = surgery,
  dx_year        = as.integer(dx_year),
  time_mo        = time_mo,
  vital_status   = vital_status,
  cause_of_death = cause_of_death,
  stringsAsFactors = FALSE
)

## -- 6. the properties the rest of the package relies on --------------------
# If a change to the model above breaks one of these, the examples, tests and
# vignette break too -- so fail here rather than there.
tab <- table(sim_gist$site, sim_gist$age_grp)
stopifnot(
  "n must be 1200"                = nrow(sim_gist) == 1200L,
  "no missing outcome"            = !anyNA(sim_gist$time_mo) &&
                                    !anyNA(sim_gist$vital_status) &&
                                    !anyNA(sim_gist$cause_of_death),
  "follow-up non-negative"        = all(sim_gist$time_mo >= 0),
  "every site x age cell used"    = all(tab >= 15L),
  "deaths between 35% and 55%"    = {
    d <- mean(sim_gist$vital_status == "Dead"); d > 0.35 && d < 0.55 },
  "both causes in every site"     = {
    dead <- sim_gist$vital_status == "Dead"
    gist <- dead & sim_gist$cause_of_death == COD_DEAD_CANCER
    oth  <- dead & sim_gist$cause_of_death != COD_DEAD_CANCER
    all(tapply(gist, sim_gist$site, any)) && all(tapply(oth, sim_gist$site, any)) &&
      all(tapply(gist, sim_gist$age_grp, any)) && all(tapply(oth, sim_gist$age_grp, any))
  },
  "some unknown causes of death"  = sum(sim_gist$cause_of_death == COD_DEAD_UNKNOWN) > 0L,
  "tau = 60 admissible everywhere" =
    min(tapply(sim_gist$time_mo, sim_gist$site, max)) > 60 &&
    min(tapply(sim_gist$time_mo, sim_gist$age_grp, max)) > 60,
  "size_grp has missing values"   = anyNA(sim_gist$size_grp),
  "surgery has missing values"    = anyNA(sim_gist$surgery),
  "enough events per coefficient" = sum(sim_gist$vital_status == "Dead") > 17L * 10L
)

cat("sim_gist:", nrow(sim_gist), "rows,", ncol(sim_gist), "columns\n")
cat("  deaths          :", sum(sim_gist$vital_status == "Dead"),
    sprintf("(%.1f%%)\n", 100 * mean(sim_gist$vital_status == "Dead")))
cat("  of the cancer   :", sum(sim_gist$cause_of_death == COD_DEAD_CANCER), "\n")
cat("  unknown cause   :", sum(sim_gist$cause_of_death == COD_DEAD_UNKNOWN), "\n")
cat("  size_grp missing:", sum(is.na(sim_gist$size_grp)), "\n")
cat("  surgery missing :", sum(is.na(sim_gist$surgery)), "\n")
cat("  follow-up (mo)  :", min(sim_gist$time_mo), "-", max(sim_gist$time_mo), "\n")
print(tab)

usethis::use_data(sim_gist, overwrite = TRUE)
