test_that("the counts in the cells add up to the group sizes", {
  d   <- prep_sim()
  tab <- tab_baseline(d, by = "site",
                      vars = c("age", "age_grp", "sex", "stage"))

  # cells read "144 (19.5%)"; pull the count and the percentage back out, so
  # that the formatting is tested along with the arithmetic
  n_of   <- function(x) as.numeric(sub("^([0-9]+).*$", "\\1", x))
  pct_of <- function(x) as.numeric(sub("^.*\\(([0-9.]+)%\\)$", "\\1", x))

  for (v in c("age_grp", "sex", "stage")) {
    rows <- tab[tab$variable == v & tab$row_type == "level", , drop = FALSE]
    expect_equal(nrow(rows), nlevels(d[[v]]))
    for (g in c("Overall", levels(d$site))) {
      n_g <- if (g == "Overall") nrow(d) else sum(d$site == g)
      expect_equal(sum(n_of(rows[[g]])), n_g)
      expect_equal(sum(pct_of(rows[[g]])), 100, tolerance = 0.05)
    }
  }

  # and one cell checked against the data directly
  r <- tab[tab$variable == "sex" & tab$label == "Female", ]
  expect_equal(n_of(r$Stomach), sum(d$site == "Stomach" & d$sex == "Female"))
})

test_that("the SMD matches the formula it documents", {
  d   <- prep_sim()
  two <- d[d$site %in% c("Stomach", "Small intestine"), , drop = FALSE]
  two$site <- droplevels(two$site)
  tab <- tab_baseline(two, by = "site", vars = c("age", "sex"))

  g1 <- two[two$site == "Stomach", ]
  g2 <- two[two$site == "Small intestine", ]

  # continuous: (m1 - m2) / sqrt((s1^2 + s2^2) / 2)
  want <- (mean(g1$age) - mean(g2$age)) /
    sqrt((stats::var(g1$age) + stats::var(g2$age)) / 2)
  got <- tab$smd[tab$variable == "age" & tab$row_type == "label"]
  expect_equal(abs(got), abs(want), tolerance = 1e-8)

  # binary: |p1 - p2| / sqrt((p1(1-p1) + p2(1-p2)) / 2)
  p1 <- mean(g1$sex == "Male"); p2 <- mean(g2$sex == "Male")
  want <- abs(p1 - p2) / sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
  got <- tab$smd[tab$variable == "sex" & tab$row_type == "label"]
  expect_equal(abs(got), want, tolerance = 1e-8)

  # with three groups the reported value is the largest over pairs, and
  # smd_pair has to name that pair rather than an arbitrary one
  tab3 <- tab_baseline(d, by = "site", vars = "sex")
  pairs <- utils::combn(levels(d$site), 2, simplify = FALSE)
  each <- vapply(pairs, function(p) {
    q1 <- mean(d$sex[d$site == p[1]] == "Male")
    q2 <- mean(d$sex[d$site == p[2]] == "Male")
    abs(q1 - q2) / sqrt((q1 * (1 - q1) + q2 * (1 - q2)) / 2)
  }, numeric(1))
  got3 <- tab3$smd[tab3$variable == "sex" & tab3$row_type == "label"]
  expect_equal(abs(got3), max(each), tolerance = 1e-8)
  expect_true(grepl(pairs[[which.max(each)]][1], got3_pair <-
                      tab3$smd_pair[tab3$variable == "sex" &
                                      tab3$row_type == "label"], fixed = TRUE))
  expect_true(grepl(pairs[[which.max(each)]][2], got3_pair, fixed = TRUE))
})

test_that("a missing grouping value is refused unless told what to do with it", {
  d <- prep_sim()
  d$site[1:5] <- NA

  expect_error(tab_baseline(d, by = "site", vars = "sex"), "missing")

  # the column is labelled "(Missing)" -- parenthesised, so that it cannot be
  # mistaken for a real level of the grouping variable
  kept <- tab_baseline(d, by = "site", vars = "sex", by_missing = "level",
                       missing_text = "Missing")
  expect_true("(Missing)" %in% names(kept))
  n_of <- function(x) as.numeric(sub("^([0-9]+).*$", "\\1", x))
  expect_equal(sum(n_of(kept[["(Missing)"]][kept$variable == "sex" &
                                              kept$row_type == "level"])),
               sum(is.na(d$site)))

  # and missing_text is what gets parenthesised
  renamed <- tab_baseline(d, by = "site", vars = "sex", by_missing = "level",
                          missing_text = "Not recorded")
  expect_true("(Not recorded)" %in% names(renamed))

  dropped <- tab_baseline(d, by = "site", vars = "sex", by_missing = "drop")
  expect_false(any(grepl("Missing", names(dropped))))
  n_of <- function(x) as.numeric(sub("^([0-9]+).*$", "\\1", x))
  rows <- dropped[dropped$variable == "sex" & dropped$row_type == "level", ]
  expect_equal(sum(n_of(rows$Overall)), sum(!is.na(d$site)))
})
