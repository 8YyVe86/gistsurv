test_that("event_os and event_cr cannot contradict each other", {
  d <- prep_sim()

  # the invariant the function promises, checked on the returned columns
  expect_true(all(d$event_cr[d$event_os == 0L] == 0L))
  expect_true(all(d$event_cr[d$event_os == 1L] %in% 1:2))
  expect_equal(sum(d$event_os == 1L),
               sum(d$event_cr == 1L) + sum(d$event_cr == 2L))

  # and the cross-tabulation it hands back must agree with those columns
  tab <- attr(d, "os_cr_table")
  expect_equal(as.vector(tab), as.vector(table(d$event_os, d$event_cr)))
  expect_equal(sum(tab), nrow(d))

  cnt <- attr(d, "prep_surv")$counts
  expect_equal(as.integer(cnt[["deaths"]]), sum(d$event_os == 1L))
  expect_equal(as.integer(cnt[["cr_disease"]]), sum(d$event_cr == 1L))
  expect_equal(as.integer(cnt[["cr_other"]]), sum(d$event_cr == 2L))
})

test_that("an undeclared or contradictory cause stops the call", {
  # a code nobody declared, on a dead subject: the whole reason for closed
  # coding. Open coding would fold this into "death from other causes".
  bad <- sim_gist
  bad$cause_of_death[which(bad$vital_status == "Dead")[1]] <- "Dead (new code)"
  expect_error(prep_sim(bad), "cause")

  # a living subject carrying a declared cause of death
  bad <- sim_gist
  bad$cause_of_death[which(bad$vital_status == "Alive")[1]] <- GIST_VALUE
  expect_error(prep_sim(bad), "contradict")

  # a dead subject carrying the censoring placeholder
  bad <- sim_gist
  bad$cause_of_death[which(bad$vital_status == "Dead")[1]] <- NA
  expect_error(prep_sim(bad), "contradict")

  # censoring a dead subject in event_cr is refused outright
  expect_error(prep_sim(unknown_to = 0L), "0")
})

test_that("an unknown cause of death goes where it is told, and is counted", {
  to2 <- prep_sim(unknown_to = 2L)
  to1 <- prep_sim(unknown_to = 1L)

  n_unk <- sum(sim_gist$cause_of_death %in% UNK_VALUE)
  expect_gt(n_unk, 0L)   # the fixture must contain some, or this proves nothing

  expect_equal(as.integer(attr(to2, "prep_surv")$counts[["cause_unknown_recoded"]]),
               n_unk)

  # moving them from 2 to 1 must move exactly that many subjects, and nobody
  # else: overall survival is untouched by the choice.
  expect_equal(sum(to1$event_cr == 1L), sum(to2$event_cr == 1L) + n_unk)
  expect_equal(sum(to1$event_cr == 2L), sum(to2$event_cr == 2L) - n_unk)
  expect_equal(to1$event_os, to2$event_os)
})
