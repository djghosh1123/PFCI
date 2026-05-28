test_that("simulate_pfci_toy returns correct structure", {
  skip_if_not_installed("pcalg")
  skip_if_not_installed("graph")
  sim <- simulate_pfci_toy(p = 20, n = 50, edge_prob = 0.05, seed = 1)
  expect_true(is.list(sim))
  expect_true(is.matrix(sim$X))
  expect_equal(nrow(sim$X), 50)
  expect_equal(ncol(sim$X), 20)
  expect_true(is.matrix(sim$truth$skel))
})

test_that("pfci_fit returns pfci_fit object", {
  skip_if_not_installed("pcalg")
  sim <- simulate_pfci_toy(p = 20, n = 50, edge_prob = 0.05, seed = 1)
  fit <- pfci_fit(sim$X, alpha = 0.05)
  expect_s3_class(fit, "pfci_fit")
  expect_true(is.matrix(fit$amat))
  expect_equal(nrow(fit$amat), 20)
  expect_true(is.numeric(fit$time$total))
})

test_that("pfci_metrics returns correct fields", {
  skip_if_not_installed("pcalg")
  sim <- simulate_pfci_toy(p = 20, n = 50, edge_prob = 0.05, seed = 1)
  fit <- pfci_fit(sim$X, alpha = 0.05)
  met <- pfci_metrics(sim, fit)
  expect_true(is.list(met))
  expect_true(all(c("SHD", "F1_total", "MCC",
    "Precision", "Recall") %in% names(met)))
  expect_true(met$F1_total >= 0 && met$F1_total <= 1)
  expect_true(met$MCC >= -1 && met$MCC <= 1)
})

test_that("pfci_metrics compute_marks works", {
  skip_if_not_installed("pcalg")
  sim <- simulate_pfci_toy(p = 20, n = 50, edge_prob = 0.05, seed = 1)
  fit <- pfci_fit(sim$X, alpha = 0.05)
  met <- pfci_metrics(sim, fit, compute_marks = TRUE)
  expect_true(all(c("F1_dir", "F1_bidir",
    "F1_arrow", "F1_tail") %in% names(met)))
})

test_that("simulate_with_latent returns correct structure", {
  skip_if_not_installed("pcalg")
  sim <- simulate_with_latent(p_obs = 20, gamma = 0.05,
    n = 50, seed_graph = 1, seed_data = 2)
  expect_true(is.list(sim))
  expect_true(is.matrix(sim$X))
  expect_equal(ncol(sim$X), 20)
  expect_true(is.matrix(sim$truth$skel))
})

test_that("metrics_with_latent returns correct fields", {
  skip_if_not_installed("pcalg")
  sim <- simulate_with_latent(p_obs = 20, gamma = 0.05,
    n = 50, seed_graph = 1, seed_data = 2)
  fit <- pfci_fit(sim$X, alpha = 0.05)
  met <- metrics_with_latent(sim, fit)
  expect_true(is.list(met))
  expect_true(all(c("SHD", "F1_total", "MCC", "Time") %in% names(met)))
})

test_that("pfci_fit rejects bad input", {
  expect_error(pfci_fit("not a matrix"))
})
