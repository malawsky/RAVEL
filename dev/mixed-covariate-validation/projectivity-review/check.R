# Run from the repository root. No fitted data or external files are needed.
for (file in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(file)
output <- "dev/mixed-covariate-validation/projectivity-review"
q <- 2L
triangle <- which(lower.tri(matrix(0, q, q), diag = TRUE), arr.ind = TRUE)
Qself <- matrix(c(1, .1, .1, 1.2), q)
Qsib <- matrix(c(.2, .03, .03, .3), q)
parameter_names <- setNames(lapply(c("mother", "father", "child_self", "child_sibling"),
  function(role) paste0("D_", role, "_", seq_len(q))),
  c("mother", "father", "child_self", "child_sibling"))
cv <- list(q = q, eta = c(0, 0, Qself[triangle], Qsib[triangle]),
  triangle = triangle, index = list(mean = 1:2, self = 3:5, sibling = 6:8),
  has_siblings = TRUE, parameter_names = parameter_names)
par <- c(state_rho = atanh(.3), state_variance = log(.7),
  loading1_mother = log(.7), loading1_father = log(.8), loading1_child = log(.6),
  loading2_mother = log(.5), loading2_father = log(.65), loading2_child = log(.55),
  noise1 = log(.8), noise2 = log(1.1),
  setNames(c(.03, -.04, -.02, .05, .08, -.06, -.01, .02),
    unlist(parameter_names, use.names = FALSE)))
# The free sibling D deliberately differs from half the parental D sum.
stopifnot(any(par[parameter_names$child_sibling] !=
  .5 * (unname(par[parameter_names$mother]) + unname(par[parameter_names$father]))))
records <- list()
for (mode in c("estimated_state", "equilibrium", "n_gen"))
  for (pi in c(0, .4, .8)) for (rho in c(-.2, 0, .3)) {
    spec <- list(mode = mode, pi = pi, generations = 3L, covariate = cv)
    current <- par; current["state_rho"] <- atanh(rho)
    vp <- switch(mode, estimated_state = 1.2, equilibrium = 1 / (1 - rho),
      n_gen = (1 - rho * ((1 + rho) / 2)^2) / (1 - rho))
    vc <- .5 * vp * (1 + rho) + .5
    Tp <- vp * matrix(c(1, rho, rho, 1), 2L)
    Ap <- (1 - pi) * Tp + pi * diag(2L)
    Rp_raw <- sqrt(1 - pi) * Tp %*% solve(Ap)
    one <- .ravel_mixed_covariate_blocks(current, spec, 1L)
    two <- .ravel_mixed_covariate_blocks(current, spec, 2L)
    for (n in 1:6) {
      b <- .ravel_mixed_covariate_blocks(current, spec, n)
      analytic <- matrix(0, n + 2L, n + 2L)
      analytic[1:2, 1:2] <- Rp_raw / sqrt(vp)
      for (i in seq_len(n)) {
        analytic[i + 2L, 1:2] <- (.5 * colSums(Rp_raw) - .5 * sqrt(1 - pi)) / sqrt(vc)
        analytic[i + 2L, i + 2L] <- sqrt(1 - pi) / sqrt(vc)
      }
      extra <- if (n == 1L) 0 else max(abs(b$R[1:3, 4:(n + 2L), drop = FALSE]))
      e <- b$R %*% b$D
      common <- max(abs(e[1:3, seq_len(q), drop = FALSE] - one$R %*% one$D))
      pair <- if (n == 1L) 0 else
        max(abs(e[1:4, seq_len(2L * q), drop = FALSE] - two$R %*% two$D))
      records[[length(records) + 1L]] <- data.frame(mode = mode, pi = pi, rho = rho,
        n_children = n, analytic_R_max_error = max(abs(b$R - analytic)),
        additional_sibling_R_max = extra, common_trio_RD_max_error = common,
        common_two_child_RD_max_error = pair,
        conditional_visible_min_eigenvalue = min(eigen(b$conditional_visible,
          symmetric = TRUE, only.values = TRUE)$values))
    }
  }
results <- do.call(rbind, records)
stopifnot(max(as.matrix(results[c("analytic_R_max_error", "additional_sibling_R_max",
  "common_trio_RD_max_error", "common_two_child_RD_max_error")])) < 1e-12)
utils::write.csv(results, file.path(output, "checks.csv"), row.names = FALSE)
cat(nrow(results), "family/mode/pi/rho checks passed; maximum algebra error",
  max(results$analytic_R_max_error), "\n")
