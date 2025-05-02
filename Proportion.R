# Confidence Interval and Sample Size for One Proportion v0.0.1
# Author: ZangKail 2025-04-19

# Reference: R source: power.binom.test{MESS}, BinomCI{clinfun}
# Reference: Handbook: JMP, G-power, Medcalc, Analyse-it, Stata, SAS

# Confidence Interval for One Proportion
ci_prop <- function(p, n, alpha = 0.05, sides = "two.sided", method = "wilson") {

  if (length(p) != 1 ||
    !is.numeric(p) ||
    p >= 1 ||
    p <= 0)
    stop("p should be between 0 and 1")
  if (length(n) != 1 || !is.numeric(n) || n < 1)
    stop("n should be an integer greater than 0")
  if (length(alpha) != 1 ||
    !is.numeric(alpha) ||
    alpha >= 0.5 ||
    alpha <= 0)
    stop("alpha should be between 0 and 0.5")
  if (length(method) != 1 ||
    !is.character(method) ||
    !method %in% c("wilson", "wald", "waldcc", "agresti-coull", "jeffreys", "wilsoncc", "clopper-pearson"))
    stop("method should be one of the strings: wilson, wald, waldcc, agresti-coull, jeffreys, wilsoncc, clopper-pearson")
  if (length(sides) != 1 ||
    !is.character(sides) ||
    !sides %in% c("two.sided", "left", "right"))
    stop("sides should be one of the strings: two.sided, left, right")

  if (sides != "two.sided")
    alpha <- 2 * alpha

  n <- floor(n)
  x <- round(p * n)
  q <- 1 - p
  Z <- qnorm(1 - alpha / 2)

  ci.2sided <- switch(method,
                      "wald" = {
                        term <- Z * sqrt(p * q / n)
                        CI.lower <- p - term
                        CI.upper <- p + term
                        c(CI.lower, CI.upper)
                      },
                      "waldcc" = {
                        term <- Z * sqrt(p * q / n) + 1 / (2 * n)
                        CI.lower <- p - term
                        CI.upper <- p + term
                        c(CI.lower, CI.upper)
                      },
                      "wilson" = {
                        term1 <- (x + Z^2 / 2) / (n + Z^2)
                        term2 <- Z * sqrt(n) / (n + Z^2) * sqrt(p * q + Z^2 / (4 * n))
                        CI.lower <- term1 - term2
                        CI.upper <- term1 + term2
                        c(CI.lower, CI.upper)
                      },
                      "wilsoncc" = {
                        CI.lower <- (2 * x + Z^2 -
                          1 -
                          Z * sqrt(Z^2 - 2 - 1 / n + 4 * p * (n * q + 1))) / (2 * (n + Z^2))
                        CI.upper <- (2 * x +
                          Z^2 +
                          1 +
                          Z * sqrt(Z^2 + 2 - 1 / n + 4 * p * (n * q - 1))) / (2 * (n + Z^2))
                        c(CI.lower, CI.upper)
                      },
                      "agresti-coull" = {
                        x.tilde <- x + Z^2 / 2
                        n.tilde <- n + Z^2
                        p.tilde <- x.tilde / n.tilde
                        q.tilde <- 1 - p.tilde
                        term <- Z * sqrt(p.tilde * q.tilde) / sqrt(n.tilde)
                        CI.lower <- p.tilde - term
                        CI.upper <- p.tilde + term
                        c(CI.lower, CI.upper)
                      },
                      "jeffreys" = {
                        CI.lower <- qbeta(alpha / 2, x + 0.5, n - x + 0.5)
                        CI.upper <- qbeta(1 - alpha / 2, x + 0.5, n - x + 0.5)
                        c(CI.lower, CI.upper)
                      },
                      "clopper-pearson" = {
                        CI.lower <- qbeta(alpha / 2, x, n - x + 1)
                        CI.upper <- qbeta(1 - alpha / 2, x + 1, n - x)
                        c(CI.lower, CI.upper)
                      }
  )

  if (sides == "left") {
    ci <- max(0, ci.2sided[1])
  }
  else if (sides == "right") {
    ci <- min(1, ci.2sided[2])
  }
  else {
    ci <- c(max(0, ci.2sided[1]), min(1, ci.2sided[2]))
  }

  cat("\n-----------------------------------------------------",
      "\nConfidence Interval for One Proportion",
      "\n\nParameters:",
      "\n  proportion = ", p,
      "\n  sample size = ", n,
      "\n  alpha = ", alpha,
      "\n  sides = ", sides,
      "\n  method = ", method,
      "\n\nEstimated:",
      "\n  confidence interval = ", round(ci, 4),
      "\n-----------------------------------------------------\n"
  )

  return(invisible(ci))

}

# ci_prop(0.95, 300, alpha = 0.05, sides = "two.sided", method = "wilson")

# Sample Size: Confidence Interval for One Proportion
size_prop_ci <- function(p, width, alpha = 0.05, sides = "two.sided", method = "wilson", interval = c(1, 100000)) {

  if (length(width) != 1 ||
    !is.numeric(width) ||
    width >= 1 ||
    width <= 0) stop("width should be between 0 and 1")
  if (length(interval) != 2 ||
    !is.numeric(interval) ||
    interval[1] < 1 ||
    interval[1] > interval[2])
    stop("interval should be integer and left <= right")

  f <- function(n) {
    if (sides == "two.sided") { abs(diff(ci_prop(p = p, n = n, alpha = alpha, sides = sides, method = method))) - width }
    else { abs(p - ci_prop(p = p, n = n, alpha = alpha, sides = sides, method = method)) - width }
  }

  capture.output(n <- try(ceiling(uniroot(f, interval = interval)$root), silent = TRUE))
  if (inherits(n, "try-error")) stop("failed to find solution in the given interval")

  cat("\n-----------------------------------------------------",
      "\nSample Size: Confidence Interval for One Proportion",
      "\n\nParameters:",
      "\n  proportion = ", p,
      "\n  width = ", width,
      "\n  alpha = ", alpha,
      "\n  sides = ", sides,
      "\n  method = ", method,
      "\n\nEstimated:",
      "\n  sample size = ", n,
      "\n-----------------------------------------------------\n"
  )

  return(invisible(n))

}

# size_prop_ci(0.95, 0.05, alpha = 0.05, sides = "two.sided", method = "wilson")

# Power for One Proportion
power_prop <- function(p0, pa, n, alpha = 0.05, method = "wald", alternative = "greater") {

  if (length(p0) != 1 ||
    !is.numeric(p0) ||
    p0 >= 1 ||
    p0 <= 0)
    stop("p0 should be between 0 and 1")
  if (length(pa) != 1 ||
    !is.numeric(pa) ||
    pa >= 1 ||
    pa <= 0)
    stop("pa should be between 0 and 1")
  if (length(n) != 1 || !is.numeric(n) || n < 1)
    stop("n should be an integer greater than 0")
  if (length(alpha) != 1 ||
    !is.numeric(alpha) ||
    alpha >= 0.5 ||
    alpha <= 0)
    stop("alpha should be between 0 and 0.5")
  if (length(method) != 1 ||
    !is.character(method) ||
    !method %in% c("wilson", "wald", "waldcc", "wilsoncc", "clopper-pearson"))
    stop("method should be one of the strings: wilson, wald, waldcc, wilsoncc, clopper-pearson")
  if (length(alternative) != 1 ||
    !is.character(alternative) ||
    !alternative %in% c("less", "greater", "two.sided"))
    stop("sides should be one of the strings: two.sided, less, greater")
  if (alternative == "greater" && pa <= p0)
    stop("pa should be greater than p0 for 'greater' alternative")
  if (alternative == "less" && pa >= p0)
    stop("pa should be less than p0 for 'less' alternative")

  n <- floor(n)

  power <- switch(method,
                  "clopper-pearson" = {
                    switch(alternative,
                           "less" = {
                             pbinom(qbinom(1 - alpha, size = n, prob = p0, lower.tail = FALSE) - 1, size = n, prob = pa)
                           },
                           "greater" = {
                             pbinom(qbinom(1 - alpha, size = n, prob = p0), size = n, prob = pa, lower.tail = FALSE)
                           },
                           "two.sided" = {
                             lx <- qbinom(alpha / 2, size = n, prob = p0)
                             ux <- qbinom(1 - alpha / 2, size = n, prob = p0)
                             pbinom(lx, size = n, prob = pa) + pbinom(ux, size = n, prob = pa, lower.tail = FALSE)
                           }
                    )
                  },
                  "wald" = {
                    x <- sqrt(n) * (pa - p0)
                    y <- sqrt(pa - pa^2)
                    switch(alternative,
                           "less" = { pnorm(-x / y - qnorm(1 - alpha)) },
                           "greater" = { pnorm(x / y - qnorm(1 - alpha)) },
                           "two.sided" = { pnorm(x / y - qnorm(1 - alpha / 2)) + pnorm(-x / y - qnorm(1 - alpha / 2))
                           }
                    )
                  },
                  "waldcc" = {
                    x <- sqrt(n) * (pa - p0) - 1 / (2 * sqrt(n))
                    c <- 1 / (2 * sqrt(n))
                    y <- sqrt(pa - pa^2)
                    switch(alternative,
                           "less" = { pnorm((-x - c) / y - qnorm(1 - alpha)) },
                           "greater" = { pnorm((x - c) / y - qnorm(1 - alpha)) },
                           "two.sided" = {
                             pnorm((-x - c) / y - qnorm(1 - alpha / 2)) + pnorm((x - c) / y - qnorm(1 - alpha / 2))
                           }
                    )
                  },
                  "wilson" = {
                    x <- sqrt(n) * (pa - p0)
                    y <- sqrt(pa - pa^2)
                    z <- sqrt((p0 - p0^2) / (pa - pa^2))
                    switch(alternative,
                           "less" = { pnorm(-x / y - z * qnorm(1 - alpha)) },
                           "greater" = { pnorm(x / y - z * qnorm(1 - alpha)) },
                           "two.sided" = {
                             pnorm(x / y - z * qnorm(1 - alpha / 2)) + pnorm(-x / y - z * qnorm(1 - alpha / 2))
                           }
                    )
                  },
                  "wilsoncc" = {
                    x <- sqrt(n) * (pa - p0) - 1 / (2 * sqrt(n))
                    c <- 1 / (2 * sqrt(n))
                    y <- sqrt(pa - pa^2)
                    z <- sqrt((p0 - p0^2) / (pa - pa^2))
                    switch(alternative,
                           "less" = { pnorm((-x - c) / y - z * qnorm(1 - alpha)) },
                           "greater" = { pnorm((x - c) / y - z * qnorm(1 - alpha)) },
                           "two.sided" = {
                             pnorm((-x - c) / y - z * qnorm(1 - alpha / 2)) +
                               pnorm((x - c) / y - z * qnorm(1 - alpha / 2))
                           }
                    )
                  }
  )

  cat("\n-----------------------------------------------------",
      "\nPower for One Proportion",
      "\n\nParameters:",
      "\n  unacceptable value p0 = ", p0,
      "\n  expected value pa = ", pa,
      "\n  sample size = ", n,
      "\n  alpha = ", alpha,
      "\n  alternative = ", alternative,
      "\n  method = ", method,
      "\n\nEstimated:",
      "\n  power = ", round(power, 4),
      "\n-----------------------------------------------------\n"
  )

  return(invisible(power))

}

# power_prop(0.9, 0.95, 300)

# Sample Size for One Propotion
size_prop <- function(p0, pa, alpha = 0.05, beta = 0.02, alternative = "greater", method = "wald",
                      interval = c(1, 100000)) {

  if (length(beta) != 1 ||
    !is.numeric(beta) ||
    beta >= 0.5 ||
    beta <= 0)
    stop("beta should be between 0 and 0.5")
  if (length(interval) != 2 ||
    !is.numeric(interval) ||
    interval[1] < 1 ||
    interval[1] > interval[2])
    stop("interval should be integer and left <= right")

  f <- function(n)
    power_prop(p0 = p0, pa = pa, alpha = alpha, n = n, method = method, alternative = alternative) - (1 - beta)

  capture.output(n <- try(ceiling(uniroot(f, interval = interval)$root), silent = TRUE))
  if (inherits(n, "try-error")) stop("failed to find solution in the given interval")

  cat("\n-----------------------------------------------------",
      "\nSample Size for One Proportion",
      "\n\nParameters:",
      "\n  unacceptable value p0 = ", p0,
      "\n  expected value pa = ", pa,
      "\n  alpha = ", alpha,
      "\n  beta = ", beta,
      "\n  alternative = ", alternative,
      "\n  method = ", method,
      "\n\nEstimated:",
      "\n  sample size = ", n,
      "\n-----------------------------------------------------\n"
  )

  return(invisible(n))

}

# size_prop(0.65, 0.8, alpha = 0.1, beta = 0.1, alternative = "greater", method = "clopper-pearson")
