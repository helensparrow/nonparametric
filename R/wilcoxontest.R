
wilcoxontest <- function() {
  myWilcoxon <- function(X, Y, alpha = 0.05, alternative = "two.sided", digits = 4) {
    alpha <- ifelse(alternative == "two.sided", alpha/2, alpha)

    # Wilcoxon Two-Sample Test
    m <- length(X)
    n <- length(Y)

    # test statistics
    W.pre <- data.frame(Values = c(X, Y),
                        Data = c(rep("X",length(X)), rep("Y",length(Y))),
                        Ranks = rank(c(X, Y)))
    W = sum(W.pre$Ranks[W.pre$Data =="Y"])
    U <- sum(outer(X, Y, FUN = "<"))

    # Exact
    # critical values
    walpha <- cv.upper <- qwilcox(1-alpha, m, n) + 1
    cv.lower <- n*m - walpha
    cv.wilcox <- c(ifelse(alternative == "greater", NA, cv.lower),
                   ifelse(alternative == "less", NA, cv.upper))

    # Asymptotic
    # critical values
    eu <- m * n / 2
    varu <- m * n * (m + n + 1) / 12

    cv.wilcox.a <- qnorm(c(ifelse(alternative == "greater", NA, alpha),
                           ifelse(alternative == "less", NA, 1 - alpha)),
                         eu, sqrt(varu))


    if(alternative == "two.sided") {
      actAlpha <- pwilcox(cv.lower, m, n) + pwilcox(cv.upper-1, m, n, lower.tail=F) # P(X<= cv.l) + P(X>= cv.u)
      p.val <- ifelse(U < m*n/2, 2*pwilcox(U, m, n), 2*pwilcox(U-1, m, n, lower.tail=F))

      actAlpha.a <- alpha * 2
      p.val.a <- 2*pnorm(U, eu, sqrt(varu))
    }
    if(alternative == "greater") {
      actAlpha <- pwilcox(cv.upper-1, m, n, lower.tail = FALSE) # P(X>= cv.u)
      p.val <- pwilcox(U-1, m, n, lower.tail = FALSE)

      actAlpha.a <- alpha
      p.val.a <- pnorm(U, eu, sqrt(varu), lower.tail = FALSE)
    }
    if(alternative == "less") {
      actAlpha <- pwilcox(cv.lower, m, n) # P(X<= cv.l)
      p.val <- pwilcox(U, m, n)

      actAlpha.a <- alpha
      p.val.a <- pnorm(U, eu, sqrt(varu))
    }


    # Estimate delta
    ## HL Estimator of delta
    D <- sort(c(outer(Y, X, FUN = "-")))
    delta.hat <- median(D)

    # confidence intervals (Moses)
    # Moses CI
    w.alpha <- qwilcox(1-alpha, m, n) + 1
    c.alpha <- n * m + 1 - w.alpha
    CI.moses <- D[c(ifelse(alternative == "less", NA, c.alpha),
                    ifelse(alternative == "greater", NA, w.alpha))]

    # Asymptotic Version of Moses
    c.alpha.a <- floor(eu - qnorm(1-alpha) * sqrt(varu))
    w.alpha.a <- n * m + 1 - c.alpha.a
    CI.moses.a <- D[c(ifelse(alternative == "less", NA, c.alpha.a),
                      ifelse(alternative == "greater", NA, w.alpha.a))]


    results <- data.frame('Exact' = c(W, U, cv.wilcox, p.val, actAlpha, delta.hat, CI.moses),
                          'Asymptotic' = c(NA, NA, cv.wilcox.a, p.val.a, actAlpha.a, NA, CI.moses.a))

    results <- data.frame('Quantity' = c("Test Statistic: W", "Test Statistic: U", "CV Lower",
                                         "CV Upper", "P-Value", "Actual Alpha",
                                         "HL Delta Estimate", "Moses CI Lower", "Moses CI Upper"),
                          round(results, digits))

    return(kable(results, digits = digits))
  }
}
