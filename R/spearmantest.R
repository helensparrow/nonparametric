
spearmantest.exact <- function() {
  Si <- rank(data$dpp)
  Ri <- rank(data$pct)

  rs <- 1 - (6*sum((Si-Ri)^2))/(n*(n^2-1))

  rsalpha <- qSpearman(p = 0.95, r = n)

  # P-Value
  pval <- (1-pSpearman(rs, n))

  results <- data.frame('Quantity' = c("Test Statistic (rs)","Critical Value", "P-Value"),
                        'By-Hand' = c(rs, rsalpha, pval),
                        'Built-In Func' = c(rsa$estimate, NA, rsa$p.value))
}


