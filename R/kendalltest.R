
kendalltest.exact <- function() {
  n <- nrow(data)

  Q <- function(pi, pj) { ifelse((pj[1]-pi[1])*(pj[2]-pi[2]) < 0, -1, 1) }
  Qstar<- function(pi, pj) { ifelse((pj[1]-pi[1])*(pj[2]-pi[2]) < 0, -1,
                                    ifelse((pj[1]-pi[1])*(pj[2]-pi[2]) > 0, 1, 0)) }

  Qij <- numeric(length = n*(n-1)/2)
  counter <- 1
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      Qij[counter] <- Qstar(data[i,], data[j,])
      counter <- counter + 1
    }
  }

  K <- sum(Qij)
  Kbar <- K/(n*(n-1)/2)

  # Critical Value
  kalpha <- qKendall(p = 0.95, N = n)

  # P-Value
  p.val <- 1 - pKendall(Kbar, N = n)

  results <- data.frame('Quantity' = c("Test Statistic (K)", "Test Statistic (KBar)","Critical Value",
                                       "P-Value"),
                        'By-Hand' = c(K, Kbar, kalpha, p.val),
                        'Built-In Func' = c(NA, tt$estimate, NA, tt$p.value))
}



kendalltest.asymptotic <- function() {
  EK <- 0
  varK <- (n*(n-1)*(2*n+5))/18

  kalpha.a <- qnorm(0.95, EK, sqrt(varK))

  p.val <- 1-pnorm(K, EK, sqrt(varK))
  # Reject HO

  results <- data.frame('Quantity' = c("Test Statistic (K)", "Critical Value", "P-Value"),
                        'Asymptotic' = c(K, kalpha.a, p.val))
}


kendalltest.tauCI <- function() {
  tau.hat <- Kbar

  Ci <- numeric(length = n)
  counter <- 1
  for(i in 1:n) {
    Qi <- 0
    for(t in 1:n) {
      if(t == i) next
      Qi <- Qi + Q(data[i,], data[t,])
    }
    Ci[counter] <- Qi
    counter <- counter + 1
  }

  Cbar <- mean(Ci)
  sig2hat <- (2/(n*(n-1)))*(((2*(n-2))/(n*(n-1)^2))*sum((Ci-Cbar)^2)+1-tau.hat^2)

  # two tail test
  CI.l <- tau.hat - qnorm(1-0.025)*sqrt(sig2hat)

  # Built-in
  library(NSM3)
  kci <- kendall.ci(x = data$dpp, y = data$pct, alpha = 0.05, type = "l")

  results <- data.frame('Quantity' = c("Test Statistic (Kbar)","CI Lower", "CI Upper"),
                        'By-Hand' = c(Kbar, CI.l, 1),
                        'Built-In Func' = c(NA, -0.035, 1))
}

