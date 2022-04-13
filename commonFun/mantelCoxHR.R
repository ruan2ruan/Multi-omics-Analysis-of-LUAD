mantelCoxHR <- function(time, event, group, conf.level = 0.95){
  
  ## Compute HR according to Mantel-Cox formula
  ## See Kirkwood & Sterne, p. 283.
  ##
  ## Kaspar Rufibach, November 2008
  
  alpha <- 1 - conf.level
  group <- as.factor(group)
  s.obj <- survfit(Surv(time, event) ~ group, conf.type = "none") 
  s <- summary(s.obj)
  dat <- data.frame(s$time, s$strata, s$n.risk, s$n.event) 
  dat <- dat[order(s$time), ] 
  lev <- levels(dat$s.strata)
  res <- matrix(NA, ncol = 4, nrow = length(dat[, 1]))
  
  for (i in 1:length(dat[, 1])){
    tmp <- dat[i, ]
    if (tmp[2] == lev[1]){res[i, 1:2] <- c(tmp[[3]], tmp[[4]])}
    if (tmp[2] == lev[2]){res[i, 3:4] <- c(tmp[[3]], tmp[[4]])}
  }
  
  ## res2 is the table on p. 285 in Kirkwood & Sterne
  res2 <- res 
  n2 <- length(res2[, 1])
  for (i in 1:n2){
    g1 <- res2[i, 1:2]
    if (is.na(g1[1]) == TRUE){res2[i, 1:2] <- c(sum(time[group == levels(group)[1]] >= dat$s.time[i], na.rm = TRUE), 0)}
    
    g2 <- res2[i, 3:4]
    if (is.na(g2[1]) == TRUE){res2[i, 3:4] <- c(sum(time[group == levels(group)[2]] >= dat$s.time[i], na.rm = TRUE), 0)}
  } 
  dimnames(res2)[[2]] <- c("n0i", "d0i", "n1i", "d1i")
  
  ## compute Mantel-Cox HR
  h0i <- res2[, "n0i"] - res2[, "d0i"] 
  h1i <- res2[, "n1i"] - res2[, "d1i"]
  ni <- res2[, "n0i"] + res2[, "n1i"]
  di <- res2[, "d0i"] + res2[, "d1i"]
  hi <- h0i + h1i
  Q <- sum(res2[, "d1i"] * h0i / ni)
  R <- sum(res2[, "d0i"] * h1i / ni)
  U <- sum(res2[, "d1i"] - di * res2[, "n1i"] / ni)
  V <- sum(di * hi * res2[, "n1i"] * res2[, "n0i"] / (ni ^ 2 * (ni - 1)), na.rm = TRUE)
  hr <- Q / R
  
  ## compute confidence interval
  se.log.hr <- sqrt(V / (Q * R))
  ci.hr <- exp(log(hr) + qnorm(c(alpha / 2, 1 - alpha / 2)) * se.log.hr)
  
  ## p-value logrank test
  chi2 <- U ^ 2 / V
  p.val <- pchisq(chi2, df = 1, lower.tail = FALSE)
  
  ## compare to HR computed from Cox-regression 
  obj2 <- Surv(time, event) 
  fit2 <- coxph(obj2 ~ group)
  cfit2 <- as.numeric(coef(fit2))
  hr.cox <- exp(cfit2)
  
  return(list("mantelCox.hr" = hr, "ci.hr" = ci.hr, "p.val.logrank" = p.val, "coxph.hr" = hr.cox))
}