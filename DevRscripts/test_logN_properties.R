

res <- matrix(data = NA, nrow = 100, ncol = 10)
for(repeats in 1:100)
{
  for(s in 1:10)
  {
    X <- c(rlnorm(n = 10000, meanlog = -(s^2)/2, s), rlnorm(n = 10000, meanlog = -((2*s)^2)/2, 2*s))
    res[repeats, s] <- mean(X)
  }
}

layout(mat = matrix(1:10, ncol=2))
for(i in 1:10)
{
  xlim <- range(log10(res[,i]))
  xlim[2] <- ifelse(max(xlim) < 0, 0, max(xlim))
  xlim[1] <- ifelse(min(xlim) > 0, 0, min(xlim))
  
  hist(log10(res[,i]), xlab=expression("log"[10]~"(arth. mean)"), main=i, nclass = 100, xlim=xlim)
  abline(v = log10(mean(res[,i])), col="red")
  abline(v = log10(median(res[,i])), col="red", lty=2)
  q <- log10(quantile(res[,i], probs = c(0.025, 0.975)))
  abline(v = q, col="blue")
}
