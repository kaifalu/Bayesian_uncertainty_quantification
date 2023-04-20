## ----libraries, message=FALSE, warning=FALSE, cache=FALSE----------------
library("plyr")
library("dplyr")
library("ggplot2")
library("gridExtra")

## ----set_seed------------------------------------------------------------
set.seed(2)

## ----Problem 2 posterior-----------------------------------------------------
a = b = 1
e = data.frame(year=1:4, y=c(36,64,67,64), n=c(95,150,171,152))
for (i in 1:nrow(e)) {
  e$y_new[i] <- sum(e$y[1:i])
  e$n_new[i] <- sum(e$n[1:i])
}
e$a = a + e$y_new
e$b = b + e$n_new - e$y_new
plot(0, 0, type="n", main="Andre Dawkins's 3-point percentage", xlab=expression(theta), ylab="Posterior",
     xlim=c(0,1), ylim=c(0,max(dbeta(e$y_new/e$n_new,e$a,e$b))))
for (i in 1:nrow(e)) curve(dbeta(x, e$a[i], e$b[i]), col=i, lwd=2, lty=i, add=TRUE)
legend("topright", c("Season 1", "Season 1+2", "Season 1+2+3", "Season 1+2+3+4"), col=1:4, lwd=2, lty=1:4)


## ----Chapter 3 joint_posterior-----------------------------------------------------
a = b = 1
d = data.frame(year=1:4, y=c(36,64,67,64), n=c(95,150,171,152))
d$a = a + d$y
d$b = b + d$n - d$y
plot(0, 0, type="n", main="Andre Dawkins's 3-point percentage", xlab=expression(theta), ylab="Posterior",
     xlim=c(0,1), ylim=c(0,max(dbeta(d$y/d$n,d$a,d$b))))
for (i in 1:nrow(d)) curve(dbeta(x, d$a[i], d$b[i]), col=i, lwd=2, lty=i, add=TRUE)
legend("topright", paste("Season", 1:nrow(d)), col=1:4, lwd=2, lty=1:4)