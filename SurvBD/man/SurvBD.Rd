\name{SurvBD}

\alias{SurvBD}

\title{Two-Sample Survival Ball Test for Survival Analysis}

\description{
Performs the Ball Divergence to compare the distribution of survival times in two samples of right-censored data.
}

\usage{
SurvBD(time, delta, group, method = 'IPW', accelerate = TRUE, perm = 499)
}

\arguments{
  \item{time}{the follow up time.}
  \item{delta}{the status indicator is 0=right censored, 1=event.}
  \item{group}{a vector indicating to which group each observation belongs, may contain two values only.}
  \item{method}{the weight of the statistic. Possible values are "IPW" for Inverse Probability Weight and "KMW" for Kaplan-Meier Weight. Default: \code{method = 'IPW'}.}
  \item{accelerate}{whether to use acceleration technology. When the sample size is large, this can help to get the result more quickly. Default: \code{accelerate = TRUE}.}
  \item{perm}{the number of permutation replications. Default: \code{perm = 499}.}
}


\value{
  \item{obs.sta}{the observed statistic.}
  \item{pvalue}{the \eqn{p}-value calculated by permutation.}
}

\references{
Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. \emph{ Ball Divergence: Nonparametric two sample test. Annals of statistics}. 46 (2018), no. 3, 1109–1137. doi:10.1214/17-AOS1579. https://projecteuclid.org/euclid.aos/1525313077
}

\author{Na You, Xueyi He}


\examples{
# generate the data from weibull model
Time1 <- (1/9)*log(1- (9/0.6)*log(runif(100)))
Time2 <- (1/2)*log(1- (2/2)*log(runif(100)))
Censor1 <- runif(100, 0, 1)
Censor2 <- runif(100, 0, 1)
delta1 = ifelse(Censor1 <= Time1, 0, 1)
delta2 = ifelse(Censor2 <= Time2, 0, 1)
data <- data.frame(time = c(Time1, Time2),
                  delta = c(delta1, delta2),
                  group = c(rep(1, 100), rep(2, 100)))
library(survival)
library(survminer)
# plot
fit <- survfit(Surv(time, delta) ~ group, data=data)
p <- ggsurvplot(fit, data = data, linetype = "strata")
print(p)
# calculate the observed statistic and p-value
SurvBD(data$time, data$delta, data$group,
      method = 'IPW', accelerate = TRUE, perm = 199)

# real data
library(KMsurv)
data(alloauto)
SurvBD(alloauto$time, alloauto$delta, alloauto$type,
       method = 'IPW', accelerate = TRUE, perm = 199)
SurvBD(alloauto$time, alloauto$delta, alloauto$type,
       method = 'IPW', accelerate = FALSE, perm = 199)
SurvBD(alloauto$time, alloauto$delta, alloauto$type,
       method = 'KMW', perm = 199)
}
\keyword{survival, Ball Divergence}

