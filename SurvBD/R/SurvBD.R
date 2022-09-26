#' Title Survival Ball
#'
#' @param time
#' @param delta
#' @param group
#' @param method
#' @param accelerate
#' @param perm
#'
#' @return statistic and p-value for SurvBD
#' @import survival
#' @export

SurvBD <- function(time, delta, group, method = 'IPW', accelerate = TRUE, perm = 499){

  observed.sta <- SurvBD_sta(time, delta, group, method = method, accelerate = accelerate)

  perm.sta <- sapply(1:perm, function(x){
    group.b = sample(group)
    statistic <- SurvBD_sta(time, delta, group.b, method = method, accelerate = accelerate)
    statistic
  })

  pvalue <- (sum(perm.sta >= observed.sta) + 1) /(perm + 1)

  return(c('obs.sta' = observed.sta, 'pvalue' = pvalue))
}


SurvBD_sta = function(time, delta, group, method = 'IPW', accelerate = TRUE) {
  if (length(time)!=length(delta) | length(time)!=length(group)) {
    stop("lengths of 'time', 'delta' and 'group' differ")
  }
  grp <- names(table(group))
  if(length(grp) != 2){
    stop("'group' requires only 2 categories.")
  }
  time1 <- time[group == grp[1]]
  time2 <- time[group == grp[2]]
  delta1 <- delta[group == grp[1]]
  delta2 <- delta[group == grp[2]]
  n1 <- length(time1)
  n2 <- length(time2)
  #
  delta1 = delta1[order(time1)]
  time1 = sort(time1)
  delta2 = delta2[order(time2)]
  time2 = sort(time2)
  #
  fit1 <- survival::survfit(survival::Surv(time1,1-delta1)~1, type='kaplan-meier')
  fit2 <- survival::survfit(survival::Surv(time2,1-delta2)~1, type='kaplan-meier')
  sc1 = sc_mapping(fit1$time, fit1$surv, time1, length(fit1$time), length(time1))
  sc2 = sc_mapping(fit2$time, fit2$surv, time2, length(fit2$time), length(time2))
  #
  if(method == 'IPW'){
    if(accelerate == FALSE){
      sc1AtSamp2 = sc_mapping(fit1$time, fit1$surv, time2, length(fit1$time), length(time2))
      sc2AtSamp1 = sc_mapping(fit2$time, fit2$surv, time1, length(fit2$time), length(time1))

      SurvBD_IPW(time1 = time1, time2 = time2,
                 delta1 = delta1, delta2 = delta2,
                 sc1 = sc1, sc2 = sc2,
                 sc1AtSamp2 = sc1AtSamp2, sc2AtSamp1 = sc2AtSamp1,
                 n1 = n1, n2 = n2)
    } else{
      # distance + rank-based --> accelerate from n^3 to n^2log(n)
      # (2*X_i-X_j)_(m*m), (2*Y_i-Y_j)_(n*n)
      XX = outer(2*time1, time1, FUN = "-")
      YY = outer(2*time2, time2, FUN = "-")
      X = as.numeric(XX)
      Y = as.numeric(YY)
      # rank
      rank_X <- rank(time1, ties.method = 'max')
      rank_Y <- rank(time2, ties.method = 'max')
      rank_XY <- rank(c(time1, time2), ties.method = 'max')

      rank_XX <- rank(X, ties.method = 'max')  # 2*X_i-X_j
      rank_YY <- rank(Y, ties.method = 'max') # 2*Y_i-Y_j

      rank_1_1 <- rank(c(time1, X), ties.method = 'max')
      rank_1_2 <- rank(c(time1, Y), ties.method = 'max')
      rank_2_1 <- rank(c(time2, X), ties.method = 'max')
      rank_2_2 <- rank(c(time2, Y), ties.method = 'max')

      SurvBD_IPW_accelerated(time1 = time1, time2 = time2,
                             delta1 = delta1, delta2 = delta2,
                             sc1 = sc1, sc2 = sc2,
                             rank_X = rank_X, rank_Y = rank_Y,
                             rank_XY = rank_XY,
                             rank_XX = rank_XX, rank_YY = rank_YY,
                             rank_1_1 = rank_1_1, rank_1_2 = rank_1_2,
                             rank_2_1 = rank_2_1, rank_2_2 = rank_2_2,
                             n1 = n1, n2 = n2)
    }
  } else{
    bracket_x = ((n1 - 1:n1) / (n1 - 1:n1 + 1))^delta1
    cumprod_x = c(1, cumprod(bracket_x))
    w_x = delta1 / (n1 - 1:n1 + 1) * cumprod_x[1:n1]
    w_x = w_x / sum(w_x)
    bracket_y = ((n2 - 1:n2) / (n2 - 1:n2 + 1))^delta2
    cumprod_y = c(1, cumprod(bracket_y))
    w_y = delta2 / (n2 - 1:n2 + 1) * cumprod_y[1:n2]
    w_y = w_y / sum(w_y)
    #
    SurvBD_KMW(time1 = time1, time2 = time2,
               delta1 = delta1, delta2 = delta2,
               weight1 = w_x, weight2 = w_y,
               n1 = n1, n2 = n2)
  }


}

