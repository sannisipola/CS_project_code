# Computational Statistics, presentation code


data <- read.delim("C:\\R\\table1.txt", sep="")
data

transform_grouped_data_to_non_grouped <- function(data) {
  # input data must be 4 columns: A, B, cases and controls
  a = b = response = numeric()
  for(i in 1:nrow(data)){
    response <- c( rep(1, data[i,3]) , rep(0, data[i,4]) )
    a <- rep(data[i,1], length(response)) # A of this iteration
    b <- rep(data[i,2], length(response)) # B of this iteration
    A <- append(A, a)
    B <- append(B, b)
    RESPONSE <- append(RESPONSE, response)
  }
  return(cbind(A,B,RESPONSE))
}

oral_cancer <- transform_grouped_data_to_non_grouped(data[,c(1,2,7,8)])

oral_cancer <- list(  y11=c(rep(1, data[1,7]), rep(0,data[1,8])), 
                      y10=c(rep(1, data[2,7]), rep(0, data[2,8])), 
                      y01=c(rep(1, data[3,7]), rep(0, data[3,8])), 
                      y00=c(rep(1, data[4,7]), rep(0, data[4,8])) 
                     ) 
                       
rep(1,data[1,7])
c( rep(1, data[1,7]) , rep(0,data[1,8]) )
#---------------------------------------------------------------
# Table 4 / Simulation 2

# setting to make mimic the oral cancer example:
n00 <- 23
n01 <- 26
n10 <- 18
n11 <- 391
p00 <- 0.130
p01 <- 0.308
p10 <- 0.333
p11 <- 0.575

trueRERI <- 3.73 # given in paper

(p11/(1-p11)) * ((1-p00)/p00) - (p10/(1-p10))*((1-p00)/p00) - (p01/(1-p01))*((1-p00)/p00) +1 # = 3.734504

#------------------------------------------------------------------
# 1. Parametric bootstrap with continuity correction
#------------------------------------------------------------------

# I could not get the same results as in the paper, I tried three following formulas for RERI:


# RERI based on relative risk, with continuity correction (page 554)
calc_reri_rr_based_cc <- function(y00, y01, y10, y11) {
  return ((y11/(n11 + 0.5)) * (n00/(y00 + 0.5)) 
        - (y10/(n10 + 0.5)) * (n00/(y00 + 0.5)) 
        - (y01/(n01 +  1 )) * (n00/(y00 + 0.5))
          +1)
}       

# Equation (2) (page 553)
calc_reri_equation2 <- function(y00, y01, y10, y11) {
  return ((y11/(n11 - y11)) * ((n00 - y00)/y00) 
        - (y10/(n10 - y10)) * ((n00 - y00)/y00) * (y01/(n01 - y01)) * ((n00 - y00)/y00) 
          + 1)
}       

# Equation (3) (page 553)
calc_reri_equation3 <- function(y00, y01, y10, y11) {
  return ((y11/(n11-y11 + .5)) * ((n00-y00)/(y00 + .5)) 
        - (y10/(n10-y10 + .5)) * ((n00-y00)/(y00 + .5))
        - (y01/(n01-y01 + .5)) * ((n00-y00)/(y00 + .5)) 
        + 1)
}       
# Equation (2) trial
newfunction <- function(y00, y01, y10, y11) {
  return ((y11/(n11 - y11)) * ((n00 - y00)/y00) 
        - (y10/(n10 - y10)) * ((n00 - y00)/y00) * (y01/(n01 - y01)) * ((n00 - y00)/y00) 
          + 1)
}
#----------------------------------------------------------------
# Monte carlo simulation
SIM <- 500

set.seed(1)
MC_CIs <- matrix(nrow=SIM, ncol=2) # col 1 is for lower bounds and col 2 is for upper bounds
I <- numeric(SIM)
RERI_simulations <- numeric(SIM)
for(i in 1:SIM){
  # parametric bootstrap samples with continuity correction
  k <- 1000 # sample size, 1000 in the paper
  y00 <- rbinom(k, n00, p00)
  y01 <- rbinom(k, n01, p01)
  y10 <- rbinom(k, n10, p10)
  y11 <- rbinom(k, n11, p11)

  # k bootstrap estimates for RERI
  RERI_bootstraps <- calc_reri_equation3(y00, y01, y10, y11) # plug in here function for RERI
  RERI_simulations[i] <- mean(RERI_bootstraps)
  # CI calculated from k bootstrap samples
  CI <- quantile(RERI_bootstraps, probs=c(0.025,0.975), type=4)
  
  I[i] <- CI[[1]] <= trueRERI & trueRERI <= CI[[2]] # indicator of coverage rate
  MC_CIs[i,] <- c(CI[[1]],CI[[2]]) # bootstrap CI to each row of the matrix
}

MC_RERI_estimate <- mean(RERI_simulations)
MC_coverage_rate <- sum(I)/SIM
MC_CI <- c(mean(MC_CIs[,1]) , mean(MC_CIs[,2]))
cat("mean RERI: ",MC_RERI_estimate, "MC_CI: ", MC_CI, "   MC-coverage_rate: ", MC_coverage_rate)


# true RERI given in paper 3.73,  CI (-7.25 , 20.77)
# 1. gives: mean RERI:  0.567687   MC_CI:  -3.134205 3.280999    MC-coverage_rate:  0.03
# 2. gives: mean RERI:  NaN        MC_CI:  NA NA                 MC-coverage_rate:  0
# 3. gives: mean RERI:  4.195587   MC_CI:  -2.359156 16.97303    MC-coverage_rate:  1 ----- SAME AS IN TABLE 2 !



#------------------------------------------------------------------
# 2. Stratified non-parametric bootstrap with continuity correction
#------------------------------------------------------------------

oral_cancer <- list(  Y11=c(rep(1, data[1,7]), rep(0,data[1,8])), 
                      Y10=c(rep(1, data[2,7]), rep(0, data[2,8])), 
                      Y01=c(rep(1, data[3,7]), rep(0, data[3,8])), 
                      Y00=c(rep(1, data[4,7]), rep(0, data[4,8])) 
) 

n00 <- length(oral_cancer$Y00)
n01 <- length(oral_cancer$Y01)
n10 <- length(oral_cancer$Y10)
n11 <- length(oral_cancer$Y11)
y_function <- function(sample){sum(sample)}

k <- 1000 # sample size, 1000 in the paper
y00 <- bootstrap(oral_cancer$Y00, k, y_function)$thetastar
y01 <- bootstrap(oral_cancer$Y01, k, y_function)$thetastar
y10 <- bootstrap(oral_cancer$Y10, k, y_function)$thetastar
y11 <- bootstrap(oral_cancer$Y11, k, y_function)$thetastar

k_bootstrap_RERIs <- calc_reri_equation3(y00, y01, y10, y11)
RERI_bootstraps <- k_bootstrap_RERIs

# mean RERI:  4.321775 MC_CI:  -1.779268 16.44699    MC-coverage_rate:  1
# paper: -2.22 to 16.74 

#------------------------------------------------------------------
# Non parametric bootstrap with continuity correction
#------------------------------------------------------------------

n_all <- n00+n01+n10+n11
attach(oral_cancer)
Y_all <- c(Y00, Y01, Y10, Y11)
y_all <- bootstrap(Y_all, k, y_function)$thetastar

# how do I split my y_all to y00,y01,y10,y11 ????????????????????????????

k_bootstrap_RERIs <- calc_reri_equation3(y00, y01, y10, y11)
RERI_bootstraps <- k_bootstrap_RERIs

#-------------------------------------------------------------------

data

calc_reri_equation3 <- function(y00, y01, y10, y11) {
    return ((y11/(n11-y11 + .5)) * ((n00-y00)/(y00 + .5)) 
          - (y10/(n10-y10 + .5)) * ((n00-y00)/(y00 + .5))
          - (y01/(n01-y01 + .5)) * ((n00-y00)/(y00 + .5)) 
          + 1)
}

boot(data$OC_Cases, R=k, statistic=calc_reri_equation3)

     