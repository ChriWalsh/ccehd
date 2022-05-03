#' Simulate a data set with interactive fixed effects
#'
#' \code{simulate_data} simulates data as in simulation section of Vogt, M.,
#' Walsh, C., and Linton, 0. (2022) "Estimation of High-Dimensional Panel Data
#' Models with Interactive Fixed Effects". The example provides the code used
#' to create the data set \bold{data_example.rda} of the package.
#'
#' @param obs_N number of cross-section units.
#' @param obs_T number of time periods.
#' @param MU vector of mean individual loadings.
#' @param RHO pairwise correlation coefficient of the regressors.
#' @param p1 number of additional regressors that only depend on first factor
#' @param p2 number of additional regressors that only depend on second factor
#' @param p3 number of additional regressors that only depend on third factor
#' @param K the number of unobserved factors. This has only been tested with
#'          the default value of K=3.
#' @return List containing the balanced panel data with data$y containing
#'   the dependent variables and data$x the regressors. Both are sorted such
#'   that first T observations are those for unit 1, followed by the T
#'   observations for unit 2, etc. Finally, data$f contains the unobserved
#'   factors.
#'
#'@examples
#'\dontrun{
#'# Simulate an example data set as in Scenario B of the paper with
#'# N = 50, T = 50, rho = 0.25, p = 900
#'#----------------------------------------------------------------
#'
#'# Set the parameters
#'#-------------------
#'
#'# Set the number of cross-sections, time periods
#'obs_N  <- 50
#'obs_T  <- 50
#'# Specify the number of additional "irrelevant" regressors
#'p1 <- p2 <- p3 <- 299
#'# Specify the pairwise correlation among the regressors
#'rho <- 0.25
#'# Specify the mean for the factor loadings
#'MU <- c(rep(1,3),c(diag(0.5,nrow=3)),rep(1,(p1+p2+p3)))
#'
#'# Simulate the data
#'set.seed(1048)
#'data_example <- simulate_data(obs_N=obs_N, obs_T=obs_T, MU=MU, RHO=rho, p1=p1, p2=p2, p3=p3)
#'
#'# Rearrange the order of the additional "irrelevant" regressors
#'# for the example data of the package.
#'#----------------------------------------------------------------
#'
#'# Rearrange the "additional" columns so that they are not in blocks
#'indices <- c(1:3)
#'
#'for (s in 1:p1){
#'  indices <- c(indices,3+s,3+p1+s,3+2*p1+s)
#'}
#'
#'data_example$x <- data_example$x[,indices]
#'
#'save(data_example, file="data_example.rda",compress=TRUE)
#'}




simulate_data <- function(obs_N,obs_T,MU,RHO,p1,p2,p3,K=3){

#==============================================================================#
#  Inputs:
#  - obs_N -- number of units n
#  - obs_T -- number of time periods T
#  - MU    -- mean vector of factor loadings
#  - RHO   -- correlation between factor loadings
#  - p1    -- number of additional regressors of type 1
#  - p2    -- number of additional regressors of type 2
#  - p3    -- number of additional regressors of type 3
#  - K     -- number of common factors
#
#  Output:
#  - X_all -- matrix with nT rows and p=3+p1+p2+p3 columns,
#             each column contains the data of one regressor
#  - Y_all -- vector of length nT containing data of dependent variable
#==============================================================================#

#==============================================================================#
# (1)  parameter definitions
#==============================================================================#

p <- 3+p1+p2+p3

#==============================================================================#
# (2)  common unobserved factors
#      [(T x K) matrix]
#
#      The common factors are generated as independent zero mean
#      unit variance AR(1) processes with AR coefficient 0.5 and
#      Gaussian innovations.
#==============================================================================#

# Set a "burn-in" phase
BURN_IN <- 50

# Generate storage variables for factors
F_all <- matrix(NA,ncol=K,nrow=(obs_T+BURN_IN))

# Set the autoregressive coefficient and the standard deviation of
# the innovation to ensure unit variance factors
COEF_F <- 0.5
SD_F <- sqrt(1-COEF_F^2)

# Set the initial values for the AR processes
INIT_F <- 0

# Begin the AR(1) recursion
F_all[1,] <- rep(INIT_F,K)

for(t in 2:(BURN_IN+obs_T)){
   F_all[t,] <- COEF_F * F_all[t-1,] + stats::rnorm(n=K,mean=0,sd=SD_F)
}

# Drop the "burn-in" phase
F_all <- F_all[-(1:BURN_IN),]

#==============================================================================#
# (3)  idiosyncratic part of the regressors for unit i
#      [(T x p) matrix]
#
#      The idiosyncratic parts of the regressors are
#      standard Gaussian white noise.
#==============================================================================#

reg_idio <- function(){
  Z_i <- matrix(stats::rnorm(n=p*obs_T,mean=0,sd=1),ncol=p,nrow=obs_T)
  return(Z_i)
}

#==============================================================================#
# (4)  idiosyncratic errors for unit i
#      [vector of length T]
#
#      The idiosyncratic errors for unit i are standard Gaussian white noise.
#==============================================================================#

error_idio <- function(){
  varepsilon_i  <- stats::rnorm(n=obs_T,mean=0,sd=1)
  return(varepsilon_i)
}

#==============================================================================#
# (5)  factor loadings for unit i
#==============================================================================#

loadings <- function(){

  GAMMA_Y <- rep(NA,K)
  GAMMA_X <- matrix(0,nrow=p,ncol=K)

  SIGMA <- matrix(RHO,ncol=(4*K+p1+p2+p3),nrow=(4*K+p1+p2+p3))
  diag(SIGMA) <- 1
  temp <- MASS::mvrnorm(n=1,mu=MU,Sigma=SIGMA)

  GAMMA_Y <- temp[1:3]
  GAMMA_X[1:3,] <- matrix(temp[4:12],ncol=3)
  if(p1 != 0)
    GAMMA_X[4:(3+p1),1] <- temp[13:(12+p1)]
  if(p2 != 0)
    GAMMA_X[(4+p1):(3+p1+p2),2] <- temp[(13+p1):(12+p1+p2)]
  if(p3 != 0)
    GAMMA_X[(4+p1+p2):(3+p1+p2+p3),3] <- temp[(13+p1+p2):(12+p1+p2+p3)]

  return(list(x=GAMMA_X,y=GAMMA_Y))
}

#==============================================================================#
# (6)  dependent variable and regressors
#==============================================================================#

# Generate storage matrix for regressors and the dependent variable
X_all <- matrix(NA,ncol=p,nrow=obs_T*obs_N)
Y_all <- rep(NA, length=obs_T*obs_N)

# Set the coefficient vector
BETA <- c(1,1,1,rep(0,p1+p2+p3))

# Construct the scale matrix for idiosyncratic part of regressors
if(p1+p2+p3==0)
  LAMBDA <- rep(1,3)
if(p1+p2+p3!=0)
  LAMBDA <- c(rep(1,3),sqrt(3.25-(MU[-(1:12)])^2))

# Iterate over the units
for(i in 1:obs_N){

  # indices of unit i
  indices <- ((i-1)*obs_T + 1):(i*obs_T)

  # loadings of unit i
  temp <- loadings()
  gamma_vec_i <- temp$y
  GAMMA_i <- temp$x

  # regressors of unit i
  X_all[indices,] <- F_all %*% t(GAMMA_i) + matrix(rep(LAMBDA,obs_T),nrow=obs_T,byrow=TRUE) * reg_idio()

  # dependent variable of unit i
  Y_all[indices] <- X_all[indices,] %*% t(t(BETA)) + F_all %*% t(t(gamma_vec_i)) + error_idio()
}

return(list(y=Y_all,x=X_all,f=F_all))

}



