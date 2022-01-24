# File containing extensions to the NIMBLE library for more efficient sampling of variance-covariance matrices

## 0. ------ DEREGISTER ANY PREVIOUSLY REGISTERED DISTRIBUTIONS ------
# In some versions of NIMBLE the distribution registeration proceedures can return an error if the distributions
# have already been previously registered.  This little section of code deregisters the functions to ensure no
# errors are thrown if this file is sourced more than once.
if(exists("distributions", nimbleUserNamespace)) {
  # List of distributions defined in this source file
  distributionNames <- c(
    "dinvWishMNormMargSigma")
  # Find those distributions that are defined in this source file and see if they are already registered
  isDefinedDist <- distributionNames %in% nimbleUserNamespace$distributions$namesVector
  if(any(isDefinedDist)) {
    # If the distributions are registered then deregister them
    deregisterDistributions(distributionNames[isDefinedDist])
  }
}

## 1. ------ DEFINE HELPER FUNCTIONS ------

### 1.1. ==== Define the multivariate log-gamma function ====
lmgamma <- nimbleFunction(
  run = function(
    p = double(0),    # The dimensionality of the multivariate gamma distribution (currently only defined for integer values for p)
    a = double(0)     # The value to evaluate the function at
  ) {
    # Return type declaration
    returnType(double(0))
    # Sanity check the inputs
    if(p < 1.0) {
      stop("invalid value for the number of dimensions")
    }
    # Initialise a set of inputs to pass to the log-gamma function
    gamInputs <- a + (1.0 - seq(1, floor(p), 1.0)) / 2.0
    # Calculate the multivariate gamma function based on the sum of the log-gamma calls (plus a scaling factor)
    outValue <- (p * (p - 1.0) / 4.0) * log(pi) + sum(lgamma(gamInputs))
    return(outValue)
  }
)

## 2. ------ DEFINE DISTRIBUTION FOR INVERSE WISHART-MULTIVARIATE NORMAL DISTRIBUTION WITH KNOWN MEAN ------

### 2.1. ==== Define the density function ====
dinvWishMNormMargSigma <- nimbleFunction(
  run= function(
    x = double(2),
    wishPsi = double(2),
    wishV = double(0),
    mnormMu = double(1),
    N = double(0),
    log = integer(0, default = 0)
  ) {
    ## 2.1.1. ---- Declare the return type ----
    # Return type declaration
    returnType(double(0))
    ## 2.1.2. ---- Sanity test the parameters ----
    # Test for consistency between the dimensions
    numDims <- dim(x)[1]
    numData <- dim(x)[2]
    if(numDims <= 0) {
      stop("invalid dimension structure for the response variable")
    }
    if(numDims != dim(wishPsi)[1] | numDims != dim(wishPsi)[2]) {
      stop("dimension of inverse-Wishart parameter does not match diemsnsion structure of response variable")
    }
    if(numDims != length(mnormMu)) {
      stop("dimension of mean for the multivariate normal distribution does not match dimension structure of response variable")
    }
    if(numData != N) {
      stop("invalid entry for the number of draws made from the distribution")
    }
    if(numData <= 0) {
      # If there are only zero columns then likelihood is one
      if(log) {
        return(0.0)
      } else {
        return(1.0)
      }
    }
    if(wishV <= numDims - 1.0) {
      # If the degrees of freedom has an invalid value then return a zero-likelihood
      if(log) {
        return(-Inf)
      } else {
        return(0.0)
      }
    }
    ## 2.1.3. ---- Calculate the log density ----
    # Calculate the matrix A = x x^T
    diffMat <- x - matrix(rep(mnormMu, numData), nrow = numDims, ncol = numData)
    aMatrix <- diffMat %*% t(diffMat)
    # Use the formula derived here to calculate the log density: https://en.wikipedia.org/wiki/Inverse-Wishart_distribution#Conjugate_distribution
    logDens <- (wishV / 2.0) * logdet(wishPsi) + lmgamma(numDims, (wishV + numData) / 2.0) - (numData * numDims / 2.0) * log(pi) - ((wishV + numData) / 2.0) * logdet(wishPsi + aMatrix) - lmgamma(numDims, wishV / 2.0)
    if(log == 0) {
      # Export the output probability on the natural scale if requested
      logDens <- exp(logDens)
    }
    return(logDens)
  }
)

### 2.2 ==== Define the sampling function ====
rinvWishMNormMargSigma <- nimbleFunction(
  run = function(
    n = integer(0),
    wishPsi = double(2),
    wishV = double(0),
    mnormMu = double(1),
    N = double(0)
  ) {
    ## 2.2.1. ---- Declare the return type ----
    # Return type declaration
    returnType(double(2))
    ## 2.2.2. ---- Sanity test the parameters ----
    # Ensure that only one sample is requested
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else {
      cat("rinvWishMNormMargSigma only allows n = 1; using n = 1\n")
    }
    # Test for consistency between the dimensions
    numDims <- length(mnormMu)
    if(numDims <= 0) {
      stop("invalid dimension structure for the response variable")
    }
    if(numDims != dim(wishPsi)[1] | numDims != dim(wishPsi)[2]) {
      stop("dimension of inverse-Wishart parameter does not match diemsnsion structure of response variable")
    }
    if(wishV <= numDims - 1.0) {
      stop("invalid degrees of freedom set for the inverse-Wishart distribution")
    }
    if(N <= 0) {
      stop("invalid entry for the number of draws made from the distribution")
    }
    ## 2.2.3. ---- Sample from the inverse-Wishart multivariate Normal distribution ----
    # Calculate the Cholesky factor of the psi parameter
    psiChol <- chol(wishPsi)
    outVals <- matrix(0, nrow = numDims, ncol = N)
    for(dataIter in 1:n) {
      # Sample the precision matrix from the Wishart distribution
      precChol <- chol(rwish_chol(1, psiChol, wishV, FALSE))
      # Sample from the multivariate normal distribution using the sampled precision matrix
      outVals[1:numDims, dataIter] <- rmnorm_chol(1, mnormMu, precChol, TRUE)
    }
    return(outVals)
  }
)

## 3. ------ REGISTER THE DISTRIBUTIONS ------
# Register the distribution with NIMBLE
registerDistributions(list(
  ### 3.1. ==== Register rinvWishMNormSigma distribution ====
  dinvWishMNormMargSigma = list(
    ## 3.1.1. ---- Define the BUGS code to call the distribution ----
    BUGSdist = "dinvWishMNormMargSigma(wishPsi, wishV, mnormMu, N)",
    ## 3.1.2. ---- Set the input and output types and dimension structure ----
    types = c(
      "value = double(2)", "wishPsi = double(2)", "wishV = double(0)", "mnormMu = double(1)", "N = double(0)"
    ),
    ## 3.1.3. ---- Define the cumulative probability and quantile function availability
    pqAvail = FALSE
  )
))
