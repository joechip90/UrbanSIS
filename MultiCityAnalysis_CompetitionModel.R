library(parallel)
library(coda)
library(nimble)
library(abind)

# Initialise appropriate parameters for MCMC analysis
mcmcSamples <- ifelse(exists("mcmcSamples"), mcmcSamples, 150000)
mcmcChains <- ifelse(exists("mcmcChains"), mcmcChains, 4)
mcmcBurnIn <- ifelse(exists("mcmcBurnIn"), mcmcBurnIn, 10000)
numVarsWanted <- ifelse(exists("numVarsWanted"), numVarsWanted, 4000)
numPredsWanted <- ifelse(exists("numPredsWanted"), numPredsWanted, 4000)
monitorOneThin <- ifelse(exists("monitorOneThin"), monitorOneThin, floor(mcmcSamples * mcmcChains / numVarsWanted))
monitorTwoThin <- ifelse(exists("monitorTwoThin"), monitorTwoThin, floor(mcmcSamples * mcmcChains / numPredsWanted))

cat("Running ", mcmcChains, " chains with a burn-in of ", mcmcChains, " samples and ", mcmcSamples, " samples collected after burn-in\n", sep = "")
cat("Thinning interval is set at ", monitorOneThin, " for variables (for a total of ", numVarsWanted, " samples across all chains) and ", monitorTwoThin, " for predictions (for a total of ", numPredsWanted, " samples across all chains)\n", sep = "")

# Set the location of the pan trap and pollinator collection data
urbanSISRepository <- file.path(Sys.getenv("WORKSPACE_URBANSIS"), "WP2 - Mapping - 15885002")
outputLocation <- file.path(urbanSISRepository, "MultiCityAnalysis_Output")
#outputLocation <- file.path(Sys.getenv("WORKSPACE_URBANSIS_ANALYSIS"), "MultiCityAnalysis_Output")
# Import the processed pantrap data
inputData <- readRDS(file.path(urbanSISRepository, "MultiCity_ProcessedPantrapData.rds"))
pantrapData <- inputData$pantrapData

# Remove the output directory if it already exists
if(dir.exists(outputLocation)) {
  unlink(outputLocation, recursive = TRUE)
}
dir.create(outputLocation)

# Set the variables to monitor
varsToMonitorOne <- c("specCoef", "specVar")            # Monitors the log of mean species counts at each city
varsToMonitorTwo <- c("transPred")                      # Monitors the predictions made at each city

# Define a function to setup and run the model
runModelMCMC <- function(curChainNumber, seedArray, speciesNames, pantrapData, varsToMonitor, predsToMonitor, mcmcSamples, mcmcBurnIn, mcmcChains, monitorOneThin, monitorTwoThin, outputLocation) {
  # Create a log file for the current chain
  sink(file = file.path(outputLocation, paste("MCMCLogFile_chain", curChainNumber, ".txt", sep = "")))
  # Load the nimble and coda libraries
  library(coda)
  library(nimble)
  source("https://raw.githubusercontent.com/joechip90/UrbanSIS/master/NIMBLEExtensions.R")
  # Retrieve the current seed from the seed array
  curSeed <- seedArray[curChainNumber]
  # Retrieve the data for the current city
  inData <- list(
    trapCounts = as.matrix(pantrapData[, gsub(" ", ".", speciesNames, fixed = TRUE)])
  )
  # Setup the constants used in the model
  inConstants <- list(
    numSpecies = length(speciesNames),
    numData = nrow(inData$trapCounts),
    numCities = length(levels(pantrapData$cityID)),
    cityID = as.integer(pantrapData$cityID)
  )
  # Set initial values
  specMeans <- log(apply(X = inData$trapCounts, FUN = mean, MARGIN = 2))
  specVar <- apply(X = inData$trapCounts, FUN = var, MARGIN = 2)
  initialValues <- list(
    specVar = specVar,
    # interactionMatrix = diag(specVar),
    covPred = matrix(0, nrow = inConstants$numSpecies, ncol = inConstants$numData) + diag(specVar),
    # meanPred = matrix(rep(specMeans, inConstants$numData), nrow = inConstants$numData, ncol = inConstants$numSpecies, byrow = TRUE),
    specCoef = matrix(rep(specMeans, inConstants$numCities), nrow = inConstants$numCities, ncol = inConstants$numSpecies, byrow = TRUE)
  )
  # Specify the NIMBLE model
  modelSpecCode <- nimbleCode({
    ###### SET PRIORS ######
    for(priorSpecIter in 1:numSpecies) {
      # Set a prior for the variance in the number of individuals for each species
      specVar[priorSpecIter] ~ dgamma(0.01, 0.01)
      for(priorCityIter in 1:numCities) {
        # Each species has a different mean count in each city
        specCoef[priorCityIter, priorSpecIter] ~ dnorm(0.0, 0.01)
      }
    }
    # Initialise an \"uninformative prior\" for the matrix of species interactions (covariance-variance matrix)
    # However see http://dx.doi.org/10.1080/00273171.2015.1065398 for explanation that this might be fairly informative
    # in situations where the variance is small
    specVarMat[1:numSpecies, 1:numSpecies] <- diag(specVar[1:numSpecies])
    # interactionMatrix[1:numSpecies, 1:numSpecies] ~ dinvwish(specVarMat[1:numSpecies, 1:numSpecies], numSpecies)
    zeroVec[1:numSpecies] <- rep(0.0, numSpecies)
    covPred[1:numSpecies, 1:numData] ~ dinvWishMNormMargSigma(specVarMat[1:numSpecies, 1:numSpecies], numSpecies, zeroVec[1:numSpecies], numData)
    ###### SET LIKELIHOOD ######
    for(dataIter in 1:numData) {
      # Define the distribution of mean predictors at each location (once species interactions are taken into account)
      # meanPred[dataIter, 1:numSpecies] ~ dmnorm(specCoef[cityID[dataIter], 1:numSpecies], cov = interactionMatrix[1:numSpecies, 1:numSpecies])
      meanPred[dataIter, 1:numSpecies] <- specCoef[cityID[dataIter], 1:numSpecies] + covPred[1:numSpecies, dataIter]
      for(specIter in 1:numSpecies) {
        # Use the inverse link function
        transPred[dataIter, specIter] <- exp(meanPred[dataIter, specIter])
        # Relate the output to the number of individuals found in the traps (through the Poison distribution)
        trapCounts[dataIter, specIter] ~ dpois(transPred[dataIter, specIter])
      }
    }
  })
  # Setup the model
  uncompiledModel <- nimbleModel(code = modelSpecCode, constants = inConstants, data = inData, inits = initialValues)
  compiledModel <- compileNimble(uncompiledModel)
  # Define the MCMC object
  uncompiledMCMC <- buildMCMC(uncompiledModel, monitors = varsToMonitor, monitors2 = predsToMonitor, thin = monitorOneThin, thin2 = monitorTwoThin, enableWAIC = FALSE)
  compiledMCMC <- compileNimble(uncompiledMCMC, project = uncompiledModel)
  # Run the MCMC
  mcmcOutput <- runMCMC(compiledMCMC, niter = mcmcSamples + mcmcBurnIn, nburnin = mcmcBurnIn,
                        thin = monitorOneThin, thin2 = monitorTwoThin, nchains = 1,
                        samplesAsCodaMCMC = TRUE, WAIC = FALSE, summary = FALSE, setSeed = curSeed)
  # Turn-off the log file for the current chain
  sink()
  mcmcOutput
}

# debug(runModelMCMC)
#mcmcOutput <- lapply(X = 1:mcmcChains,  FUN = runModelMCMC,
#                     seedArray = floor(runif(mcmcChains, 0.0, .Machine$integer.max)),
#                     speciesNames = inputData$speciesNames, pantrapData = pantrapData,
#                     varsToMonitor = varsToMonitorOne, predsToMonitor = varsToMonitorTwo,
#                     mcmcSamples = mcmcSamples, mcmcBurnIn = mcmcBurnIn, mcmcChains = mcmcChains,
#                     monitorOneThin = monitorOneThin, monitorTwoThin = monitorTwoThin, outputLocation = outputLocation)

# Initialise a cluster with a process for each chain
chainCluster <- makeCluster(mcmcChains)
# Run the model
mcmcOutput <- parLapply(cl = chainCluster, X = 1:mcmcChains, fun = runModelMCMC,
                        seedArray = floor(runif(mcmcChains, 0.0, .Machine$integer.max)),
                        speciesNames = inputData$speciesNames, pantrapData = pantrapData,
                        varsToMonitor = varsToMonitorOne, predsToMonitor = varsToMonitorTwo,
                        mcmcSamples = mcmcSamples, mcmcBurnIn = mcmcBurnIn, mcmcChains = mcmcChains,
                        monitorOneThin = monitorOneThin, monitorTwoThin = monitorTwoThin, outputLocation = outputLocation)
# Stop the cluster once the processing is complete
stopCluster(chainCluster)

# Retrieve the estimated means for each of the species
sampledSpeciesMeans <- do.call(abind, lapply(X = 1:length(mcmcOutput), FUN = function(curIndex, mcmcOutput, speciesNames, cityNames, mcmcBurnIn, monitorOneThin) {
  # Retrieve the current output
  curOutput <- mcmcOutput[[curIndex]]
  # Retrieve all the samples of the species means
  allSamples <- as.data.frame(curOutput$samples)
  allSamples <- allSamples[, grepl("^specCoef", colnames(allSamples), perl = TRUE)]
  # Reorder the samples into one multi-dimensional array
  outValues <- apply(X = as.matrix(allSamples), FUN = function(curRow, speciesNames, cityNames) {
    outMat <- matrix(curRow, nrow = length(cityNames), ncol = length(speciesNames), dimnames = list(cityNames, speciesNames), byrow = TRUE)
    outMat
  }, MARGIN = 1, speciesNames = speciesNames, cityNames = cityNames)
  dim(outValues) <- c(length(cityNames), length(speciesNames), nrow(allSamples))
  dimnames(outValues) <- list(cityNames, speciesNames, paste("chain", curIndex, "_mcmcSample", mcmcBurnIn + 1:nrow(allSamples) * monitorOneThin, sep = ""))
  outValues
}, mcmcOutput = mcmcOutput, speciesNames = gsub(" ", ".", inputData$speciesNames, fixed = TRUE), cityNames = levels(inputData$pantrapData$cityID), mcmcBurnIn = mcmcBurnIn, monitorOneThin = monitorOneThin))

# Retrieve the predicted values for each of the species at each of the data points
sampledPredictions <- do.call(abind, lapply(X = 1:length(mcmcOutput), FUN = function(curIndex, mcmcOutput, speciesNames, sampleID, mcmcBurnIn, monitorTwoThin) {
  # Retrieve the current output
  curOutput <- mcmcOutput[[curIndex]]
  # Retrieve all the samples of the species means
  allSamples <- as.data.frame(curOutput$samples2)
  # Reorder the samples into one multi-dimensional array
  outValues <- apply(X = as.matrix(allSamples), FUN = function(curRow, speciesNames, sampleID) {
    outMat <- matrix(curRow, nrow = length(sampleID), ncol = length(speciesNames), dimnames = list(sampleID, speciesNames), byrow = TRUE)
    outMat
  }, MARGIN = 1, speciesNames = speciesNames, sampleID = sampleID)
  dim(outValues) <- c(length(sampleID), length(speciesNames), nrow(allSamples))
  dimnames(outValues) <- list(sampleID, speciesNames, paste("chain", curIndex, "_mcmcSample", mcmcBurnIn + 1:nrow(allSamples) * monitorTwoThin, sep = ""))
  outValues
}, mcmcOutput = mcmcOutput, speciesNames = gsub(" ", ".", inputData$speciesNames, fixed = TRUE), sampleID = as.character(inputData$pantrapData$sampleID), mcmcBurnIn = mcmcBurnIn, monitorTwoThin = monitorTwoThin))

# Retrieve the sampled variances for the different species
sampledVariances <- do.call(abind, lapply(X = 1:length(mcmcOutput), FUN = function(curIndex, mcmcOutput, speciesNames, mcmcBurnIn, monitorOneThin) {
  # Retrieve the current output
  curOutput <- mcmcOutput[[curIndex]]
  # Retrieve all the samples of the species means
  allSamples <- as.data.frame(curOutput$samples)
  allSamples <- allSamples[, grepl("^specVar", colnames(allSamples), perl = TRUE)]
  # Reorder the samples into one multi-dimensional array
  outValues <- apply(X = as.matrix(allSamples), FUN = function(curRow, speciesNames) {
    setNames(curRow, speciesNames)
  }, MARGIN = 1, speciesNames = speciesNames)
  dim(outValues) <- c(length(speciesNames), nrow(allSamples))
  dimnames(outValues) <- list(speciesNames, paste("chain", curIndex, "_mcmcSample", mcmcBurnIn + 1:nrow(allSamples) * monitorOneThin, sep = ""))
  outValues
}, mcmcOutput = mcmcOutput, speciesNames = gsub(" ", ".", inputData$speciesNames, fixed = TRUE), mcmcBurnIn = mcmcBurnIn, monitorOneThin = monitorOneThin))

# Retrieve the interaction matrix from the mcmc samples
sampledInteractionMatrix <- do.call(abind, lapply(X = (dimnames(sampledPredictions)[[3]])[dimnames(sampledPredictions)[[3]] %in% dimnames(sampledSpeciesMeans)[[3]]], FUN = function(curRowName, sampledSpeciesMeans, sampledPredictions, sampledVariances, cityID) {
  # Retrieve the current species means
  curSpeciesMean <- sapply(X = cityID, FUN = function(curCity, specMeans) {
    specMeans[curCity, ]
  }, specMeans = sampledSpeciesMeans[, , curRowName])
  # Retrieve the current prediction (minus the predicted mean to centre the distribution around zero)
  curPrediction <- log(t(sampledPredictions[, , curRowName])) - curSpeciesMean
  # Calculate the A matrix
  aMatrix <- curPrediction %*% t(curPrediction)
  outMat <- solve(rinvwish_chol(1, chol(aMatrix + diag(sampledVariances[, curRowName])), nrow(sampledVariances) + ncol(curPrediction), TRUE))
  dim(outMat) <- c(dim(outMat), 1)
  dimnames(outMat) <- list(rownames(curSpeciesMean), rownames(curSpeciesMean), curRowName)
  outMat
}, sampledSpeciesMeans = sampledSpeciesMeans, sampledPredictions = sampledPredictions, sampledVariances = sampledVariances, cityID = as.integer(pantrapData$cityID)))

## Retrieve the interaction matrix from the mcmc samples
#sampledInteractionMatrix <- do.call(abind, lapply(X = 1:length(mcmcOutput), FUN = function(curIndex, mcmcOutput, speciesNames, sampledPredictions) {
#  # Retrieve the current output
#  curOutput <- mcmcOutput[[curIndex]]
#  # Retrieve the sampled elements of the interaction matrix
#  allSamples <- as.data.frame(curOutput$samples)
#  allSamples <- allSamples[, grepl("^interactionMatrix", colnames(allSamples), perl = TRUE)]
#  # Reorder the samples into one multi-dimensional array
#  outValues <- apply(X = as.matrix(allSamples), FUN = function(curRow, speciesNames) {
#    matrix(curRow, nrow = length(speciesNames), ncol = length(speciesNames), dimnames = list(speciesNames, speciesNames))
#  }, MARGIN = 1, speciesNames = speciesNames)
#  dim(outValues) <- c(length(speciesNames), length(speciesNames), nrow(allSamples))
#  dimnames(outValues) <- list(speciesNames, speciesNames, paste("chain", curIndex, "_mcmcSample", 1:nrow(allSamples), sep = ""))
#  outValues
#}, mcmcOutput = mcmcOutput, speciesNames = gsub(" ", ".", inputData$speciesNames, fixed = TRUE), sampledPredictions = sampledPredictions))

# Merge the completed MCMC objects into single coda MCMC lists
variablesMCMC <- do.call(mcmc.list, lapply(X = mcmcOutput, FUN = function(curMCMC) {
  curMCMC$samples
}))
predictionsMCMC <- do.call(mcmc.list, lapply(X = mcmcOutput, FUN = function(curMCMC) {
  curMCMC$samples2
}))

# Save the model output to the analysis output directory
saveRDS(list(
  sampledInteractionMatrix = sampledInteractionMatrix,
  sampledSpeciesMeans = sampledSpeciesMeans,
  sampledPredictions = sampledPredictions,
  sampledVariances = sampledVariances,
  variablesMCMC = variablesMCMC,
  predictionsMCMC = predictionsMCMC
), file = file.path(outputLocation, "mcmcOutput.rds"))
