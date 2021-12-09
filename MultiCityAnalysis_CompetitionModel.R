library(parallel)
library(coda)
library(nimble)

# Initialise appropriate parameters for MCMC analysis
mcmcSamples <- 150000
mcmcChains <- 4
mcmcBurnIn <- 10000
numVarsWanted <- 500
numPredsWanted <- 500
monitorOneThin <- floor(mcmcSamples * mcmcChains / numVarsWanted)
monitorTwoThin <- floor(mcmcSamples * mcmcChains / numPredsWanted)

# Set the location of the pan trap and pollinator collection data
urbanSISRepository <- file.path(Sys.getenv("WORKSPACE_URBANSIS"), "WP2 - Mapping - 15885002")
outputLocation <- file.path(urbanSISRepository, "MultiCityAnalysis_Output")
# outputLocation <- file.path(Sys.getenv("WORKSPACE_URBANSIS_ANALYSIS"), "MultiCityAnalysis_Output")
# Import the processed pantrap data
inputData <- readRDS(file.path(urbanSISRepository, "MultiCity_ProcessedPantrapData.rds"))
pantrapData <- inputData$pantrapData

if(dir.exists(outputLocation)) {
  unlink(outputLocation, recursive = TRUE)
}
dir.create(outputLocation)
# Add trapping duration information to the pantrap data
pantrapData <- cbind(pantrapData, data.frame(
  # Get an estimate for the number of trapping days
  estTrappingDays = abs(as.numeric(difftime(
    as.Date(paste(
      pantrapData$yearBegin,
      pantrapData$monthBegin,
      ifelse(is.na(pantrapData$dayBegin), 1, pantrapData$dayBegin),
      sep = "-")),
    as.Date(paste(
      ifelse(is.na(pantrapData$dayEnd) & pantrapData$monthEnd >= 12, 1, 0) + pantrapData$yearEnd,
      ifelse(is.na(pantrapData$dayEnd), ifelse(pantrapData$monthEnd >= 12, 1, pantrapData$monthEnd + 1), pantrapData$monthEnd),
      ifelse(is.na(pantrapData$dayEnd), 1, pantrapData$dayEnd),
      sep = "-"))
  )))
))
# Set the variables to monitor
varsToMonitorOne <- c("interactionMatrix", "specCoef")
varsToMonitorTwo <- c("transPred")
# Define a function to setup and run the model
runModelMCMC <- function(curIndex, seedArray, speciesNames, pantrapData, varsToMonitor, predsToMonitor, mcmcSamples, mcmcBurnIn, mcmcChains, monitorOneThin, monitorTwoThin, outputLocation) {
  library(nimble)
  # Retrieve the current city
  curCity <- levels(pantrapData$cityID)[curIndex]
  isCurCity <- as.character(pantrapData$cityID) == curCity
  sink(file = file.path(outputLocation, paste("MCMCLogFile_", curCity, ".txt", sep = "")))
  # Retrieve the current seed from the seed array
  curSeed <- seedArray[curIndex]
  # Retrieve those species that are present in the city
  specMatrix <- as.matrix(pantrapData[isCurCity, gsub(" ", ".", speciesNames, fixed = TRUE)])
  isPresentSpecies <- apply(X = specMatrix, FUN = sum, MARGIN = 2) > 0
  curSpeciesNames <- speciesNames[isPresentSpecies]
  # Retrieve the months that are present in the data set
  monthsPresent <- sort(unique(pantrapData$monthBegin))
  # Retrieve the data for the current city
  inData <- list(
    trapCounts = specMatrix[, isPresentSpecies],
    estTrappingDays = pantrapData$estTrappingDays[isCurCity]
  )
  # Specify the NIMBLE model
  modelSpecCode <- nimbleCode({
    ###### SET PRIORS ######
    for(priorSpecIter in 1:numSpecies) {
      # Set a prior for the variance in the number of individuals for each species
      specVar[priorSpecIter] ~ dgamma(0.01, 0.01)
      for(priorMonthIter in 1:numMonths) {
        # Set a prior for the background mean number of individuals for each species caught per trapping day
        specCoef[priorMonthIter, priorSpecIter] ~ dnorm(0.0, 0.01)
      }
    }
    # Initialise an \"uninformative prior\" for the matrix of species interactions (covariance-variance matrix)
    # However see http://dx.doi.org/10.1080/00273171.2015.1065398 for explanation that this might be fairly informative
    # in situations where the variance is small
    specVarMat[1:numSpecies, 1:numSpecies] <- diag(specVar[1:numSpecies])
    interactionMatrix[1:numSpecies, 1:numSpecies] ~ dinvwish(specVarMat[1:numSpecies, 1:numSpecies], numSpecies)
    ###### SET LIKELIHOOD ######
    for(dataIter in 1:numData) {
      # Define the distribution of mean predictors at each location (once species interactions are taken into account)
      expTrapMean[dataIter, 1:numSpecies] <- specCoef[indexedMonth[dataIter], 1:numSpecies]
      meanPred[dataIter, 1:numSpecies] ~ dmnorm(expTrapMean[dataIter, 1:numSpecies], cov = interactionMatrix[1:numSpecies, 1:numSpecies])
      for(specIter in 1:numSpecies) {
        # Use the inverse link function
        transPred[dataIter, specIter] <- exp(meanPred[dataIter, specIter])
        transPredTrap[dataIter, specIter] <- transPred[dataIter, specIter] * estTrappingDays[dataIter]
        trapCounts[dataIter, specIter] ~ dpois(transPredTrap[dataIter, specIter])
      }
    }
  })
  # Setup the constants used in the model
  inConstants <- list(
    numSpecies = length(curSpeciesNames),
    numData = nrow(inData$trapCounts),
    numMonths = length(monthsPresent),
    indexedMonth = as.integer(factor(pantrapData$monthBegin, as.character(monthsPresent)))
  )
  # Set initial values
  specMeans <- log(apply(X = inData$trapCounts, FUN = mean, MARGIN = 2))
  specVar <- apply(X = inData$trapCounts, FUN = var, MARGIN = 2)
  initialValues <- list(
    specVar = specVar,
    interactionMatrix = diag(specVar),
    meanPred = matrix(rep(specMeans, inConstants$numData), nrow = inConstants$numData, ncol = inConstants$numSpecies),
    specCoef = matrix(rep(specMeans / mean(inData$estTrappingDays), inConstants$numMonths), byrow = TRUE, nrow = inConstants$numMonths, ncol = inConstants$numSpecies)
  )
  # Setup the model
  uncompiledModel <- nimbleModel(code = modelSpecCode, constants = inConstants, data = inData, inits = initialValues)
  compiledModel <- compileNimble(uncompiledModel)
  # Define the MCMC object
  uncompiledMCMC <- buildMCMC(uncompiledModel, monitors = varsToMonitor, monitors2 = predsToMonitor, thin = monitorOneThin, thin2 = monitorTwoThin, enableWAIC = FALSE)
  compiledMCMC <- compileNimble(uncompiledMCMC, project = uncompiledModel)
  # Run the MCMC
  mcmcOutput <- runMCMC(compiledMCMC, niter = mcmcSamples + mcmcBurnIn, nburnin = mcmcBurnIn,
    thin = monitorOneThin, thin2 = monitorTwoThin, nchains = mcmcChains,
    samplesAsCodaMCMC = TRUE, WAIC = FALSE, summary = FALSE, setSeed = curSeed)
  sink()
  list(
    mcmcOutput = mcmcOutput,
    specMatrix = inData$trapCounts)
}
# Retrieve the number of cities in the dataset
numCities <- length(levels(pantrapData$cityID))
# Set the number of cores to use
numCores <- min(numCities, detectCores())
# Initialise a cluster
chainCluster <- makeCluster(numCores)
# Run the parallelised version of NIMBLE
allChainsOutput <- parLapply(cl = chainCluster, X = 1:numCities, fun = runModelMCMC,
  seedArray = floor(runif(numCities, 0.0, .Machine$integer.max)),
  speciesNames = inputData$speciesNames, pantrapData = pantrapData,
  varsToMonitor = varsToMonitorOne,
  predsToMonitor = varsToMonitorTwo, mcmcSamples = mcmcSamples, mcmcBurnIn = mcmcBurnIn, mcmcChains = mcmcChains,
  monitorOneThin = monitorOneThin, monitorTwoThin = monitorTwoThin, outputLocation = outputLocation)
names(allChainsOutput) <- levels(pantrapData$cityID)
# Close the cluster
stopCluster(chainCluster)
# Save the model output to the analysis output directory
saveRDS(allChainsOutput, file = file.path(outputLocation, "mcmcOutput.rds"))
