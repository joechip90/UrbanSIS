library(sf)
library(lubridate)
library(httr)
library(jsonlite)
library(nimble)

source(paste(Sys.getenv("WORKSPACE_URBANSIS_GITHUB"), "RetrieveWeatherData.R", sep = "/"))

# Set the workspace to perform the analysis in
workspaceLoc <- paste(Sys.getenv("WORKSPACE_URBANSIS_ANALYSIS"), "OsloBeeAnalysis", sep = "/")
# Import the transect and pantrap data
processedData <- readRDS(paste(workspaceLoc, "ProcessedPantrapAndTransectData.rds", sep = "/"))

# Weather elements to retrieve from met data
weatherElements <- c("air_temperature", "wind_speed")
# Create a border region for plotting and retrieval of weather station data
locCoords <- st_coordinates(processedData$locations)
locXRange <- range(locCoords[, "X"], na.rm = TRUE)
locYRange <- range(locCoords[, "Y"], na.rm = TRUE)
paddingPercent <- 10 / 100
borderPolygon <- st_sfc(st_polygon(list(matrix(c(
  locXRange[1] - diff(locXRange) * paddingPercent, locYRange[1] - diff(locYRange) * paddingPercent,
  locXRange[2] + diff(locXRange) * paddingPercent, locYRange[1] - diff(locYRange) * paddingPercent,
  locXRange[1] - diff(locXRange) * paddingPercent, locYRange[2] + diff(locYRange) * paddingPercent,
  locXRange[2] + diff(locXRange) * paddingPercent, locYRange[2] + diff(locYRange) * paddingPercent,
  locXRange[1] - diff(locXRange) * paddingPercent, locYRange[1] - diff(locYRange) * paddingPercent
), ncol = 2, byrow = TRUE))), crs = st_crs(processedData$locations))

# Retrieve the weather station data in the region for the time periods in question
timePeriod <- range(as.POSIXct(unlist(strsplit(as.character(processedData$samples$samplingTime), "/", fixed = TRUE))), na.rm = TRUE)
timePeriod[1] <- as.POSIXct(timePeriod[1]) - hours(24)
timePeriod[2] <- as.POSIXct(timePeriod[2]) + hours(24)
metData <- retrieveNorwegianWeatherStationData(borderPolygon, timePeriod, weatherElements,
  Sys.getenv("METEOROLOGISK_INSTITUTT_FROSTID"), Sys.getenv("METEOROLOGISK_INSTITUTT_FROSTSECRET"))

# Process the weather data from the meterological institute and link them to the relevant samples
processWeatherData <- function(elementText, pantrapLocs, pantrapSamples, metLocs, metSamples) {
  # Create a distance matrix between each pantrap location and the weather stations
  wsDistMat <- st_distance(pantrapLocs, metLocs)
  # Find the lowest level at which the element is monitored
  lowestLevel <- min(metSamples[as.character(metSamples$elementId) == elementText, "levelValue"], na.rm = TRUE)
  sapply(X = 1:nrow(pantrapSamples), FUN = function(curIndex, pantrapLocs, pantrapSamples, metLocs, wsDistMat, elementSamples) {
    outVal <- NA
    if(nrow(elementSamples) > 0) {
      # Retrieve the current trap number and sampling date range at the current index
      curTrapNum <- pantrapSamples$trapID[curIndex]
      curDateRange <- as.POSIXct(strsplit(as.character(pantrapSamples$samplingTime[curIndex]), "/", fixed = TRUE)[[1]])
      trapIndex <- which(as.character(curTrapNum) == rownames(pantrapLocs))
      # Retrieve the element samples within the current time range
      curElements <- elementSamples[elementSamples$timeStamp >= curDateRange[1] & elementSamples$timeStamp <= curDateRange[2], ]
      if(nrow(curElements) == 0) {
        # In situations where there are no samples in the current time period then instead find the nearest time period
        timeVals <- apply(X = cbind(
          abs(difftime(elementSamples$timeStamp, curDateRange[1], units = "secs")),
          abs(difftime(elementSamples$timeStamp, curDateRange[2], units = "secs"))
        ), FUN = min, MARGIN = 1, na.rm = TRUE)
        curElements <- elementSamples[timeVals == min(timeVals, na.rm = TRUE), ]
      }
      # Retrieve the mean values of the relevant element at each weather station
      meanElement <- sapply(X = rownames(metLocs)[order(wsDistMat[trapIndex, ])], FUN = function(curSourceId, elementSamples) {
        mean(elementSamples$value[as.character(elementSamples$id) == curSourceId], na.rm = TRUE)
      }, elementSamples = curElements)
      outVal <- meanElement[!is.na(meanElement)]
      outVal <- ifelse(length(outVal) > 0, outVal[1], NA)
    }
    outVal
  }, pantrapLocs = pantrapLocs, pantrapSamples = pantrapSamples, metLocs = metLocs, wsDistMat = wsDistMat,
    # Filter out those meterological samples that only have the current element being tested
    elementSamples = metSamples[as.character(metSamples$elementId) == elementText & metSamples$levelValue == lowestLevel, ])
}
processedData$samples <- cbind(processedData$samples, as.data.frame(setNames(lapply(X = weatherElements, FUN = processWeatherData, pantrapLocs = processedData$locations,
  pantrapSamples = processedData$samples, metLocs = metData$locationData, metSamples = metData$observationData), weatherElements)))

# Merge the processed sample data with the information in the location frame
mergedSampleData <- merge(processedData$samples, processedData$locations, by = "trapID")

# Retrieve the resource density of the flowers (that have some variation in the dataset)
isResourceColumn <- grepl("^resourceDensity_", colnames(mergedSampleData), perl = TRUE)
suitableResCols <- colnames(mergedSampleData)[isResourceColumn][sapply(X = colnames(mergedSampleData)[isResourceColumn], FUN = function(curCol, mergedSampleData) {
  colExtents <- range(mergedSampleData[, curCol], na.rm = TRUE)
  outVal <- FALSE
  if(all(is.finite(colExtents))) {
    outVal <- ifelse(colExtents[1] == colExtents[2], FALSE, TRUE)
  }
  outVal
}, mergedSampleData = mergedSampleData)]
predictorVariables <- c(suitableResCols, "air_temperature", "wind_speed", "effortMeasure", "sampleType")
# Filter out those samples that do not have full data attached to them (quite a few transects didn't have corresponding food resource data)
mergedSampleData <- mergedSampleData[!apply(X = as.matrix(mergedSampleData[, predictorVariables]), FUN = anyNA, MARGIN = 1), ]
# Filter out those bee species that have variation in the filtered dataset
usableSpecies <- processedData$pollinatorSpecies[apply(X = mergedSampleData[, processedData$pollinatorSpecies], FUN = function(curColumn) {
  length(unique(curColumn[!is.na(curColumn)])) >= 2
}, MARGIN = 2)]

# Initialise appropriate paramaeters for MCMC analysis
mcmcSamples <- 100000
mcmcChains <- 4
mcmcBurnIn <- 10000
numVarsWanted <- 1000
numPredsWanted <- 100
monitorOneThin <- floor(mcmcSamples * mcmcChains / numVarsWanted)
monitorTwoThin <- floor(mcmcSamples * mcmcChains / numPredsWanted)
# Specify the NIMBLE code to run the model
modelSpecCode <- nimbleCode({
  ##### SET PRIORS #####
  for(priorSpecIter in 1:numSpecies) {
    # Set hyperprior for variances in species counts (in log-scale)
    specVar[priorSpecIter] ~ dgamma(0.01, 0.01)
    # Set the LASSO regularisation hyperprior and the respective coefficient as a random effect
    lassoRate[priorSpecIter] ~ dgamma(0.01, 0.01)
    for(priorCoeffIter in 1:numCoeffs) {
      regCoeffs[priorCoeffIter, priorSpecIter] ~ ddexp(0.0, lassoRate[priorSpecIter])
    }
    # Set the intercept regresion component
    interceptCoeff[priorSpecIter] ~ dnorm(0.0, 0.01)
    # Set the trap methodology coefficient
    pantrapCoeff[priorSpecIter] ~ dnorm(0.0, 0.01)
  }
  # Initialise an "uninformative prior" for the matrix of species interactions (covariance-variance matrix)
  # However see http://dx.doi.org/10.1080/00273171.2015.1065398 for explanation that this might be fairly informative
  # in situations where the variance is small
  specVarMat[1:numSpecies, 1:numSpecies] <- diag(specVar[1:numSpecies])
  interactionMatrix[1:numSpecies, 1:numSpecies] ~ dinvwish(specVarMat[1:numSpecies, 1:numSpecies], numSpecies)
  ##### SET LIKELIHOOD #####
  for(dataIter in 1:numData) {
    for(specIter in 1:numSpecies) {
      # Poisson error model for the observation of counts of each species (including correction for sampling effort)
      obsCounts[dataIter, specIter] ~ dpois(exp(meanPred[dataIter, specIter] + effortMeasure[dataIter]))
      # Calculate the linear predictor for each species
      linPred[dataIter, specIter] <- inprod(covValues[dataIter, 1:numCoeffs], regCoeffs[1:numCoeffs, specIter]) + pantrapCoeff[specIter] * isPantrap[dataIter] + interceptCoeff[specIter]
    }
    # Define the distribution of mean predictors at each location (once species interactions are taken into account)
    meanPred[dataIter, 1:numSpecies] ~ dmnorm(linPred[dataIter, 1:numSpecies], cov = interactionMatrix[1:numSpecies, 1:numSpecies])
  }
})
# Produce a matrix of regression covariates and centre and scale them
regCovariates <- apply(X = as.matrix(cbind(as.data.frame(mergedSampleData[, suitableResCols]), data.frame(
  airTemperature = mergedSampleData[, "air_temperature"],
  windSpeed = mergedSampleData[, "wind_speed"]
))), FUN = function(curCol) {
  (curCol - mean(curCol, na.rm = TRUE)) / sd(curCol, na.rm = TRUE)
}, MARGIN = 2)
# Define the constants to be used in the model
inConstants <- list(
  numCoeffs = ncol(regCovariates),
  numData = nrow(regCovariates),
  numSpecies = length(usableSpecies),
  isPantrap = ifelse(as.character(mergedSampleData$sampleType) == "pantrap", 1.0, 0.0),
  effortMeasure = mergedSampleData$effortMeasure
)
# Define the data to be used in the model
inData <- list(
  covValues = regCovariates,
  obsCounts = as.matrix(mergedSampleData[, usableSpecies])
)
# Define the variables to be monitored
varsToMonitor <- c("regCoeffs", "pantrapCoeff", "interceptCoeff", "interactionMatrix")
predsToMonitor <- c("linPred", "meanPred")
# Initialise the parameters
initialValues <- list(
  specVar = rep(1.0, inConstants$numSpecies),
  lassoRate = rep(1.0, inConstants$numSpecies),
  interceptCoeff = log(apply(X = inData$obsCounts, FUN = mean, MARGIN = 2)),
  pantrapCoeff = rep(0.0, inConstants$numSpecies),
  regCoeffs = matrix(0.0, nrow = inConstants$numCoeffs, ncol = inConstants$numSpecies),
  interactionMatrix = diag(inConstants$numSpecies),
  meanPred = matrix(rep(log(apply(X = inData$obsCounts, FUN = mean, MARGIN = 2)), rep(inConstants$numData, inConstants$numSpecies)), ncol = inConstants$numSpecies)
)
# Setup the model 
uncompiledModel <- nimbleModel(code = modelSpecCode, constants = inConstants, data = inData, inits = initialValues)
compiledModel <- compileNimble(uncompiledModel)
# Define the MCMC object
uncompiledMCMC <- buildMCMC(uncompiledModel, monitors = varsToMonitor, monitors2 = predsToMonitor, thin = monitorOneThin, thin2 = monitorTwoThin, enableWAIC = TRUE)
compiledMCMC <- compileNimble(uncompiledMCMC, project = uncompiledModel)
# Run the MCMC
mcmcOutput <- runMCMC(compiledMCMC, niter = mcmcSamples + mcmcBurnIn, nburnin = mcmcBurnIn,
  thin = monitorOneThin, thin2 = monitorTwoThin, nchains = mcmcChains,
  samplesAsCodaMCMC = TRUE, WAIC = TRUE, summary = TRUE)

# Save the model output to the analysis output directory
saveRDS(list(
  mcmcOutput = mcmcOutput,
  pollinatorSpecies = usableSpecies,
  covariateNames = colnames(regCovariates),
  analysisData = mergedSampleData
), file = paste(workspaceLoc, "NetworkAnalysisOutput.rds", sep = "/"))
