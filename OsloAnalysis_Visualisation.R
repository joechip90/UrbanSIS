library(coda)
library(nimble)
library(sf)
library(ggplot2)
library(DHARMa)

# Set the workspace location
workspaceLocation <- file.path(Sys.getenv("WORKSPACE_URBANSIS"), "WP2 - Mapping - 15885002", "OsloAnalysis_Output")
# The location of the analysis output
outLocation <- paste(workspaceLocation, "mcmcOutput.rds", sep = "/")
# Import the analysis outputs
analysisOutput <- readRDS(outLocation)

# Import the coordinates of the bee hives
beehiveLocation <- paste(Sys.getenv("WORKSPACE_ESTIMAP"), "ESTIMAP Pollination", "Pollination potential - beehives", "DATA", sep = "/")
beehivePoints <- st_transform(st_read(beehiveLocation, "pollination_potential_beehives_20181219"), st_crs(st_sf(analysisOutput$analysisData)))

# Download the coastline data
testFolder <- tempdir()
download.file("https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip",
  paste(testFolder, "coastLine.zip", sep = "/"))
unzip(paste(testFolder, "coastLine.zip", sep = "/"), exdir = paste(testFolder, "coastline", sep = "/"))
coastlineData <- st_transform(st_read(paste(testFolder, "coastline", "GSHHS_shp", "f", sep = "/"), layer = "GSHHS_f_L1"), st_crs(st_sf(analysisOutput$analysisData)))

# Retrieve the interaction matrix from the mcmc samples
retrieveInteractionMatrix <- function(mcmcOutput, speciesNames) {
  # Import all the samples from the different mcmc chains into a single data.frame
  allSamples <- do.call(rbind, lapply(X = mcmcOutput$samples, FUN = as.data.frame))
  allSamples <- allSamples[, grepl("^interactionMatrix", colnames(allSamples), perl = TRUE)]
  # Reorder the samples into one multi-dimensional array
  outValues <- apply(X = as.matrix(allSamples), FUN = function(curRow, speciesNames) {
    matrix(curRow, nrow = length(speciesNames), ncol = length(speciesNames), dimnames = list(speciesNames, speciesNames))
  }, MARGIN = 1, speciesNames = speciesNames)
  dim(outValues) <- c(length(speciesNames), length(speciesNames), nrow(allSamples))
  dimnames(outValues) <- list(speciesNames, speciesNames, paste("mcmcSample", 1:nrow(allSamples), sep = ""))
  outValues
}
sampledInteractionMatrix <- retrieveInteractionMatrix(analysisOutput$mcmcOutput, analysisOutput$pollinatorSpecies)

# Retrieve the matrix of fixed-effect coefficient from the mcmc samples
retrieveCoefficientMatrix <- function(mcmcOutput, speciesNames, covariateNames) {
  # Import all the samples from the different mcmc chains into a single data.frame
  allSamples <- do.call(rbind, lapply(X = mcmcOutput$samples, FUN = as.data.frame))
  allSamples <- allSamples[, grepl("^regCoeffs", colnames(allSamples), perl = TRUE)]
  # Reorder the samples into one multi-dimensional array
  outValues <- apply(X = as.matrix(allSamples), FUN = function(curRow, speciesNames, covariateNames) {
    matrix(curRow, nrow = length(covariateNames), ncol = length(speciesNames), dimnames = list(covariateNames, speciesNames))
  }, MARGIN = 1, speciesNames = speciesNames, covariateNames = covariateNames)
  dim(outValues) <- c(length(covariateNames), length(speciesNames), nrow(allSamples))
  dimnames(outValues) <- list(covariateNames, speciesNames, paste("mcmcSample", 1:nrow(allSamples), sep = ""))
  outValues
}
sampledCoefficientMatrix <- retrieveCoefficientMatrix(analysisOutput$mcmcOutput, analysisOutput$pollinatorSpecies, analysisOutput$covariateNames)

# Function to create a plot of the regression coefficients
createCoefficientPlot <- function(sampledCoefficientMatrix, italXLabs = TRUE) {
  # Create a data frame that summarises the MCMC samples
  summaryFrame <- do.call(rbind, lapply(X = dimnames(sampledCoefficientMatrix)[[1]], FUN = function(curRowIndex, sampledCoefficientMatrix) {
    do.call(rbind, lapply(X = dimnames(sampledCoefficientMatrix)[[2]], FUN = function(curColIndex, curRowIndex, sampledCoefficientMatrix) {
      data.frame(
        covariateName = factor(curRowIndex, dimnames(sampledCoefficientMatrix)[[1]]),
        speciesName = factor(curColIndex, dimnames(sampledCoefficientMatrix)[[2]]),
        meanEffect = mean(sampledCoefficientMatrix[curRowIndex, curColIndex, ]),
        probPos = mean(c(
          sum(sampledCoefficientMatrix[curRowIndex, curColIndex, ] > 0.0),
          sum(sampledCoefficientMatrix[curRowIndex, curColIndex, ] >= 0.0))) / dim(sampledCoefficientMatrix)[3]
      )
    }, curRowIndex = curRowIndex, sampledCoefficientMatrix = sampledCoefficientMatrix))
  }, sampledCoefficientMatrix = sampledCoefficientMatrix))
  # Create plot friendly labels for the covariates
  speciesLabels <- gsub("_", " ", dimnames(sampledCoefficientMatrix)[[2]], fixed = TRUE)
  coeffLabels <- gsub("_", " ", gsub("^resourceDensity_", "", dimnames(sampledCoefficientMatrix)[[1]], perl = TRUE), fixed = TRUE)
  coeffLabels[coeffLabels == "airTemperature"] <- "Air temperature"
  coeffLabels[coeffLabels == "windSpeed"] <- "Wind speed"
  summaryFrame$speciesName <- factor(as.character(summaryFrame$speciesName), levels = dimnames(sampledCoefficientMatrix)[[2]], labels = speciesLabels)
  summaryFrame$covariateName <- factor(as.character(summaryFrame$covariateName), levels = dimnames(sampledCoefficientMatrix)[[1]], labels = coeffLabels)
  summaryFrame <- cbind(summaryFrame, data.frame(absMean = abs(summaryFrame$meanEffect), isSignif = summaryFrame$probPos > 0.95 | summaryFrame$probPos < 0.05))
  # Create the covariate plot
  ggplot(summaryFrame, aes(x = covariateName, y = speciesName)) + geom_point(aes(fill = probPos, size = absMean), shape = 21, colour = "white") +
    geom_tile(aes(colour = isSignif), fill = NA) + scale_color_manual(values = c(NA, "grey"), guide = "none") +
    scale_fill_gradient2(midpoint = 0.5, na.value = "white") + xlab("") + ylab("") + theme_minimal() +
    guides(size = guide_legend(override.aes = list(colour = "black"))) +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, face = ifelse(italXLabs, "italic", "plain")), axis.text.y = element_text(face = "italic"))
}

# Function to create an interaction plot between the species
createInteractionPlot <- function(sampledInteractionMatrix) {
  # Function to create a summary data frame from the sampled interaction matrix
  createSummaryFrame <- function(sampledInteractionMatrix, apFun) {
    do.call(rbind, lapply(X = dimnames(sampledInteractionMatrix)[[1]], FUN = function(curRowIndex, sampledInteractionMatrix, apFun) {
      do.call(rbind, lapply(X = dimnames(sampledInteractionMatrix)[[2]], FUN = function(curColIndex, curRowIndex, sampledInteractionMatrix, apFun) {
        sendHoneyBeeToLast <- function(inNames) {
          c(sort(inNames[inNames != "Apis_mellifera"]), "Apis_mellifera")
        }
        fromLevelNames <- sendHoneyBeeToLast(dimnames(sampledInteractionMatrix)[[1]])
        toLevelNames <- sendHoneyBeeToLast(dimnames(sampledInteractionMatrix)[[2]])
        data.frame(
          fromSpec = factor(curRowIndex, fromLevelNames, labels = gsub("_", " ", fromLevelNames, fixed = TRUE)),
          toSpec = factor(curColIndex, toLevelNames, labels = gsub("_", " ", toLevelNames, fixed = TRUE)),
          value = apFun(sampledInteractionMatrix[curRowIndex, curColIndex, ]),
          isSignif = NA
        )
      }, curRowIndex = curRowIndex, sampledInteractionMatrix = sampledInteractionMatrix, apFun = apFun))
    }, sampledInteractionMatrix = sampledInteractionMatrix, apFun = apFun))
  }
  # Create a data frame of interation strength effects
  effectSizeInteraction <- createSummaryFrame(sampledInteractionMatrix, mean)
  # Create a data frame of probability of positive interaction
  probPosInteraction <- createSummaryFrame(sampledInteractionMatrix, function(inVals) {
    sum(ifelse(inVals > 0, 1.0, 0.0)) / length(inVals)
  })
  # Define the "significant" elements in the interaction matrix
  probPosInteraction$isSignif <- ifelse(probPosInteraction$value > 0.95 | probPosInteraction$value < 0.05, TRUE, FALSE)
  effectSizeInteraction$isSignif <- probPosInteraction$isSignif
  # Remove those elements on the opposite diagonal
  effectSizeInteraction$value <- ifelse(as.integer(effectSizeInteraction$fromSpec) > as.integer(effectSizeInteraction$toSpec), effectSizeInteraction$value, NA)
  effectSizeInteraction$isSignif <- ifelse(is.na(effectSizeInteraction$value), FALSE, effectSizeInteraction$isSignif)
  probPosInteraction$value <- ifelse(as.integer(probPosInteraction$fromSpec) < as.integer(probPosInteraction$toSpec), probPosInteraction$value, NA)
  probPosInteraction$isSignif <- ifelse(is.na(probPosInteraction$value), FALSE, probPosInteraction$isSignif)
  # Add in auxilliary information for the effect size information
  effectSizeInteraction <- cbind(effectSizeInteraction, data.frame(
    absValue = abs(effectSizeInteraction$value),
    isPositive = as.factor(ifelse(effectSizeInteraction$value > 0.0, "Positive", "Negative")),
    posVal = ifelse(effectSizeInteraction$isSignif,
      ifelse(effectSizeInteraction$value > 0.0, max(probPosInteraction$value - 0.5, na.rm = TRUE) * 0.8, min(probPosInteraction$value - 0.8, na.rm = TRUE) * 0.2) + 0.5,
      0.5)
  ))
  # Create a plot of the probability of positive interaction
  ggplot(probPosInteraction, aes(x = fromSpec, y = toSpec)) + geom_tile(aes(fill = value), colour = "white") +
    geom_tile(aes(colour = isSignif), data = effectSizeInteraction[!is.na(effectSizeInteraction$value), ], fill = NA) +
    geom_point(aes(size = absValue, fill = posVal), data = effectSizeInteraction[!is.na(effectSizeInteraction$value), ], shape = 21) +
    scale_fill_gradient2(midpoint = 0.5, na.value = "white") + scale_color_manual(values = c(NA, "grey"), guide = "none") +
    xlab("") + ylab("") + theme_minimal() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), axis.text.y = element_text(face = "italic"))
}

# Function to create a violin plot of the MCMC samples of regression coefficients
createCoefficientViolinPlot <- function(sampledCoefficientMatrix) {
  # Create a data frame that reorders the MCMC samples
  reorderFrame <- do.call(rbind, lapply(X = dimnames(sampledCoefficientMatrix)[[1]], FUN = function(curRowIndex, sampledCoefficientMatrix) {
    do.call(rbind, lapply(X = dimnames(sampledCoefficientMatrix)[[2]], FUN = function(curColIndex, curRowIndex, sampledCoefficientMatrix) {
      covNames <- dimnames(sampledCoefficientMatrix)[[1]]
      covLabels <- sapply(X = strsplit(gsub("([A-Z])", " \\1", covNames, perl = TRUE), " ", fixed = TRUE), FUN = function(curVals) {
        paste(tolower(curVals), collapse = " ")
      })
      covLabels <- paste(toupper(substring(covLabels, 1, 1)), substring(covLabels, 2), sep = "")
      data.frame(
        covariateName = factor(rep(curRowIndex, length(sampledCoefficientMatrix[curRowIndex, curColIndex, ])), covNames, labels = covLabels),
        speciesName = factor(rep(curColIndex, length(sampledCoefficientMatrix[curRowIndex, curColIndex, ])), dimnames(sampledCoefficientMatrix)[[2]],
          labels = gsub("_", " ", dimnames(sampledCoefficientMatrix)[[2]], fixed = TRUE)),
        values = sampledCoefficientMatrix[curRowIndex, curColIndex, ],
        probPos = rep(mean(c(
          sum(sampledCoefficientMatrix[curRowIndex, curColIndex, ] > 0.0),
          sum(sampledCoefficientMatrix[curRowIndex, curColIndex, ] >= 0.0))) / dim(sampledCoefficientMatrix)[3],
          dim(sampledCoefficientMatrix)[3])
      )
    }, curRowIndex = curRowIndex, sampledCoefficientMatrix = sampledCoefficientMatrix))
  }, sampledCoefficientMatrix = sampledCoefficientMatrix))
  ggplot(reorderFrame, aes(x = speciesName, y = values)) + geom_violin(aes(fill = probPos), scale = "width") +
    scale_fill_gradient2(midpoint = 0.5, na.value = "white") + theme_classic() +
    theme(legend.title = element_blank(), axis.text.y = element_text(face = "italic")) +
    xlab("") + ylab("") + facet_wrap(~ covariateName) + coord_flip() + geom_hline(yintercept = 0.0, linetype = "dashed", colour = "black")
}

ggsave(paste(workspaceLocation, "vegCovariatesAll.pdf", sep = "/"), plot = createCoefficientPlot(sampledCoefficientMatrix[1:124, , ]), width = 36, height = 8)
ggsave(paste(workspaceLocation, "vegCovariates1.pdf", sep = "/"), plot = createCoefficientPlot(sampledCoefficientMatrix[1:31, , ]), width = 11, height = 8)
ggsave(paste(workspaceLocation, "vegCovariates2.pdf", sep = "/"), plot = createCoefficientPlot(sampledCoefficientMatrix[32:62, , ]), width = 11, height = 8)
ggsave(paste(workspaceLocation, "vegCovariates3.pdf", sep = "/"), plot = createCoefficientPlot(sampledCoefficientMatrix[63:93, , ]), width = 11, height = 8)
ggsave(paste(workspaceLocation, "vegCovariates4.pdf", sep = "/"), plot = createCoefficientPlot(sampledCoefficientMatrix[94:124, , ]), width = 11, height = 8)
ggsave(paste(workspaceLocation, "interactionMatrix.pdf", sep = "/"), plot = createInteractionPlot(sampledInteractionMatrix), width = 11, height = 8)
ggsave(paste(workspaceLocation, "weatherCovariates.pdf", sep = "/"), plot = createCoefficientViolinPlot(sampledCoefficientMatrix[125:126, , ]), width = 8, height = 11)

# Retrieve the linear predictions and full effect predictions (and the difference between them)
predFrame <- do.call(rbind, lapply(X = analysisOutput$mcmcOutput$samples2, FUN = as.data.frame))
fullPred <- as.matrix(predFrame[, grepl("^meanPred", colnames(predFrame), perl = TRUE)])
# Create a prediction matrix for use in DHARMa plots (posterior predictive distribution)
predMatrix <- do.call(rbind, lapply(X = 1:length(analysisOutput$pollinatorSpecies), FUN = function(curSpecIndex, fullPred, speciesNames) {
  specMat <- fullPred[, grepl(paste("\\s+", curSpecIndex, "\\]$", sep = ""), colnames(fullPred), perl = TRUE)]
  outMat <- t(apply(X = specMat, FUN = function(curRow) {
    rpois(length(curRow), exp(curRow))
  }, MARGIN = 2))
  rownames(outMat) <- paste(speciesNames[curSpecIndex], gsub("^meanPred\\[", "", gsub(",\\s+\\d+\\]$", "", colnames(specMat), perl = TRUE), perl = TRUE), sep = "_")
  outMat
}, fullPred = fullPred, speciesNames = analysisOutput$pollinatorSpecies))
# Create a vector of mean predictions
meanPreds <- as.double(sapply(X = 1:length(analysisOutput$pollinatorSpecies), FUN = function(curSpecIndex, fullPred, speciesNames) {
  specMat <- fullPred[, grepl(paste("\\s+", curSpecIndex, "\\]$", sep = ""), colnames(fullPred), perl = TRUE)]
  outMat <- exp(apply(X = specMat, FUN = median, MARGIN = 2))
  names(outMat) <- paste("location", gsub("^meanPred\\[", "", gsub(",\\s+\\d+\\]$", "", colnames(specMat), perl = TRUE), perl = TRUE), sep = "_")
  outMat
}, fullPred = fullPred, speciesNames = analysisOutput$pollinatorSpecies))
# Create a vector of species counts
speciesCounts <- as.integer(as.matrix(analysisOutput$analysisData[, analysisOutput$pollinatorSpecies]))
# Create the DHARMa object
dharmOb <- createDHARMa(predMatrix, speciesCounts, meanPreds, integerResponse = T)
# Plot the DHARMa object
pdf(file = paste(workspaceLocation, "DHARMaPlot.pdf", sep = "/"), width = 14, height = 8)
plot(dharmOb)
dev.off()

# Function to create an MCMC plot of the regression coefficients for a species
createMCMCPlot <- function(speciesIndex, speciesNames, covNames, mcmcOutput, workspaceLocation) {
  curSpecies <- speciesNames[speciesIndex]
  # Combine the samples into one giant data frame
  allSamplesFrame <- cbind(do.call(rbind, lapply(X = as.list(mcmcOutput), FUN = function(curElement) {
    as.data.frame(curElement)
  })), data.frame(
    mcmcChain = rep(1:length(mcmcOutput), sapply(X = mcmcOutput, FUN = nrow)),
    mcmcTime = unlist(lapply(X = mcmcOutput, FUN = function(curSamples) { 1:nrow(curSamples) }))
  ))
  # Create a set of coefficient labels
  coeffLabels <- gsub("_", " ", gsub("resourceDensity_", "Resource density ", covNames, fixed = TRUE), fixed = TRUE)
  coeffLabels[coeffLabels == "airTemperature"] <- "Air temperature"
  coeffLabels[coeffLabels == "windSpeed"] <- "Wind speed"
  # Retrieve a data frame of coefficient values
  coeffFrame <- do.call(rbind, lapply(X = 1:length(covNames), FUN = function(covIndex, speciesIndex, allSamplesFrame, covNames, coeffLabels) {
    data.frame(
      coefficientType = factor(rep(covNames[covIndex], nrow(allSamplesFrame)), levels = covNames, labels = coeffLabels),
      coefficientValues = allSamplesFrame[, paste("regCoeffs[", covIndex, ", ", speciesIndex, "]", sep = "")],
      mcmcChain = allSamplesFrame$mcmcChain,
      mcmcTime = allSamplesFrame$mcmcTime
    )
  }, speciesIndex = speciesIndex, allSamplesFrame = allSamplesFrame, covNames = covNames, coeffLabels = coeffLabels))
  # Retrieve a data frame of interaction coefficients
  interactionFrame <- do.call(rbind, lapply(X = 1:length(speciesNames), FUN = function(compIndex, speciesIndex, allSamplesFrame, speciesNames) {
    outFrame <- NULL
    if(compIndex == speciesIndex) {
      outFrame <- data.frame(
        coefficientType = factor(rep("Variance", nrow(allSamplesFrame)), levels = "Variance", labels = "Observation variance"),
        coefficientValues = allSamplesFrame[, paste("interactionMatrix[", compIndex, ", ", speciesIndex, "]", sep = "")],
        mcmcChain = allSamplesFrame$mcmcChain,
        mcmcTime = allSamplesFrame$mcmcTime
      )
    } else {
      outFrame <- data.frame(
        coefficientType = factor(rep(paste("Covariance", speciesNames, sep = "_"), nrow(allSamplesFrame)), levels = paste("Covariance", speciesNames, sep = "_"),
          labels = paste("Interaction ", gsub("_", " ", speciesNames, fixed = TRUE), sep = " ")),
        coefficientValues = allSamplesFrame[, paste("interactionMatrix[", compIndex, ", ", speciesIndex, "]", sep = "")],
        mcmcChain = allSamplesFrame$mcmcChain,
        mcmcTime = allSamplesFrame$mcmcTime
      )
    }
    outFrame
  }, speciesIndex = speciesIndex, allSamplesFrame = allSamplesFrame, speciesNames = speciesNames))
  # Integrate the entire data frame of coefficients
  fullCoeffFrame <- rbind(coeffFrame, interactionFrame, data.frame(
    coefficientType = factor(rep(c("Intercept", "PantrapEffect"), rep(nrow(allSamplesFrame), 2)), levels = c("Intercept", "PantrapEffect"), labels = c("Intercept", "Pantrap effect")),
    coefficientValues = c(allSamplesFrame[, paste("interceptCoeff[", speciesIndex, "]", sep = "")], allSamplesFrame[, paste("pantrapCoeff[", speciesIndex, "]", sep = "")]),
    mcmcChain = allSamplesFrame$mcmcChain,
    mcmcTime = allSamplesFrame$mcmcTime
  ))
  fullCoeffFrame$mcmcChain <- as.factor(paste("Chain", fullCoeffFrame$mcmcChain, sep = " "))
  outPlot <- ggplot(fullCoeffFrame, aes(mcmcTime, coefficientValues)) + geom_line(aes(colour = mcmcChain)) + facet_wrap(coefficientType ~ ., scales = "free", nrow = length(levels(fullCoeffFrame$coefficientType))) +
    theme_classic() + xlab("Time") + ylab("Parameter value") + theme(legend.title = element_blank())
  ggsave(paste(workspaceLocation, "/mcmcOutput_", curSpecies, ".pdf", sep = ""), plot = outPlot, width = 11, height = 2 * length(levels(fullCoeffFrame$coefficientType)), limitsize = FALSE)
}
# Produce an MCMC plot for each of the pollinator species
lapply(X = 1:length(analysisOutput$pollinatorSpecies), FUN = createMCMCPlot, speciesNames = analysisOutput$pollinatorSpecies,
  covNames = analysisOutput$covariateNames, mcmcOutput = analysisOutput$mcmcOutput$samples, workspaceLocation = workspaceLocation)

# Make a padded bounding box around the point data
locCoords <- st_coordinates(st_sf(analysisOutput$analysisData))
locXRange <- range(locCoords[, "X"], na.rm = TRUE)
locYRange <- range(locCoords[, "Y"], na.rm = TRUE)
paddingPercent <- 10 / 100
borderPolygon <- st_sfc(st_polygon(list(matrix(c(
  locXRange[1] - diff(locXRange) * paddingPercent, locYRange[1] - diff(locYRange) * paddingPercent,
  locXRange[1] - diff(locXRange) * paddingPercent, locYRange[2] + diff(locYRange) * paddingPercent,
  locXRange[2] + diff(locXRange) * paddingPercent, locYRange[2] + diff(locYRange) * paddingPercent,
  locXRange[2] + diff(locXRange) * paddingPercent, locYRange[1] - diff(locYRange) * paddingPercent,
  locXRange[1] - diff(locXRange) * paddingPercent, locYRange[1] - diff(locYRange) * paddingPercent
), ncol = 2, byrow = TRUE))), crs = st_crs(st_sf(analysisOutput$analysisData)))
# Intersect to create an Oslo fjord coastline data
osloFjordPolygon <- st_intersection(borderPolygon, coastlineData)
beehivePointsLocal <- beehivePoints[borderPolygon, ]
makeCompetitionMap <- function(curSpecies, compSpecies, speciesList, analysisFrame, fullPred, mcmcOutput, referencePolygons, beehivePoints, workspaceLocation) {
  # Retrieve the indeces in the species names vector that represent the target species and the competative species
  targetIndex <- which(speciesList == curSpecies)
  compIndex <- which(speciesList == compSpecies)
  # Retrieve the pantrap fixed effects
  expPantrapTarget <- mean(unlist(lapply(X = mcmcOutput, FUN = function(curElement, targetIndex) {
    curElement[, paste("pantrapCoeff[", targetIndex, "]", sep = "")]
  }, targetIndex = targetIndex)))
  expPantrapComp <- mean(unlist(lapply(X = mcmcOutput, FUN = function(curElement, compIndex) {
    curElement[, paste("pantrapCoeff[", compIndex, "]", sep = "")]
  }, compIndex = compIndex)))
  # Retrieve the expected density of the target species (correcting for the pantrap fixed effect)
  expDensityTarget <- apply(X = as.matrix(fullPred[, grepl(paste("^meanPred\\[\\d+, ", targetIndex, "\\]$", sep = ""), colnames(fullPred), perl = TRUE)]), FUN = mean, MARGIN = 2) -
    ifelse(as.character(analysisFrame$sampleType) == "pantrap", expPantrapTarget, 0.0)
  # Retrieve the expected density of the competative species (correcting for the pantrap fixed effect)
  expDensityComp <- apply(X = as.matrix(fullPred[, grepl(paste("^meanPred\\[\\d+, ", compIndex, "\\]$", sep = ""), colnames(fullPred), perl = TRUE)]), FUN = mean, MARGIN = 2) -
    ifelse(as.character(analysisFrame$sampleType) == "pantrap", expPantrapComp, 0.0)
  # Retrieve the expected interaction strength
  expInteractionStrength <- mean(unlist(lapply(X = mcmcOutput, FUN = function(curElement, targetIndex, compIndex) {
    curElement[, paste("interactionMatrix[", targetIndex, ", ", compIndex, "]", sep = "")]
  }, targetIndex = targetIndex, compIndex = compIndex)))
  # Aggregate the competition predictions at each location
  compPredFrame <- st_sf(do.call(rbind, lapply(X = unique(analysisFrame$trapID), FUN = function(curTrap, curSpecies, predsFrame) {
    outFrame <- data.frame(
      trapID = curTrap,
      expComp = mean(predsFrame[predsFrame$trapID == curTrap, "expComp"]),
      geometry = predsFrame[predsFrame$trapID == curTrap, "geometry"][1]
    )
    outFrame
  }, predsFrame = data.frame(
    trapID = analysisFrame$trapID,
    expComp = expDensityComp * expInteractionStrength,
    geometry = analysisFrame$geometry
  ))))
  # Aggregate the target density prediction at each location
  targetPredFrame <- st_sf(do.call(rbind, lapply(X = unique(analysisFrame$trapID), FUN = function(curTrap, curSpecies, predsFrame) {
    outFrame <- data.frame(
      trapID = curTrap,
      expTarget = mean(predsFrame[predsFrame$trapID == curTrap, "expTarget"]),
      geometry = predsFrame[predsFrame$trapID == curTrap, "geometry"][1]
    )
    outFrame
  }, predsFrame = data.frame(
    trapID = analysisFrame$trapID,
    expTarget = expDensityTarget,
    geometry = analysisFrame$geometry
  ))))
  # Create a plot of the competition
  compPlot <- ggplot(compPredFrame) + geom_sf(data = osloFjordPolygon) + geom_sf(aes(size = abs(expComp), fill = expComp), shape = 21) +
    geom_sf(data = beehivePoints, shape = 24, size = 2, fill = "yellow", colour = "black") +
    scale_fill_gradient(low = ifelse(expInteractionStrength > 0.0, "white", "red"), high = ifelse(expInteractionStrength > 0.0, "blue", "white"), na.value = "white") +
    guides(size = "none") +
    xlab("") + ylab("") + theme_minimal() + theme(legend.title = element_blank())
  # Create a plot of the density of the target species
  targetPlot <- ggplot(targetPredFrame) + geom_sf(data = osloFjordPolygon) + geom_sf(aes(size = expTarget, fill = exp(expTarget)), shape = 21) +
    geom_sf(data = beehivePoints, shape = 24, size = 2, fill = "yellow", colour = "black") +
    scale_fill_gradient(low = "white", high = "purple", na.value = "white") +
    guides(size = "none") +
    xlab("") + ylab("") + theme_minimal() + theme(legend.title = element_blank())
  ggsave(paste(workspaceLocation, "/competitionMap_", curSpecies, ".pdf", sep = ""), plot = compPlot, width = 11, height = 8)
  ggsave(paste(workspaceLocation, "/densityMap_", curSpecies, ".pdf", sep = ""), plot = targetPlot, width = 11, height = 8)
}
# Produce competition maps for each of the species
lapply(X = analysisOutput$pollinatorSpecies, FUN = makeCompetitionMap, compSpecies = "Apis_mellifera", speciesList = analysisOutput$pollinatorSpecies, fullPred = fullPred,
  analysisFrame = analysisOutput$analysisData, mcmcOutput = analysisOutput$mcmcOutput$samples, referencePolygons = osloFjordPolygon, beehivePoints = beehivePointsLocal, workspaceLocation = workspaceLocation)
