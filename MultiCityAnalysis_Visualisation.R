library(coda)
library(nimble)
library(ggplot2)
library(DHARMa)
library(condMVNorm)
library(parallel)
library(MASS)

# The location of the analysis output
urbanSISRepository <- file.path(Sys.getenv("WORKSPACE_URBANSIS"), "WP2 - Mapping - 15885002")
# outLocation <- file.path(urbanSISRepository, "MultiCityAnalysis_Output", "mcmcOutput.rds")
outputLocation <- file.path("C:/Users/joseph.chipperfield/OneDrive - NINA/Work/UrbanBees", "MultiCityAnalysis_Output")
# Import the analysis outputs
analysisOutput <- readRDS(file.path(outputLocation, "mcmcOutput.rds"))
# Import the processed pantrap data
inputData <- readRDS(file.path(urbanSISRepository, "MultiCity_ProcessedPantrapData.rds"))

# Function to create an interaction plot between the species
createInteractionPlot <- function(sampledInteractionMatrix) {
  # Function to create a summary data frame from the sampled interaction matrix
  createSummaryFrame <- function(sampledInteractionMatrix, apFun) {
    do.call(rbind, lapply(X = dimnames(sampledInteractionMatrix)[[1]], FUN = function(curRowIndex, sampledInteractionMatrix, apFun) {
      do.call(rbind, lapply(X = dimnames(sampledInteractionMatrix)[[2]], FUN = function(curColIndex, curRowIndex, sampledInteractionMatrix, apFun) {
        sendHoneyBeeToLast <- function(inNames) {
          c(sort(inNames[inNames != "Apis.mellifera"]), "Apis.mellifera")
        }
        fromLevelNames <- sendHoneyBeeToLast(dimnames(sampledInteractionMatrix)[[1]])
        toLevelNames <- sendHoneyBeeToLast(dimnames(sampledInteractionMatrix)[[2]])
        data.frame(
          fromSpec = factor(curRowIndex, fromLevelNames, labels = gsub(".", " ", fromLevelNames, fixed = TRUE)),
          toSpec = factor(curColIndex, toLevelNames, labels = gsub(".", " ", toLevelNames, fixed = TRUE)),
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
  probPosInteraction$isSignif <- ifelse(probPosInteraction$value > 0.60 | probPosInteraction$value < 0.40, TRUE, FALSE)
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
    geom_point(aes(size = absValue, fill = posVal), data = effectSizeInteraction[!is.na(effectSizeInteraction$value), ], shape = 21) +
    scale_fill_gradient2(midpoint = 0.5, na.value = "white") +
    xlab("") + ylab("") + theme_minimal() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), axis.text.y = element_text(face = "italic"))
}
interactionPlot <- createInteractionPlot(analysisOutput$sampledInteractionMatrix)
ggsave(file = file.path(outputLocation, "interactionMatrix.pdf"), plot = interactionPlot, width = 16, height = 11)

# Order the bee species with honeybee first
firstSpecies <- c("Apis mellifera")
specOrder <- c(firstSpecies, inputData$speciesNames[!(inputData$speciesNames %in% firstSpecies)])
specOrder <- specOrder[gsub(" ", ".", specOrder, fixed = TRUE) %in% colnames(analysisOutput$sampledSpeciesMeans)]
specOrder <- specOrder[1:4]

# Function to plot diversity estimates with changes in the numbers of different bee species
diversityEstimation <- function(xSpecies, sampledSpeciesMeans, sampledInteractionMatrix, pantrapData, outputLocation) {
  cat("\tProcessing ", xSpecies, "...", sep = "")
  # Calculate evenness according to the Evar metric described in https://www.jstor.org/stable/pdf/3545749.pdf
  eVarCalc <- function(inVals) {
    # Remove zero-valued entries
    aboveZero <- inVals[inVals > 0.0]
    outVal <- NA
    if(length(aboveZero) > 0) {
      # Use the arctan formula described in the paper
      outVal <- 1.0 - (2.0 / pi) * atan(sum(sapply(X = aboveZero, FUN = function(curVal, aboveZero) {
        (log(curVal) - sum(log(aboveZero)) / length(aboveZero))^2
      }, aboveZero = aboveZero)) / length(aboveZero))
    }
    outVal
  }
  # Calculate the overall means of the species
  specMeans <- apply(X = sampledSpeciesMeans, FUN = function(curMat, cityTotals) {
    cityWeightAv <- apply(X = curMat, FUN = function(curSpec, cityTotals) {
      sum(curSpec * cityTotals) / sum(cityTotals)
    }, MARGIN = 2, cityTotals = cityTotals)
    mean(cityWeightAv)
  }, MARGIN = 2, cityTotals = table(pantrapData$cityID))
  # Convert the species names to ones that can be used to extract interaction terms
  xSpeciesConv <- gsub(" ", ".", xSpecies, fixed = TRUE)
  # Get the range of values to test for the species
  xSpeciesVals <- seq(log(0.5), log(max(pantrapData[, xSpeciesConv] * 2.0) + 0.5), length.out = 30)
  # Function to sample from conditional distribution of species abundances given a certain value for a species
  sampCondDist <- function(curXVal, xSpeciesConv, specMeans, sampledInteractionMatrix) {
    # Retrieve the index of the relvant species
    specIndex <- which(names(specMeans) == xSpeciesConv)[1]
    # Calculate a mean vector for the conditional distribution
    newMeans <- specMeans[-specIndex] + (sampledInteractionMatrix[-specIndex, specIndex] / sampledInteractionMatrix[specIndex, specIndex]) * (curXVal - specMeans[specIndex])
    # Calculate a new variance-covariance matrix for the conditional distribution
    newVarCov <- sampledInteractionMatrix[-specIndex, -specIndex] - ((sampledInteractionMatrix[-specIndex, specIndex, drop = FALSE] / sampledInteractionMatrix[specIndex, specIndex]) %*% sampledInteractionMatrix[specIndex, -specIndex, drop = FALSE])
    # Sample from the conditional distribution
    outExp <- rep(exp(curXVal), length(specMeans))
    outExp[-specIndex] <- exp(mvrnorm(1, newMeans, newVarCov))
    setNames(outExp, names(specMeans))
  }
  cat(" processing richness...")
  # Iterate over each of the x-values and calculate
  richValues <- do.call(rbind, lapply(X = xSpeciesVals, FUN = function(curXVal, xSpeciesConv, specMeans, sampledInteractionMatrix) {
    # Create species means for each of the interaction matrix slices
    sampledSpecMeans <- apply(X = sampledInteractionMatrix, FUN = function(curMat, curXVal, xSpeciesConv, specMeans) {
      sampCondDist(curXVal, xSpeciesConv, specMeans, curMat)
    }, curXVal = curXVal, xSpeciesConv = xSpeciesConv, specMeans = specMeans, MARGIN = 3)
    # Create the occurrence probabilities for each species
    occProbs <- apply(X = apply(X = sampledSpecMeans, FUN = function(curSamp) {
      1.0 - dpois(0, curSamp)
    }, MARGIN = 2), MARGIN = 1, FUN = mean)
    # Remove the target species from the occurrence probability calculation
    occProbs <- occProbs[names(specMeans) != xSpeciesConv]
    # Calculate the richness probabilities
    richProbs <- rep(1.0, length(occProbs) + 1)
    for(curIndex in 1:length(occProbs)) {
      curProb <- occProbs[curIndex]
      failProbs <- rep(0.0, curIndex + 1)
      passProbs <- rep(0.0, curIndex + 1)
      for(inIndex in 1:curIndex) {
        failProbs[inIndex] <- richProbs[inIndex] * (1.0 - curProb)
        passProbs[inIndex + 1] <- richProbs[inIndex] * curProb
      }
      richProbs[1:(curIndex + 1)] <- failProbs + passProbs
    }
    # Create an output data frame for plotting
    data.frame(
      richnessProbs = richProbs,
      richnessVals = 0:length(occProbs),
      logSpecAbund = rep(curXVal, length(occProbs) + 1),
      expRichness = rep(sum(0:length(occProbs) * richProbs), length(occProbs) + 1)
    )
  }, specMeans = specMeans, xSpeciesConv = xSpeciesConv, sampledInteractionMatrix = sampledInteractionMatrix))
  cat(" processing evenness...")
  # Create samples of evenness from stuff
  evenValues <- do.call(rbind, lapply(X = xSpeciesVals, FUN = function(curXVal, xSpeciesConv, specMeans, sampledInteractionMatrix) {
    # Create species means for each of the interaction matrix slices
    sampledEvenness <- as.double(apply(X = sampledInteractionMatrix, FUN = function(curMat, curXVal, xSpeciesConv, specMeans) {
      outExp <- sampCondDist(curXVal, xSpeciesConv, specMeans, curMat)
      replicate(1000, eVarCalc(rpois(length(outExp) - 1, outExp[names(outExp) != xSpeciesConv])), simplify = "array")
    }, curXVal = curXVal, xSpeciesConv = xSpeciesConv, specMeans = specMeans, MARGIN = 3))
    # Create a data frame for output plotting
    data.frame(
      evennessVals = sampledEvenness,
      logSpecAbund = rep(curXVal, length(sampledEvenness)),
      expEvenness = rep(mean(sampledEvenness, na.rm = TRUE), length(sampledEvenness))
    )
  }, xSpeciesConv = xSpeciesConv, specMeans = specMeans, sampledInteractionMatrix = sampledInteractionMatrix))
  # Create a data frame of pantrap richness/evenness
  pantrapRichEven <- t(apply(X = as.matrix(pantrapData[, names(specMeans)[names(specMeans) != xSpeciesConv]]), FUN = function(curRow) {
    outVec <- setNames(c(sum(curRow > 0), eVarCalc(curRow)), c("richness", "evenness"))
    outVec
  }, MARGIN = 1))
  # Append the abundance and city information
  pantrapRichEven <- cbind(data.frame(
    specAbund = pantrapData[, xSpeciesConv],
    city = pantrapData$cityID
  ), as.data.frame(pantrapRichEven))
  # Find richness values that constrain the plotting the region to only the more probable values
  richLims <- c(-0.5, max(richValues[richValues$richnessProbs > 0.0001, "richnessVals"], pantrapRichEven$richness) + 0.5)
  # Function to condense the data frames for neater plotting
  condenseDataFrame <- function(inData) {
    # Get the condensed data frame
    condensedFrame <- unique(inData)
    # Count the number of rows that have the same values
    rowCounts <- apply(X = as.matrix(condensedFrame), FUN = function(curRow, inData) {
      sum(apply(X = as.matrix(inData), FUN = function(compRow, curRow) {
        all(compRow == curRow)
      }, curRow = curRow, MARGIN = 1))
    }, inData = inData, MARGIN = 1)
    cbind(condensedFrame, data.frame(rowCount = rowCounts))
  }
  # Condense the richness and evenness data frames for easier plotting
  condensedRichness <- condenseDataFrame(pantrapRichEven[, names(pantrapRichEven) != "evenness"])
  condensedEvenness <- condenseDataFrame(pantrapRichEven[, names(pantrapRichEven) != "richness"])
  # Set some variables controlling the plotting offsets
  zeroLine <- min(richValues$logSpecAbund) - 0.5 # Position to put the zero count data
  cityRichOffset <- 0.3   # The gap between plotting the city symbols (richness)
  cityEvenOffset <- 0.1   # The gap between plotting the city symbols (evenness)
  textOffset <- 0.08  # The horizontal gap between the city symbols and the text
  # Create a plot for the expected richness and the richness probabilities
  richPlot <- ggplot(data = richValues, mapping = aes(x = logSpecAbund, y = richnessVals)) + geom_tile(mapping = aes(fill = richnessProbs)) +
    geom_line(mapping = aes(x = logSpecAbund, y = expRichness)) + ylim(richLims) +
    geom_point(mapping = aes(x = ifelse(specAbund == 0, zeroLine, log(specAbund)), y = richness + (as.integer(city) - 2) * cityRichOffset, shape = city), data = condensedRichness) + 
    geom_text(aes(x = ifelse(specAbund == 0, zeroLine, log(specAbund)) + textOffset, y = richness + (as.integer(city) - 2) * cityRichOffset, label = as.character(rowCount)), data = condensedRichness, vjust = "center", hjust = "left", size = rel(2)) + 
    geom_vline(xintercept = zeroLine, linetype = "dotted") +
    xlab("ln(Sampled abundance)") + ylab("Richness") + labs(fill = "Probability", shape = "City") + theme_classic() +
    scale_fill_gradient(low = rgb(255, 255, 255, maxColorValue = 255), high = rgb(255, 130, 171, maxColorValue = 255))
  # Create a plot for the expected evenness 
  evenPlot <- ggplot(data = evenValues[!is.na(evenValues$evennessVals), ], mapping = aes(x = logSpecAbund, y = evennessVals)) +
    # stat_summary(fun.min = function(inX) {quantile(inX, probs = 0.025)}, fun.max = function(inX) {quantile(inX, probs = 0.975)}, fun = mean, colour = "palegreen4") +
    # stat_summary(fun.min = min, fun.max = max, geom = "ribbon", fill = "palegreen4", alpha = 0.2) +
    geom_bin_2d(mapping = aes(fill = log2(..count..))) +
    geom_line(mapping = aes(x = logSpecAbund, y = expEvenness)) +
    geom_count(mapping = aes(x = ifelse(specAbund == 0, zeroLine, log(specAbund)) + (as.integer(city) - 2) * cityEvenOffset, y = evenness, shape = city, group = city), data = pantrapRichEven) +
    scale_size(range = c(rel(0.9), rel(3)), breaks = 2^(0:floor(log2(max(condensedEvenness$rowCount, na.rm = TRUE))))) +
    # geom_point(mapping = aes(x = ifelse(specAbund == 0, zeroLine, log(specAbund)), y = evenness + (as.integer(city) - 2) * cityEvenOffset, shape = city), data = condensedEvenness) + 
    # geom_text(aes(x = ifelse(specAbund == 0, zeroLine, log(specAbund)) + textOffset, y = evenness + (as.integer(city) - 2) * cityEvenOffset, label = as.character(rowCount)), data = condensedEvenness, vjust = "center", hjust = "left", size = rel(2)) + 
    geom_vline(xintercept = zeroLine, linetype = "dotted") +
    xlab("ln(Sampled abundance)") + ylab("Evenness") + labs(shape = "City", size = "Count", fill = expression(log[2]*"(Prediction Density)")) + theme_classic() +
    scale_fill_gradient(low = rgb(240, 255, 240, maxColorValue = 255), high = rgb(69, 139, 116, maxColorValue = 255))
  # Make sub-folders for the richness and evenness plots
  richnessFolder <- file.path(outputLocation, "richnessPlots")
  evennessFolder <- file.path(outputLocation, "evennessPlots")
  if(!dir.exists(richnessFolder)) {
    dir.create(richnessFolder)
  }
  if(!dir.exists(evennessFolder)) {
    dir.create(evennessFolder)
  }
  ggsave(file.path(richnessFolder, paste("richness_", xSpeciesConv, ".pdf", sep = "")), richPlot, width = 8, height = 6)
  ggsave(file.path(evennessFolder, paste("evenness_", xSpeciesConv, ".pdf", sep = "")), evenPlot, width = 8, height = 6)
  cat("\n")
  list(
    richnessPlot = richPlot,
    evennessPlot = evenPlot
  )
}
# Run the diversity plots
diversityPlots <- setNames(lapply(X = specOrder, FUN = diversityEstimation, sampledSpeciesMeans = analysisOutput$sampledSpeciesMeans,
  sampledInteractionMatrix = analysisOutput$sampledInteractionMatrix, pantrapData = inputData$pantrapData, outputLocation = outputLocation), specOrder)

# Function to plot the Copula estimation
copulaEstimation <- function(xSpecies, ySpecies, sampledSpeciesMeans, sampledInteractionMatrix, pantrapData, outputLocation) {
  # Calculate the overall means of the species
  specMeans <- apply(X = sampledSpeciesMeans, FUN = function(curMat, cityTotals) {
    cityWeightAv <- apply(X = curMat, FUN = function(curSpec, cityTotals) {
      sum(curSpec * cityTotals) / sum(cityTotals)
    }, MARGIN = 2, cityTotals = cityTotals)
    mean(cityWeightAv)
  }, MARGIN = 2, cityTotals = table(pantrapData$cityID))
  # Convert the species names to ones that can be used to extract interaction terms
  xSpeciesConv <- gsub(" ", ".", xSpecies, fixed = TRUE)
  ySpeciesConv <- gsub(" ", ".", ySpecies, fixed = TRUE)
  # Get the range of values to test for the species
  xSpeciesVals <- seq(log(0.5), log(max(pantrapData[, xSpeciesConv]) * 1.5), length.out = 30)
  ySpeciesVals <- seq(log(0.5), log(max(pantrapData[, ySpeciesConv]) * 1.5), length.out = 50)
  # Iterate over each combination of the species values
  copulaDensity <- sapply(X = xSpeciesVals, FUN = function(xCurVal, ySpeciesVals, specMeans, xSpeciesConv, ySpeciesConv, sampledInteractionMatrix) {
    sapply(X = ySpeciesVals, FUN = function(yCurVal, xCurVal, specMeans, xSpeciesConv, ySpeciesConv, sampledInteractionMatrix) {
      # Calculate the indeces of the species being calculated
      xIndex <- which(names(specMeans) == xSpeciesConv)[1]
      yIndex <- which(names(specMeans) == ySpeciesConv)[1]
      # Calculate the conditional density of the species value pair
      condDens <- mean(apply(X = sampledInteractionMatrix, FUN = function(tempCovMat, specMeans, yIndex, xIndex, xCurVal, yCurVal) {
        # When calculating marginal densities of multivariate distributions you can simply calculate the density from a reduced multivariate
        # normal distribution with the non-normal components removed.  Proof is here: https://web.archive.org/web/20100117200722/http://fourier.eng.hmc.edu/e161/lectures/gaussianprocess/node7.html
        dcmvnorm(yCurVal, mean = specMeans[c(xIndex, yIndex)], sigma = tempCovMat[c(xIndex, yIndex), c(xIndex, yIndex)], dependent.ind = 2, given.ind = 1, X.given = xCurVal, check.sigma = FALSE)
        # Calculate the conditioning variables
        # tempGiven <- specMeans
        # tempGiven[xIndex] <- xCurVal
        # Calculate the conditional density of the multivariate normal distribution
        # dcmvnorm(yCurVal, mean = specMeans, sigma = tempCovMat, dependent.ind = yIndex, given.ind = (1:length(specMeans))[-yIndex], X.given = tempGiven[-yIndex], log = TRUE, check.sigma = FALSE)
      }, specMeans = specMeans, yIndex = yIndex, xIndex = xIndex, xCurVal = xCurVal, yCurVal = yCurVal, MARGIN = 3))
      condDens
    }, xCurVal = xCurVal, specMeans = specMeans, xSpeciesConv = xSpeciesConv, ySpeciesConv = ySpeciesConv, sampledInteractionMatrix = sampledInteractionMatrix)
  }, ySpeciesVals = ySpeciesVals, specMeans = specMeans, xSpeciesConv = xSpeciesConv, ySpeciesConv = ySpeciesConv,
    sampledInteractionMatrix = sampledInteractionMatrix)
  # Reorder the density frame for plotting
  densFrame <- cbind(expand.grid(yVals = ySpeciesVals, xVals = xSpeciesVals), data.frame(density = log10(as.numeric(copulaDensity))))
  xZeroVal <- min(xSpeciesVals) - diff(range(xSpeciesVals)) * 0.06
  yZeroVal <- min(ySpeciesVals) - diff(range(ySpeciesVals)) * 0.10
  # Create a data frame containing the pantrap data for the species
  pantrapFrame <- data.frame(
    xSpecies = ifelse(pantrapData[, xSpeciesConv] <= 0, xZeroVal, log(pantrapData[, xSpeciesConv])),
    ySpecies = ifelse(pantrapData[, ySpeciesConv] <= 0, yZeroVal, log(pantrapData[, ySpeciesConv])),
    cityID = pantrapData$cityID
  )
  # Condense the pantrap frame so there is not too much overplotting
  condensedPantrap <- unique(pantrapFrame)
  condensedPantrap <- cbind(condensedPantrap, data.frame(
    count = sapply(X = 1:nrow(condensedPantrap), FUN = function(curIndex, condensedPantrap, pantrapFrame) {
      sum(condensedPantrap[curIndex, 1] == pantrapFrame[, 1] & condensedPantrap[curIndex, 2] == pantrapFrame[, 2] & condensedPantrap[curIndex, 3] == pantrapFrame[, 3])
    }, condensedPantrap = condensedPantrap, pantrapFrame = pantrapFrame)
  ))
  # Plot the copula along with the pantrap data
  outplot <- ggplot(densFrame, aes(x = xVals, y = yVals)) + geom_tile(aes(fill = density)) +
    geom_point(aes(x = xSpecies, y = ySpecies + (as.integer(cityID) - 2) * 0.025 * diff(range(ySpeciesVals)), shape = cityID), data = condensedPantrap, size = 2.5) + 
    geom_text(aes(x = xSpecies + 0.05, y = ySpecies + (as.integer(cityID) - 2) * 0.025 * diff(range(ySpeciesVals)), label = as.character(count)), data = condensedPantrap, vjust = "center", hjust = "left", size = rel(2)) + 
    geom_vline(xintercept = xZeroVal, linetype = "dotted") + geom_hline(yintercept = yZeroVal, linetype = "dotted") +
    xlab(bquote(paste("ln(", italic(.(xSpecies)), ")", sep = ""))) + ylab(bquote(paste("ln(", italic(.(ySpecies)), ")", sep = ""))) + labs(fill = expression(log[10]*"(Density)"), shape = "City") + theme_classic() +
    scale_fill_gradient(low = rgb(255, 255, 255, maxColorValue = 255), high = rgb(255, 130, 171, maxColorValue = 255), na.value = rgb(255, 255, 255, maxColorValue = 255))
  ggsave(file = file.path(outputLocation, paste("copulaPlot_", xSpeciesConv, "_", ySpeciesConv, ".pdf", sep = "")), plot = outplot, width = 11, height = 7)
  outplot
}
# Create a folder to store the copula plots
copulaLocation <- file.path(outputLocation, "copulaPlots")
if(!dir.exists(copulaLocation)) {
  dir.create(copulaLocation)
}
# Go over every combination of species interaction and calculate the copula plots
copulaCluster <- makeCluster(detectCores())
copulaPlots <- setNames(parLapply(cl = copulaCluster, X = specOrder, fun = function(xSpecies, speciesNames, sampledSpeciesMeans, sampledInteractionMatrix, pantrapData, outputLocation, copulaEstimation) {
  library(condMVNorm)
  library(ggplot2)
  setNames(lapply(X = speciesNames[speciesNames != xSpecies], FUN = function(ySpecies, xSpecies, sampledSpeciesMeans, sampledInteractionMatrix, pantrapData, outputLocation) {
    cat("Processing copula between ", xSpecies, " and ", ySpecies, "...\n", sep = "")
    copulaEstimation(xSpecies, ySpecies, sampledSpeciesMeans, sampledInteractionMatrix, pantrapData, outputLocation)
  }, xSpecies = xSpecies, sampledSpeciesMeans = sampledSpeciesMeans, sampledInteractionMatrix = sampledInteractionMatrix, pantrapData = pantrapData, outputLocation = outputLocation),
    gsub(" ", ".", speciesNames[speciesNames != xSpecies], fixed = TRUE))
}, speciesNames = specOrder, sampledSpeciesMeans = analysisOutput$sampledSpeciesMeans, sampledInteractionMatrix = analysisOutput$sampledInteractionMatrix,
  pantrapData = inputData$pantrapData, outputLocation = copulaLocation, copulaEstimation = copulaEstimation), gsub(" ", ".", specOrder, fixed = TRUE))
stopCluster(copulaCluster)

saveRDS(
  list(
    interactionPlot = interactionPlot,
    richnessPlots = setNames(lapply(X = diversityPlots, FUN = function(curPlot) {curPlot$richnessPlot}), gsub(" ", ".", specOrder, fixed = TRUE)),
    evennessPlots = setNames(lapply(X = diversityPlots, FUN = function(curPlot) {curPlot$evennessPlot}), gsub(" ", ".", specOrder, fixed = TRUE)),
    copulaPlots = copulaPlots
  ),
  file = file.path(outputLocation, "plotObjects.rds"))
