library(openxlsx)
library(sf)

# Set the location of the pan trap and pollinator collection data
urbanSISRepository <- paste(Sys.getenv("WORKSPACE_URBANSIS"), "WP2 - Mapping - 15885002/Field DATA final", sep = "/")
pantrapDataLocation <- paste(urbanSISRepository, "Pantrap_data.xlsx", sep = "/")
pollinatorTransectDataLocation <- paste(urbanSISRepository, "Pollinators_Transect.xlsx", sep = "/")
pollinatorPantrapDataLocation <- paste(urbanSISRepository, "Pantrap_bee_ID_FO_for_analysis.xlsx", sep = "/")
hoverflyPantrapDataLocation <- paste(urbanSISRepository, "Pantrap_hoverflies_Jon_for_analysis.xlsx", sep = "/")
foodResourceDataLocation <- paste(urbanSISRepository, "Food_resources.xlsx", sep = "/")

# Set the output location for the processed data
outputLocation <- paste(urbanSISRepository, "ProcessedPantrapAndTransectData.rds", sep = "/")

# Retrieve the raw data from the respective workbooks
rawPantrapData <- readWorkbook(pantrapDataLocation, sheet = 1, startRow = 2)
rawPollinatorTransectData <- readWorkbook(pollinatorTransectDataLocation, sheet = 1)
rawPollinatorPantrapData <- as.data.frame(t(readWorkbook(pollinatorPantrapDataLocation, sheet = 1, rowNames = TRUE, colNames = FALSE)))
rawHoverflyPantrapData <- readWorkbook(hoverflyPantrapDataLocation, sheet = 1)
rawFoodResourceData <- readWorkbook(foodResourceDataLocation, sheet = 1)

# Function to tidy the species names of the transect and pantrap data
tidyPollinatorNames <- function(inSpecNames) {
  outNames <- gsub("\\s+male$", "", inSpecNames, perl = TRUE)
  outNames <- gsub("\\s+female$", "", outNames, perl = TRUE)
  outNames <- gsub("\\s+Coll\\.$", "", outNames, perl = TRUE)
  outNames <- gsub("[\\.\\s]+", "_", outNames, perl = TRUE)
  outNames <- gsub("_+$", "", outNames, perl = TRUE)
  outNames <- paste(toupper(substring(outNames, 1, 1)), substring(outNames, 2), sep = "")
  # There are a couple of spelling mistakes of species names in the dataset
  outNames <- ifelse(outNames == "Bombus_pascourum", "Bombus_pascuorum", outNames)
  outNames <- ifelse(outNames == "Apis_melifera", "Apis_mellifera", outNames)
  # Remove some unidentified samples from the dataset
  outNames <- ifelse(outNames == "Solitary_bee_not_collected", NA, outNames)
  ifelse(is.na(inSpecNames), NA, outNames)
}

# Process the pantrap location data
processedPantrapData <- data.frame(
  trapID = rawPantrapData[["Trap_number"]],
  description = rawPantrapData[["Location.description"]],
  trapXCoord = rawPantrapData[["Trap_DD_x"]],
  trapYCoord = rawPantrapData[["Trap_DD_y"]],
  firstDeploymentDate = convertToDateTime(rawPantrapData[["1_Date"]] + rawPantrapData[["1_Time"]]),
  firstDeploymentTemp = rawPantrapData[["1_Temperature"]],
  firstDeploymentSkycover = rawPantrapData[["1_Skycover"]],
  firstDeploymentWind = rawPantrapData[["1_Wind"]],
  firstDeploymentComment = rawPantrapData[["1_Comments"]],
  firstRetrievalDate = convertToDateTime(rawPantrapData[["2_Date"]] + rawPantrapData[["2_Time"]]),
  firstRetrievalTemp = rawPantrapData[["2_Temperature"]],
  firstRetrievalSkycover = rawPantrapData[["2_Skycover"]],
  firstRetrievalWind = rawPantrapData[["2_Wind"]],
  firstRetrievalComment = rawPantrapData[["2_Comments"]],
  secondDeploymentDate = convertToDateTime(rawPantrapData[["3_Date"]] + rawPantrapData[["3_Time"]]),
  secondDeploymentTemp = rawPantrapData[["3_Temperature"]],
  secondDeploymentSkycover = rawPantrapData[["3_Skycover"]],
  secondDeploymentWind = rawPantrapData[["3_Wind"]],
  secondDeploymentComment = rawPantrapData[["3_Comments"]],
  secondRetrievalDate = convertToDate(rawPantrapData[["4_Date"]])
)
rownames(processedPantrapData) <- processedPantrapData$trapID

# Retrieve the unique species names for the flower resources
rawFoodResourceData$Flowering_species_latin <- gsub("-", " ", rawFoodResourceData$Flowering_species_latin, perl = TRUE)
floweringSpeciesList <- unique(rawFoodResourceData$Flowering_species_latin)
floweringSpeciesList <- sort(floweringSpeciesList[!is.na(as.character(floweringSpeciesList)) & as.character(floweringSpeciesList) != "unknown"])
# Process the food resources data
processedFoodResourceData <- do.call(rbind, lapply(X = unique(rawFoodResourceData$Trap_number), FUN = function(curID, rawFoodResourceData, speciesList) {
  # Function to tidy some of the cover records
  coverTidy <- function(curVal) {
    as.double(ifelse(grepl("^\\s*missing\\s*$", as.character(curVal), perl = TRUE), rep(NA, length(curVal)), curVal))
  }
  # Retrieve the resource data associated with the current trap
  curResourceData <- rawFoodResourceData[rawFoodResourceData$Trap_number == curID, ]
  # Create an output frame with the general transect information
  outFrame <- cbind(data.frame(
    trapID = curID,
    transectDate = convertToDate(curResourceData$F_Date)[1],
    transectObserverOne = as.character(curResourceData$F_Observer)[1],
    transectObserverTwo = as.character(curResourceData$F_Observer2)[1],
    grassCover = mean(coverTidy(curResourceData$Grass), na.rm = TRUE)[1],
    shrubCover = mean(coverTidy(curResourceData$Shrub), na.rm = TRUE)[1],
    treeCover = mean(coverTidy(curResourceData$Tree), na.rm = TRUE)[1],
    urbanCover = mean(coverTidy(curResourceData$Urban), na.rm = TRUE)[1],
    waterCover = mean(coverTidy(curResourceData$Water_edge), na.rm = TRUE)[1],
    transectWidth = mean(curResourceData$Transect_width_cm, na.rm = TRUE)[1] / 100.0,
    transectLength = mean(curResourceData$Transect_length_m, na.rm = TRUE)[1],
    transectComments = as.character(curResourceData$kommentar)[1]
  ), as.data.frame(
    # Add columns for the count type
    setNames(as.list(rep(NA, length(speciesList))), paste("countType", gsub("\\s+", "_", speciesList, perl = TRUE), sep = "_"))
  ),as.data.frame(
    # Add columns for the total resource availability
    setNames(as.list(rep(0.0, length(speciesList))), paste("resourceTotal", gsub("\\s+", "_", speciesList, perl = TRUE), sep = "_"))
  ), as.data.frame(
    # Add columns for the resource density
    setNames(as.list(rep(0.0, length(speciesList))), paste("resourceDensity", gsub("\\s+", "_", speciesList, perl = TRUE), sep = "_"))
  ))
  # Create a data frame of composite food resources
  curSpeciesList <- unique(as.character(curResourceData$Flowering_species_latin))
  curSpeciesList <- curSpeciesList[!is.na(curSpeciesList) & curSpeciesList %in% speciesList]
  if(length(curSpeciesList) > 0) {
    # Create a data frame with one row for each species
    totalSpeciesResource <- do.call(rbind, lapply(X = curSpeciesList, FUN = function(curSpecies, curResourceData, curArea) {
      speciesLgl <- as.character(curResourceData$Flowering_species_latin) == curSpecies
      data.frame(
        species = curSpecies,
        countType = paste(curResourceData[speciesLgl, "F_number_dont_use"], "(", curResourceData[speciesLgl, "F_Counted"], ")", sep = "", collapse = " + "),
        resourceTotal = sum(curResourceData[speciesLgl, "F_Final"], na.rm = TRUE),
        resourceDensity = sum(curResourceData[speciesLgl, "F_Final"], na.rm = TRUE) / curArea
      )
    }, curResourceData = curResourceData, curArea = outFrame$transectWidth * outFrame$transectLength))
    # Import the occurrence records and the count type
    outFrame[1, 
             paste("countType", gsub("\\s+", "_", totalSpeciesResource$species, perl = TRUE), sep = "_")
             ] <- totalSpeciesResource$countType
    outFrame[1, c(
      paste("resourceTotal", gsub("\\s+", "_", totalSpeciesResource$species, perl = TRUE), sep = "_"),
      paste("resourceDensity", gsub("\\s+", "_", totalSpeciesResource$species, perl = TRUE), sep = "_")
    )] <- c(
      totalSpeciesResource$resourceTotal,
      totalSpeciesResource$resourceDensity
    )
  }
  outFrame
}, rawFoodResourceData = rawFoodResourceData, speciesList = floweringSpeciesList))
rownames(processedFoodResourceData) <- processedFoodResourceData$trapID
# Merge the location-level data into one big data frame
locationLevelData <- merge(processedPantrapData, processedFoodResourceData, by = "trapID", all = TRUE)
rownames(locationLevelData) <- locationLevelData$trapID

# Fix a couple of errors in the pantrap dataset
# On first deployment, pantrap 162 was forgotten about and not collected
locationLevelData[locationLevelData$trapID == 162, "firstRetrievalDate"] <- locationLevelData[locationLevelData$trapID == 162, "secondDeploymentDate"]

# Retrieve the species names in the data sets
pollinatorPantrapSpeciesList <- tidyPollinatorNames(sort(colnames(rawPollinatorPantrapData)[!(colnames(rawPollinatorPantrapData) %in% c("Trap_number", "Art\\Prøve_nr_old", "Dato"))]))
pollinatorTransectSpeciesList <- unique(as.character(rawPollinatorTransectData$Pollinator_species_latin))
pollinatorTransectSpeciesList <- sort(unique(tidyPollinatorNames(pollinatorTransectSpeciesList[
  !is.na(pollinatorTransectSpeciesList) &
  !(pollinatorTransectSpeciesList %in% c("Solitary bee not collected"))])))
pollinatorSpeciesList <- sort(unique(c(pollinatorPantrapSpeciesList, pollinatorTransectSpeciesList, "Syrphidae")))
# Retrieve the trap IDs in the sample-level data
pollinatorTrapIDs <- sort(unique(c(
  # Retrieve all the trap IDs assocciated 
  as.integer(gsub(" ", "", rawPollinatorPantrapData[, "Trap_number"], fixed = TRUE)),
  rawPollinatorTransectData$Trap_number,
  rawHoverflyPantrapData$Trap_number)))

# Process the pollinator occurrence data (from the pollinator pantrap and transect files - supplemented by the hoverfly datafile)
sampleLevelData <- do.call(rbind, lapply(X = pollinatorTrapIDs, FUN = function(curTrapNum, rawPollinatorPantrapData, rawPollinatorTransectData, rawHoverflyPantrapData, pollinatorSpeciesList, locationLevelData) {
  # Function to tidy the columns in the pantrap data and properly convert them to numeric values
  convColumnToNum <- function(inValues) {
    as.numeric(gsub(" ", "", as.character(inValues), fixed = TRUE))
  }
  # Function to lookup the pantrap collection times from a collection dates
  getPantrapTimes <- function(inDates, curTrapNum, locationLevelData) {
    sapply(X = convertToDate(inDates), FUN = function(curDate, curTrapNum, locationLevelData) {
      outVal <- NA
      curDateFormat <- format(curDate, "%Y-%m-%d")
      # Test whether the patrap was collected on the first or second collection
      isCurTrap <- curTrapNum == locationLevelData$trapID
      isFirstCollection <- isCurTrap & !is.na(locationLevelData$firstRetrievalDate) & curDateFormat == format(locationLevelData$firstRetrievalDate, "%Y-%m-%d")
      isSecondCollection <- isCurTrap & !is.na(locationLevelData$secondRetrievalDate) & curDateFormat == format(locationLevelData$secondRetrievalDate, "%Y-%m-%d")
      if(any(isFirstCollection)) {
        # If this sample represents the first collection then retrieve those dates and times
        outVal <- paste(
          format(locationLevelData[which(isFirstCollection)[1], c("firstDeploymentDate", "firstRetrievalDate")], "%Y-%m-%d %H:%M:%S %Z"),
        collapse = "/")
      } else if(any(isSecondCollection)) {
        # If this sample represents the second collection then retrieve thoose dates and times
        outVal <- paste(
          paste(
            format(locationLevelData[which(isSecondCollection)[1], "secondDeploymentDate"], "%Y-%m-%d %H:%M:%S %Z"),
            format(locationLevelData[which(isSecondCollection)[1], "secondRetrievalDate"], "%Y-%m-%d"),
          sep = "/"),
          # Need to add a time to the second retrieval dates because only the date was recorded
          "12:00:00", "CEST", sep = " ")
      }
      outVal
    }, curTrapNum = curTrapNum, locationLevelData = locationLevelData)
  }
  # Retrieve the pollinator data in the current trap ID
  curPantrapData <- rawPollinatorPantrapData[as.integer(gsub(" ", "", rawPollinatorPantrapData[, "Trap_number"], fixed = TRUE)) == curTrapNum, , drop = FALSE]
  curTransectData <- rawPollinatorTransectData[rawPollinatorTransectData$Trap_number == curTrapNum, ]
  curHoverflyData <- rawHoverflyPantrapData[rawHoverflyPantrapData$Trap_number == curTrapNum, ]
  processedPantrapData <- NULL
  processedTransectData <- NULL
  # Import the pantrap data
  if(nrow(curPantrapData) > 0) {
    processedPantrapData <- cbind(data.frame(
      trapID = convColumnToNum(curPantrapData[, "Trap_number"]),
      samplingTime = getPantrapTimes(convColumnToNum(curPantrapData[, "Dato"]), curTrapNum, locationLevelData),
      sampleType = factor(rep("pantrap", nrow(curPantrapData)), levels = c("pantrap", "transect")),
      samplingTimeInterval = rep(NA, nrow(curPantrapData)),
      effortMeasure = rep(NA, nrow(curPantrapData))
    ),
      as.data.frame(matrix(0, nrow = nrow(curPantrapData), ncol = length(pollinatorSpeciesList), dimnames = list(NULL, pollinatorSpeciesList)))
    )
    # Retrieve the pantrap data
    pantrapSpecData <- apply(X = curPantrapData[, colnames(curPantrapData)[tidyPollinatorNames(colnames(curPantrapData)) %in% pollinatorSpeciesList]], FUN = function(curCol) {
      outVal <- convColumnToNum(curCol)
      ifelse(is.na(outVal), 0, outVal)
    }, MARGIN = 2)
    if(is.null(dim(pantrapSpecData)) || length(dim(pantrapSpecData) < 2)) {
      dim(pantrapSpecData) <- c(1, length(pantrapSpecData))
    }
    colnames(pantrapSpecData) <- tidyPollinatorNames(colnames(pantrapSpecData))
    processedPantrapData[, colnames(pantrapSpecData)] <- pantrapSpecData
    rownames(processedPantrapData) <- paste(processedPantrapData$trapID, processedPantrapData$samplingTime, sep = " ")
  }
  # Import the hoverfly data
  if(nrow(curHoverflyData) > 0) {
    hoverflySamplingDates <- getPantrapTimes(curHoverflyData$Trap_retrival.date, curTrapNum, locationLevelData)
    isInPantrapData <- hoverflySamplingDates %in% processedPantrapData$samplingTime
    # Update the values in the hoverfly column if the pantrap sampling date already exists in the output data frame
    processedPantrapData[
      paste(curTrapNum, hoverflySamplingDates[isInPantrapData], sep = " "),
      "Syrphidae"] <- curHoverflyData$Hoverflies[isInPantrapData]
    # Add a row if the sample date is non in the output data frame
    if(any(!isInPantrapData)) {
      processedPantrapData <- rbind(processedPantrapData, cbind(
        data.frame(
          trapID = rep(curTrapNum, sum(!isInPantrapData)),
          samplingTime = hoverflySamplingDates[isInPantrapData],
          sampleType = factor(rep("pantrap", sum(!isInPantrapData)), levels = c("pantrap", "transect")),
          samplingTimeInterval = rep(NA, sum(!isInPantrapData)),
          effortMeasure = rep(NA, sum(!isInPantrapData))
        ),
        as.data.frame(matrix(
          as.numeric(sapply(X = curHoverflyData$Hoverflies[!isInPantrapData], FUN = function(curVal, isHoverfly) {
            curVal * ifelse(isHoverfly, 1, 0)
          }, isHoverfly = pollinatorSpeciesList == "Syrphidae")), byrow = TRUE, nrow = sum(!isInPantrapData),
            ncol = length(pollinatorSpeciesList), dimnames = list(NULL, pollinatorSpeciesList)
        ))
      ))
    }
  }
  # Import the transect data
  processedTransectData <- do.call(rbind, lapply(X = unique(curTransectData$I_Date), FUN = function(curDate, curTransectData, curTrapNum, pollinatorSpeciesList) {
    # Retrieve the transect data sampled on the current date
    curDateTransect <- curTransectData[curTransectData$I_Date == curDate, ]
    # Get the range of date/time values that
    dateTimeRange <- convertToDateTime(range(c(curDate + as.double(curDateTransect$I_Time_start), curDate + as.double(curDateTransect$I_Time_end)), na.rm = TRUE))
    # The species found on the transect
    foundSpecies <- table(unlist(apply(X = as.matrix(curDateTransect[, c("Pollinator_species_latin", "Flying_nr", "Sitting_nr")]), FUN = function(curRow) {
      rep(tidyPollinatorNames(curRow[1]), sum(as.integer(curRow[2:3]), na.rm = TRUE))
    }, MARGIN = 1)))
    # Initialise an output data frame
    outFrame <- cbind(data.frame(
      trapID = curTrapNum,
      samplingTime = paste(format(dateTimeRange, "%Y-%m-%d %H:%M:%S %Z"), collapse = "/"),
      sampleType = factor("transect", levels = c("pantrap", "transect")),
      samplingTimeInterval = NA,
      effortMeasure = NA
    ), as.data.frame(setNames(as.list(rep(0, length(pollinatorSpeciesList))), pollinatorSpeciesList)))
    # Update the species counts
    outFrame[, names(foundSpecies)] <- foundSpecies
    outFrame
  }, curTransectData = curTransectData, curTrapNum = curTrapNum, pollinatorSpeciesList = pollinatorSpeciesList))
  outFrame <- rbind(processedPantrapData, processedTransectData)
  rownames(outFrame) <- paste(outFrame$trapID, as.character(outFrame$sampleType), "_", sapply(
    X = strsplit(as.character(outFrame$samplingTime), "/", fixed = TRUE), FUN = function(curVal) {
      paste(as.integer(as.POSIXct(curVal[1])), "/", as.integer(as.POSIXct(curVal[2])), sep = "")
  }), sep = "")
  outFrame
}, rawPollinatorPantrapData = rawPollinatorPantrapData, rawPollinatorTransectData = rawPollinatorTransectData,
  rawHoverflyPantrapData = rawHoverflyPantrapData, pollinatorSpeciesList = pollinatorSpeciesList, locationLevelData = locationLevelData))
# Calculate the number of days of the sampling interval (used for the pantrap data)
sampleLevelData$samplingTimeInterval <- sapply(X = strsplit(as.character(sampleLevelData$samplingTime), "/", fixed = TRUE), FUN = function(curSplit) {
    outVal <- NA
    if(length(curSplit) >= 2) {
      outVal <- abs(difftime(as.POSIXct(curSplit[1]), as.POSIXct(curSplit[2]), units = "days"))
    }
    outVal
  })
sampleLevelData$effortMeasure <- log(ifelse(as.character(sampleLevelData$sampleType) == "pantrap",
  sampleLevelData$samplingTimeInterval,
  # Calculate the size of the transect for the transect-based analyses
  apply(X = as.matrix(locationLevelData[as.character(sampleLevelData$trapID), c("transectWidth", "transectLength")]), FUN = prod, MARGIN = 1)))

# Write the location and sample level data
write.csv2(locationLevelData, file = paste(urbanSISRepository, "locationData.csv", sep = "/"), row.names = FALSE, na = "")
write.csv2(sampleLevelData, file = paste(urbanSISRepository, "sampleData.csv", sep = "/"), row.names = FALSE, na = "")
# Save the processed data files
saveRDS(
  list(
    locations = st_as_sf(
      locationLevelData,
      coords = c("trapXCoord", "trapYCoord"),
      crs = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"),
    samples = sampleLevelData,
    pollinatorSpecies = pollinatorSpeciesList,
    flowerSpecies = colnames(locationLevelData)[grepl("^resourceDensity_", colnames(locationLevelData), perl = TRUE)]
  ),
  file = outputLocation)
