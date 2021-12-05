library(openxlsx)

# Set the location of the pan trap and pollinator collection data
urbanSISRepository <- file.path(Sys.getenv("WORKSPACE_URBANSIS"), "WP2 - Mapping - 15885002")

# Retrieve the Aarhus bee pantrap data
rawAarhusData <- read.xlsx(
  file.path(urbanSISRepository, "Denmark_Data", "Aarhus_PanTrap_2016.xlsx"),
  "rekonstructed data")
# Retrieve the Copenhagen bee pantrap data
rawCopenhagenData <- read.xlsx(
  file.path(urbanSISRepository, "Denmark_Data", "Copenhagen_PanTrap_2020.xlsx"),
  "bee data")
# Retrieve the Oslo bee pantrap data
rawOsloData <- read.xlsx(
  file.path(urbanSISRepository, "Field DATA final", "Pantrap_bee_ID_FO_for_analysis.xlsx"),
  "Sheet1")
# Retrieve the Oslo pantrap information
rawOsloPantrapInfo <- read.xlsx(
  file.path(urbanSISRepository, "Field DATA final", "Pantrap_data.xlsx"),
  "Sheet1", startRow = 2
)

# Set the output location for the processed data
outputLocation <- file.path(urbanSISRepository, "MultiCity_ProcessedPantrapData.rds")

# Some species have multiple names/misspellings in the datasets.  Here is a list of synonyms for species that have those
speciesAliases <- list(
  "Apis mellifera" = c("Apis mellifera", "Apis melifera", "Apis melliferra"),
  "Bombus lucorum Coll." = c("Bombus lucorum", "Bombus lucorum Coll.", "Bombus terrestris-kompleks", "Bombus terrestris kompleks", "Bombus terrestris"),
  "Colletes daviesanus" = c("Colletes davesianus", "Colletes daviesanus"),
  "Halictus rubicundus" = c("Halictus rubicundus", "Halictus rubicundus "),
  "Sphecodes geoffrellus" = c("Sphecodes geoffrellus", "Sphecodes geofrellus"),
  "Bombus pascuorum" = c("Bombus_pascourum", "Bombus_pascuorum")
)

# Process the bee data from each of the sites
processedBeeData <- rbind(
  # Process the Aarhus data
  data.frame(
    species = paste(rawAarhusData$Slægt, rawAarhusData$Art, sep = " "),
    count = rawAarhusData$Count,
    siteName = rawAarhusData$Lokalitet,
    cityID = rep("Aarhus", nrow(rawAarhusData)),
    siteID = paste("Aarhus", gsub(" ", "", rawAarhusData$Lokalitet, fixed = TRUE), sep = "_"),
    sampleID = paste("Aarhus", gsub(" ", "", rawAarhusData$Lokalitet, fixed = TRUE), rawAarhusData$Måned, 2016, sep = "_"),
    dayBegin = rep(NA, nrow(rawAarhusData)),
    monthBegin = rawAarhusData$Måned,
    yearBegin = rep(2016, nrow(rawAarhusData)),
    dayEnd = rep(NA, nrow(rawAarhusData)),
    monthEnd = rawAarhusData$Måned,
    yearEnd = rep(2016, nrow(rawAarhusData)),
    xCoord = rep(NA, nrow(rawAarhusData)),
    yCoord = rep(NA, nrow(rawAarhusData))
  ),
  # Process the Copenhagen data
  data.frame(
    species = paste(rawCopenhagenData$Slægt, rawCopenhagenData$Art, sep = " "),
    count = rawCopenhagenData$Antal,
    siteName = paste("Copenhagen", "Garden", rawCopenhagenData$Have, sep = " "),
    cityID = rep("Copenhagen", nrow(rawCopenhagenData)),
    siteID = paste("Copenhagen", rawCopenhagenData$Have, sep = "_"),
    sampleID = paste("Copenhagen", rawCopenhagenData$Have,
      rawCopenhagenData[, which(colnames(rawCopenhagenData) == "Dag")[1]],
      rawCopenhagenData[, which(colnames(rawCopenhagenData) == "Måned")[1]],
      rawCopenhagenData[, which(colnames(rawCopenhagenData) == "År")[1]],
      rawCopenhagenData[, which(colnames(rawCopenhagenData) == "Dag")[2]],
      rawCopenhagenData[, which(colnames(rawCopenhagenData) == "Måned")[2]],
      rawCopenhagenData[, which(colnames(rawCopenhagenData) == "År")[2]], sep = "_"),
    dayBegin = rawCopenhagenData[, which(colnames(rawCopenhagenData) == "Dag")[1]],
    monthBegin = rawCopenhagenData[, which(colnames(rawCopenhagenData) == "Måned")[1]],
    yearBegin = rawCopenhagenData[, which(colnames(rawCopenhagenData) == "År")[1]],
    dayEnd = rawCopenhagenData[, which(colnames(rawCopenhagenData) == "Dag")[2]],
    monthEnd = rawCopenhagenData[, which(colnames(rawCopenhagenData) == "Måned")[2]],
    yearEnd = rawCopenhagenData[, which(colnames(rawCopenhagenData) == "År")[2]],
    xCoord = rep(NA, nrow(rawCopenhagenData)),
    yCoord = rep(NA, nrow(rawCopenhagenData))
  ),
  # Process the Oslo data
  do.call(rbind, apply(X = as.matrix(rawOsloData[3:nrow(rawOsloData), ]), FUN = function(curRow, trapNumbers) {
    data.frame(
      species = rep(curRow[1], length(trapNumbers)),
      count = ifelse(is.na(curRow[2:length(curRow)]), 0, as.integer(curRow[2:length(curRow)])),
      siteName = paste("Oslo Trap", trapNumbers, sep = " "),
      cityID = rep("Oslo", length(trapNumbers)),
      siteID = paste("Oslo", trapNumbers, sep = "_"),
      sampleID = paste("Oslo", trapNumbers, sep = "_"),
      dayBegin = rep(NA, length(trapNumbers)),
      monthBegin = rep(NA, length(trapNumbers)),
      yearBegin = rep(NA, length(trapNumbers)),
      dayEnd = rep(NA, length(trapNumbers)),
      monthEnd = rep(NA, length(trapNumbers)),
      yearEnd = rep(NA, length(trapNumbers)),
      xCoord = rep(NA, length(trapNumbers)),
      yCoord = rep(NA, length(trapNumbers))
    )
  }, MARGIN = 1, trapNumbers = as.integer(colnames(rawOsloData)[2:ncol(rawOsloData)])))
)
processedBeeData$cityID <- factor(processedBeeData$cityID, levels = c("Aarhus", "Copenhagen", "Oslo"))
processedBeeData$siteID <- factor(processedBeeData$siteID)
processedBeeData$sampleID <- factor(processedBeeData$sampleID)
# Harmonise the species names in the dataset
processedBeeData$species <- sapply(X = as.character(processedBeeData$species), FUN = function(curSpecies, speciesAliases) {
  # Test to see if the species name is in the alias list
  isAnAlias <- sapply(X = speciesAliases, FUN = function(curAlias, curSpecies) {
    curSpecies %in% curAlias
  }, curSpecies = curSpecies)
  outName <- curSpecies
  if(any(isAnAlias)) {
    outName <- names(speciesAliases)[isAnAlias]
  }
  outName
}, speciesAliases = speciesAliases)
processedBeeData$species <- factor(processedBeeData$species, levels = sort(unique(processedBeeData$species)))
speciesNames <- levels(processedBeeData$species)
# Reorder the processed data so that each species appears as a column and add in the extra information for the Oslo pan traps
processedBeeData <- do.call(rbind, lapply(X = unique(processedBeeData$sampleID), FUN = function(curSample, processedBeeData, rawOsloPantrapInfo) {
  speciesList <- levels(processedBeeData$species)
  # Get the elements of processed data that are in the currently being processed sample
  curData <- processedBeeData[processedBeeData$sampleID == curSample, ]
  # Get the total number of each species present in the current sample
  specCounts <- setNames(lapply(X = speciesList, FUN = function(curSpecies, curData) {
    # Retrieve the rows of the data set that belong to the current species
    specData <- curData[curData$species == curSpecies, ]
    outVal <- 0
    if(nrow(specData) > 0) {
      outVal <- sum(specData$count, na.rm = TRUE)
    }
    outVal
  }, curData = curData), speciesList)
  # If the trap is from Oslo than append the collection and spatial data
  isOsloTrap <- paste("Oslo", rawOsloPantrapInfo$Trap_number, sep = "_") == curSample
  # Initialise some information about the trap
  trapInfo <- curData[1, c("siteName", "cityID", "siteID", "sampleID", "dayBegin", "monthBegin", "yearBegin", "dayEnd", "monthEnd", "yearEnd", "xCoord", "yCoord")]
  if(any(isOsloTrap)) {
    # Retrieve the collection dates
    trapInfo[1, c("yearBegin", "monthBegin", "dayBegin")] <- as.integer(strsplit(as.character(convertToDate(rawOsloPantrapInfo[isOsloTrap, "1_Date"])), "-", fixed = TRUE)[[1]])
    trapInfo[1, c("yearEnd", "monthEnd", "dayEnd")] <- as.integer(strsplit(as.character(convertToDate(rawOsloPantrapInfo[isOsloTrap, "2_Date"])), "-", fixed = TRUE)[[1]])
    # Retrieve the coordinates of the traps
    trapInfo[1, "xCoord"] <- rawOsloPantrapInfo[isOsloTrap, "Trap_DD_x"]
    trapInfo[1, "yCoord"] <- rawOsloPantrapInfo[isOsloTrap, "Trap_DD_y"]
  }
  cbind(as.data.frame(specCounts), trapInfo)
}, processedBeeData = processedBeeData, rawOsloPantrapInfo = rawOsloPantrapInfo))
# Remove the species with zero counts in any location
speciesToKeep <- sapply(X = gsub(" ", ".", speciesNames, fixed = TRUE), FUN = function(curSpecies, processedBeeData) {
  any(processedBeeData[, curSpecies] > 0)
}, processedBeeData = processedBeeData)
# Save the processed data to the output location
saveRDS(list(
  pantrapData = processedBeeData[, c(speciesToKeep, rep(TRUE, ncol(processedBeeData) - length(speciesToKeep)))],
  speciesNames = speciesNames[speciesToKeep]
), file = outputLocation)
