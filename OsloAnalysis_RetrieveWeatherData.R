retrieveNorwegianWeatherStationData <- function(spatialExtent, timeRange, weatherElements, frostUserID, frostUserSecret) {
  # Test to ensure that the spatial extent is a correctly initialised sf object
  if(!any(class(spatialExtent) == "sfc_POLYGON")) {
    stop("spatial extent must be a feature geometry column of polygon type")
  }
  # Import the time range data
  inTimeRange <- tryCatch(range(as.POSIXct(timeRange), na.rm = TRUE), error = function(err) {
    stop("error importing input time range: ", err)
  })
  # Import the weather elements
  inWeatherElements <- tryCatch(as.character(weatherElements), error = function(err) {
    stop("error importing the weather elements: ", err)
  })
  inWeatherElements <- paste(inWeatherElements, collapse = "%2C")
  # Import the user ID
  inFrostID <- tryCatch(as.character(frostUserID), error = function(err) {
    stop("error importing the FROST user ID: ", err)
  })
  # Import the secret
  inFrostSecret <- tryCatch(as.character(frostUserSecret), error = function(err) {
    stop("error importing the FROST secret: ", err)
  })
  # Transform the input spatial features to the same coordinate reference system as used by the FROST API
  frostCRS <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  transBoundary <- st_transform(spatialExtent, crs = frostCRS)
  boundaryText <- gsub("\\s+", "\\%20", st_as_text(transBoundary), perl = TRUE)
  # Create the time range string
  timeString <- gsub("\\s+", "%20", paste(format(inTimeRange, "%Y-%m-%dT%H:%M:%S"), collapse = "/"), perl = TRUE)
  # Initialise the authentication
  authApp <- oauth_app("FROST", key = inFrostID, secret = inFrostSecret)
  authEnd <- oauth_endpoint(access = "https://frost.met.no/auth/accessToken", authorize = NULL)
  authToken <- oauth2.0_token(endpoint = authEnd, app = authApp, client_credentials = TRUE)
  # Retrieve the source data from the FROST API (only getting those sources with information from some part of the date range in
  # the sampled data)
  sourceData <- GET(url = paste("https://frost.met.no/sources/v0.jsonld?",
    "geometry=", boundaryText, "&",
    "types=SensorSystem&",
    "validtime=", gsub("T\\d\\d:\\d\\d:\\d\\d", "", timeString, perl = TRUE), # Source API has slightly different time range interface than observations
  sep = ""),  config(token = authToken))
  # Initialise an output list
  outList <- list(locationData = NULL, observationData = NULL)
  if(sourceData$status_code == 200) {
    # Import the location data
    outList$locationData <- fromJSON(content(sourceData, as = "text", content = "application/json", encoding = "UTF-8"))$data
    colnames(outList$locationData) <- gsub("@", "", colnames(outList$locationData), fixed = TRUE)
    # Reorder the geometry information
    outList$locationData <- st_as_sf(cbind(outList$locationData[, colnames(outList$locationData) != "geometry"],
      as.data.frame(matrix(unlist(outList$locationData$geometry$coordinates), byrow = TRUE, ncol = 2, dimnames = list(NULL, c("geomX", "geomY"))))
    ), coords = c("geomX", "geomY"), crs = frostCRS)
    # For some reason the inclusion of certain source ID causes an internal error on the FROST API (even though they appear to be
    # valid IDs in other tables).  As a result, we define a function that goes through the observation data one source ID at a time
    # and then removes those results that give an error
    processSource <- function(curSource, timeStr, elementStr, authToken) {
      message("Processing observation data for source ", curSource, "...", appendLF = FALSE)
      # Retrieve the observation data from the FROST API for the current source
      obsData <- GET(url = paste("https://frost.met.no/observations/v0.jsonld?",
        "sources=", curSource, "&",
        "referencetime=", timeStr, "&",
        "elements=", elementStr, sep = ""), config(token = authToken))
      # Put a delay on the processing to avoid overwhelming the API and getting locked out
      Sys.sleep(5)
      outFrame <- NULL
      if(obsData$status_code == 200) {
        # Reformat the information from the source into a table
        sourceTable <- fromJSON(content(obsData, as = "text", type = "application/json", encoding = "UTF-8"))$data
        outFrame <- do.call(rbind, lapply(X = 1:nrow(sourceTable), FUN = function(curIndex, sourceTable) {
          # Retrieve the current observation data
          curObservations <- sourceTable$observations[[curIndex]]
          # Rearrange the data into a data frame
          data.frame(
            id = rep(strsplit(sourceTable$sourceId[curIndex], ":", fixed = TRUE)[[1]][1], nrow(curObservations)),
            measurementId = 1:nrow(curObservations),
            arrayNum = rep(strsplit(sourceTable$sourceId[curIndex], ":", fixed = TRUE)[[1]][2], nrow(curObservations)),
            timeStamp = rep(strptime(
              gsub("\\..*$", "", gsub("T", " ", sourceTable$referenceTime[curIndex], fixed = TRUE), perl = TRUE),
              "%Y-%m-%d %H:%M:%S"), nrow(curObservations)),
            elementId = curObservations$elementId,
            value = curObservations$value,
            unit = curObservations$unit,
            levelType = curObservations$level$levelType,
            levelUnit = curObservations$level$unit,
            levelValue = curObservations$level$value,
            timeOffset = curObservations$timeOffset,
            timeResolution = curObservations$timeResolution,
            timeSeriesId = curObservations$timeSeriesId,
            performanceCategory = curObservations$performanceCategory,
            exposureCategory = curObservations$exposureCategory,
            qualityCode = curObservations$qualityCode
          )
        }, sourceTable = sourceTable))
        message(" completed")
      } else if(obsData$status_code == 404 || obsData$status_code == 412) {
        message(" completed (no data found)")
      } else {
        message(" failed with exit status code ", obsData$status_code)
      }
      outFrame
    }
    # Create the observation data frame
    outList$observationData <- do.call(rbind, lapply(X = outList$locationData$id, FUN = processSource,
      timeStr = gsub("/", "%2F", gsub(":", "%3A", timeString, fixed = TRUE), fixed = TRUE),
      elementStr = inWeatherElements, authToken = authToken))
    # Retrieve those IDs that are actually attached to observation data
    hasObsData <- unique(outList$observationData$id)
    # Remove those element from the location data that don't have observation data attached
    outList$locationData <- outList$locationData[outList$locationData$id %in% hasObsData, ]
    # Convert the location data to the coordinate reference system in the boundary geometry
    outList$locationData <- st_transform(outList$locationData, crs = st_crs(spatialExtent))
    rownames(outList$locationData) <- outList$locationData$id
    rownames(outList$observationData) <- paste(outList$observationData$id, outList$observationData$arrayNum,
      as.integer(outList$observationData$timeStamp), outList$observationData$measurementId, sep = "_")
  } else if(sourceData$staus_code == 404) {
    warning("no location data found that matches query")
  } else {
    stop("error status (code ", sourceData$staus_code, ") revealed when returning location data")
  }
  outList
}