library(sf)
library(sp)
library(tidyverse)

#' Calculate what rows of input layer fall within every stops buffer area
#'
#' @param stops an point layer
#' @param input_layer an polygon layer
#' @param proj numeric ESPG number of projection
#' @param buffer_length an single numeric value (note: the unit is based purely on the projection)
#'
#' @return a point layer based on stops with one extra column contains the results
#' @export
#'
#' @examples
what_within_each_stops <- function(stops, input_layer, proj, buffer_length){
  stops <- st_transform(stops, crs = proj) # reproject stops
  input_layer <- st_transform(input_layer, crs = proj) # reproject input layers
  buffer <- st_buffer(stops, dist = buffer_length) # create buffer
  inter <- st_intersects(buffer,input_layer,sparse = TRUE) # intersect buffer with input layer
  holder <- data.frame(1:dim(inter)[1]) # create holder to store results
  for (i in 1:dim(inter)[1]) { # write results into holder (loop over intersect result)
    if (length(unlist(inter[i]))>0) {
      holder[i,1] <- paste(as.character(unlist(inter[i])), collapse = ",") # gsub(" ", "", substr(paste(inter[i]),start = 3, stop = str_length(paste(inter[i]))-1), fixed = TRUE) # concatenate order number of input layer rows into one string separate with comma and remove white space
    } else {
      holder[i,1] <- NA
    }
  }
  output <- as.data.frame(stops) # create output variable based on stops
  output[,ncol(output) + 1] <- holder[1] # insert holder as one column to the output
  colnames(output)[ncol(output)] <- "input_rows_within" # rename the inserted column
  output <- st_as_sf(output, sf_column_name = "geometry", crs = proj) # project and transform output as sf object
  return(output) # return output
}

#' Convert "input_rows_within" from string to numeric
#'
#' @param row a single numeric value of the stop row
#' @param column the character of the column name that stores the intersect result; if using the result from the previous function, then it should be "input_rows_within"
#' @param dataframe the name of the dataframe that store the result
#'
#' @return the numeric row number
#' @export
#'
#' @examples
str_to_num <- function(row, column, dataframe){
  output <- as.numeric(unlist(str_split(st_drop_geometry(dataframe[row,c(column)]),","))) # convert string to numeric by splitting with comma
  return(output) # return output
}

#' Retrieve rows from input layer with what_within_each_stops result
#'
#' @param row a single numeric value of the stop row
#' @param column the character of the column name that stores the intersect result; if using the result from the previous function, then it should be "input_rows_within"
#' @param inter_stops the name of the sf object and dataframe that store the intersect result
#' @param input_layer an polygon layer which you intersect with stops layer
#'
#' @return rows from input_layer based on the numeric row number just generated
#' @export
#'
#' @examples
return_input_layer_rows <- function(row, column, inter_stops, input_layer){
  output <- as.numeric(unlist(str_split(st_drop_geometry(inter_stops[row,c(column)]),","))) # convert string to numeric by splitting with comma
  return(input_layer[str_to_list(row,column,inter_stops),]) # return input_layer rows based on the numeric row number just generated
}

#' Calculate median demographic value that are nearest to stops
#'
#' @param stops a point layer
#' @param input_layer an polygon layer
#' @param column the character name of the column that you want to calculate in input_layer
#' @param proj numeric ESPG number of projection
#'
#' @return a point layer based on stops with one extra column contains the results
#' @export
#'
#' @examples
nearest_median_value <- function(stops, input_layer, column, proj){
  stops <- st_transform(stops, crs = proj) # reproject stops layer
  input_layer <- st_transform(input_layer, crs = proj) # reproject input layer
  cent <- st_centroid(input_layer) # get centroid of input layer polygons
  dist <- as.data.frame(st_distance(cent,stops)) # calculate the distance of centroids to stops and store them into a dataframe
  output <- as.data.frame(stops) # create output based on stops
  initialcol <- ncol(output) # store the number of output columns
  for (i in 1:nrow(stops)) { # loop over the stops
    colnumber <- paste("V",i,sep = "") # generate column name for each stop
    output[i, initialcol + 1] <- st_drop_geometry(cent[which(dist[,c(colnumber)]==min(dist[,c(colnumber)])),c(column)])[1,1]  # grasp the stop's smallest distance row and create a new column to store it (use [1,1] to ensure return 1 median value when have two centroid have the same distance to one stop)
  }
  colnames(output)[initialcol + 1] <- column <- column # rename the variable to the name in the input layer
  output <- st_as_sf(output, sf_column_name = "geometry", crs = proj) # project and transform output as sf object
  return(output) # return output
}

#' Calculate average demographic value that are within certain buffer of stops
#'
#' @param stops a point layer
#' @param input_layer an polygon layer
#' @param proj numeric ESPG number of projection
#' @param buffer_length an single numeric value (note: the unit is based purely on the projection)
#' @param column the character name of the column that you want to calculate in input_layer
#' @param NA_omit TRUE or FALSE; whether omit NA value in the column that you want to calculate
#'
#' @return a point layer based on stops with one extra column contains the results
#' @export
#'
#' @examples
average_value_in_buffer <- function(stops, input_layer, proj, buffer_length, column, NA_omit = FALSE){
  stops <- st_transform(stops, crs = proj) # reproject stops layer
  input_layer <- st_transform(input_layer, crs = proj) # reproject input layer
  input_layer <- as.data.frame(input_layer) # coerce sf object into dataframe only
  input_layer[,c("input_layer_ID")] <- c(1:nrow(input_layer)) # insert ID column
  input_layer <- st_as_sf(input_layer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  input_layer[,c("input_layer_area")] <- as.numeric(st_area(input_layer)) # calculate polygon areas
  
  stops_buffer <- st_buffer(stops, dist = buffer_length) # create buffer of stops layer
  stops_buffer <- as.data.frame(stops_buffer) # coerce sf object into dataframe only
  stops_buffer[,c("stops_buffer_ID")] <- c(1:nrow(stops_buffer)) # insert ID column
  stops_buffer <- st_as_sf(stops_buffer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  stops_buffer[,c("stops_buffer_area")] <- as.numeric(st_area(stops_buffer)) # calculate buffer areas (should be the same value for every row)
  
  inter <- st_intersection(stops_buffer,input_layer) # intersect the buffer and input_layer
  inter[,c("inter_area")] <- as.numeric(st_area(inter)) # calculate intersect areas
  
  output <- as.data.frame(stops) # create output based on stops layer
  colnumber <- ncol(output) # calculate the number of columns of output
  frame_inter <- as.data.frame(inter) # coerce sf object into dataframe only
  if (NA_omit) { # if NA_omit equals TRUE than omit column with NA values
    frame_inter <- frame_inter[which(!is.na(frame_inter[,c(column)])),]
  }
  for (i in 1:nrow(stops)) { # loop over each stop row
    total_value_inside <- 0 # create holder for total value inside each buffer
    total_area_within <- 0 # create holder for total intersect area within each buffer
    for (x in 1:dim(inter[which(inter$stops_buffer_ID==i),0])[1]) { # loop over each stop row's intersected area
      average_value <- frame_inter[which(frame_inter$stops_buffer_ID==i),c(column)][x] # retrieve average value in the intersected area
      area_within <- frame_inter[which(frame_inter$stops_buffer_ID==i),c("inter_area")][x] # retrieve area size in the intersected area
      total_value_inside <- total_value_inside + average_value*area_within # add this intersected area's total value to total value holder
      total_area_within <- total_area_within + area_within # add this intersected area's size to total area within each buffer holder
    }
    output[i, colnumber + 1] <- total_value_inside/total_area_within # calculate the average value in the buffer area (only uses intersect area in case buffer is at edges that are not fully covered by input layer)
  }
  colnames(output)[colnumber + 1] <- column # rename the result column to the original name in the input layer
  output <- st_as_sf(output, sf_column_name = "geometry", crs = proj) # project and transform output as sf object
  return(output) # return output
}

#' Calculate total demographic value that are within certain buffer of stops
#'
#' @param stops a point layer
#' @param input_layer an polygon layer
#' @param proj numeric ESPG number of projection
#' @param buffer_length an single numeric value (note: the unit is based purely on the projection)
#' @param column the character name of the column that you want to calculate in input_layer
#' @param NA_omit TRUE or FALSE; whether omit NA value in the column that you want to calculate
#'
#' @return a polygon layer based on stops buffer with extra columns contain the results
#' @export
#'
#' @examples
total_value_in_buffer <- function(stops, input_layer, proj, buffer_length, column, NA_omit = FALSE){
  stops <- st_transform(stops, crs = proj) # reproject stops layer
  input_layer <- st_transform(input_layer, crs = proj) # reproject input layer
  input_layer <- as.data.frame(input_layer) # coerce sf object into dataframe only
  input_layer[,c("input_layer_ID")] <- c(1:nrow(input_layer)) # insert ID column
  input_layer <- st_as_sf(input_layer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  input_layer[,c("input_layer_area")] <- as.numeric(st_area(input_layer)) # calculate polygon areas
  
  stops_buffer <- st_buffer(stops, dist = buffer_length) # create buffer of stops layer
  stops_buffer <- as.data.frame(stops_buffer) # coerce sf object into dataframe only
  stops_buffer[,c("stops_buffer_ID")] <- c(1:nrow(stops_buffer)) # insert ID column
  stops_buffer <- st_as_sf(stops_buffer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  stops_buffer[,c("stops_buffer_area")] <- as.numeric(st_area(stops_buffer)) # calculate buffer areas (should be the same value for every row)
  
  inter <- st_intersection(stops_buffer,input_layer) # intersect the buffer and input_layer
  inter[,c("inter_area")] <- as.numeric(st_area(inter)) # calculate intersect areas
  
  output <- as.data.frame(stops_buffer) # create output based on stops layer
  colnumber <- ncol(output) # calculate the number of columns of output
  frame_inter <- as.data.frame(inter) # coerce sf object into dataframe only
  if (NA_omit) { # if NA_omit equals TRUE than omit column with NA values
    frame_inter <- frame_inter[which(!is.na(frame_inter[,c(column)])),]
  }
  for (i in 1:nrow(stops)) { # loop over each stop row
    total_value_inside <- 0 # create holder for total value inside each buffer
    for (x in 1:dim(inter[which(inter$stops_buffer_ID==i),0])[1]) { # loop over each stop row's intersected area
      total_value <- frame_inter[which(frame_inter$stops_buffer_ID==i),c(column)][x] # retrieve the total value in the intersected area
      area_within <- frame_inter[which(frame_inter$stops_buffer_ID==i),c("inter_area")][x] # retrieve the size of the intersected area
      input_layer_area <- frame_inter[which(frame_inter$stops_buffer_ID==i),c("input_layer_area")][x] # retrieve the size of the input layer polygon
      total_value_inside <- total_value_inside + total_value*(area_within/input_layer_area) # calculate this stop buffer's value based on the area size and add this value to the total value (if the buffer is bigger than the input layer than the total value is underestimated)
    }
    output[i, colnumber + 1] <- total_value_inside # write each stop's total value to output
  }
  colnames(output)[colnumber + 1] <- column # rename the result column to the original name in the input layer
  output <- st_as_sf(output, sf_column_name = "geometry", crs = proj) # project and transform output as sf object
  return(output) # return output
}

#' Calculate percentage of area type that are within certain buffer of stops
#'
#' @param stops a point layer
#' @param input_layer an polygon layer
#' @param proj numeric ESPG number of projection
#' @param buffer_length an single numeric value (note: the unit is based purely on the projection)
#' @param column the character name of the column that you want to calculate in input_layer
#' @param type the character name of the type in the column that you want to calculate in input_layer
#' @param NA_omit TRUE or FALSE; whether omit NA value in the column that you want to calculate
#'
#' @return a polygon layer based on stops buffer with extra columns contain the results
#' @export
#'
#' @examples
area_pct_in_buffer <- function(stops, input_layer, proj, buffer_length, column, type, NA_omit = FALSE){
  stops <- st_transform(stops, crs = proj) # reproject stops layer
  input_layer <- st_transform(input_layer, crs = proj) # reproject input layer
  input_layer <- as.data.frame(input_layer) # coerce sf object into dataframe only
  input_layer[,c("input_layer_ID")] <- c(1:nrow(input_layer)) # insert ID column
  input_layer <- st_as_sf(input_layer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  input_layer[,c("input_layer_area")] <- as.numeric(st_area(input_layer)) # calculate polygon areas
  
  stops_buffer <- st_buffer(stops, dist = buffer_length) # create buffer of stops layer
  stops_buffer <- as.data.frame(stops_buffer) # coerce sf object into dataframe only
  stops_buffer[,c("stops_buffer_ID")] <- c(1:nrow(stops_buffer)) # insert ID column
  stops_buffer <- st_as_sf(stops_buffer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  stops_buffer[,c("stops_buffer_area")] <- as.numeric(st_area(stops_buffer)) # calculate buffer areas (should be the same value for every row)
  
  inter <- st_intersection(stops_buffer,input_layer) # intersect the buffer and input_layer
  inter[,c("inter_area")] <- as.numeric(st_area(inter)) # calculate intersect areas
  
  output <- as.data.frame(stops_buffer) # create output based on stops layer
  colnumber <- ncol(output) # calculate the number of columns of output
  frame_inter <- as.data.frame(inter) # coerce sf object into dataframe only
  if (NA_omit) { # if NA_omit equals TRUE than omit column with NA values
    frame_inter <- frame_inter[which(!is.na(frame_inter[,c(column)])),]
  }
  for (i in 1:nrow(stops)) { # loop over each stop row
    total_area_inside <- 0 # create holder for total area column_type inside each buffer
    if (dim(inter[which(frame_inter$stops_buffer_ID==i & frame_inter[,c(column)] == type),0])[1]!=0){ # check whether there are column type area in the buffer
      for (x in 1:dim(inter[which(frame_inter$stops_buffer_ID == i & frame_inter[,c(column)] == type),0])[1]) { # loop over each stop row's intersected area
        area_within <- frame_inter[which(frame_inter$stops_buffer_ID == i & frame_inter[,c(column)] == type),c("inter_area")][x] # retrieve the size of the intersected area with the input column type
        total_area_inside <- total_area_inside + area_within # record to total area column_type inside each buffer
      }
      buffer_area <- frame_inter[which(frame_inter$stops_buffer_ID==i),c("stops_buffer_area")][1] # retrieve the size of the input layer polygon (all are the same, therefore only retrieve the first one)
      output[i, colnumber + 1] <- total_area_inside*100/buffer_area # write each stop's total value to output
    } else {
      output[i, colnumber + 1] <- 0 # if there is no column type area inside, record 0 to it
    }
  }
  colnames(output)[colnumber + 1] <- column # rename the result column to the original name in the input layer
  output <- st_as_sf(output, sf_column_name = "geometry", crs = proj) # project and transform output as sf object
  return(output) # return output
}

#' Calculate area size of certain type that are within certain buffer of stops
#'
#' @param stops a point layer
#' @param input_layer an polygon layer
#' @param proj numeric ESPG number of projection
#' @param buffer_length an single numeric value (note: the unit is based purely on the projection)
#' @param column the character name of the column that you want to calculate in input_layer
#' @param type the character name of the type in the column that you want to calculate in input_layer
#' @param NA_omit TRUE or FALSE; whether omit NA value in the column that you want to calculate
#'
#' @return a polygon layer based on stops buffer with extra columns contain the results
#' @export
#'
#' @examples
area_in_buffer <- function(stops, input_layer, proj, buffer_length, column, type, NA_omit = FALSE){
  stops <- st_transform(stops, crs = proj) # reproject stops layer
  input_layer <- st_transform(input_layer, crs = proj) # reproject input layer
  input_layer <- as.data.frame(input_layer) # coerce sf object into dataframe only
  input_layer[,c("input_layer_ID")] <- c(1:nrow(input_layer)) # insert ID column
  input_layer <- st_as_sf(input_layer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  input_layer[,c("input_layer_area")] <- as.numeric(st_area(input_layer)) # calculate polygon areas
  
  stops_buffer <- st_buffer(stops, dist = buffer_length) # create buffer of stops layer
  stops_buffer <- as.data.frame(stops_buffer) # coerce sf object into dataframe only
  stops_buffer[,c("stops_buffer_ID")] <- c(1:nrow(stops_buffer)) # insert ID column
  stops_buffer <- st_as_sf(stops_buffer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  stops_buffer[,c("stops_buffer_area")] <- as.numeric(st_area(stops_buffer)) # calculate buffer areas (should be the same value for every row)
  
  inter <- st_intersection(stops_buffer,input_layer) # intersect the buffer and input_layer
  inter[,c("inter_area")] <- as.numeric(st_area(inter)) # calculate intersect areas
  
  output <- as.data.frame(stops_buffer) # create output based on stops layer
  colnumber <- ncol(output) # calculate the number of columns of output
  frame_inter <- as.data.frame(inter) # coerce sf object into dataframe only
  if (NA_omit) { # if NA_omit equals TRUE than omit column with NA values
    frame_inter <- frame_inter[which(!is.na(frame_inter[,c(column)])),]
  }
  for (i in 1:nrow(stops)) { # loop over each stop row
    total_area_inside <- 0 # create holder for total area column_type inside each buffer
    if (dim(inter[which(frame_inter$stops_buffer_ID==i & frame_inter[,c(column)] == type),0])[1]!=0){ # check whether there are column type area in the buffer
      for (x in 1:dim(inter[which(frame_inter$stops_buffer_ID == i & frame_inter[,c(column)] == type),0])[1]) { # loop over each stop row's intersected area
        area_within <- frame_inter[which(frame_inter$stops_buffer_ID == i & frame_inter[,c(column)] == type),c("inter_area")][x] # retrieve the size of the intersected area with the input column type
        total_area_inside <- total_area_inside + area_within # record to total area column_type inside each buffer
      }
      output[i, colnumber + 1] <- total_area_inside # write each stop's total value to output
    } else {
      output[i, colnumber + 1] <- 0 # if there is no column type area inside, record 0 to it
    }
  }
  colnames(output)[colnumber + 1] <- column # rename the result column to the original name in the input layer
  output <- st_as_sf(output, sf_column_name = "geometry", crs = proj) # project and transform output as sf object
  return(output) # return output
}

#' Calculate length of input layer that are within certain buffer of stops
#'
#' @param stops a point layer
#' @param input_layer an polyline layer
#' @param proj numeric ESPG number of projection
#' @param buffer_length an single numeric value (note: the unit is based purely on the projection)
#'
#' @return a polygon layer based on stops buffer with extra columns contain the results
#' @export
#'
#' @examples
length_in_buffer <- function(stops, input_layer, proj, buffer_length){
  stops <- st_transform(stops, crs = proj) # reproject stops layer
  input_layer <- st_transform(input_layer, crs = proj) # reproject input layer
  input_layer <- as.data.frame(input_layer) # coerce sf object into dataframe only
  input_layer[,c("input_layer_ID")] <- c(1:nrow(input_layer)) # insert ID column
  input_layer <- st_as_sf(input_layer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  input_layer[,c("input_layer_length")] <- as.numeric(st_length(input_layer)) # calculate polyline length
  
  stops_buffer <- st_buffer(stops, dist = buffer_length) # create buffer of stops layer
  stops_buffer <- as.data.frame(stops_buffer) # coerce sf object into dataframe only
  stops_buffer[,c("stops_buffer_ID")] <- c(1:nrow(stops_buffer)) # insert ID column
  stops_buffer <- st_as_sf(stops_buffer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  stops_buffer[,c("stops_buffer_area")] <- as.numeric(st_area(stops_buffer)) # calculate buffer areas (should be the same value for every row)
  
  inter <- st_intersection(stops_buffer,input_layer) # intersect the buffer and input_layer
  inter[,c("inter_length")] <- as.numeric(st_length(inter)) # calculate intersect polyline length
  
  output <- as.data.frame(stops_buffer) # create output based on stops layer
  colnumber <- ncol(output) # calculate the number of columns of output
  frame_inter <- as.data.frame(inter) # coerce sf object into dataframe only
  for (i in 1:nrow(stops)) { # loop over each stop row
    total_length_inside <- 0 # create holder for total length inside each buffer
    if (dim(inter[which(frame_inter$stops_buffer_ID==i),0])[1]!=0){ # check whether there are intersected polyline in the buffer
      for (x in 1:dim(inter[which(frame_inter$stops_buffer_ID == i),0])[1]) { # loop over each stop row's intersected area
        length_within <- frame_inter[which(frame_inter$stops_buffer_ID == i),c("inter_length")][x] # retrieve the size of the intersected length
        total_length_inside <- total_length_inside + length_within # record to total length inside each buffer
      }
      output[i, colnumber + 1] <- total_length_inside # write each stop's total value to output
    } else {
      output[i, colnumber + 1] <- 0 # if there is no length inside, record 0 to it
    }
  }
  colnames(output)[colnumber + 1] <- "total_length_within_buffer"
  output <- st_as_sf(output, sf_column_name = "geometry", crs = proj) # project and transform output as sf object
  return(output) # return output
}

#' Calculate land use entropy for a certain buffer area
#' Below are the formula used in this function to calculate entropy value
#' $$Entropy = -\sum_{k=1}^nP_k*\frac{ln(P_k)}{ln(n)}$$
#'
#' @param stops a point layer
#' @param input_layer an polygon layer
#' @param proj numeric ESPG number of projection
#' @param buffer_length an single numeric value (note: the unit is based purely on the projection)
#' @param column character column name that contains the land use categories
#' @param exclude_intermediate whether to exclude intermediate land use area percentage results in the output
#'
#' @return a point layer based on stops containing the entropy value for each stop
#' @export
#'
#' @examples
calculate_entropy <- function(stops, input_layer, proj, buffer_length, column, exclude_intermediate = TRUE){
  stops <- st_transform(stops, crs = proj) # reproject stops layer
  input_layer <- st_transform(input_layer, crs = proj) # reproject input layer
  holder <- stops # create a holder to contain intermediate results
  holder <- as.data.frame(holder) # convert sf object to dataframe only
  counter <- 1 # create a counter to track number of categories
  colnumber <- ncol(holder) # record the original column number of holder
  types <- rownames(as.data.frame(summary(as.factor(as.data.frame(input_layer)[,c(column)])))) # record the number of categories in the input layer
  for (i in types) { # loop over each category
    holder[,colnumber+counter] <- as.data.frame(area_pct_in_buffer(stops, input_layer, proj, buffer_length, column,i, FALSE))[,c(column)] # calculate the area percentage for each category
    colnames(holder)[colnumber+counter] <- paste(column,i,sep = "_") # rename the output column name to distinguish each category's result
    counter <- counter + 1 # increment counter value
  }
  holder$entropy <- NA # create new column to store entropy value
  for (x in 1:nrow(holder)) { # loop over every row of holder
    entropy1 <- 0 # create holder to store each rows individual entropy value
    for (i in types) { # loop over each category
      P_k <- holder[x, c(paste(column,i,sep = "_"))]/100 # calculate the P_k value
      if (P_k != 0) { # omit P_k = 0 case to avoid ln(0) exist in the calculation
        entropy1 <- entropy1 + P_k*log(P_k)/log(length(types)) # add each summation term to the entropy value
        holder$entropy[x] <- -entropy1 # record the individual entropy value to the holder column
      }
    }
  }
  output <- holder # create output dataframe based on holder
  if (exclude_intermediate) { # check whether to exclude intermediate results or not
    cols.to.remove <- NULL # create holder to contain the names of column that want to be removed
    for (i in types){ # loop over each category
      cols.to.remove <- c(cols.to.remove, paste(column,i,sep = "_")) # add each category's column name to the name holder
    }
    output <- output[, ! names(output) %in% cols.to.remove, drop = F] # remove all the columns have name included in the name holder
  }
  output <- st_as_sf(output, sf_column_name = "geometry", crs = proj) # project and transform output as sf object
  return(output) # return output
}

#' Calculate total traffic volume within stop buffer areas
#'
#' @param stops a point layer
#' @param input_layer an polyline layer
#' @param proj numeric ESPG number of projection
#' @param buffer_length an single numeric value (note: the unit is based purely on the projection)
#' @param column character column name that contains traffic volume
#'
#' @return a polygon layer based on stops buffer with extra columns contain the results
#' @export
#'
#' @examples
total_volume_in_buffer <- function(stops, input_layer, proj, buffer_length, column){
  stops <- st_transform(stops, crs = proj) # reproject stops layer
  input_layer <- st_transform(input_layer, crs = proj) # reproject input layer
  input_layer <- as.data.frame(input_layer) # coerce sf object into dataframe only
  input_layer[,c("input_layer_ID")] <- c(1:nrow(input_layer)) # insert ID column
  input_layer <- st_as_sf(input_layer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  
  stops_buffer <- st_buffer(stops, dist = buffer_length) # create buffer of stops layer
  stops_buffer <- as.data.frame(stops_buffer) # coerce sf object into dataframe only
  stops_buffer[,c("stops_buffer_ID")] <- c(1:nrow(stops_buffer)) # insert ID column
  stops_buffer <- st_as_sf(stops_buffer, sf_column_name = "geometry", crs = proj) # coerce dataframe back into sf object
  
  inter <- st_intersection(stops_buffer,input_layer) # intersect the buffer and input_layer
  
  output <- as.data.frame(stops_buffer) # create output based on stops layer
  colnumber <- ncol(output) # calculate the number of columns of output
  frame_inter <- as.data.frame(inter) # coerce sf object into dataframe only
  for (i in 1:nrow(stops)) { # loop over each stop row
    total_volume_inside <- 0 # create holder for total volume inside each buffer
    if (dim(inter[which(frame_inter$stops_buffer_ID==i),0])[1]!=0){ # check whether there are intersected polyline in the buffer
      for (x in 1:dim(inter[which(frame_inter$stops_buffer_ID == i),0])[1]) { # loop over each stop row's intersected area
        volume_within <- frame_inter[which(frame_inter$stops_buffer_ID == i),column][x] # retrieve the size of the intersected volume
        if (!is.na(volume_within)) { # get rid of NA in the data source
          total_volume_inside <- total_volume_inside + volume_within # record to total volume inside each buffer
        }
      }
      output[i, colnumber + 1] <- total_volume_inside # write each stop's total value to output
    } else {
      output[i, colnumber + 1] <- 0 # if there is no volume inside, record 0 to it
    }
  }
  colnames(output)[colnumber + 1] <- "total_volume_within_buffer"
  output <- st_as_sf(output, sf_column_name = "geometry", crs = proj) # project and transform output as sf object
  return(output) # return output
}
