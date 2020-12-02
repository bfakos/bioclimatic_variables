# R code for calculating bioclimatic variables (BCVs)
# Updated version (02 Dec 2020) of Supporting Information S4 of Bede-Fazekas, Á., & Somodi, I. (2020). The way bioclimatic variables are calculated has impact on potential distribution models. Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.13488


# databases: a named list, each element is a data.frame (rows for points, columns for monthly climate data)
#            all list elements have tha same structure (i.e. number of rows, order of the rows)
#            column name of the list elements are P1, P2, ..., Tmean1, ..., Tmin1, ...Tmax1, ..., Tmax12
#            the first list element is the reference database
# returned value: a named list, each element is a named list of one (for the reference database) or two (for future databases) data.frames called "static" and "dynamic"
#                 all data.framed have the same structure (i.e. number of rows, order of the rows), which is similar to that of the input databases
#                 column name of the data.frames are bcv1, bcv2, ..., bcv19
calculate_BCVs <- function(databases){
	database_names <- names(databases) # name of the databases, first one is the reference
	P_columns <- paste0("P", as.character(1:12)) # name of the precipitation columns
	Tmean_columns <- paste0("Tmean", as.character(1:12)) # name of the mean temperature columns
	Tmin_columns <- paste0("Tmin", as.character(1:12)) # name of the minimum temperature columns
	Tmax_columns <- paste0("Tmax", as.character(1:12)) # name of the maximum temperature columns
	number_of_columns_before_P <- 0 # the number of the columns befor the first precipitation column
	number_of_columns_before_Tmean <- length(P_columns) # the number of the columns before the first mean temperature column
	number_of_columns_before_Tmin <- length(P_columns) + length(Tmean_columns) # the number of the columns before the first inimum temperature column
	number_of_columns_before_Tmax <- length(P_columns) + length(Tmean_columns) + length(Tmin_columns) # the number of the columns before the first maximum temperature column
	row_number <- nrow(databases[[1]]) # the number of rows (i.e. points) in the databases
	specific_period_names <- c("wettestQ", "driestQ", "warmestQ", "coldestQ", "wettestM", "driestM", "warmestM", "coldestM") # name of the specific months and quarters
	months_and_quarters <- list() # create an initially empty list
	for (database_name in database_names) { # iterate through the databases
		months_and_quarters[[database_name]] <- data.frame(matrix(data = NA, ncol = length(specific_period_names), nrow = row_number)) # add an empty dta frame as a new list element
		colnames(months_and_quarters[[database_name]]) <- specific_period_names # set column names
		Tmean_values_of_quarters <- rowMeans(databases[[database_name]][, Tmean_columns[1:3]]) # calculate the JFM mean temperature values
		P_values_of_quarters <- rowSums(databases[[database_name]][, P_columns[1:3]]) # alculate the JFM precipitation values
		for (month in 2:12){ # iterate throuh months
			Tmean_values_of_quarters <- cbind(Tmean_values_of_quarters, rowMeans(databases[[database_name]][, Tmean_columns[c(month, month %% 12 + 1, (month + 1) %% 12 + 1)]])) # append the mean temperature values of the next quarter as a new column
			P_values_of_quarters <- cbind(P_values_of_quarters, rowSums(databases[[database_name]][, P_columns[c(month, month %% 12 + 1, (month + 1) %% 12 + 1)]])) # append the precipitation values of the next quarter as a new column
		} # for month
		months_and_quarters[[database_name]][, "wettestM"] <- apply(X = databases[[database_name]][, P_columns], MARGIN = 1, FUN = which.max) # calculate the indices of the wettest month
		months_and_quarters[[database_name]][, "driestM"] <- apply(X = databases[[database_name]][, P_columns], MARGIN = 1, FUN = which.min) # calculate the indices of the driest month
		months_and_quarters[[database_name]][, "warmestM"] <- apply(X = databases[[database_name]][, Tmean_columns], MARGIN = 1, FUN = which.max) # calculate the indices of the warmest month
		months_and_quarters[[database_name]][, "coldestM"] <- apply(X = databases[[database_name]][, Tmean_columns], MARGIN = 1, FUN = which.min) # calculate the indices of the coldest month
		months_and_quarters[[database_name]][, "wettestQ"] <- apply(X = P_values_of_quarters, MARGIN = 1, FUN = which.max) # calculate the indices of the wettest quarter
		months_and_quarters[[database_name]][, "driestQ"] <- apply(X = P_values_of_quarters, MARGIN = 1, FUN = which.min) # calculate the indices of the driest quarter
		months_and_quarters[[database_name]][, "warmestQ"] <- apply(X = Tmean_values_of_quarters, MARGIN = 1, FUN = which.max) # calculate the indices of the warmest quarter
		months_and_quarters[[database_name]][, "coldestQ"] <- apply(X = Tmean_values_of_quarters, MARGIN = 1, FUN = which.min) # calculate the indices of the coldest quarter
	} # for database_name
	output <- list() # create an empty list
	for (database_name in database_names) { # iterate through the databases
		if (database_name == database_names[1]) { # if reference database is processed
			temporal_approaches <- "static" # there is no sense to use the dynamic approach
		} else { # if future database is processed
			temporal_approaches <- c("static", "dynamic") # the two temporal approaches
		} # else
		output[[database_name]] <- list() # create an empty list as a new list element
		for (temporal_approach in temporal_approaches) { # iterate through the one or two temporal approaches
			if (temporal_approach == "static") { # is static approach is used
				database_name_for_selection <- database_names[1] # use the reference database for selecting the specific month/quarter
			} else { # if dynamic approach is used
				database_name_for_selection <- database_name # use the database under processing for selecting the specific month/quarter
			} # else
			database <- databases[[database_name]] # create a copy of the original database of monthly climate values
			database$bcv1 <- rowMeans(database[, Tmean_columns]) # BCV1 = Annual Mean Temperature
			database$bcv2 <- (rowSums(database[, Tmax_columns]) - rowSums(x = database[, Tmin_columns])) / 12 # BCV2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
			database$bcv4 <- sqrt(rowMeans((database[, Tmean_columns] - rowMeans(database[, Tmean_columns])) ^ 2)) * 100 # BCV4 = Temperature Seasonality (standard deviation * 100)
			database$bcv5 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmax + months_and_quarters[[database_name_for_selection]][, "warmestM"])]) # BCV5 = Max Temperature of Warmest Month
			database$bcv6 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmin + months_and_quarters[[database_name_for_selection]][, "coldestM"])]) # BCV6 = Min Temperature of Coldest Month
			database$bcv7 <- database$bcv5 - database$bcv6 # BCV7 = Temperature Annual Range (BCV5-BCV6)
			database$bcv3 <- (database$bcv2 / database$bcv7) * 100 # BCV3 = Isothermality (BCV2/BCV7) (* 100)
			indices_of_first_month <- months_and_quarters[[database_name_for_selection]][, "wettestQ"] # get the indices of the first month of the quarter
			values1 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + indices_of_first_month)]) # values of the first month of the quarter
			values2 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + indices_of_first_month %% 12 + 1)]) # values of the second month of the quarter
			values3 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + (indices_of_first_month + 1) %% 12 + 1)]) # values of the third month of the quarter
			database$bcv8 <- rowMeans(cbind(values1, values2, values3)) # BCV8 = Mean Temperature of Wettest Quarter
			indices_of_first_month <- months_and_quarters[[database_name_for_selection]][, "driestQ"] # get the indices of the first month of the quarter
			values1 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + indices_of_first_month)]) # values of the first month of the quarter
			values2 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + indices_of_first_month %% 12 + 1)]) # values of the second month of the quarter
			values3 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + (indices_of_first_month + 1) %% 12 + 1)]) # values of the third month of the quarter
			database$bcv9 <- rowMeans(cbind(values1, values2, values3)) # BCV9 = Mean Temperature of Driest Quarter
			indices_of_first_month <- months_and_quarters[[database_name_for_selection]][, "warmestQ"] # get the indices of the first month of the quarter
			values1 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + indices_of_first_month)]) # values of the first month of the quarter
			values2 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + indices_of_first_month %% 12 + 1)]) # values of the second month of the quarter
			values3 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + (indices_of_first_month + 1) %% 12 + 1)]) # values of the third month of the quarter
			database$bcv10 <- rowMeans(cbind(values1, values2, values3)) # BCV10 = Mean Temperature of Warmest Quarter
			indices_of_first_month <- months_and_quarters[[database_name_for_selection]][, "coldestQ"] # get the indices of the first month of the quarter
			values1 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + indices_of_first_month)]) # values of the first month of the quarter
			values2 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + indices_of_first_month %% 12 + 1)]) # values of the second month of the quarter
			values3 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_Tmean + (indices_of_first_month + 1) %% 12 + 1)]) # values of the third month of the quarter
			database$bcv11 <- rowMeans(cbind(values1, values2, values3)) # BCV11 = Mean Temperature of Coldest Quarter
			database$bcv12 <- rowSums(database[, P_columns]) # BCV12 = Annual Precipitation
			database$bcv13 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + months_and_quarters[[database_name_for_selection]][, "wettestM"])]) # BCV13 = Precipitation of Wettest Month
			database$bcv14 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + months_and_quarters[[database_name_for_selection]][, "driestM"])]) # BCV14 = Precipitation of Driest Month
			database$bcv15 <- sqrt(rowMeans((database[, P_columns] - rowMeans(database[, P_columns])) ^ 2)) / rowMeans(database[, P_columns]) # BCV15 = Precipitation Seasonality (Coefficient of Variation)
			indices_of_first_month <- months_and_quarters[[database_name_for_selection]][, "wettestQ"] # get the indices of the first month of the quarter
			values1 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + indices_of_first_month)]) # values of the first month of the quarter
			values2 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + indices_of_first_month %% 12 + 1)]) # values of the second month of the quarter
			values3 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + (indices_of_first_month + 1) %% 12 + 1)]) # values of the third month of the quarter
			database$bcv16 <- rowSums(cbind(values1, values2, values3)) # BCV16 = Precipitation of Wettest Quarter
			indices_of_first_month <- months_and_quarters[[database_name_for_selection]][, "driestQ"] # get the indices of the first month of the quarter
			values1 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + indices_of_first_month)]) # values of the first month of the quarter
			values2 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + indices_of_first_month %% 12 + 1)]) # values of the second month of the quarter
			values3 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + (indices_of_first_month + 1) %% 12 + 1)]) # values of the third month of the quarter
			database$bcv17 <- rowSums(cbind(values1, values2, values3)) # BCV17 = Precipitation of Driest Quarter
			indices_of_first_month <- months_and_quarters[[database_name_for_selection]][, "warmestQ"] # get the indices of the first month of the quarter
			values1 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + indices_of_first_month)]) # values of the first month of the quarter
			values2 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + indices_of_first_month %% 12 + 1)]) # values of the second month of the quarter
			values3 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + (indices_of_first_month + 1) %% 12 + 1)]) # values of the third month of the quarter
			database$bcv18 <- rowSums(cbind(values1, values2, values3)) # BCV18 = Precipitation of Warmest Quarter
			indices_of_first_month <- months_and_quarters[[database_name_for_selection]][, "coldestQ"] # get the indices of the first month of the quarter
			values1 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + indices_of_first_month)]) # values of the first month of the quarter
			values2 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + indices_of_first_month %% 12 + 1)]) # values of the second month of the quarter
			values3 <- as.numeric(database[cbind(c(1:row_number), number_of_columns_before_P + (indices_of_first_month + 1) %% 12 + 1)]) # values of the third month of the quarter
			database$bcv19 <- rowSums(cbind(values1, values2, values3)) # BCV19 = Precipitation of Coldest Quarter
			output[[database_name]][[temporal_approach]] <- database[, paste0("bcv", as.character(1:19))] # append a new list element to the output list
		} # for temporal_approach
	} # for database_names
	return(output) # return the output list
} # calculate_BCVs()
