################################################################################
################################################################################
#### Hydrological regime model for rock pool dynamics                       ####
#### Author: Neil Collier & Rylan Shearn                                    ####
#### Revision History:                                                      ####
#### Date started: September, 2014                                          ####
#### Date completed: January, 2015                                          ####
#### Cut down to minimal working version (RS): 2016.09.14                   ####
####                                                                        ####
#### License information:                                                   ####
#### https://github.com/rylanshearn/makeItRain/blob/master/LICENSE          ####
################################################################################
################################################################################
########## - This code works with BOM datasets extracted from database #########

rm(list=ls()) # Clear workspace
setwd("~/bin/makeItRain")
set.seed(7210)

### Import datasets
frogRock <- read.csv("frogRock.csv", header=TRUE, stringsAsFactors=FALSE)
boyaginRock <- read.csv("boyaginRock.csv", header=TRUE, stringsAsFactors=FALSE)
gorgeRock <- read.csv("gorgeRock.csv", header=TRUE, stringsAsFactors=FALSE)
mtChudalup <- read.csv("mtChudalup.csv", header=TRUE, stringsAsFactors=FALSE)
puntapinRock <- read.csv("puntapinRock.csv", header=TRUE, stringsAsFactors=FALSE)
wanarraRock <- read.csv("wanarraRock.csv", header=TRUE, stringsAsFactors=FALSE)

### Import monthly pan evaporation data.
panMonth <- read.csv("panMonth.csv")

### Rename columns and melt data.
library(reshape2)
colnames(panMonth)[2:13] <- c("Jan","Feb","March","April","May","June",
    "July","August","Sep","Oct","Nov","Dec")
panMonthL <- melt(panMonth, id.vars="Pool", variable.name="Month",
    value.name="Pan")

### Remove duplicate rows in the dataframe.
panUnique <- unique(panMonthL)
panUnique # Subset this dataframe to select pan values for simulation.

### Plot the pan data
library(ggplot2) # Load ggplot2
panPlot <- ggplot(panUnique, aes(x=Month, y=Pan)) +
   geom_bar(stat="identity", fill="grey75",
       colour="white") +
   theme_bw() +
   ylab("Mean pan evaporation (mm)") +
   facet_wrap(~Pool)+
   theme(legend.position="") +
   scale_x_discrete(labels=c("J","F","M","A","M","J","J","A","S","O","N","D"))

### Save plot to file
ggsave("pan_month.png", panPlot, 
    path="~/bin/makeItRain")
################################################################################
########### - Part 1: Function to format data and fit tweedie model - ##########
################################################################################
formFit <- function(data, site){
    ## Format dataframe
    colnames(data)[6] <- "Rainfall" # Rename rainfall variable
    data$Month <- factor(data$Month,
        labels=c("January", "February", "March", "April", "May", "June", "July", 
            "August", "September", "October", "November", "December"), 
            levels=1:12)
    data$Site <- site
    simData <- data[ ,c("Site", "Year", "Month", "Day", "Rainfall")]
    
    ## Tweedie fitting function
	rainFit <- function(data){
    require(plyr)
    require(tweedie)
    require(statmod)
    fit <- tweedie.profile(Rainfall~1, xi.vec=seq(1.1, 1.9, length=9),
                           do.plot=FALSE, data=data)
    ## Variables to keep in a dataframe
    meanRain = mean(data$Rainfall, na.rm = TRUE) # Mean rainfall per event
    xi <- fit$xi.max # xi is synonomous with power
    phi <- fit$phi.max # Dispersion parameter
    return(data.frame(Site=data[ ,1][1], MeanRainfall=meanRain, xi=xi, phi=phi))
    }
    
    ## Apply 'rainPar' function to data
    rainPar <- function(data){
    	    require(plyr)
    	    ddply(data, .(Site, Month), rainFit)
    	    }
d_fit <- rainPar(simData)
d_fit
}
################################################################################
################# - Part 2: Simulate the rainfall dataset - ####################
################################################################################
#### data = Dataframe returned from formFit().
#### sims = Number of simulations (years).

mir <- function(data, sims){
	
	# Warnings and stops
	if(is.data.frame(data) == FALSE)
	stop("Rainfall simulation function requires a dataframe of model parameters.
	Please use the dataframe returned by the formFit() function (Part 1).")
	
    # Stop is number of simmulations argument is missing
    if (missing(sims))
    stop("Number of simulations not specified")
   
	rainSim <- function(data){ # add data and specify number of years for simulation
    require(reshape2)
    monthLength <- function(data){   # Set the length of each month's rainfall vector
        if (data$Month == "February")
            return(28)
        if (data$Month == "September" | data$Month == "April" | data$Month == "June" |
                data$Month == "November")
                return(30)
        else
            return(31)
    }

    Month <- rep(data$Month, monthLength(data))
    Day <- 1:monthLength(data) # Day sequence
    d <- data.frame(Month, Day)
    nsims <- sims
    dMatrix <- matrix(0, nrow=monthLength(data), ncol=nsims) # Empty matrix to store rainfall.
    # for loop to simulate multiple years
    for (i in 1:nsims){
        dMatrix[ ,i] <- rtweedie(n = monthLength(data), mu = data$MeanRainfall,
                             power = data$xi, phi = data$phi)
    }
    
    d <- cbind(d,dMatrix)
    dmelt <- melt(d, id.vars = c("Month", "Day"), variable.name = "Year",
                  value.name = "Rainfall") # Melt the data to long form
    dmelt # Return the melted dataframe
    }
	
	## Apply the 'rainSim' function using the estimated parameters from 'rainFit'
	rainSimData <-function(data){
    require(plyr)
    ddply(data, .(Site, Month), rainSim)
    }
    rainSimResults <- rainSimData(data) # Return the simulated rainfall dataframe
}
################################################################################
################ - Part 3: Simulate the pool volume dynamics - #################
################################################################################
#### Part 3a: Add a vector of daily pan evaporation estimates for each month
evap <- function(data, panVec){ ## days argument is the number of days in the month

	# Stop simulation if 'panVec' length != (12)
	if (length(panVec) != 12)
    stop("Pan input data doesn't have 12 months of data. Check the input vector:
    it must have twelve values.")
    
    for (i in 1:nrow(data)){

	if (data$Month[i] == "January")
		    data$panEvap[i] <- panVec[1]/31
		if (data$Month[i] == "February")
		    data$panEvap[i] <- panVec[2]/28
		if (data$Month[i] == "March")
		    data$panEvap[i] <- panVec[3]/31
		if (data$Month[i] == "April")
		    data$panEvap[i] <- panVec[4]/30
		if (data$Month[i] == "May")
		    data$panEvap[i] <- panVec[5]/31
		if (data$Month[i] == "June")
		    data$panEvap[i] <- panVec[6]/30
		if (data$Month[i] == "July")
		    data$panEvap[i] <- panVec[7]/31
		if (data$Month[i] == "August")
		    data$panEvap[i] <- panVec[8]/31
		if (data$Month[i] == "September")
		    data$panEvap[i] <- panVec[9]/30
		if (data$Month[i] == "October")
		    data$panEvap[i] <- panVec[10]/31
		if (data$Month[i] == "November")
		    data$panEvap[i] <- panVec[11]/30
		else
		    panVec[12]/31
		}
		data
}
################################################################################
########## - Part 3b: Function to simulate the pool volume dynamics ############
################################################################################
funcVol <- function(data, maxDepth){
	
	# Stop if maxDepth argument is missing
    if (missing(maxDepth))
    stop("Pool depth not specified")
	
	poolVol <- function(data){
        data$poolVol <- 0
        data$DayOfYear <- seq(1:365) # Add a vector of day of the year
        ## Populate the 1st time point
        data$poolVol[1] <- data$Rainfall[1] - data$panEvap[1]
        if (data$poolVol[1] < 0)
            data$poolVol[1] <- 0
        if (data$poolVol[1] > maxDepth)
            data$poolVol[1] <- maxDepth
        else
            data$poolVol[1]
        for (i in 2:length(data$poolVol)){
        ## Remaining vector
        data$poolVol[i] <- data$poolVol[i-1] + data$Rainfall[i] - data$panEvap[i]
        if (data$poolVol[i] < 0)
            data$poolVol[i] <- 0
        if (data$poolVol[i] >= maxDepth)
            data$poolVol[i] <- maxDepth
        else
            data$poolVol[i]
    }
    data
}

	poolSim <- function(data){
    d <- ddply(data, .(Year), poolVol)
    return(d)
    }
    poolSim(data)
}
################################################################################
####### - Part 4: Extracting eco/evo metrics relevant to the study #############
################################################################################
# Conditional summing of sequences in a vector. Locate sequences of positive
# numbers in the vector, perform calculations on these sequences and return
# them in a data object.
# Metrics needed:
# 1. Number of dessication events per pool per year (Dessication frequency).
# 2. Sum the events when no water is present in the rock pool.
# Found some code on stackoverflow. This code is a re-engineered version. The
# original code can be found here:
# http://stackoverflow.com/questions/8537783/conditional-cumulative-sum

# Calculate the wet events and mean rain per event (per day)
wetEvents <- function(data){
	event <- data$poolVol
	y <- ifelse(data$poolVol > 0, 1, 0) # Vector conditional on 'x'
	r <- rle(sign(y))
	s <- diff(c(0, cumsum(event)[cumsum(r$lengths)]))[r$values==1]
	days <- r$lengths[r$values==1] # consecutive days with rain
	meanEvent <- s/r$lengths[r$values==1] # Mean rain per day of event.
	return(list(PoolVol = event, DaysWet = days, MeanRain = meanEvent))
}

# Calculate the dry events
dryEvents <- function(data){
	event <- data$poolVol
	y <- ifelse(data$poolVol == 0, 1, 0) # Vector conditional on 'x'
	r <- rle(sign(y))
	s <- diff(c(0, cumsum(event)[cumsum(r$lengths)]))[r$values==1]
	days <- r$lengths[r$values==1] # Consecutive days without rain
	return(data.frame(DaysDry = days))
}
################################################################################
########## - Function sequence for running analysis on a rock pool - ###########
################################################################################
## Run analysis for Mt Chudalup
a <- formFit(mtChudalup, "Mt Chudalup") # Format data and fit models
b <- mir(data=a, sims=50) # Simulate rainfall

# Prepare pan evaporation data
chudPan <- panUnique[panUnique$Pool=="Mt Chudalup", "Pan"]
c <- evap(data=b, panVec=chudPan) # Calculate evaporation per day
d <- funcVol(data=c, maxDepth=100) # Run the pool dynamics simulation

theme_custom <- theme(
                      panel.background=element_rect(fill="#EAEFED"),
                      panel.grid.major.y = element_blank(), 
                      panel.grid.minor.y = element_blank(),
                      panel.grid.major.x = element_blank(), 
                      panel.grid.minor.x = element_blank(),
                      legend.position="none") 
# Plot the data
f_1 <- ggplot(d, aes(x=Day, y=poolVol))+
    facet_wrap(~Month, nrow=4, ncol=3)+
    theme_custom+
    ylab("Pool depth (mm)")+
    xlab("Day")+
    geom_point(alpha=0.25, size=1.25)+
    geom_line(alpha=0.25, aes(group=Year))
f_1

# Summary statistics for dry and wet events
e <- wetEvents(d) # Calculate wet day metrics
f <- dryEvents(d) # Calculate dry day metrics

# Plot some metrics
hist(f$DaysDry, breaks=365, col="grey75")
hist(e$DaysWet, breaks=365, col="grey75")
hist(e$MeanRain, col="grey75")

################################################################################
#################### - Model fit diagnostics - #################################
# Example using Mt Chudalup data

dataMod <- function(data, site){
    ## Format dataframe
    colnames(data)[6] <- "Rainfall" # Rename rainfall variable
    data$Month <- factor(data$Month,
        labels=c("January", "February", "March", "April", "May", "June", "July", 
            "August", "September", "October", "November", "December"), 
            levels=1:12)
    data$Site <- site
    simData <- data[ ,c("Site", "Year", "Month", "Day", "Rainfall")]
}

d_obs <- dataMod(mtChudalup, "Mt Chudalup") # Observed rainfall data
head(d_obs)
pars <- formFit(mtChudalup, "Mt Chudalup") # Format data and fit models

d_sim <- mir(data=pars, sims=50)
head(d_sim)

# Put this in a function perhaps using ggplot2 instead of base graphics
diagPlot <- function(obs, sim){
	require(ggplot2)
	require(gridExtra)
	# Extract variables from dataset and rbind to new dataset
	d_1 <- obs[ ,c("Site","Year","Month","Day","Rainfall")]
	d_1$Data <- "Observed"
	d_2 <- sim[ ,c("Site","Year","Month","Day","Rainfall")]
	d_2$Data <- "Simulated"
	d_3 <- rbind(d_1, d_2)
	
	
	rainfallSum <- ddply(d_3, .(Data, Month), summarise,
                     mean_Rainfall = mean(Rainfall, na.rm=TRUE),
                     sd_Rainfall = sd(Rainfall, na.rm=TRUE) )

# Plot data
	p_1 <- ggplot(rainfallSum, aes(Month, mean_Rainfall, ymin=mean_Rainfall-sd_Rainfall,
                        ymax=mean_Rainfall+sd_Rainfall, colour=Data)) +
    geom_pointrange(position=position_dodge(width = 0.75), size=0.5) +
    theme_bw() +
    theme(legend.position="none") +
    scale_colour_manual(values = c("grey25", "red")) +
    ylab("Meand rainfall/day (mm)") +
    xlab("Month") +
    theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10)) +
    scale_x_discrete(labels=c("J","F","M","A","M","J","J","A","S","O","N","D"))

    p_2 <- ggplot(subset(d_3, Rainfall!=0), aes(x=log(Rainfall+1), group=Data)) +
        facet_wrap(~Month, nrow=4) +
        theme_bw() +
        scale_fill_manual(values = c("grey25", "red")) +
        scale_colour_manual(values = c("grey25", "red")) +
        theme(legend.position="none") +
        ylab("Frequency") +
        xlab("Rainfall (mm/day)") +
        theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10)) +
        #geom_histogram(aes(y= ..density.. , fill=Data), alpha=0.5) +
        geom_density(aes(colour=Data), alpha=0.5)
    
# Summarise data frames
    d_3_sum <- ddply(d_1, .(Year, Month), summarise,
                     mean_Rainfall = mean(Rainfall, na.rm=TRUE),
                     sd_Rainfall = sd(Rainfall, na.rm=TRUE) )
    d_4_sum <- ddply(d_2, .(Year, Month), summarise,
                     mean_Rainfall = mean(Rainfall, na.rm=TRUE),
                     sd_Rainfall = sd(Rainfall, na.rm=TRUE) )

    p_3 <- ggplot(d_3_sum, aes(x=Month, y=log(mean_Rainfall+1))) +
        geom_point(aes(group=Year), alpha=0.25, colour="grey25") +
        geom_line(aes(group=Year), alpha=0.25, colour="grey25") +
        ylab("Mean rainfall (mm/day)") +
        theme_bw()+
        theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10)) +
        scale_x_discrete(labels=c("J","F","M","A","M","J","J","A","S","O","N","D"))
    
    p_4 <- ggplot(d_4_sum, aes(x=Month, y=log(mean_Rainfall+1))) +
        geom_point(aes(group=Year), alpha=0.25, colour="red") +
        geom_line(aes(group=Year), alpha=0.25, colour="red") +
        ylab("Mean rainfall (mm/day)") +
        theme_bw()+
        theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10)) +
        scale_x_discrete(labels=c("J","F","M","A","M","J","J","A","S","O","N","D"))
    grob_1 <- arrangeGrob(p_1, p_3, p_4, nrow=3)
    g <- arrangeGrob(grob_1, p_2, ncol=2)
    ggsave("Diagnostic_plot.png", g, width=8, height=7, dpi=600)
}

diagPlot(d_obs, d_sim) # Plot the observed and simulated data.
