library(dplyr)
library(tidyr)
library(rworldmap)
library(RColorBrewer)
library(readxl)

#Plot map of countries represented in main study, and how many times
#Read in data for main study
Main_Study <- read.csv("FosterMA8_DataExtracted_Main - STUDY (1).csv")

#Remove empty rows
Main_Study <- Main_Study[!is.na(Main_Study$year),]

#Unnest countries
allCountries <- Main_Study %>%
  mutate(country = strsplit(country, ", ")) %>%
  unnest(country) 

#count country occurrence frequencies
CountryCounts <- plyr::count(allCountries$country)

#join country count data to a coarse resolution map
MapCounts <- joinCountryData2Map(CountryCounts, joinCode="NAME", nameJoinColumn="x")

#create colourblind-friendly palette
colorPalette <- brewer.pal(3,'YlGnBu')

#plot map
par(mar=c(0,0,1,0))
mapParams <- mapCountryData(MapCounts, 
                            nameColumnToPlot="freq", 
                            catMethod="categorical", 
                            addLegend=FALSE,
                            mapTitle= 'Countries in Dataset',
                            colourPalette = colorPalette)

#add legend
do.call( addMapLegendBoxes,
         c(mapParams,
           x='left',
           title="N Studies",
           cex = 0.55))

#Now for pilot data and country data combined
#Load pilot data
PilotDataStudy <- read_excel("FosterMA8_DataExtracted_Pilot_Split.xlsx", sheet = 2)

#Bind data Study sheets
#Remove VariableCategory column from MainStudy
allCountries <- subset(allCountries, select=-c(VariableCategory))
AllDataMerged <- rbind(allCountries, PilotDataStudy)

#Select unique studyID + country names
AllDataMerged <- unique(AllDataMerged[,c('study_ID','country')])

#Trim whitespace
AllDataMerged$country <- trimws(AllDataMerged$country)

#Correct type
AllDataMerged$country[AllDataMerged$country == "Armernia"] <- "Armenia"

#count country occurrence frequencies
CountryCounts <- plyr::count(AllDataMerged$country)

#join country count data to a coarse resolution map
MapCounts <- joinCountryData2Map(CountryCounts, joinCode="NAME", nameJoinColumn="x")

#create colourblind-friendly palette
colorPalette <- brewer.pal(7,'YlOrRd')

#plot map
par(mar=c(0,0,1,0))
mapParams <- mapCountryData(MapCounts, 
                            nameColumnToPlot="freq", 
                            catMethod="categorical", 
                            addLegend=FALSE,
                            mapTitle= NA,
                            colourPalette = colorPalette)

do.call( addMapLegendBoxes,
         c(mapParams,
           x='left',
           title="N Studies",
           cex = 0.55))
