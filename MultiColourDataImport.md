---
title: "MultiDataImport"
author:
- Laurel Genge
- Carlie Barnhill
- Max Berthold
- Douglas A. Campbell
date: "2021-05-27"
output:
  html_document:
    code_folding: hide
    keep_md: yes
    fig_caption: yes
    toc: TRUE
    toc_float: TRUE
bibliography: MicrobialGrowthFitting.bib
csl: plos-one.csl
---
# Run These Chunks First


# Introduction
The PSI Multicultivator is used to grow 8 x 80 ml of phytoplankton culture under a common temperature regime, with individual control of bubbling, light level, light spectral quality and photoperiod for each of the 8 culture tubes.

This .Rmd Rworkbook imports data in simple .csv long form exported from PSI Multicultivators based upon project specific values for variables set by the user.

It tidies and organizes the data.
It uses a pivot_wider and interpolation approach to get the actinic_par and OD values in line rowwise.
This requires careful 'arrange' of the rows.
It imports a metadata Catalog and merges the metadata with the imported data based upon shared values for the variables 'MC', 'Tube', and 'Filename' which should unambiguously identify a given growth trajectory measured at OD680 or OD720.

It generates preliminary data plots.
It filters the data for outliers by screening out values distant from the moving average of a window in the stream; this requires careful 'arrange' of the rows so sequential rows represent sequential time steps.

This works because the OD680 & OD720 data are only episodically, but widely, aberrant when a bubble interrupts the measurement, and if the MultiCultivator is running properly these bubble aberration events are rare.


# Load Libraries and set Project Variables

```r
# libraries; Note check actual dependencies
library(tidyverse)
library(lubridate)
library(stringr)
library(broom)
library(knitr)
library(zoo)
#library(tidyquant)
library(data.table)
library(googledrive)
library(googlesheets4)

#attach psych for harmonic.mean() for averaging rates
#library(psych)
```


```r
#"..", takes up a level in the directory path
Project <- "PICO"
RawData <- file.path("MultiColourData")
PlotsPath <- file.path("Figs")
ImportedData <- file.path("ImportedColourData")

Sensor <- "od-680"

#set the hour of day at which photoperiod starts for later generation of ToD Time of Day
StartHour <- 6
```


```r
MyWavelengths = c(405, 450, 475, 530, 615, 660, 730, "WW")
MCMIXColours = c("violet", "darkblue", "dodgerblue", "green","orange","red","purple", "black")

names(MCMIXColours) <- MyWavelengths
MCMIXColours
```

```
##          405          450          475          530          615          660 
##     "violet"   "darkblue" "dodgerblue"      "green"     "orange"        "red" 
##          730           WW 
##     "purple"      "black"
```

```r
# SensorWavebands = c(680, 720)
# SensorColours = c("red", "black")
# names(SensorColours) <- SensorWavebands
# 
# SensorColours
```

# Import MetaData

```r
# #implement read with googlesheet name instead of full url

gs4_deauth()

# imagine this is the URL or ID of a Sheet readable by anyone (with a link)

Catalog <- read_sheet("https://docs.google.com/spreadsheets/d/1ZXpwR7Gfto-uRzVdXzMpQF4frbrvMLH_IyLqonFZRSw/edit#gid=0") %>%
  # read_sheet has an annoying "feature" to set the type of columns it can't parse to a list.
  # ggplot/dplyr doesn't like working with a dataframe of lists.
  # In this case WL is set to a list since some values are numbers, some are strings, some are blank.
  # To fix this, first drop all rows missing WL, then unlist.
  # Must first drop NA rows since unlist will collapse NULL lists, then the unlisted WL is a shorter length than original WL column, which mutate doesn't like.
  drop_na(WL) %>%
  mutate(WL = unlist(WL))

as.data.frame(Catalog)

Catalog <- Catalog %>%
  mutate(ExpDate = ymd(ExpDate))
```

# List previously imported files

```r
list.files(path = ImportedData, pattern = Project, full.names = TRUE)
```

```
## [1] "ImportedColourData/20200124_PICO_MCMIX004_RUN1_TargetDataMetaFilter.Rds"
```

# List MultiCulti files
for Project that are saved in RawData

```r
MultiFiles <- list.files(path = RawData, pattern = Project, full.names = TRUE)

#check file names
MultiFiles
```

```
## [1] "MultiColourData/20200124_PICO_MCMIX004_RUN1.csv"
```

# Set Target File
Could implement function to list difference between unique(GrowthSummaryAppend$Filename) and MultiFiles if needed later; scan visually to set TargetFile for now.


```r
TargetFile <- "MultiColourData/20200124_PICO_MCMIX004_RUN1.csv"

#setting the second index number depends upon the number of "/" in the TargetFile path
TargetFileName <- str_split(string = TargetFile, "/")[[1]][2] %>%
  str_remove(pattern = ".csv")
```

## Create function using data.table::fread to skip the beginning comments and starts reading file after key word Skip

```r
fread_plus <- function(Flnm, Skip){data.table::fread(file = Flnm, skip = Skip) %>%
    mutate(Filename = Flnm, CDateTime = ymd_hms(file.info(Flnm)$ctime))
}
```

# Read in TargetFile
select(-c(V5)) only needed for some files.

```r
TargetData <- fread_plus(Flnm = TargetFile, Skip = "key") %>%
select(-c(V5))
```

```
## Warning in data.table::fread(file = Flnm, skip = Skip): Detected 5 column names
## but the data has 4 columns. Filling rows automatically. Set fill=TRUE explicitly
## to avoid this warning.
```

```r
TargetData[1:10,]
```

```
##                             key     time               abs-time  value
##  1: actinic-lights.light-1-B450 0.000000 24-Jan-2020 6:00:01 AM     NA
##  2: actinic-lights.light-1-B450 0.000278 24-Jan-2020 6:00:01 AM 0.0000
##  3: actinic-lights.light-1-B450 0.270556 24-Jan-2020 6:01:38 AM 0.5009
##  4: actinic-lights.light-1-B450 0.469167 24-Jan-2020 6:02:50 AM 1.5011
##  5: actinic-lights.light-1-B450 0.606667 24-Jan-2020 6:03:39 AM 2.5014
##  6: actinic-lights.light-1-B450 0.718889 24-Jan-2020 6:04:20 AM 3.5005
##  7: actinic-lights.light-1-B450 0.816667 24-Jan-2020 6:04:55 AM 4.5019
##  8: actinic-lights.light-1-B450 0.904444 24-Jan-2020 6:05:26 AM 5.5026
##  9: actinic-lights.light-1-B450 0.984722 24-Jan-2020 6:05:55 AM 6.5001
## 10: actinic-lights.light-1-B450 1.059722 24-Jan-2020 6:06:22 AM 7.5015
##                                            Filename           CDateTime
##  1: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08
##  2: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08
##  3: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08
##  4: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08
##  5: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08
##  6: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08
##  7: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08
##  8: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08
##  9: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08
## 10: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08
```



```r
#MultiData files are big; attempting to read multiple files in using Map overruns RAM on most personal computers.
# MultiData <- MultiFiles %>%
#   map_df(~fread_plus(Flnm = ., Skip = "key")) %>%
#   select(-c(V5))
# 
# MultiData[1:10,]

# MultiDataTarget <- MultiData %>%
   #filter(Filename == TargetFile)
```


```r
# TargetData <- TargetData %>%
#   mutate(value = if_else(grepl("actinic-lights", key), if_else(is.na(value), 0, value), value))
```

# Add ToD column, extract ExpDate, extract MC, filter superfluous rows

```r
TargetData <- TargetData %>%
  select(key, time, `abs-time`, value, Filename, CDateTime) %>%
  mutate(Tube = str_extract(key, "-[:digit:]"),
         Tube = as.numeric(str_extract(Tube, "[:digit:]")),
         abs_time = dmy_hms(`abs-time`),
         ToD = (time + StartHour) %% 24) %>%
  select(-`abs-time`)
  
 
#extract ExpDate for matching with Catalog
TargetData <- TargetData %>% 
    mutate(ExpDate = str_extract(Filename, "/202[:digit:][:digit:][:digit:][:digit:][:digit:]_"),
           ExpDate = ymd(str_extract(ExpDate, "202[:digit:][:digit:][:digit:][:digit:][:digit:]")))

#extract MC for matching with Catalog
#will break with white MC; think about it.
TargetData <- TargetData %>% 
    mutate(MC = str_extract(Filename, "MCMIX[:digit:][:digit:][:digit:]"))

#filter superfluous rows to simplify later pivot
TargetData <- TargetData %>%
  # filter(key != "thermo.temperature",
  #        key !="thermo.thermo-reg",
  #        key != "mc-airpump.airpump") %>%
  filter(str_detect(key, "od-720|od-680|actinic-lights.light")) 

#add ln(value) later

  
TargetData[1:10,]
```

```
##                             key     time  value
##  1: actinic-lights.light-1-B450 0.000000     NA
##  2: actinic-lights.light-1-B450 0.000278 0.0000
##  3: actinic-lights.light-1-B450 0.270556 0.5009
##  4: actinic-lights.light-1-B450 0.469167 1.5011
##  5: actinic-lights.light-1-B450 0.606667 2.5014
##  6: actinic-lights.light-1-B450 0.718889 3.5005
##  7: actinic-lights.light-1-B450 0.816667 4.5019
##  8: actinic-lights.light-1-B450 0.904444 5.5026
##  9: actinic-lights.light-1-B450 0.984722 6.5001
## 10: actinic-lights.light-1-B450 1.059722 7.5015
##                                            Filename           CDateTime Tube
##  1: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08    1
##  2: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08    1
##  3: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08    1
##  4: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08    1
##  5: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08    1
##  6: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08    1
##  7: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08    1
##  8: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08    1
##  9: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08    1
## 10: MultiColourData/20200124_PICO_MCMIX004_RUN1.csv 2021-05-25 15:53:08    1
##                abs_time      ToD    ExpDate       MC
##  1: 2020-01-24 06:00:01 6.000000 2020-01-24 MCMIX004
##  2: 2020-01-24 06:00:01 6.000278 2020-01-24 MCMIX004
##  3: 2020-01-24 06:01:38 6.270556 2020-01-24 MCMIX004
##  4: 2020-01-24 06:02:50 6.469167 2020-01-24 MCMIX004
##  5: 2020-01-24 06:03:39 6.606667 2020-01-24 MCMIX004
##  6: 2020-01-24 06:04:20 6.718889 2020-01-24 MCMIX004
##  7: 2020-01-24 06:04:55 6.816667 2020-01-24 MCMIX004
##  8: 2020-01-24 06:05:26 6.904444 2020-01-24 MCMIX004
##  9: 2020-01-24 06:05:55 6.984722 2020-01-24 MCMIX004
## 10: 2020-01-24 06:06:22 7.059722 2020-01-24 MCMIX004
```




## Create preliminary plot for TargetData facetted by Tube

```r
#Run <- c("MultiCulti/20200124_PICO_MCMIX004_RUN1.csv", "MultiCulti/20200124_PICO_MCMIX006_RUN2.csv")

TargetDataPlotPrelim <- TargetData %>%
  filter(grepl("od-", key)) %>%
  #filter(grepl("od-sensors-", key)) %>%
  ggplot(data = .) +
  geom_point(aes(x = time, y = value, colour = as.factor(str_detect(key, "680"))), size = 0.5) +
  scale_x_continuous(breaks=seq(0, 250, by = 125)) + 
  coord_cartesian(xlim = c(-5, 255)) +
  scale_colour_manual(values = c("black", "red")) +
  labs(y = "Optical Density (OD)", x = "Elapsed Time (h)", title = "Tubes") +
  facet_grid(cols = vars(as.factor(Tube))) +
  theme_bw()
#+ 
#    theme(legend.position="none")

TargetDataPlotPrelim
```

![](Figs/prelimplot-1.png)<!-- -->

```r
#cols = vars(Filename), 

#, fill = if_else(str_detect(key, "680"), "red", "black")

#
```

## Save Preliminary Plot to PlotsPath folder if desired
Generates .png for later presentation if needed.

```r
# ggsave(file = file.path(PlotsPath, paste(TargetFileName, "TargetDataPlotPrelim",".png",sep = "")), plot = TargetDataPlotPrelim, device = NULL, scale = 1, height=10, width= 20, units = c("cm"),dpi = 300, limitsize = TRUE)
```

# Generate par_ue column with rows aligned with OD measures
Pivot_wider to get actinic-lights data aligned with relevant sensor data.
Why are there NULL values in the OD sensor columns? Should always have a value?
Import issue? Misalignment b/t light and OD?
Need to include arrange(Filename, time, Tube) to keep things aligned!
Need to group_by and/or reorder rows appropirately; Be Careful

```r
#possible issue with data structure; there are multiple values for some of the rows of actinic light columns, so the column becomes a list.
#Can add  values_fn = 
#to arbitrarily takes the max or min etc. element of the list; but there might be a wider problem here when importing multiple files

TargetDataWide <- TargetData %>%
  pivot_wider(names_from = key, values_from = value, values_fn = list(value = max)) %>%
  arrange(Filename, MC, Tube, time)
```

Actinic light values do not align time wise with OD measures.
Interpolate NA in actinic light columns from last observation, arrange by MC & Tube
Then generate actinic_par summary column
If multiple lights are activated this chunk will give the summed par of all different colours for the tube.
If a single actinic light is activated per tube, this gives the par for that tube.
Filter rows where !is.na(actinic_par) to check for incorrect row sums.

Interpolation for Sine is not necessarily appropriate interpolation for Square photoregime; issues with propagating last actinic_par of afternoon through evening, or back-casting first actinic_par of morning.

Small glitching adding actinic_light values for tubes where actinic_light column should be 0; issue with interpolation we think.

```r
#http://publish.illinois.edu/spencer-guerrero/2014/12/11/2-dealing-with-missing-data-in-r-omit-approx-or-spline-part-1/

#https://dplyr.tidyverse.org/dev/articles/colwise.html

#Interpolation causes problems with final rows that repeat last value.

interpolate <- function(x){zoo::na.locf(x, na.rm = FALSE, fromLast = FALSE, type = "l", maxgap = Inf)}

TargetDataWide <- TargetDataWide %>%
  group_by(Tube) %>%
  arrange(Filename, MC, Tube, time) %>%
  mutate(across(.cols = starts_with("actinic-lights.light"), .fns = interpolate)) %>%
  ungroup() %>%
  mutate(actinic_par = rowSums(.[grep("actinic-lights.light", names(.))], na.rm = TRUE)) %>%
  filter(!is.na(actinic_par))

#drop original actinic light columns?
# select(!contains("actinic-lights.light")) 
```

Now that actinic_par is aligned with each row, coalesce od-sensors-X.od-720 and od-sensors-X.od-680 into 2 columns, b/c 'Tube' is already a column, so no need to identify tube X in od-sensors-X.od-680 columns.
This might cause problems later matching OD measures to actinic light colours.
Filter out rows where OD <= 0

```r
TargetDataWide <- TargetDataWide  %>%
   mutate(OD680 = rowSums(.[grep("od-680", names(.))], na.rm = TRUE),
          OD720 = rowSums(.[grep("od-720", names(.))], na.rm = TRUE)) %>%
   select(!contains("od-sensors")) %>%
  filter(OD680 > 0,
         OD720 > 0)
```


# Merge Data with meta data

```r
#This generates 'NA' values for ~1,000,000 rows of 'O2'; possibly the temperature rows?
TargetDataMeta <- left_join(x = TargetDataWide, y= Catalog, by = c("ExpDate", "MC", "Tube"), suffix  = c("_multi", "_cat"))

TargetDataMeta
```

```
## # A tibble: 9,349 x 55
##     time Filename CDateTime            Tube abs_time              ToD ExpDate   
##    <dbl> <chr>    <dttm>              <dbl> <dttm>              <dbl> <date>    
##  1  10.8 MultiCo… 2021-05-25 15:53:08     1 2020-01-24 07:05:01  16.8 2020-01-24
##  2  11   MultiCo… 2021-05-25 15:53:08     1 2020-01-24 07:06:01  17   2020-01-24
##  3  11.1 MultiCo… 2021-05-25 15:53:08     1 2020-01-24 07:06:31  17.1 2020-01-24
##  4  11.2 MultiCo… 2021-05-25 15:53:08     1 2020-01-24 07:07:01  17.2 2020-01-24
##  5  11.2 MultiCo… 2021-05-25 15:53:08     1 2020-01-24 07:07:31  17.2 2020-01-24
##  6  11.8 MultiCo… 2021-05-25 15:53:08     1 2020-01-24 07:10:31  17.8 2020-01-24
##  7  15.2 MultiCo… 2021-05-25 15:53:08     1 2020-01-24 07:31:01  21.2 2020-01-24
##  8  29.7 MultiCo… 2021-05-25 15:53:08     1 2020-01-24 08:58:01  11.7 2020-01-24
##  9  29.8 MultiCo… 2021-05-25 15:53:08     1 2020-01-24 08:58:31  11.8 2020-01-24
## 10  29.8 MultiCo… 2021-05-25 15:53:08     1 2020-01-24 08:59:01  11.8 2020-01-24
## # … with 9,339 more rows, and 48 more variables: MC <chr>,
## #   `actinic-lights.light-1-B450` <dbl>, `actinic-lights.light-2-B450` <dbl>,
## #   `actinic-lights.light-2-B470` <dbl>, `actinic-lights.light-3-B450` <dbl>,
## #   `actinic-lights.light-4-B450` <dbl>, `actinic-lights.light-5-R660` <dbl>,
## #   `actinic-lights.light-6-R660` <dbl>, `actinic-lights.light-7-R660` <dbl>,
## #   `actinic-lights.light-8-R660` <dbl>, actinic_par <dbl>, OD680 <dbl>,
## #   OD720 <dbl>, PrimaryOperator <chr>, Run <dbl>, Description <chr>,
## #   Motivation <chr>, doi <chr>, ProjectID <chr>, ID <chr>, Strain <chr>,
## #   Inoc_mL <dbl>, Media_mL <dbl>, ExpCul <dbl>, InnocDate <dbl>, Par_ue <dbl>,
## #   Photoperiod <dbl>, Temp_c <dbl>, Media <chr>, MediaDate <dbl>,
## #   Source <chr>, SourceSalinity <lgl>, Plate <lgl>, Well <lgl>, EndDate <lgl>,
## #   Salinity <lgl>, N2 <dbl>, O2 <dbl>, Ar <dbl>, CO2 <dbl>, WL <chr>,
## #   LightShape <chr>, Optode <dbl>, InocpH <lgl>, FinalpH <dbl>,
## #   Par_ueAdjusted <dbl>, DateOfAdjustment <dbl>,
## #   ElaspedHoursAtAdjustment <dbl>
```

## Create second preliminary plot for TargetData facetted by PAR, Strain & Tube
Overall trace of PAR/1000 in orange

```r
TargetDataPlotFacet <- TargetDataMeta %>% ggplot() +
  geom_point(aes(x = time, y = OD680, colour = as.factor(WL)), size = 0.1) +
  geom_point(aes(x = time, y = OD720), size = 0.01, alpha = 0.1, colour = "black") +
  geom_point(aes(x = time, y = actinic_par/1000),  colour = "orange", size = 0.0001) +
  scale_x_continuous(breaks=seq(0, 250, by = 125)) +
  coord_cartesian(xlim = c(-10, 260)) +
  scale_colour_manual(values = MCMIXColours) +
  labs(y = "Optical Density (OD)", x = "Elapsed Time (h)", subtitle = "Growth Light (µE); Strain; Tube") +
  facet_grid(rows = vars(as.factor(O2)), cols = vars(as.factor(Par_ue),Strain, as.factor(Tube))) +
  theme_bw() +  
  labs(colour = "Actinic PAR (nm)")

TargetDataPlotFacet
```

![](Figs/secondprelimplot-1.png)<!-- -->

#Run this chunk if some tubes have poor data
##This is bad practice

```r
# TargetDataMeta <- TargetDataMeta %>%
#   filter(Tube != ) %>% filter(Tube != 8)
```



## Save second preliminary plot if desired

```r
# ggsave(file = file.path(PlotsPath, paste(TargetFileName, "TargetDataPlotFacet",".png",sep = "")), plot = TargetDataPlotFacet, device = NULL, scale = 1, height=10, width= 20, units = c("cm"),dpi = 300, limitsize = TRUE)
```


# Filter OD outliers
There are packages to do this; we wrote homebuilt code.

```r
#Lag screen; too greedy with small points?
# Screen <- 10
# 
# MultiDataTargetMetaFilterTest <- MultiDataTargetMeta  %>%
#   group_by(Filename, Tube) %>%
#   arrange(Filename, Tube, time) %>%
#   mutate(IsLagOutlier = if_else(OD680 < lag(OD680)/Screen | OD680 > lag(OD680) * Screen, 1, 0))
#          
# print(paste("Outliers caught by IsLagOutlier filter:", 
#             sum(MultiDataTargetMetaFilterTest$IsLagOutlier, na.rm = TRUE)))
# 
# MultiDataTargetMetaFilterTest <- MultiDataTargetMetaFilterTest %>%
#   filter(!is.na(IsLagOutlier),
#          IsLagOutlier !=1)
# 
# MultiDataTargetMetaFilterTest  %>% 
#   ggplot() +
#   geom_point(aes(x = time, y = OD680, colour = as.factor(WL)), size = 0.1) +
#   scale_colour_manual(values = MCMIXColours) +
#   geom_point(aes(x = time, y = actinic_par/1000), colour = "orange", size = 0.1) +
#   facet_grid(cols = vars(as.factor(O2)), rows = vars(as.factor(Tube))) +
#   theme_bw()

#moving average screen
MovAvgScreen <- 1.1
MovAvgWindow <- 65


TargetDataMetaFilter <- TargetDataMeta %>%
  group_by(Filename, MC, Tube) %>%
  arrange(Filename, MC, Tube, time) %>%
  mutate(MovAvg = rollmean(OD680, MovAvgWindow, fill = "extend"),
         IsMovAvgOutlier = if_else((OD680 > MovAvgScreen*MovAvg) | (OD680 < MovAvg/MovAvgScreen), 1, 0)) %>%
  mutate(MovAvg = rollmean(OD720, MovAvgWindow, fill = "extend"),
         IsMovAvgOutlier = if_else((OD720 > MovAvgScreen*MovAvg) | (OD720 < MovAvg/MovAvgScreen), 1, 0))
  

print(paste("Outliers caught by IsMovAvgOutlier filter:", 
            sum(TargetDataMetaFilter$IsMovAvgOutlier, na.rm = TRUE)))
```

```
## [1] "Outliers caught by IsMovAvgOutlier filter: 709"
```

```r
TargetDataMetaFilter <- TargetDataMetaFilter %>%
  filter(!is.na(IsMovAvgOutlier),
         IsMovAvgOutlier  !=1)
```



```r
#manual coordinates for annotation; better to do it relative to data streams
OD680_x = 130
OD680_y = 0.2

OD720_x = 240
OD720_y = 0.22

Light_x = 230
Light_y = 0.05



TargetDataMetaFilterPlot <- TargetDataMetaFilter %>% 
  ggplot() +
  geom_point(aes(x = time, y = OD680, colour = as.factor(WL)), size = 0.1) +
  geom_point(aes(x = time, y = OD720), size = 0.01, alpha = 0.1, colour = "black") +
  geom_point(aes(x = time, y = actinic_par/1000),  colour = "orange", size = 0.0001) +
  scale_x_continuous(breaks=seq(0, 250, by = 125)) +
  coord_cartesian(xlim = c(-10, 260)) +
  scale_colour_manual(values = MCMIXColours) +
  labs(y = "Optical Density (OD)", x = "Elapsed Time (h)", subtitle = "Growth Light (µE); Strain; Tube") +
  facet_grid(rows = vars(as.factor(O2)), cols = vars(as.factor(Par_ue),Strain, as.factor(Tube))) +
  theme_bw() +  
  labs(colour = "Actinic PAR (nm)")

TargetDataMetaFilterPlot
```

![](Figs/filter data plots-1.png)<!-- -->

```r
TargetDataMetaFilterExpandPlot <- TargetDataMetaFilter %>% 
  filter(Tube == 1) %>%
  ggplot() +
  geom_point(aes(x = time, y = OD680, colour = as.factor(WL)), size = 0.1) +
  geom_point(aes(x = time, y = OD720), size = 0.1, alpha = 0.1, colour = "black") +
  scale_colour_manual(values = MCMIXColours) +
  geom_point(aes(x = time, y = actinic_par/1000), colour = "orange",  size = 0.0001) +
  scale_x_continuous(breaks=seq(0, 250, by = 125)) +
  coord_cartesian(xlim = c(0, 250)) +
  labs(y = "Optical Density (OD)", x = "Elapsed Time (h)",subtitle = "Growth Light (µE); Strain; Tube") +
  facet_grid(rows = vars(O2), cols = vars(Par_ue, Strain, Tube)) +
  annotate(geom = "text", x = OD680_x, y = OD680_y, label = "OD680", size = 5, colour = "darkblue") +
  annotate(geom = "text", x = OD720_x, y = OD720_y, label = "OD720", size = 5, colour = "black") +
  annotate(geom = "text", x = Light_x, y = Light_y, label = "Light level", size = 5, colour = "orange") +
  theme_bw() +
  labs(colour = "Actinic PAR (nm)")

TargetDataMetaFilterExpandPlot 
```

![](Figs/filter data plots-2.png)<!-- -->

```r
# TargetDataPlotFacet <- TargetDataMeta %>% ggplot() +
#   geom_point(aes(x = time, y = OD680, colour = as.factor(WL)), size = 0.1) +
#   geom_point(aes(x = time, y = actinic_par/1000),  colour = "orange", size = 0.0001) +
#   scale_x_continuous(breaks=seq(0, 250, by = 125)) +
#   scale_colour_manual(values = MCMIXColours) +
#   labs(y = "Optical Density 680nm (OD680)", x = "Elapsed Time (h)") +
#   facet_grid(rows = vars(as.factor(O2)), cols = vars(as.factor(Par_ue),Strain, as.factor(Tube))) +
#   theme_bw() +  
#   labs(colour = "Actinic PAR (nm)")
```

## Save filtered preliminary plot if desired

```r
# ggsave(file = file.path(PlotsPath, paste(TargetFileName, "TargetDataMetaFilterPlot",".png",sep = "")), plot = TargetDataMetaFilterPlot, device = NULL, scale = 1, height=10, width= 20, units = c("cm"),dpi = 300, limitsize = TRUE)
# 
# ggsave(file = file.path(PlotsPath, paste(TargetFileName, "TargetDataMetaFilterExpandPlot",".png",sep = "")), plot = TargetDataMetaFilterExpandPlot, device = NULL, scale = 1, height=10, width= 20, units = c("cm"),dpi = 300, limitsize = TRUE)
```


# Save .Rds of Imported Data from TargetFile to pass to MultiColourDataProcessLog.Rmd

```r
saveRDS(object = TargetDataMetaFilter, file = file.path(ImportedData, paste(TargetFileName, "TargetDataMetaFilter.Rds",  sep = "_")), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
```




