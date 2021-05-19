---
title: "MicroGrowth"
author: "Douglas A. Campbell, Maximilian Berthold"
date: "2021-05-19"
output:
  html_document: 
    code_folding: hide
    keep_md: yes
    fig_caption: yes
bibliography: MicrobialGrowthFitting.bib    
---
# Acknowledgements
This report uses test data accumulated by the BIOL2201 Class, 2021, Mount Allison University under the guidance of L. Barney.

# Run These Chunks First


```r
knitr::opts_chunk$set(fig.path='Figs/')
```

Set variables for project

```r
# Project <- "salinity_growth"
# FileIDSingle <- "Ab"
# Replic <- "AllSal"
# DataFolder <- "ProcessData"
# 
# #grouping variables for rows & columns, and for nesting for fitting
# RowVar <- "salinity" #expected salinity
# RowVar2 <- "CalcSal"
# ColVar <- "strain"
# NestVar <- "par_ue"
# NestVarII <- "id"
# 
# #filtering variables
# ChlAb = c("680", "720", "750")
# Wavelength_nm <- c("680min720", "680min750", "680", "720", "750")
# Par_ue = 300
# Strains <- c("CZS25K")
# StartDate <- c("20200916", "20200923", "20201013")
# TissueTreated = FALSE
# Media = "BG11+Na2CO3"
# SourceSal = "4"
# StopHours = 200
# 
# #fitting variable
# FitVar = "CorrOD"


#Set your initials code for your MicroColonyGrowth data.
MyInitials <- c("AnMa")

ColonyGoogleSheet <- "https://docs.google.com/spreadsheets/d/1DWyASxmu5DURkcegWM0Sykp9x3F3eQfMxdz0VYJ9mOI/edit?usp=sharing"
```


Additional libraries needed to run the project.
'googledrive' and 'googlesheets4' manage access to googlesheets, used for MetaData catalogs


# Introduction
Microbes grow.
We will import microbial growth data, tidy it, plot growth data vs. time and fit models of microbial growth to the data.

# Some possible models for microbial growth
Each model defines a pattern of growth.  

## Linear
Microbial population increases by a constant amount per unit time.
(Pseudo)Linear growth is sometimes observed when a constant input of a resource (ex. light) is shared among an increasing number of cells.  
The population adds a constant number of cells per unit time, but each cell is growing more slowly with each time increment as its share of the constant resource input decreases.  

Two fitted parameters; Intercept (P0) and linear slope (m)
Pt = ehours*m + P0

### Linear Example
Having figures appear from .md in GitHub requires chunk names without spaces.  
https://stackoverflow.com/questions/51531228/unable-to-view-plots-on-md-file-in-github

```r
linear_eqn <- function(x, m, Intercept){(x * m) + Intercept
}

linear_eqn_test <- function(x, m = 0.005, Intercept = 5){x*m + Intercept
}

ggplot(data.frame(x=c(1, 3000)), aes(x=x)) + 
  stat_function(fun = linear_eqn_test) +
  labs(title = "Linear example") +
  theme_bw()
```

![](Figs/linearexample-1.png)<!-- -->


## Exponential  
Microbial population increases by a constant proportion per unit time.
Exponential growth is observed when microbial cell replication is not resource limited, but is rather limited by intrinsic processes like nutrient uptake rates or DNA replication rate.

Two fitted parameters; intercept (P0) and exponential growth rate constant (mu, µ)
Pt = P0 *(exp(mu*ehours))

### Exponential Example

```r
exp_eqn <- function(x, mu, Intercept){Intercept *exp(x*mu)
}
  
exp_eqn_test <- function(x, mu = 0.005, Intercept = 5){Intercept *exp(x*mu)
}
ggplot(data.frame(x=c(1, 3000)), aes(x=x)) + 
  stat_function(fun = exp_eqn_test) +
  labs(title = "Exponential example") +
  theme_bw()
```

![](Figs/exponentialexample-1.png)<!-- -->


## Exponential rise towards a plateau
The gap between a maximum population and the current microbial population decreases by a constant proportion per unit time.
Three parameters; intercept (P0), exponential growth rate constant (mu) and plateau (Max)
Pt = Max - ((Max-Intercept)*exp(x*-mu)

### Exponential rise towards a plateau example

```r
exp_rise_plat_eqn <- function(x, mu = 0.005, Intercept){Max - ((Max-Intercept)*exp(x*-mu))
}

exp_rise_plat_eqn_test <- function(x, mu = 0.005, Intercept = 0.5, Max = 1000){Max - ((Max-Intercept)*exp(x*-mu))
}
ggplot(data.frame(x=c(1, 3000)), aes(x=x)) + 
  stat_function(fun = exp_rise_plat_eqn_test) +
  labs(title = "Exponential Rise towards Plateau example") +
  theme_bw()
```

![](Figs/risetoplateauexample-1.png)<!-- -->


## Logistic
Microbial population initially increases exponentially but slows as population approaches a plateau, Pmax.

Pmax is maximum population tenable under the prevailing conditions, sometimes termed the carrying capacity.
mu dominates the early part of the rise, transitioning towards an asymptote towards Pmax.
Three parameters; intercept (P0), exponential rate constant (mu) and maximum population Pmax.
There are multiple equivalent formulations of the logistic equation.
Pt= {Pmax e^µt }/{Pmax + P0 (e^µt- 1)}

### Logistic example

```r
logistic_eqn <- function(x, mu, Intercept, Max){(Max*Intercept*exp(mu*x))/(Max + (Intercept*(exp(mu*x)-1)))
}

logistic_eqn_test <- function(x, mu = 0.005, Intercept = 5, Max = 10000){(Max*Intercept*exp(mu*x))/(Max + (Intercept*(exp(mu*x)-1)))
}

ggplot(data.frame(x = c(1, 3000)), aes(x=x)) + 
  stat_function(fun = logistic_eqn_test) +
  labs(title = "Logistic example") +
  theme_bw()
```

![](Figs/logisticexample-1.png)<!-- -->


## Modified Gompertz Fit with Lag
Gompertz equations [@zwieteringModelingBacterialGrowth1990] have shapes similar to Logistic fits but are more readily modified to include 'lag' as a parameter for an initial time period over which cell count does not (measurably) increase.
'lag' may be a true biological lag representing cellular acclimation to new growth conditions, and/or may result from cell counts below the level of instrumental resolution, giving no detectable growth.

The Gompertz equations are based on fitting the ln of (data normalized to an 'Intercept').
Therefore the plots start at '0'; ln(Intercept/Intercept) = ln(1) = 0.
The reported value for Amu is comparable to mu from a Logistic or Exponential fit.
But the reported value for Amax = ln(Max/Intercept).
To retrieve an estimate for Max we need to take an antilog and multiply through by the Intercept.
Chosing what value to use as an intercept for the Gompertz equations is the main challenge.  
The first value or a single minimum value may be aberrant.
An average of kth(minimum values) is possible.
The advantage of the normalization is that the Intercept parameter is no longer fit, which improves statistical power for the fitting of additional parameters such as Lag.


Pt = ModGompertzEqn <- (Amax*(exp(-exp((Amu*exp(1))/Amax*(Lag-x)+1))))
Pt = ModGompertzEqn <- ((ln(Max/Intercept)*(exp(-exp((Amu*exp(1))/ln(Max/Intercept)*(Lag-x)+1))))

ModGompertzEqn <- function(Amax,Amu,Lag,x){(Amax*(exp(-exp((Amu*exp(1))/Amax*(Lag-x)+1))))}

### Modified Gompertz Fit with Lag example

```r
#natural log(RFU/min(RFU)); works, but confusing
mod_gomp_lag_eqn <-  function(x, Amax,Amu,Lag){(Amax*(exp(-exp((Amu*exp(1))/Amax*(Lag-x)+1))))
  }

mod_gomp_lag_eqn_test <-  function(x, Amax = 1000/10, Amu  = 0.5, Lag = 1000 ){(Amax*(exp(-exp((Amu*exp(1))/Amax*(Lag-x)+1))))
  }

ggplot(data.frame(x = c(1, 3000)), aes(x=x)) + 
  stat_function(fun = mod_gomp_lag_eqn_test) +
  labs(title = "Modified Gompertz with Lag example") +
  theme_bw()
```

![](Figs/modgompertzexample-1.png)<!-- -->


# Materials and Methods 
We grow cultures. To compare results across cultures we want to fit the growth data with models to get fitted parameters.


## Import data. 
The import function may open a Google log in to get an authorization code in a separate browser window. Copy/paste the authorization code into the 'Console' window below.
See next chunk for manual alternative if needed; download .csv from GoogleSheet, then upload the .csv to RStudio.cloud.


```r
gs4_deauth()
ColonyData <- read_sheet(ColonyGoogleSheet)

ColonyData
```

```
## # A tibble: 1,087 x 10
##    Initials_4letter YYYYMMDDHHMM Substrate Treatment col1mm2 col2mm2 col3mm2
##    <chr>            <list>       <chr>     <chr>       <dbl>   <dbl>   <dbl>
##  1 JuDe             <dbl [1]>    WhiteRoll Control      7.28   NA      NA   
##  2 JuDe             <dbl [1]>    WhiteRoll Control     13.0    NA      NA   
##  3 JuDe             <dbl [1]>    WhiteRoll Control     14.4    NA      NA   
##  4 JuDe             <dbl [1]>    WhiteRoll Control     18.9    NA      NA   
##  5 JuDe             <dbl [1]>    WhiteRoll Control     26.1    NA      NA   
##  6 JuDe             <dbl [1]>    WhiteRoll Control     30.0    NA      NA   
##  7 JuDe             <dbl [1]>    WhiteRoll Control     40.9    NA      NA   
##  8 JuDe             <dbl [1]>    WhiteRoll Control     44.4    NA      NA   
##  9 JuDe             <dbl [1]>    WhiteRoll Control     59.5     1.72    7.25
## 10 JuDe             <dbl [1]>    WhiteRoll Control     64.6     1.90   12.4 
## # … with 1,077 more rows, and 3 more variables: Temp_C <dbl>, ...9 <lgl>,
## #   ...10 <dbl>
```



```r
# ImportGrowth <- readRDS(file = file.path(DataFolder,paste(Project,FileIDSingle, "GrowthLong.Rds", sep  = "_"),fsep = .Platform$file.sep))
# 
# #https://stackoverflow.com/questions/64111558/passing-variables-into-the-names-glue-parameter-of-tidyrpivot-wider
# 
# GrowthLong <- ImportGrowth %>%
#   filter(Wavelength %in% ChlAb) %>%
#   pivot_wider(id_cols = -c(OD, AvBlankOD), names_from = Wavelength, values_from = CorrOD, names_glue = "nm_{Wavelength}") %>% 
#   mutate("680min720" = nm_680 - nm_720, "680min750" = nm_680 - nm_750) %>%
#   pivot_longer(cols = c(nm_680, nm_720, nm_750, "680min720", "680min750"), values_to = "CorrOD", names_to = "Wavelength") %>%
#   mutate(Wavelength = str_remove(Wavelength, "nm_")) %>%
#   mutate(CorrOD = if_else(CorrOD < 0.001, 0.001, CorrOD)) #experimentally test to set values below 0.001 to 0.001 (~detection limit device)
# 
# 
# GrowthLong <- GrowthLong %>%
#   group_by(!!sym(NestVar), !!sym(RowVar2), !!sym(ColVar), !!sym(NestVarII), Wavelength, exp_date) %>%
#   #mutate(OD = OD - 0.08) #rough approach to correct for media blank, replace with true media-blank once script runs
#   mutate(logODminOD = log(CorrOD/min(CorrOD))) %>% #problem with log leads to NaN, when value is very small
#   ungroup()
```



Alternative
Doug: Manually Download Data from Googlesheets as .csv.
Doug: Manually Upload .csv to folder containing .Rproj
Read in .csv into ColonyData object.min

```r
# ColonyData <- read_csv("MicroColonyData2021.csv", col_names = TRUE)
# ColonyData
```


Example Error correction chunk if some data columns have converted to 'list' format upon import by 'plucking' the first entry in each cell in each column and converting to numeric format.
Could make this more sophisticated.

```r
# ColonyData <- ColonyData %>%
#   mutate (YYYYMMDDHHMM = as.numeric(as.character(pluck(col1mm2, 1))),
#           col1mm2_test = as.numeric(as.character(pluck(col1mm2, 1))),
#           col2mm2_test = as.numeric(as.character(pluck(col2mm2, 1))),
#           col3mm2_test = as.numeric(as.character(pluck(col3mm2, 1))))
# 
# ColonyData
```

Filter the class data to include only 'MyData' (defined by 'MyInitials' set above). 
Convert the data to long format for analyses.
Convert the YYYYMMDDHHMM data column to a formatted datetime column,
Generate a numeric ehours column for elapsed time.
Add a ln(normalized area column) for Gompertz fitting

```r
#filter data,long format for analyses
MyData <- ColonyData %>%
  filter(Initials_4letter == MyInitials) %>%
  pivot_longer(cols = colnames(ColonyData[5:7]), names_to = 'Colonies', values_to = "Areamm2") %>%
  mutate(Areamm2 = as.numeric(Areamm2))

#convert YYYYMMDDHHMM to datetime format
MyData <- MyData %>%
  group_by(Substrate, Treatment) %>%
  mutate(YYYYMMDDHHMM = ymd_hm(YYYYMMDDHHMM)) %>%
  mutate(ehours = as.numeric(((YYYYMMDDHHMM - min(YYYYMMDDHHMM, na.rm = TRUE)))/3600)) %>%
  ungroup()

#ln(Area/minArea)
MyData <- MyData %>%
  group_by(Substrate, Treatment) %>%
  mutate(lnAreaNorm = log(Areamm2/min(Areamm2, na.rm = TRUE))) %>%
  filter(!is.na(Areamm2)) %>%
  ungroup()
```


# Results  
Upload representative colony images to your RStudio.cloud.
Add representative images of your colony growth over time to your .Rmd.
Will work with .JPG and .png; not sure about other graphic formats.

```r
knitr::include_graphics(file.path("202101240946.jpg"))
```

<img src="202101240946.jpg" width="50%" height="70%" />


## Example of Mould Colony on Naan, 202101240946.JPG

## Plot the data with a separate facet for each colony.
Or, plot the data with all colonies under a common combination of substrate and treatment in a single panel.

```r
#set Y axis limit range
#YLIM = c(0,max(MyData$Areamm2, na.rm = TRUE))

MyData %>%
  ggplot() +
  geom_point(aes(x = ehours, y = Areamm2, colour = Colonies)) +
  facet_grid(cols = vars(Treatment), rows = vars(Substrate)) +
  theme_bw() +
  labs(title = "Colony Growth",
       subtitle = MyInitials)
```

![](Figs/ColonyPlot-1.png)<!-- -->

```r
MyData %>%
  ggplot() +
  geom_point(aes(x = ehours, y = Areamm2, colour = Colonies)) +
  facet_grid(cols = vars(Colonies), rows = vars(Substrate, Treatment)) +
  theme_bw() +
  labs(title = "Colony Growth",
       subtitle = MyInitials)
```

![](Figs/ColonyPlot-2.png)<!-- -->



```r
  # filter(strain %in% Strains) %>%
  # filter(exp_date %in% StartDate) %>%
  # filter(Wavelength %in% Wavelength_nm) %>%
  # #filter(tissue_treated_plate == TissueTreated) %>%
  # #filter(media == Media) %>%
  # filter(SourceSal == source_salinity) %>%
  # ggplot()+
  # geom_point(aes(x = E_hours, y = !!sym(FitVar), colour = as.factor(Wavelength))) +
  # coord_cartesian(ylim = c(0, 2)) +
  # facet_grid(cols = vars(!!sym(RowVar2)), rows = vars(!!sym(NestVar))) +
  #   theme_bw() +
  # labs(caption = paste(Wavelength_nm, "nm; PAR", Par_ue, "uE"))
```


```r
#set filtering variables at top
#pick filtering option as either a simple threshold or a factor change between OD and OD_lead
# 
# Screen_OD <- 2
# 
# FitData <- as_tibble(GrowthLong) %>%
#   filter(strain %in% Strains) %>%
#   filter(exp_date %in% StartDate) %>% #if multiple dates are applied, use %in% instead of ==
#   filter(Wavelength %in% Wavelength_nm) %>%
#   #filter(tissue_treated_plate == TissueTreated) %>%
#   #filter(media == Media) %>%
#   filter(source_salinity == SourceSal) %>%
#   group_by(well) %>%
#   arrange(well, E_hours) %>%
#   mutate(OD_lead = lead(CorrOD)) %>%
#   #filter(CorrOD < (OD_lead / Screen_OD)) %>%
#   filter(CorrOD < Screen_OD) %>%
#   filter(E_hours < StopHours)
# 
# test <- FitData %>%
#   group_by(!!sym(NestVar), !!sym(RowVar2), !!sym(RowVar), !!sym(ColVar), Wavelength, exp_date) %>%
#   mutate(LagSeed = E_hours[which.min(CorrOD)])
```




```r
# FitData %>%
#   filter(strain %in% Strains) %>%
#   filter(exp_date %in% StartDate) %>%
#   filter(Wavelength %in% Wavelength_nm) %>%
#   #filter(tissue_treated_plate == TissueTreated) %>%
#   filter(media == Media) %>%
#   ggplot()+
#   geom_point(aes(x = E_hours, y = !!sym(FitVar), colour = as.factor(Wavelength))) +
#   coord_cartesian(ylim = c(0, 2)) +
#   facet_grid(cols = vars(!!sym(RowVar2)), rows = vars(!!sym(NestVar))) +
#     theme_bw() +
#   labs(caption = paste(Wavelength_nm, "nm; PAR", Par_ue, "uE"))
```


# Colony Size vs. Hours Elapsed since First Measure.

Plot the ln(Areamm2) to check whether growth follows an exponential pattern
Default log in R is log base E, not log base 10.

```r
MyData %>%
  ggplot() +
  geom_point(aes(x = ehours, y = log(Areamm2), colour = Colonies)) +
  facet_grid(cols = vars(Colonies), rows = vars(Substrate, Treatment)) +
  theme_bw() +
  labs(title = "Colony Growth",
       subtitle = MyInitials)
```

![](Figs/ln data plot-1.png)<!-- -->


## Let R find the best fit parameters for different models of the data.

Each user measured multiple colonies, possibly on multiple substrates with multiple treatments.
We create a 'nested' data frame where all data in columns YYYYMMDDHHMM, Areamm2, ehours, from a single user, a single colony, on a single substrate, with a single treatment are 'nested' together for later fitting of growth curves.

```r
MyData_nest <- as_tibble(MyData) %>%
  #filter(Treatment != "Hot") %>%
  nest(data = c(YYYYMMDDHHMM, Areamm2, lnAreaNorm, ehours, Temp_C))
```

## Fit and plot treatment colony specific linear growth trajectories using nest purrr:map & broom::augment
This chunk uses complicated code from the 'Tidyverse' package.
R will iteratively vary the parameters of the fitting equations to minimize the residuals (discrepancies) between the data points (nested) and the points predicted by the model with a given set of parameters.

```r
#define starting, lower, and upper bounds for fit parameters to constrain linear equations
linear_eqn_start <- list(m = 1, Intercept = min(MyData$Areamm2, na.rm = TRUE))
linear_eqn_lower <- c(0,0)

linear_eqn_upper <- c(max(MyData$Areamm2, na.rm = TRUE), max(MyData$Areamm2, na.rm = TRUE))

#, lower = linear_eqn_lower, upper = linear_eqn_upper

colony_lin  <- MyData_nest %>% 
  mutate(
  fit = map(data, possibly(~nlsLM(Areamm2 ~ linear_eqn(x = ehours, m, Intercept), data = .x, start = linear_eqn_start), otherwise = NULL)),
  predict = map(fit, possibly(augment, otherwise = NULL)),
  tidied =  map(fit, possibly(tidy, otherwise = NULL)),
  param = map(fit, possibly(glance, otherwise = NULL))
  )

colony_lin %>%
  unnest(predict)  %>% 
  ggplot() +  
  geom_line(aes(x = ehours, y = .fitted), size = 0.5) +
  geom_line(aes(x = ehours, y = .resid, alpha = 0.1), linetype = "dashed", colour = "red", show.legend = FALSE) +
  facet_grid(cols = vars(Colonies), rows = vars(Treatment, Substrate)) +
  geom_point(data = MyData, aes(x = ehours, y = Areamm2), colour = "green") +
  theme_bw() +
  labs(y = ".fitted Areamm^2",
       title = "Colony Growth, Linear Fits",
       subtitle = MyInitials,
       caption = "Linear line, residuals dashed red line")
```

![](Figs/colonylineargrowth-1.png)<!-- -->

```r
#parameters of colony specific linear fits
colony_lin %>% 
  unnest(tidied) %>%
  select(-c(data, fit, predict, param)) %>%
  mutate_if(is.numeric, round, digits = 2)
```

```
## # A tibble: 10 x 11
##    Initials_4letter Substrate Treatment ...4   ...5 Colonies term  estimate
##    <chr>            <chr>     <chr>     <lgl> <dbl> <chr>    <chr>    <dbl>
##  1 AnMa             Cranewoo… Control   NA       NA col1mm2  m         4.77
##  2 AnMa             Cranewoo… Control   NA       NA col1mm2  Inte…   -99.1 
##  3 AnMa             Cranewoo… Control   NA       NA col2mm2  m         2.91
##  4 AnMa             Cranewoo… Control   NA       NA col2mm2  Inte…  -123.  
##  5 AnMa             Cranewoo… Dark      NA       NA col1mm2  m         1.57
##  6 AnMa             Cranewoo… Dark      NA       NA col1mm2  Inte…    13.6 
##  7 AnMa             Cranewoo… Dark      NA       NA col2mm2  m         1.85
##  8 AnMa             Cranewoo… Dark      NA       NA col2mm2  Inte…     8.5 
##  9 AnMa             Cranewoo… Dark      NA       NA col3mm2  m         1.39
## 10 AnMa             Cranewoo… Dark      NA       NA col3mm2  Inte…    50.4 
## # … with 3 more variables: std.error <dbl>, statistic <dbl>, p.value <dbl>
```

## Fit and plot treatment colony specific exponential growth trajectories using nest purrr:map & broom::augment

```r
#define starting, lower, and upper bounds for fit parameters to constrain exponential fit
exp_eqn_start <- list(mu = 0.01, Intercept = min(MyData$Areamm2, na.rm = TRUE))
exp_eqn_lower <- c(0,0)

exp_eqn_upper <- c(1, max(MyData$Areamm2, na.rm = TRUE))
#, lower = exp_eqn_lower, upper = exp_eqn_upper

colony_exp <- MyData_nest %>% 
  mutate(
  fit = map(data, possibly(~nlsLM(Areamm2 ~ exp_eqn(x = ehours, mu, Intercept), data = .x, start = exp_eqn_start), otherwise = NULL)),
  predict = map(fit, possibly(augment, otherwise = NULL)),
  tidied =  map(fit, possibly(tidy, otherwise = NULL)),
  param = map(fit,possibly(glance, otherwise = NULL))
  )

colony_exp %>%
  unnest(predict)  %>% 
  ggplot() +  
  geom_line(aes(x = ehours, y = .fitted), size = 0.5) +
   geom_line(aes(x = ehours, y = .resid, alpha = 0.1), linetype = "dashed", colour = "red", show.legend = FALSE) +
  facet_grid(cols = vars(Colonies), rows = vars(Treatment, Substrate)) +
  geom_point(data = MyData, aes(x = ehours, y = Areamm2), colour = "green") +
  theme_bw() +
  labs(y = ".fitted Areamm^2",
       title = "Colony Growth, Exponential Fits",
       subtitle = MyInitials,
       caption = "Exponential line, residuals dashed red line")
```

![](Figs/colonyexponentialgrowth-1.png)<!-- -->

```r
#parameters of colony specific exponential fits
colony_exp %>% 
  unnest(tidied) %>%
  select(-c(data, fit, predict, param)) %>%
  mutate_if(is.numeric, round, digits = 2)
```

```
## # A tibble: 10 x 11
##    Initials_4letter Substrate Treatment ...4   ...5 Colonies term  estimate
##    <chr>            <chr>     <chr>     <lgl> <dbl> <chr>    <chr>    <dbl>
##  1 AnMa             Cranewoo… Control   NA       NA col1mm2  mu        0.03
##  2 AnMa             Cranewoo… Control   NA       NA col1mm2  Inte…    26.2 
##  3 AnMa             Cranewoo… Control   NA       NA col2mm2  mu        0.03
##  4 AnMa             Cranewoo… Control   NA       NA col2mm2  Inte…     6.85
##  5 AnMa             Cranewoo… Dark      NA       NA col1mm2  mu        0.03
##  6 AnMa             Cranewoo… Dark      NA       NA col1mm2  Inte…    21.5 
##  7 AnMa             Cranewoo… Dark      NA       NA col2mm2  mu        0.04
##  8 AnMa             Cranewoo… Dark      NA       NA col2mm2  Inte…    18.4 
##  9 AnMa             Cranewoo… Dark      NA       NA col3mm2  mu        0.02
## 10 AnMa             Cranewoo… Dark      NA       NA col3mm2  Inte…    54.7 
## # … with 3 more variables: std.error <dbl>, statistic <dbl>, p.value <dbl>
```

## Fit and plot treatment colony specific logistic growth trajectories using nest purrr:map & broom::augment

```r
#define starting, lower, and upper bounds for fit parameters to constrain logistic fit
logistic_eqn_start<-list(mu = 0.05, Intercept = min(MyData$Areamm2, na.rm = TRUE), Max = max(MyData$Areamm2, na.rm = TRUE))

logistic_eqn_lower <- c(0.001, (min(MyData$Areamm2, na.rm = TRUE) * 0.1), (max(MyData$Areamm2, na.rm = TRUE) * 0.5))

logistic_eqn_upper <- c(1, (min(MyData$Areamm2, na.rm = TRUE) * 4), (max(MyData$Areamm2, na.rm = TRUE) * 2))

#, lower = logistic_eqn_lower, upper = logistic_eqn_upper
colony_log <- MyData_nest  %>% 
  #filter(Colonies == "col1mm2") %>%
  mutate(
  fit = map(data, possibly(~nlsLM(Areamm2 ~ logistic_eqn(x = ehours, mu, Intercept, Max), data = .x, start = logistic_eqn_start), otherwise = NULL)),
  predict = map(fit, possibly(augment, otherwise = NULL)),
  tidied =  map(fit, possibly(tidy, otherwise = NULL)),
  param = map(fit,possibly(glance, otherwise = NULL))
  )

colony_log %>% 
  unnest(predict) %>%
  ggplot() +  
  geom_line(aes(x = ehours, y = .fitted), size = 0.5) +
   geom_line(aes(x = ehours, y = .resid, alpha = 0.1), linetype = "dashed", colour = "red", show.legend = FALSE) +
  geom_point(data = MyData, aes(x = ehours, y = Areamm2), colour = "green") +
  facet_grid(cols = vars(Colonies), rows = vars(Treatment, Substrate)) +
  theme_bw() +
  labs(y = ".fitted Areamm^2",
       title = "Colony Growth, Logistic Fits",
       subtitle = MyInitials,
       caption = "Logistic line, residuals dashed red line")
```

![](Figs/colonylogisticgrowth-1.png)<!-- -->

```r
#parameters of colony specific logistic fits
colony_log %>%
  unnest(tidied) %>%
  select(-c(data, fit, predict, param)) %>%
  mutate_if(is.numeric, round, digits = 2)
```

```
## # A tibble: 15 x 11
##    Initials_4letter Substrate Treatment ...4   ...5 Colonies term  estimate
##    <chr>            <chr>     <chr>     <lgl> <dbl> <chr>    <chr>    <dbl>
##  1 AnMa             Cranewoo… Control   NA       NA col1mm2  mu     1.00e-1
##  2 AnMa             Cranewoo… Control   NA       NA col1mm2  Inte…  1.90e-1
##  3 AnMa             Cranewoo… Control   NA       NA col1mm2  Max    4.84e+2
##  4 AnMa             Cranewoo… Control   NA       NA col2mm2  mu     1.80e-1
##  5 AnMa             Cranewoo… Control   NA       NA col2mm2  Inte…  0.     
##  6 AnMa             Cranewoo… Control   NA       NA col2mm2  Max    2.04e+2
##  7 AnMa             Cranewoo… Dark      NA       NA col1mm2  mu     3.00e-2
##  8 AnMa             Cranewoo… Dark      NA       NA col1mm2  Inte…  2.15e+1
##  9 AnMa             Cranewoo… Dark      NA       NA col1mm2  Max   -5.04e+8
## 10 AnMa             Cranewoo… Dark      NA       NA col2mm2  mu     4.00e-2
## 11 AnMa             Cranewoo… Dark      NA       NA col2mm2  Inte…  1.84e+1
## 12 AnMa             Cranewoo… Dark      NA       NA col2mm2  Max   -1.02e+9
## 13 AnMa             Cranewoo… Dark      NA       NA col3mm2  mu     2.00e-2
## 14 AnMa             Cranewoo… Dark      NA       NA col3mm2  Inte…  5.47e+1
## 15 AnMa             Cranewoo… Dark      NA       NA col3mm2  Max   -3.88e+8
## # … with 3 more variables: std.error <dbl>, statistic <dbl>, p.value <dbl>
```

# In progress; add modified Gompertz with lag

## Fit and plot treatment colony specific modified Gompertz with lag growth trajectories using nest purrr:map & broom::augment.  
Note fit is run on lnAreaNorm

```r
mod_gomp_lag_eqn_start <- list(Amax= max(MyData$lnAreaNorm, na.rm = TRUE), Amu = 0.005, Lag = 10)

mod_gomp_lag_eqn_lower <- c((max(MyData$Areamm2, na.rm = TRUE) * 0.25), 0.001, 0)
mod_gomp_lag_eqn_upper <- c((max(MyData$Areamm2, na.rm = TRUE) * 2), 0.1, max(MyData$ehours, na.rm = TRUE))

logistic_eqn_upper <- c(1, (min(MyData$Areamm2, na.rm = TRUE) * 4), (max(MyData$Areamm2, na.rm = TRUE) * 2))

colony_modgompertzlag <- MyData_nest  %>% 
  #filter(Colonies == "col1mm2") %>%
  mutate(
  fit = map(data, possibly(~nlsLM(lnAreaNorm ~ mod_gomp_lag_eqn(x = ehours, Amax, Amu, Lag), data = .x, start =  mod_gomp_lag_eqn_start), otherwise = NULL)),
  predict = map(fit, possibly(augment, otherwise = NULL)),
  tidied =  map(fit, possibly(tidy, otherwise = NULL)),
  param = map(fit,possibly(glance, otherwise = NULL))
  )

colony_modgompertzlag %>% 
  unnest(predict) %>%
  ggplot() +  
  geom_line(aes(x = ehours, y = .fitted), size = 0.5) +
   geom_line(aes(x = ehours, y = .resid, alpha = 0.1), linetype = "dashed", colour = "red", show.legend = FALSE) +
  geom_point(data = MyData, aes(x = ehours, y = lnAreaNorm), colour = "green") +
  facet_grid(cols = vars(Colonies), rows = vars(Treatment, Substrate)) +
  theme_bw() +
  labs(y = ".fitted Areamm^2",
       title = "Colony Growth, Modified Gompertz with Lag",
       subtitle = MyInitials,
       caption = "Modified Gompertz with Lag, residuals dashed red line")
```

![](Figs/modgompertzlaggrowth-1.png)<!-- -->

```r
#parameters of colony specific logistic fits
colony_modgompertzlag %>%
  unnest(tidied) %>%
  select(-c(data, fit, predict, param)) %>%
  mutate_if(is.numeric, round, digits = 2)
```

```
## # A tibble: 15 x 11
##    Initials_4letter Substrate Treatment ...4   ...5 Colonies term  estimate
##    <chr>            <chr>     <chr>     <lgl> <dbl> <chr>    <chr>    <dbl>
##  1 AnMa             Cranewoo… Control   NA       NA col1mm2  Amax      6.04
##  2 AnMa             Cranewoo… Control   NA       NA col1mm2  Amu       0.05
##  3 AnMa             Cranewoo… Control   NA       NA col1mm2  Lag      -4.17
##  4 AnMa             Cranewoo… Control   NA       NA col2mm2  Amax      4.14
##  5 AnMa             Cranewoo… Control   NA       NA col2mm2  Amu       0.22
##  6 AnMa             Cranewoo… Control   NA       NA col2mm2  Lag      60.0 
##  7 AnMa             Cranewoo… Dark      NA       NA col1mm2  Amax      2.09
##  8 AnMa             Cranewoo… Dark      NA       NA col1mm2  Amu       0.12
##  9 AnMa             Cranewoo… Dark      NA       NA col1mm2  Lag       1.41
## 10 AnMa             Cranewoo… Dark      NA       NA col2mm2  Amax      2.21
## 11 AnMa             Cranewoo… Dark      NA       NA col2mm2  Amu       0.2 
## 12 AnMa             Cranewoo… Dark      NA       NA col2mm2  Lag       8.11
## 13 AnMa             Cranewoo… Dark      NA       NA col3mm2  Amax      2.53
## 14 AnMa             Cranewoo… Dark      NA       NA col3mm2  Amu       0.05
## 15 AnMa             Cranewoo… Dark      NA       NA col3mm2  Lag     -28.7 
## # … with 3 more variables: std.error <dbl>, statistic <dbl>, p.value <dbl>
```


## Compare Linear, Exponential, and Logistic Models using ANOVA
Doug to implement if feasible;  too slow, too long



## Organize all class data

```r
ClassData <- ColonyData %>%
  pivot_longer(cols = colnames(ColonyData[5:7]), names_to = 'Colonies', values_to = "Areamm2") %>%
  filter(Areamm2 != "NULL") %>%
  filter(!is.na(Areamm2)) %>%
  mutate(Areamm2 = as.numeric(Areamm2)) %>%
  mutate(YYYYMMDDHHMM = ymd_hm(YYYYMMDDHHMM)) %>%
  group_by(Initials_4letter) %>%
  mutate(ehours = as.numeric((YYYYMMDDHHMM - min(YYYYMMDDHHMM, na.rm = TRUE))/3600)) %>% 
  filter(ehours < 320)

ClassData <- ClassData %>%
  group_by(Initials_4letter, Colonies, Substrate, Treatment) %>%
  mutate(lnAreaNorm = log(Areamm2/min(Areamm2, na.rm = TRUE))) %>%
  filter(!is.na(Areamm2)) %>%
  ungroup()
```

## Plot class data facet by Substrate and Treatment

```r
ClassDataPlot <- ClassData %>%
  ggplot() +
  geom_point(aes(x = ehours, y = Areamm2, colour = Initials_4letter)) +
  facet_grid(rows = vars(Substrate), cols = vars(Treatment)) +
  theme_bw() +
  labs(title = "Class Colony Growth")

 ClassDataPlot 
```

![](Figs/classdataplot-1.png)<!-- -->

# Exploratory plot of all combinations of substrate and treatment run by the class.


```r
Class_nest <- as_tibble(ClassData) %>%
  nest(data = c(Temp_C, Initials_4letter, YYYYMMDDHHMM, Areamm2, lnAreaNorm, ehours, Colonies))
```

## Fit logistic model to pooled class data, nested by Substrate and Treatment
Some data 'nests' have too little data to fit or too little variation in the data to support a valid fit.  The range of minimum and maximum values across different Substrate/Treatment 'nests' also makes it difficult to assign starting, minimum and maximum values for parameters.
Doug added 'possibly' as an error catch that allows the fit of nests to proceed even if one nest is unable to be fit.

```r
#define starting values for fit parameters to constrain logistic fit
#logistic_eqn_start<-list(Max = max(ClassData$Areamm2, na.rm = TRUE), mu = 0.05, Intercept = min(ClassData$Areamm2, na.rm = TRUE))

logistic_eqn_start<-list(mu = 0.05, Intercept = 1, Max = 300)

#define lower, and upper bounds for fit parameters to constrain logistic fit
#logistic_eqn_lower<-c(50,0.001,1)

#logistic_eqn_upper<-c(1000,1,20)

#fit a pooled model using the logistic equation nested by temp_C

Class_log <- Class_nest %>%
  mutate(
   fit = map(data, possibly(~nlsLM(Areamm2 ~ logistic_eqn(x = ehours, mu, Intercept, Max), data = .x, start = logistic_eqn_start), otherwise = NULL)),
  predict = map(fit, possibly(augment, otherwise = NULL)),
  tidied =  map(fit, possibly(tidy, otherwise = NULL)),
  param = map(fit, possibly(glance, otherwise = NULL))
  )

# , lower = logistic_eqn_lower, upper = logistic_eqn_upper

Class_log %>%
  unnest(predict) %>%
  ggplot() +
  geom_point(aes(x = ehours, y = Areamm2)) +
  geom_line(aes(x = ehours, y = .fitted), size = 0.5) +
geom_line(aes(x = ehours, y = .resid, alpha = 0.1), linetype = "dashed", colour = "red", show.legend = FALSE) +
  facet_grid(cols = vars(Substrate), rows = vars(Treatment)) +
  theme_bw() +
  labs(y = "Areamm^2",
       title = "Colony Growth, Logistic Fits",
       caption = "Logistic line, residuals dashed red line, points data")
```

![](Figs/logmodelclass-1.png)<!-- -->

```r
#parameters of colony specific logistic fits
Class_log %>%
  unnest(tidied) %>%
  select(-c(data, fit, predict, param)) %>%
  mutate_if(is.numeric, round, digits = 2)
```

```
## # A tibble: 33 x 9
##    Substrate   Treatment ...3   ...4 term   estimate std.error statistic p.value
##    <chr>       <chr>     <lgl> <dbl> <chr>     <dbl>     <dbl>     <dbl>   <dbl>
##  1 WhiteRoll   Control   NA       NA mu         0.05      0.12      0.46    0.65
##  2 WhiteRoll   Control   NA       NA Inter…     6.82     16.3       0.42    0.68
##  3 WhiteRoll   Control   NA       NA Max       30.4       6.59      4.61    0   
##  4 WhiteRoll   Wet       NA       NA mu         0.09      0.13      0.73    0.48
##  5 WhiteRoll   Wet       NA       NA Inter…     0         0.06      0.06    0.95
##  6 WhiteRoll   Wet       NA       NA Max     1010.      705.        1.43    0.17
##  7 Naan        Dark      NA       NA mu         0.05      0        11.7     0   
##  8 Naan        Dark      NA       NA Inter…     0.81      0.39      2.09    0.04
##  9 Naan        Dark      NA       NA Max      835.       74.8      11.2     0   
## 10 SobeysBake… Control   NA       NA mu         0.06      0.03      2.34    0.03
## # … with 23 more rows
```


# Bibliography



