library(readxl)
library(tidyverse)
library(rcarbon)


#Producing Calibrated Dates
CD<-read_excel("data/MainData.xlsx")|> 
  filter(!is.na(`D Uncal`),
         !is.na(`Error`))

Caldates = calibrate(x= CD$`D Uncal`,
                     errors = CD$Error,
                     calCurves = "intcal20")

#sub-setting radiocarbon dates, calibrated probability greater than or equal to .5
which.CalDates(Caldates,
               BP<=6000&BP>=0,
               p=0.5)


#creating function to plot SPDs
SPDSplot <-
  function(x, y) {
    
    plot(x)
    plot(x,
         runm=200,
         add = TRUE,
         type="simple",
         col="darkorange",
         lwd=1.5,
         lty=2)
    title(main = y)
  }


#creating function to plot Monte-Carlo results
Expontentialplot <-
  function(x, y) {
    
    plot(x)
    title(main = y)
    legend( x = 
              "topright",
            box.col = "black",
            box.lwd = 1 ,  
            legend=c(
              "Observed SPD",
              "Simulation Envelope",
              "Positive Deviation",
              "Negative Deviation"),  
            fill = c(
              "black",
              "grey",
              "lightcoral",
              "cornflowerblue"))
    
  }

#creating function to plot CKDE results 
CKDEplot <- 
  function(x,y) {
    
    plot(x)
    title(main = y)
    legend( x = 
              "topleft",
            box.col = "black",
            box.lwd = 1 ,  
            legend=c(
              "Observed SPD",
              "Simulation Envelope"),  
            fill = c("black", "grey"))
  }






#Running analyses by culture 
analysis_by_culture <- 
  #Reading Data Into R
  read_excel("data/MainData.xlsx")|> 
  filter(!is.na(`D Uncal`),
         !is.na(`Error`)) |>
  #grouping by culture
  dplyr::group_by(
    `I or S`) |>
  #Calibrating Dates, Binning, and Adding DCal/bins Column
  dplyr::summarize( 
    #Visual binning value sensitivity test
    `binsensitivity` = 
      list(
        binsense(
          x = 
            rcarbon::calibrate(x = `D Uncal`,
                               errors = Error,
                               calCurves = "intcal20"),
          y = `Site`,
          h = seq(0,500,200),
          timeRange=c(6000,1000))
      ),
    #Binning Sites
    `bins` = 
      list(
        binPrep(
          sites= `Site`,
          ages= `D Uncal`,
          h=100)
      ),
    #Calibrating Dates
    `DCal`= 
      list(
        calibrate(
          x=`D Uncal`,
          errors=Error,
          calCurves='intcal20')
      ),
    #adding errors and site name columns - for later functions
    `errors` =
      list(Error),
    `Site Name` = 
      list(`Site Name`)
  )|>
  #Calculating SPDs with and without bins, and Adding SPDs Column
  dplyr::rowwise() |>
  dplyr::mutate(
    SPDs = list(
      spd(`DCal`,
          bins = `bins`,
          timeRange=c(6000,500)
      )
    ),
    #sensitivity test
    SPDsunbinned = 
      list(
        spd(
          `DCal`,
          timeRange=c(6000,500)
        )
      ),
    #testing SPDs against null hyp of exponential growth
    exptest = 
      list(
        modelTest(x = DCal,
                  errors=errors,
                  bins=bins,
                  nsim=1000,
                  timeRange=c(6000,500),
                  model="exponential",
                  runm=NA)
      ),
    # smaller timeframe - bias test
    exptest_timeframe = 
      list(
        modelTest(x = DCal,
                  errors=errors,
                  bins=bins,
                  nsim=1000,
                  timeRange=c(5000,500),
                  model="exponential",
                  runm=NA)
      ),
    #  random sampling dates for CKDE
    ckdesamples =
      list(
        sampleDates(DCal,
                    bins=bins,
                    nsim=1000,
                    verbose=FALSE,
                    boot = TRUE)
      ),
    #computing CKDE - bw=50 - small sample size 
    ckde = 
      list(
        ckde(
          ckdesamples,
          timeRange = c(6000,500),
          bw = 50)
      ),
    #CKDE with sampling by site-bias test 
    ckdes_samp_bins = 
      list(
        sampleDates(
          DCal,
          bins = `Site Name`,
          nsim = 1000,
          verbose = FALSE,
          boot = TRUE
        )
      ),
    ckde_samp_sites = 
      list(
        ckde(
          ckdes_samp_bins,
          timeRange = c(6000,500),
          bw = 50
        )
      )
  )


#plotting SPDS
testSPDplots =
  list (
    purrr::map2(
      .x = analysis_by_culture$SPDs,
      .y = analysis_by_culture$`I or S`,
      .f = SPDSplot)
  )
#plotting exp test
purrr::map2(
  .x=analysis_by_culture$exptest,
  .y = analysis_by_culture$`I or S`,
  .f = Expontentialplot
)
#plotting exp test - time range - sensitivity test
purrr::map2(
  .x=analysis_by_culture$exptest_timeframe,
  .y = analysis_by_culture$`I or S`,
  .f = Expontentialplot
)
#plotting CKDE test      
purrr::map2(
  .x = analysis_by_culture$ckde,
  .y = analysis_by_culture$`I or S`,
  .f = CKDEplot
)
#plotting CKDE - bias test      
purrr::map2(
  .x = analysis_by_culture$ckde_samp_sites,
  .y = analysis_by_culture$`I or S`,
  .f = CKDEplot
)




analysis_by_Region <- 
  #Reading Data Into R
  read_excel("data/MainData.xlsx")|> 
  filter(!is.na(`D Uncal`),
         !is.na(`Error`)) |>
  #grouping by region
  dplyr::group_by(
    `Region`) |>
  #Calibrating Dates, Binning, and Adding DCal/bins Column
  dplyr::summarize( 
    #visual bin value sensitivity test
    binsensitivity = 
      list(
        binsense(
          x = 
            rcarbon::calibrate(x = `D Uncal`,
                               errors = Error,
                               calCurves = "intcal20"),
          y = `Site`,
          h = seq(0,500,200),
          timeRange=c(6000,1000))
      ),
    #binning
    `bins` = 
      list(
        binPrep(
          sites= `Site`,
          ages= `D Uncal`,
          h=100)
      ),
    #Calibrating Dates
    `DCal`= 
      list(
        calibrate(
          x=`D Uncal`,
          errors=Error,
          calCurves='intcal20')
      ),
    #adding site name column - for later functions
    `Site Name` = 
      list(`Site Name`)
  )|>
  #Calculating SPDs with and without bins, and Adding SPDs Column
  dplyr::rowwise() |>
  dplyr::mutate(
    SPDs = list(
      spd(`DCal`,
          bins = `bins`,
          timeRange=c(6000,500)
      )
    ),
    #sensitivity test
    SPDsunbinned = 
      list(
        spd(
          `DCal`,
          timeRange=c(6000,500)
        )
      ),
    #random sampling for CKDE
    ckdesamples =
      list(
        sampleDates(DCal,
                    bins=bins,
                    nsim=1000,
                    verbose=FALSE,
                    boot = TRUE
        )
      ),
    #computing CKDE- BW of 50-small sample size 
    ckde = 
      list(
        ckde(
          ckdesamples,
          timeRange = c(6000,500),
          bw = 50)
      ),
    #CKDE samp sites -bias test 
    ckdes_samp_sites = 
      list(
        sampleDates(
          DCal,
          bins = `Site Name`,
          nsim = 1000,
          verbose = FALSE,
          boot = TRUE
        )
      ),
    ckde_region = 
      list(
        ckde(
          ckdes_samp_sites,
          timeRange = c(6000,500),
          bw = 50
        )
      )
  )

#plotting SPDS
purrr::map2(
  .x = analysis_by_Region$SPDs,
  .y = analysis_by_Region$`Region`,
  .f = SPDSplot
)
#plotting CKDE test 
purrr::map2(
  .x = analysis_by_Region$ckde,
  .y = analysis_by_Region$`Region`,
  .f = CKDEplot
)
#plotting CKDE - bias test 
purrr::map2(
  .x = analysis_by_Region$ckde_region,
  .y = analysis_by_Region$`Region`,
  .f = CKDEplot
)




analysis_by_site <- 
  #Reading Data Into R
  read_excel("data/MainData.xlsx")|> 
  filter(!is.na(`D Uncal`),
         !is.na(`Error`)) |>
  #grouping by site
  dplyr::group_by(
    `Site Name`) |>
  #Calibrating Dates, Binning, and Adding DCal/bins Column
  dplyr::summarize(
    #visual bin value sensitivity test
    `binsensitivity` = 
      list(
        binsense(
          x = 
            rcarbon::calibrate(x = `D Uncal`,
                               errors = Error,
                               calCurves = "intcal20"),
          y = `Site Name`,
          h = seq(0,500,200),
          timeRange=c(6000,1000))
      ),
    #binning
    `bins` = 
      list(
        binPrep(
          sites= `Site Name`,
          ages= `D Uncal`,
          h=100)
      ),
    #calibrating dates
    `DCal`= 
      list(
        calibrate(
          x=`D Uncal`,
          errors=Error,
          calCurves='intcal20')
      )
  )|>
  #Calculating SPDs and Adding SPDs Column
  dplyr::rowwise() |>
  dplyr::mutate(
    SPDs = list(
      spd(`DCal`,
          bins = `bins`,
          timeRange=c(6000,500)
      )
    )
  )

#plotting SPDS
purrr::map2(
  .x = analysis_by_site$SPDs,
  .y = analysis_by_site$`Site Name`,
  .f = SPDSplot
)


#Finding the Max PrDense for plotting pop. movement-By Site
Densitiesbysite <-
  analysis_by_site |>
  dplyr::rowwise() |>
  dplyr::mutate(
    prdensdate =
      list(
        SPDs$grid |>
          tibble::as_tibble() |>
          filter(PrDens == max(PrDens))),
    #Finding interquartile range
    prdensquat = 
      list(
        quantile(
          sample(
            SPDs$grid$calBP,
            size = 99999,
            replace =TRUE,
            prob = SPDs$grid$PrDens))))

#printing density quartiles
Densitiesbysite$prdensquat



#Finding the Max PrDense for plotting pop. movement-by region
Densitiesbyregion <-
  analysis_by_Region |>
  dplyr::rowwise() |>
  dplyr::mutate(
    prdensdate =
      list(
        SPDs$grid |>
          tibble::as_tibble() |>
          filter(PrDens == max(PrDens))),
    #Finding interquartile range
    prdensquat = 
      list(
        quantile(
          sample(
            SPDs$grid$calBP,
            size = 99999,
            replace =TRUE,
            prob = SPDs$grid$PrDens))))

#Printing Quartiles
Densitiesbyregion$prdensquat 





#Supplemental plots and tests

#Running analyses by culture- binning parameter 50
Suptestsbycult50 <- 
  #Reading Data Into R
  read_excel("data/MainData.xlsx")|> 
  filter(!is.na(`D Uncal`),
         !is.na(`Error`)) |>
  #grouping by culture
  dplyr::group_by(
    `I or S`) |>
  #Calibrating Dates, Binning, and Adding DCal/bins Column
  dplyr::summarize( 
    `bins` = list(
      binPrep(
        sites= `Site`,
        ages= `D Uncal`,
        h=50)
    ),
    `DCal`= 
      list(
        calibrate(
          x=`D Uncal`,
          errors=Error,
          calCurves='intcal20')
      ),
    `errors` =
      list(Error))|>
  #Calculating SPDs and Adding SPDs Column
  dplyr::rowwise() |>
  dplyr::mutate(
    SPDs = list(
      spd(`DCal`,
          bins = `bins`,
          timeRange=c(6000,500)
      )
    ),
    #testing SPDs against null hyp of exponential growth
    exptest50 = 
      list(
        modelTest(x = DCal,
                  errors=errors,
                  bins=bins,
                  nsim=1000,
                  timeRange=c(6000,500),
                  model="exponential",
                  runm=NA)
      )
  )

Suptestsbycult200 <- 
  #Reading Data Into R
  read_excel("data/MainData.xlsx")|> 
  filter(!is.na(`D Uncal`),
         !is.na(`Error`)) |>
  #grouping by culture
  dplyr::group_by(
    `I or S`) |>
  #Calibrating Dates, Binning, and Adding DCal/bins Column
  dplyr::summarize( 
    `bins` = list(
      binPrep(
        sites= `Site`,
        ages= `D Uncal`,
        h=200)
    ),
    `DCal`= 
      list(
        calibrate(
          x=`D Uncal`,
          errors=Error,
          calCurves='intcal20')
      ),
    `errors` =
      list(Error))|>
  #Calculating SPDs and Adding SPDs Column
  dplyr::rowwise() |>
  dplyr::mutate(
    SPDs = list(
      spd(`DCal`,
          bins = `bins`,
          timeRange=c(6000,500)
      )
    ),
    #testing SPDs against null hyp of exponential growth
    exptest200 = 
      list(
        modelTest(x = DCal,
                  errors=errors,
                  bins=bins,
                  nsim=1000,
                  timeRange=c(6000,500),
                  model="exponential",
                  runm=NA)
      )
  )

#Running analyses by culture- binning parameter 300
Suptestsbycult300 <- 
  #Reading Data Into R
  read_excel("data/MainData.xlsx")|> 
  filter(!is.na(`D Uncal`),
         !is.na(`Error`)) |>
  #grouping by culture
  dplyr::group_by(
    `I or S`) |>
  #Calibrating Dates, Binning, and Adding DCal/bins Column
  dplyr::summarize( 
    `bins` = list(
      binPrep(
        sites= `Site`,
        ages= `D Uncal`,
        h=300)
    ),
    `DCal`= 
      list(
        calibrate(
          x=`D Uncal`,
          errors=Error,
          calCurves='intcal20')
      ),
    `errors` =
      list(Error))|>
  #Calculating SPDs and Adding SPDs Column
  dplyr::rowwise() |>
  dplyr::mutate(
    SPDs = list(
      spd(`DCal`,
          bins = `bins`,
          timeRange=c(6000,500)
      )
    ),
    #testing SPDs against null hyp of exponential growth
    exptest300 = 
      list(
        modelTest(x = DCal,
                  errors=errors,
                  bins=bins,
                  nsim=1000,
                  timeRange=c(6000,500),
                  model="exponential",
                  runm=NA)
      )
  )

#Running analyses by culture- binning parameter 500
Suptestsbycult500 <- 
  #Reading Data Into R
  read_excel("data/MainData.xlsx")|> 
  filter(!is.na(`D Uncal`),
         !is.na(`Error`)) |>
  #grouping by culture
  dplyr::group_by(
    `I or S`) |>
  #Calibrating Dates, Binning, and Adding DCal/bins Column
  dplyr::summarize( 
    `bins` = list(
      binPrep(
        sites= `Site`,
        ages= `D Uncal`,
        h=500)
    ),
    `DCal`= 
      list(
        calibrate(
          x=`D Uncal`,
          errors=Error,
          calCurves='intcal20')
      ),
    `errors` =
      list(Error))|>
  #Calculating SPDs Adding SPDs Column
  dplyr::rowwise() |>
  dplyr::mutate(
    SPDs = list(
      spd(`DCal`,
          bins = `bins`,
          timeRange=c(6000,500)
      )
    ),
    #testing SPDs against null hyp of exponential growth
    exptest500 = 
      list(
        modelTest(x = DCal,
                  errors=errors,
                  bins=bins,
                  nsim=1000,
                  timeRange=c(6000,500),
                  model="exponential",
                  runm=NA)
      )
  )



#plotting exptest - binning value 50     
purrr::map2(
  .x = Suptestsbycult50$exptest50,
  .y = suptestsbycult$`I or S`,
  .f = Expontentialplot
)
#plotting exptest - binning value 200     
purrr::map2(
  .x = Suptestsbycult200$exptest200,
  .y = suptestsbycult$`I or S`,
  .f = Expontentialplot
)
#plotting exptest - binning value 300     
purrr::map2(
  .x = Suptestsbycult300$exptest300,
  .y = suptestsbycult$`I or S`,
  .f = Expontentialplot
)
#plotting exptest - binning value 500     
purrr::map2(
  .x = Suptestsbycult500$exptest500,
  .y = suptestsbycult$`I or S`,
  .f = Expontentialplot
)


