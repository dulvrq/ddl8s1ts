
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ddl8s1ts

<!-- badges: start -->

<!-- badges: end -->

`ddl8s1ts` is an R package for detecting and mapping forest disturbance
using Landsat 8 and Sentinel-1 time seires. Harmonic regression models
are applied to both Landsat 8 and Sentinel-1 time series to predict
disturbance probability. Disturbance is detected based on the time
series of disturbance probabilities predicted by Random Forests (RF).

  - The functions are tested only on Windows OS.
  - The process takes very long time to map large area and uses much
    memory and CPUs in parallel computing.

## Installation

You can install the latest version from GitHub using `devtools`:

``` r
devtools::install_github("dulvrq/ddl8s1ts")
library(ddl8s1ts)
```

## Example

Currently, `mapDisturbanceL8S1()` is an only usable function in this
package. All other functions are called by this function. Several data
are required for running `mapDisturbanceL8S1()` to detect and map
disturbance. The followings are key variables that should be prepared
beforehand.

**Landsat 8 data**

`ls_l8`: requires a file list of 6 bands images of Landsat 8 (B2, B3,
B4, B5, B6, and B7). All images must have the same extent (origin,
resolution, and crs). These should be surface reflectance products, may
be from Google Earth Engine (GEE) or other preprocessing chains. NoData
values should be `NA`, not 0 (sometimes occurred in GEE output). An
exmaple of this product is the following.

``` r
head(basename(ls_l8))
#> [1] "LC08_L1TP_132047_20140106_20170427_01_T1.tif"
#> [2] "LC08_L1TP_132047_20140207_20170426_01_T1.tif"
#> [3] "LC08_L1TP_132047_20140223_20170425_01_T1.tif"
#> [4] "LC08_L1TP_132047_20140311_20170425_01_T1.tif"
#> [5] "LC08_L1TP_132047_20140327_20170424_01_T1.tif"
#> [6] "LC08_L1TP_132047_20140412_20170423_01_T1.tif"
```

**Sentinel-1 data**

`ls_s1`: requires a file list of 2 bands images of Sentinel-1 (VV and
VH). All images must have the same extent (origin, resolution, and crs).
These might be from GEE or other preprocessing chains. An exmaple of
this product is the following.

``` r
head(basename(ls_s1))
#> [1] "S1A_IW_GRDH_1SDV_20141012T233136_20141012T233201_002803_003287_0803.tif"
#> [2] "S1A_IW_GRDH_1SDV_20141012T233136_20141012T233201_002803_003287_6FDB.tif"
#> [3] "S1A_IW_GRDH_1SDV_20141012T233201_20141012T233226_002803_003287_1E3D.tif"
#> [4] "S1A_IW_GRDH_1SDV_20141012T233201_20141012T233226_002803_003287_F55D.tif"
#> [5] "S1A_IW_GRDH_1SDV_20141012T233226_20141012T233251_002803_003287_3DB5.tif"
#> [6] "S1A_IW_GRDH_1SDV_20141012T233226_20141012T233251_002803_003287_9E79.tif"
```

**Julian day of Landsat 8 and Sentinel-1 data**

`l8_doys`: is a list of julian day of Landsat 8. This should be the same
order and length as `ls_l8`.  
`s1_doys`: is a list of julian day of Sentinel-1. This should be the
same order and length as `ls_s1`.

``` r
# calculate julian day from filenames
n <- 18
l8_dates <- paste(substring(basename(ls_l8), n,   n+3),
                  substring(basename(ls_l8), n+4, n+5),
                  substring(basename(ls_l8), n+6, n+7), sep = "-")
l8_doys <- as.numeric(substring(basename(ls_l8), n, n+3)) + as.numeric(strftime(l8_dates, format = "%j")) / 365
s1_dates <- paste(substring(basename(ls_s1), n,   n+3),
                  substring(basename(ls_s1), n+4, n+5),
                  substring(basename(ls_s1), n+6, n+7), sep = "-")
s1_doys <- as.numeric(substring(basename(ls_s1), n, n+3)) + as.numeric(strftime(s1_dates, format = "%j")) / 365

head(l8_doys)
#> [1] 2014.016 2014.104 2014.148 2014.192 2014.236 2014.279
head(s1_doys)
#> [1] 2014.781 2014.781 2014.781 2014.781 2014.781 2014.781
```

**Refernce data for disturbance detection**

`dt_ref`: is a dataframe of reference data for RF modeling. This should
contain column `x`, `y`, and `date`. `x` and `y` should indicate
locations of reference samples in the same crs as Landsat 8 or
Sentinel-1. `date` should indicate the timing of disturbance. If there
is no disturbance at the sample location, use `NA`.

``` r
# this is an example.
head(dt_ref)
#>        x       y       date
#> 1 783945 2179395 2016-01-31
#> 2 817905 2154075 2016-09-17
#> 3 838485 2147385 2016-03-23
#> 4 829215 2135325 2016-03-07
#> 5 778575 2121855 2016-10-07
#> 6 848685 2117775 2016-11-01
```

**Others**

`ls_dem`: a filename of DEM with 2 bands (elevation and slope).
Optional. If DEM will be not used for model, use `NULL`.  
`dir_save`: a directory for save the results. NULL is acceptable.  
`VI` : a list of spectral index used for Landsat 8.  
`rf_model`: use `NULL` to build RF model. To skip RF modeling, provide a
list of two RF models.  
`startDOY`: a start date of disturbance detection.  
`endDOY`: a end date of disturbance detection.  
`mmu`: a minimum mapping unit (pixel) for final mapping. Use `NULL` to
avoid `mmu`.  
`only_rf`: If `TRUE`, only build RF models. If `FALSE`, build RF models
and map disturbance detection.

``` r
ls_dem <- NULL # do not use DEM in RF models
dir_save <- getwd()
startDOY <- 2016 # detect disturbance from 2016/01/01
endDOY   <- 2018 # detect disturbance until 2017/12/31
mmu      <- 4
VI <- c("NBR", "TCA8", "TCB8", "TCG8","TCW8") # NBR, TCA, TCB, TCG, and TCW
```

#### example 1. Full implementation of building RF models and mapping disturbance

`mapDisturbanceL8S1()` has 2 processing parts: RF model building and
disturbance mapping. A simple implementation is to process both as
follows. Please note that this function uses parallel processing and
consume much memory and CPUs. Furthermore, this process takes very long
time.

``` r
# use only_rf = F & rf_model = NULL to run entire process.
mapDisturbanceL8S1(ls_l8, ls_s1, l8_doys, s1_doys, dt_ref, ls_dem, dir_save, VI,
                             rf_model = NULL, startDOY, endDOY, mmu, only_rf = F)
```

#### example 2. Full implementation only with Landsat 8

If you would like to use either of Landsat 8 or Sentinel-1, use NULL to
either of file list.

``` r
# use ls_s1 = NULL.
mapDisturbanceL8S1(ls_l8, NULL, l8_doys, NULL, dt_ref, ls_dem, dir_save, VI,
                             rf_model = NULL, startDOY, endDOY, mmu, only_rf = F)
```

#### example 3. Implementation of RF modeling

If you would like to run only RF model building, such as in case that
samples are from large extent but would like to map only small extent,
then use `only_rf = TRUE`. This only generate RF models in a subfolder.

``` r
# use only_rf = T.
mapDisturbanceL8S1(ls_l8, ls_s1, l8_doys, s1_doys, dt_ref, ls_dem, dir_save, VI,
                             rf_model = NULL, startDOY, endDOY, mmu, only_rf = T)
```

#### example 4. Implementation of disturbance mapping

If you would like to run only disturbance mapping, such as in case of
mapping only small extent, then provide a list of two RF models (Landsat
and Sentinel-1).

``` r
# provide rf_model, a list of two RF models.
mapDisturbanceL8S1(ls_l8, ls_s1, l8_doys, s1_doys, dt_ref, ls_dem, dir_save, VI,
                             rf_model = rf_model, startDOY, endDOY, mmu, only_rf = F)
```

#### Other note:

  - It is a good idea to mask forest area before running this
    function.  
  - A sub folder is automatically generated for saving RF results and
    temporally tiles of mapping.  
  - 50% of reference samples are used for RF modeling, and the rest are
    used for evaluation.

## Reference

Shimizu, K. Ota, T. Mizoue, N. (2019) Detecting Forest Changes Using
Dense Landsat 8 and Sentinel-1 Time Series Data in Tropical Seasonal
Forests. *Remote Sensing* 11: 1899.
[doi.org/10.3390/rs11161899](https://doi.org/10.3390/rs11161899)
