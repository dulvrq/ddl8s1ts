
# package generation

library(devtools)
library(usethis)

use_package("raster")
use_package("doParallel")
use_package("randomForest")
use_package("caret")
use_package("dplyr")
use_package("gdalUtils")
use_package("foreach")
use_package("ggplot2")
use_package("data.table")

# library(raster)
# library(doParallel)
# library(randomForest)
# library(caret)
# library(dplyr)
# library(gdalUtils)

usethis::use_vignette("ddl8s1ts")
usethis::use_mit_license("Katsuto Shimizu")

usethis::use_readme_rmd()

devtools::document()
devtools::load_all()
devtools::check()
#devtools::build_vignettes()
