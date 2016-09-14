## makeItRain
An R script for basic modelling of hydroregime

### Usage
This script doesn't take input variables - for each location/pool the script should be modified accordingly.

### input data
For each rock location the following data are required:
- `.csv` file containing historical rainfall data from the nearest Bureau of Meteorology station to the pool being modelled
- `.csv` file containing mean pan evaporation data for each month for the coordinates of the pool
For each pool (on each rock) maximum pool depth is also required. Example files are available on request.

### Requirements
The following packages are required:
- reshape2
- ggplot2
- plyr
- tweedie
- statmod
- gridExtra

### Test system
Tested and running on 2016.09.14 with the following system:
```
Linux 3.19.0-26-generic x86_64
Distributor ID:	Ubuntu
Description:	Ubuntu 15.04
Release:	15.04
Codename:	vivid
R version 3.3.1 (2016-06-21) -- "Bug in Your Hair"
R Studio Version 0.99.473
```
