# <img src="mpaeu_logo.png" align="right" width="240" /> MPA Europe - Climate velocity map for European seas under recent conditions

Here you find the codes used to produce the Deliverable 2.8 of MPA Europe, "Climate velocity map for European seas under recent conditions".

The underlying data is available on Bio-ORACLE v3.0 ([publication here](https://onlinelibrary.wiley.com/doi/abs/10.1111/geb.13813) and [package here](https://github.com/bio-oracle/biooracler)).

The study area shapefile is available on the repository [iobis/mpaeu_studyarea](https://github.com/iobis/mpaeu_studyarea)

## VoCC - conversion to `terra` functions

This analysis is mainly done through the VoCC package, as described in García Molinos, J., Schoeman, D. S., Brown, C. J. and Burrows, M. T. (2019), VoCC: An R package for calculating the velocity of climate change and related climatic metrics. Methods Ecol Evol. doi:10.1111/2041-210X.13295

The original package was written using functions from the `raster` package. For better performance, we converted all package to work mainly with the [`terra` package](https://rspatial.github.io/terra/index.html). A version of this converted package is available here in the folder `VoCC-terra`. You can install it using `devtools::install_github("iobis/mpaeu_climatevelocity/VoCC-terra")`. Users should also refer to the original package in the repository [JorGarMol/VoCC](https://github.com/JorGarMol/VoCC)