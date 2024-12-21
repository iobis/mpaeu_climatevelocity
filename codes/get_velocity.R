############################## MPA Europe project ##############################
################### WP2 - Environmental layers and analyses ####################
# December 2024
# Authors: Michael Burrows, Silas C. Principe, Jorge Assis
# Contact: s.principe@unesco.org
#
########################## Calculate climate velocity ##########################

# Load packages
library(terra)
library(VoCC)
library(furrr)
library(future)

procf <- "processed-layers"
fs::dir_create(procf)


paths <- list.files("rasterLayers", full.names = T)
paths <- paths[grepl("meanSST ssp", paths)]

scenarios <- rep(paste0("ssp", c(119, 245, 370, 585)), 2)
periods <- rep(c("2040-2060", "2080-2100"), each = 4)

get_climate_velocity <- function(scenario, period, paths, comp_red = 10,
                                plot = FALSE, verbose = FALSE) {

    if (verbose) message("Scenario ", scenario, " period ", period, " started.")

    # Skip if done
    if (file.exists(file.path(procf, paste0("trajclass_", scenario, "_", period, ".tif")))) {
        return(invisible(NULL))
    }

    # Load raster layers
    if (verbose) message("Loading base layer")
    meanSST20002020 <- rast("rasterLayers/meanSST 2000-2020.tif")
    if (plot) plot(meanSST20002020)

    europext <- terra::ext(-30, 45, 25, 75)
    meanSST20002020eu <- terra::crop(meanSST20002020, europext)

    sel_path <- paths[grepl(scenario, paths)]
    sel_path <- sel_path[grepl(period, sel_path)]

    meanSSTssp <- rast(sel_path)

    meanSSTsspeu <- terra::crop(meanSSTssp, europext)

    timeint <- 2050 - 2010
    bmeanSSTsspeu <- (meanSSTsspeu - meanSST20002020eu) / timeint
    tmeanSSTsspeu <- 0.5 * (meanSSTsspeu + meanSST20002020eu)


    # Calculate local spatial gradient
    if (verbose) message("Calculating spatial gradient")
    spat_gradient <- spatGrad(tmeanSSTsspeu, th = -Inf, projected = FALSE)  

    if (plot) plot(spat_gradient)

    writeRaster(spat_gradient, file.path(procf, paste0("sgr_", scenario, "_", period, ".tif")), overwrite = T)

    # Calculate velocity 
    if (verbose) message("Calculating velocity")
    vocc_layer <- gVoCC(bmeanSSTsspeu, spat_gradient)

    if (plot) plot(vocc_layer)

    writeRaster(vocc_layer,
                file.path(procf, paste0("vocc_", scenario, "_", period, ".tif")), overwrite = T)

    # Calculate trajectories
    if (verbose) message("Calculating trajectories")
    vel <- vocc_layer[[1]]
    ang <- vocc_layer[[2]]
    mn <- tmeanSSTsspeu
    baseclim <- meanSST20002020eu

    lonlat <- na.omit(data.frame(xyFromCell(baseclim, 1:ncell(baseclim)),
                                vel[],ang[],mn[]))[,1:2] # Omits all NA values

    lonlat1 <- lonlat[seq(1,nrow(lonlat), comp_red),] # reduces computational time

    ssp_eurotraj <- voccTraj(
        lonlat = lonlat1, vel = vel, ang = ang,
        mn = baseclim, tyr = timeint, correct = T
    ) # 1970 + 70 = 2040

    if (verbose) message("Getting trajectories lines")
    ssp_eurotraj_lines <- trajLine(ssp_eurotraj)
    writeVector(vect(ssp_eurotraj_lines), 
                file.path(procf, paste0("lines_", scenario, "_", period, ".shp")), overwrite = T)

    if (verbose) message("Classifying trajectories")
    clas <- trajClas(ssp_eurotraj, vel, ang, mn,
        trajSt = 4, tyr = timeint, nmL = 20, smL = 100,
        Nend = 45, Nst = 15, NFT = 70, DateLine = FALSE
    )

    # Classify raster / build attribute table
    r <- as.factor(clas[[7]])
    rat_r <- levels(r)[[1]]
    rat_r$TrajClas <- c("N-M", "S-M", "IS", "BS", "Srce", "RS", "Cor", "Div", "Con")[sort(unique(clas[[7]][]))]
    levels(r) <- rat_r

    if (plot) {
        my_col = c('gainsboro', 'darkseagreen1', 'coral4', 'firebrick2', 'mediumblue', 'darkorange1',
                'magenta1', 'cadetblue1', 'yellow1')
        # Keep only the categories present in our raster
        my_col <- my_col[sort(unique(clas[[7]][]))]
        rasterVis::levelplot(r, col.regions = my_col, xlab = NULL, ylab = NULL, scales = list(draw=FALSE))
    }

    writeRaster(r, file.path(procf, paste0("trajclass_", scenario, "_", period, ".tif")), overwrite = TRUE)

    rm(clas)
    gc()

    if (verbose) message("Scenario ", scenario, " period ", period, " concluded.")

    return(invisible(NULL))

}

plan(multisession, workers = 2)

result <- furrr::future_map2(
    scenarios,
    periods,
    .f = get_climate_velocity,
    comp_red = 2,
    paths = paths,
    .progress = T
)

plan(sequential)
