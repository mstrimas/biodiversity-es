# Prioritizing ecosystem services and biodiversity

R project for simulataneous global prioritization of ecosystem services and biodiversity.

## Setup

To aid portability of code, avoid the use of `setwd()`, and account for the fact that the input and output data directories are too large to be under version control on GitHub:

1. **Always ensure the working directory is the project directory**. Either use RStudio and open the RStudio project file or start your R session from the project directory. Using `setwd()` will break the workflow.
2. Store the input and output files in some sensible location, which doesn't necessarily have to be within the project directory. For example, if you intend to sync the files with DropBox or similar just not the location of these directories on your hard drive. I (Matt) have these as the `data/` and `output/` sub-directories of the project directory and have set them to be ignored by git using the `.gitignore` file.
3. Specify user specific `DATA_DIR` and `OUTPUT_DIR` variables specifying the directories from 2 in the `.Rprofile` file. Look at the output of `Sys.info()[["user"]]`, then modify the if-else block using the existing code as a template to set these variables for your system. These variables default to sub-directories of the project directory, so if you go that route you can skip step 3 entirely. If this step fails, it's likely that you are not starting R from the correct directory as per 1 above.

## Data

This prioritization simulaneously prioritizes 10 ecosystem services and a broad suite of terrestrial tetrapods species: 13,144 birds, 5,118 mammals, 4,083 reptiles, and 5,978 amphibians. The ecosystem services are:

1. Timber (commercial) `commercialtimber-forest`
2. Timber (domestic) `domestictimber-forest`
3. Flood mitigation `flood-nathab`
4. Fuelwood `fuelwood-forest`
5. Livestock grazing `grazing-natnotforest`
6. Carbon `irreplaceable-carbon`
7. Nature access `natureaccess10-nathab`
8. Nitrogen retention `nitrogenretention-nathab`
9. Pollination `pollination-nathab`
10. Sediment retention `sedimentdeposition-nathab`

All biodiversity layers are derived from IUCN range polygons, rasterized at 1 km resolution, then aggregated to higher resolutions to represent the percent (0-100) of each cell that falls within each species range. Bird species are typically broken into multiple species-season combinations: breeding, non-breeding, and migration. The impact is that some bird species have each of their seasons prioritized independently.

All layers are masked to land boundary polygons prior to prioritization. Some or all of the range of some species may be lost when masking to the terrestrial boundary. Pelagic species are of particular note here, but also some species only occur on tiny islands that are not captured by the terrestrial boudnary polygon we're using. This has a couple notable consequences:

1. When using proportional targets in prioritization, the target will be a proportion of terrestrial AOH, not the total AOH.
2. Of the 29,360 biodiversity features, 1,027 have no AOH falling with the terrestrial boundary and are therefore removed from the analysis.

Raw data are all stored as GeoTIFFs within the `data/tifs/` directory in subdirectories according to data type, either `es` for ecosystem services or class level directory for biodiversity features.

## Workflow

1. `01_planning-units.R`: define grids of planning units for a variety of spatial resolutions, e.g. 3, 5, and 10 km.
2. `02_feature-prep.R`: given a single specified planning unit resolution, generate representation tables for every feature, i.e. data frames that specify the amount of each feature within each planning unit. Also, creates a master table listing all the (non-zero) features to be included in the prioritization.
3. `03_rij-matrices.R`: given a single specified planning unit resolution, stack the feature specific representation tables and create a representaiton matrix in `sparseMatrix` format.
4. `04_prioritize.R`: given a single specified planning unit resolution, perform a suite of prioritizations for 9 scenarios, consisting of all permutations of:
  - Ecosystem services targets: proportional targets of 0.3, 0.5, or 0.9
  - Budget: the maximum proportion of planning units selected in the prioritization: 0.3, 0.5, or 1
  
Prioritization results are stored in tabular and raster format in the `OUTPUT_DIR` for each scenario.