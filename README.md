# Prioritizing ecosystem services and biodiversity

R project for simultaneous global prioritization of ecosystem services and biodiversity.

## Directories

Due to size constrains, the input and output data directories aren't kept under version control. By default the input data are in `data/` and the output data are in `outputs/`, both ignored by git. These directories are defined at the top of each script, and if you wish to change them modify the `data_dir` and `output_dir` variables. This project assumes that the working directory is the project directory.

## Data

This prioritization simultaneously prioritizes 12 ecosystem services and a broad suite of terrestrial tetrapods species: 13,144 birds, 5,118 mammals, 4,083 reptiles, and 5,978 amphibians. The ecosystem services are:

1. Coastal protection `coastalprotection-norm-onshore`
2. Timber (commercial) `commercialtimber-forest`
3. Timber (domestic) `domestictimber-forest`
4. Flood mitigation `flood-nathab`
5. Fuelwood `fuelwood-forest`
6. Livestock grazing `grazing-natnotforest`
7. Carbon `vulnerable-c-total`
8. Nature access (rural) `norm-nature-access-lspop-2017-urca-rural-60`
9. Nature access (urban) `norm-nature-access-lspop-2017-urca-urban-60`
10. Nitrogen retention `nitrogenretention-attn-500km`
11. Pollination `pollination-norm-nathab`
12. Sediment retention `sedimentdeposition-attn-500km`

The analysis is conducted at a variety of resolutions: 2km, 3km, 5km, and 10km. All ecosystem services layers were provided at a native resolution of 2km and aggregated to the lower resolutions. Biodiversity layers are only used in the 5km and 10km analyses, and where provided directly at these resolutions with cell values representing the percent (0-100) of each cell that falls within each species range.

All layers are masked to land boundary polygons prior to prioritization. Some or all of the range of some species may be lost when masking to the terrestrial boundary. Pelagic species are of particular note here, but also some species only occur on tiny islands that are not captured by the terrestrial boundary polygon we're using. This has a couple notable consequences:

1. When using proportional targets in prioritization, the target will be a proportion of terrestrial AOH, not the total AOH.
2. Of the 29,360 biodiversity features, 1,027 have no AOH falling with the terrestrial boundary and are therefore removed from the analysis.

Raw data are all stored as GeoTIFFs within the `data/tifs/` directory in sub-directories according to data type, either `es/` for ecosystem services or class level directory for biodiversity features: `mammals/`, `birds/`, `reptiles/`, and `amphibians/`.

## Workflow

1. `00_aggregate.R`: aggregate the ecosystem services from the native 2km resolution to 3km, 5km, and 10km resolutions.
1. `01_planning-units.R`: define grids of planning units at 10 km, these are cells within 10 km of land that have at positive values for at least one feature.
2. `02_feature-prep.R`: generate representation tables for every feature, i.e. data frames that specify the amount of each feature within each planning unit. Also, creates a master table listing all the (non-zero) features to be included in the prioritization.
3. `03_rij-matrices.R`: stack the feature specific representation tables and create a representation matrix in `sparseMatrix` format.
4. `04_prioritize.R`: perform a suite of prioritizations for 9 scenarios, consisting of all permutations of:
  - Ecosystem services targets: proportional targets of 0.3, 0.5, or 0.9
  - Budget: the maximum proportion of planning units selected in the prioritization: 0.3, 0.5, or 1
5. `04-1_prioritize_es-only.R`: similar to 4, but ignore the biodiversity features and only prioritize for the ecosystem services features.
  
Prioritization results are stored in tabular and raster format in the `OUTPUT_DIR` for each scenario.