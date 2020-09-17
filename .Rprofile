# define user specific data and output directories
# necessary since these are not under version control
user <- Sys.info()[["user"]]
if (user == "mes335") {
  DATA_DIR <- file.path(getwd(), "data")
  OUTPUT_DIR <- file.path(getwd(), "output")
} else if (user == "") {
  DATA_DIR <- ""
  OUTPUT_DIR <- ""
} else {
  DATA_DIR <- file.path(getwd(), "data")
  OUTPUT_DIR <- file.path(getwd(), "output")
}
