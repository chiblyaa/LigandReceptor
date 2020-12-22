#' LRdatabase.
#'
#' A dataset containing 2557 ligand-receptor pairs
#'
#' @format A data frame with 2557 rows and 5 variables:
#' \describe{
#'   \item{Pair}{Gene symbols of ligand-receptor pair}
#'   \item{Ligand}{Gene symbol of Ligand}
#'   \item{Ligand.name}{Gene name of ligand}
#'   \item{Receptor}{Gene symbol of receptor}
#'   \item{Receptor.name}{Gene name of receptor}
#' }
#' @source \url{https://static-content.springer.com/esm/art%3A10.1038%2Fncomms8866/MediaObjects/41467_2015_BFncomms8866_MOESM611_ESM.xlsx}
LRdatabase <- readxl::read_excel(path = "~/Desktop/testfile.xlsx")
names(LRdatabase) <- c("Pair", "Ligand", "Ligand.name", "Receptor", "Receptor.name")
LRdatabase$Pair <- tolower(LRdatabase$Pair)
LRdatabase$Pair <- Hmisc::capitalize(LRdatabase$Pair)
LRdatabase$Ligand <- tolower(LRdatabase$Ligand)
LRdatabase$Ligand <- Hmisc::capitalize(LRdatabase$Ligand)
LRdatabase$Receptor <- tolower(LRdatabase$Receptor)
LRdatabase$Receptor <- Hmisc::capitalize(LRdatabase$Receptor)

LRdatabase <- data.frame(lapply(LRdatabase, function(x) {
  gsub("Ntf4", "Ntf5", x)
}))
"LRdatabase"
usethis::use_data(LRdatabase, overwrite = T)
