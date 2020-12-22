#' test_dataset.
#'
#' Example dataset containing gene expression values across 23 cell types
#' obtained using FindAllMarkers function from SEURAT v3.
#'
#' @format A data frame with 12814 rows and 7 variables:
#' \describe{
#'   \item{p_val}{p value}
#'   \item{avg_logFC}{log Fold Change}
#'   \item{pct.1}{}
#'   \item{pct.2}{}
#'   \item{p_val_adj}{Adjusted p value}
#'   \item{cluster}{Cluster ID}
#'   \item{gene}{Gene name}
#' }
#' @source \url{}
test_dataset <- read.csv("~/Desktop/Test dataset.csv", row.names = 1)
"test_dataset"
usethis::use_data(test_dataset)
