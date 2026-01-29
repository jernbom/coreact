# R/globals.R

# 1. Import specific functions needed globally
#' @importFrom dplyr %>% everything
#' @importFrom utils object.size
NULL

# 2. Declare global variables to silence R CMD check notes about NSE
utils::globalVariables(c(
  "idx_x", "idx_y", "size_x", "size_y", "feature_x", "feature_y",
  "idx_x_tmp", "size_x_tmp", "feat_x_tmp", "idx_y_tmp", "size_y_tmp", "feat_y_tmp",
  "p_adj", "intersection", "union", "jaccard", "overlap", "p_val" # Added these just in case
))
