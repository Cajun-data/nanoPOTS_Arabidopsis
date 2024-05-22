### Example usage:
###  test1 <- fast_cor(
###   data,
###   use = "pairwise.complete.obs",
###   method = "spearman"
###     )
### "data" is class "matrix" with samples as rows and observations as columns


fast_cor <- function(data, use = "pairwise.complete.obs", method = "spearman") 
{
  {
  n = NULL
  n <- psych::pairwiseCount(data) - 2 ## Create matrix of degrees freedom
  }
  if (method == "spearman") {
    data <- base::apply(X = data, MARGIN = 2, data.table::frankv)
  }
  r = NULL
  p1 = NULL
  p2 = NULL
  r <- coop::pcor(x = data, use = use)
  {
    
    t <- sqrt(n)*r/sqrt(1 - r^2) ## T-test
    p1 <- stats::pt(t, n)  ## One side p-value
    p2 <- stats::pt(t, n, lower.tail = FALSE) ## Unfortunately run twice...
  }
  flt.Corr.Matrix <- function(cormat, pmat1 = NULL, pmat2 = NULL,
                              df = NULL) {
    ut <- r
    ut[,] <- TRUE
    ut <- ut == 1
    flt_data <- data.frame(Gene1 = base::rownames(cormat)[base::row(cormat)[ut]], 
                           Gene2 = base::rownames(cormat)[base::col(cormat)[ut]], 
                           cor = cormat[ut])
    if (!is.null(pmat1)) 
      flt_data$p1 <- pmat1[ut]
    if (!is.null(pmat2)) 
      flt_data$p2 <- pmat2[ut]
    if (!is.null(df)) 
      flt_data$df <- df[ut]
    return(flt_data)
  }
  {
    result <- flt.Corr.Matrix(cormat = r,
                              df = n,
                              pmat1 = p1,
                              pmat2 = p2)
    result <-  base::transform(result, p = base::pmin(p1, p2)*2) %>%
      mutate(p = case_when(Gene1 == Gene2 ~ NA,
                           TRUE ~ p)) %>%
      mutate(cor = case_when(Gene1 == Gene2 ~ 1,
                             TRUE ~ cor))
    result$FDR <- stats::p.adjust(result$p, method = "BH")
    result <- base::subset(result, select = -c(p1,p2)
                           )
    result <- list(result, r)
  }
  return(result)
}

