#' @export
select.exons <- function(cohort.counts, bin.length=NULL,n.bins.reduced =10000){
  total.counts <- rowSums(cohort.counts);
  my.quantiles <- quantile(total.counts [ which(total.counts > 30) ], prob = c(0.1, 0.9), na.rm = TRUE)
  
  selected <- which(total.counts > 30 &
                      bin.length >= as.numeric(quantile(bin.length, prob = 0.05, na.rm = TRUE)) & #I remove very small exons here, because they cause instability
                      bin.length <= as.numeric(quantile(bin.length, prob = 0.95, na.rm = TRUE)) &  # no large exons
                      total.counts < my.quantiles[2])
  if ( (n.bins.reduced < length(selected)) && (n.bins.reduced > 0) ) selected <- selected[ seq(1, length(selected), length(selected) / n.bins.reduced) ]
  return(selected)
}
