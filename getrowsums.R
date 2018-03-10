#find out how many sequences per sample to figure out rarefaction
getrowsums <- function(transdataframe){
  EMF_summary <- rowSums(transdataframe)
  EMF_summarysorted <- EMF_summary [order(EMF_summary, decreasing = TRUE)]
  boxplot(EMF_summarysorted)
  print(EMF_summarysorted)
  quantile10 <- quantile(EMF_summarysorted, 0.10)
  print(quantile10)
  quantile15 <- quantile(EMF_summarysorted, 0.15)
  print(quantile15)
}