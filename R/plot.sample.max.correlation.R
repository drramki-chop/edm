#'
#'
#' This function plots the maximum pairwise correlation for each sample in the cohort.
#'
#'
#'
#'
#' @export
#' @importFrom yaml  read_yaml
#' @importFrom ggplot2 ggplot

plot.sample.max.correlation <- function(input){
  input.yaml <- yaml::read_yaml(input) 
  cnv.cor.files <- list.files(paste0(input.yaml$output.directory,"/results/individual.edm.calls/"),pattern = ".edm.stat.txt", full.names=T)
  cnv.cor.data <- do.call( rbind,lapply(1:length(cnv.cor.files), function(file) read.table(cnv.cor.files[file],header=T,stringsAsFactors=F)))
  cnv.cor.data %>% filter(test == "correlation") -> cnv.cor.data
  cnv.cor.data$sample_name <- gsub(".edm.stat.txt","",list.files(paste0(input.yaml$output.directory,"/results/individual.edm.calls/"),pattern = ".edm.stat.txt"))
  write.table(cnv.cor.data,paste0(input.yaml$output.directory,"/results/cohort.edm.correlation.txt"), row.names = F,quote=F)
  p <- ggplot(cnv.cor.data,aes(x=sample_name, y = as.numeric(value))) + 
    geom_point(size=2,alpha=0.5) + 
    theme_classic() + 
    ylim(values=c(0,1)) + 
    geom_hline(yintercept = 0.97, color="red") + 
    xlab("Sample Index") + 
    ylab("Maximum pairwise correlation") +
    ggtitle("Pairwise correlation between test and reference samples") + theme(axis.text.x = element_blank())
  pdf(paste0(input.yaml$output.directory,"/results/pairwise.cor.plot.pdf"))
  print(p)
  dev.off()
  return(NULL)
}
