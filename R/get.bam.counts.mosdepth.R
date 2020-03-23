#' @importFrom yaml  read_yaml
#' @export
get.bam.counts.mosdepth <- function(task_id,input){
# pkgPath <- find.package("ExomeDepthMod")
  input.yaml <- yaml::read_yaml(input)
  manifest <- read.table(paste0(input.yaml$output.directory,"/results/",input.yaml$manifest), header = T, sep ="\t", stringsAsFactors = F)
  names(manifest) <- c("bam","sampleID","sex")
  
  exon_path <- paste0(input.yaml$output.directory,"/results/bed_file.bed")

  # options(
  #  clustermq.scheduler = paste0(input.yaml$scheduler),
  #  clustermq.template = paste0(find.package("EDM"),"/template/",input.yaml$scheduler,"_template"),
  #  clustermq.defaults = list(conda="edm_env")
  # )
  # mosdepth_exe = "/cm/shared/apps_chop/mosdepth/0.2.1/bin/mosdepth";
  # data(exons.hg19,package = "ExomeDepth")
  # data(exons.hg19.X,package = "ExomeDepth")
  # exons <- rbind(exons.hg19,exons.hg19.X)
  # exons <- paste0(pkgPath,"/data/exons.hg19.bed")
  # input.yaml <- yaml::read_yaml("input.yaml")
  # manifest <- read.table(as.character(input.yaml$manifest),header=T,sep="\t")
  sampleID <- task_id;
  sample_name <- gsub(" ","_", manifest[sampleID,2])
  cmd = paste0("mosdepth --by ",exon_path," --mapq 20 --no-per-base --fast-mode ",input.yaml$output.directory,"/coverage/",sample_name," ",as.character(manifest[sampleID,1]));
  res <- system(cmd,intern=T)
  return(NULL)
}
