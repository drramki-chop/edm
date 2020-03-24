#' @importFrom yaml  read_yaml
#' @export
check.exons.def <- function(input){
  ## read yaml file  
  input.yaml <- yaml::read_yaml(input)
  
  ## bed-file
    if(is.null(input.yaml[["bed.file"]]) | length(input.yaml[["bed.file"]])==0 ){
      data(exons.hg19,package = "ExomeDepth")
      data(exons.hg19.X,package = "ExomeDepth")
      exons <- rbind(exons.hg19,exons.hg19.X)
      exons_chr <- exons
      exons_chr[,1] <- paste0("chr",exons_chr[,1])
      
      write.table(exons[,1:3],paste0(input.yaml$output.directory,"/results/bed_file.bed"),row.names =F,sep="\t",quote=F,col.names=F)
      exon_path <- paste0(input.yaml$output.directory,"/results/bed_file.bed")
      write.table(exons_chr[,1:3],paste0(input.yaml$output.directory,"/results/chr_bed_file.bed"),row.names =F,sep="\t",quote=F,col.names=F)
    } else {
      if(tools::file_ext(input.yaml$bed.file)=="rds"){
        exons = readRDS(input.yaml$bed.file)
        if (dim(exons)[2] == 3){
          names(exons) <- c("chromosome","start","end")
        } else {
          names(exons) <- c("chromosome","start","end","name")
          exons$chromosome <- gsub("chr","",exons$chromosome)
          
          write.table(exons[,1:3],paste0(input.yaml$output.directory,"/results/bed_file.bed"),row.names =F,sep="\t",quote=F,col.names=F)
          exon_path <- paste0(input.yaml$output.directory,"/results/bed_file.bed")
        exons_chr <- exons
         exons_chr[,1] <- paste0("chr",exons_chr[,1])
write.table(exons_chr[,1:3],paste0(input.yaml$output.directory,"/results/chr_bed_file.bed"),row.names =F,sep="\t",quote=F,col.names=F)
        }
        
      } else{
        exons <- read.table(input.yaml$bed.file,header = F)
        if (dim(exons)[2] == 3){
          names(exons) <- c("chromosome","start","end")
        } else {
          names(exons) <- c("chromosome","start","end","name")
        }
        exons$chromosome <- gsub("chr","",exons$chromosome)
        
        write.table(exons[,1:3],paste0(input.yaml$output.directory,"/results/bed_file.bed"),row.names =F,sep="\t",quote=F,col.names=F)
        exon_path <- paste0(input.yaml$output.directory,"/results/bed_file.bed")
        exons_chr <- exons
        exons_chr[,1] <- paste0("chr",exons_chr[,1])
write.table(exons_chr[,1:3],paste0(input.yaml$output.directory,"/results/chr_bed_file.bed"),row.names =F,sep="\t",quote=F,col.names=F)
      }
    }
} 
