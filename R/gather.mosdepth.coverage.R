#' @importFrom yaml read_yaml
#' @importFrom purrr map_dfc
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @export
gather.mosdepth.coverage <- function(task_id,input){
  input.yaml <- yaml::read_yaml(input)
  
 # manifest <- read.table(paste0(input.yaml$output.directory,"/results/",input.yaml$manifest), header = T, sep ="\t", stringsAsFactors = F)
  manifest <- read.table(paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,"_manifest.txt"),header=T, sep = "\t",stringsAsFactors=F)
  names(manifest) <- c("bam","sampleID","sex")
  
  # options(
  #   clustermq.scheduler = paste0(input.yaml$scheduler),
  #  clustermq.template = paste0(find.package("EDM"),"/template/",input.yaml$scheduler,"_template"),
  #  clustermq.defaults = list(conda="edm_env")
  # )
  
  coverageFolder <- paste0(input.yaml$output.directory,"/coverage/");
  cohort.name <- input.yaml$cohort.name;
  
  countFiles <- list.files(path=coverageFolder,pattern = "*.regions.bed.gz$",full.names=T)
  snames <- list.files(path=coverageFolder,pattern = "*.regions.bed.gz$")
  snames <- gsub(".regions.bed.gz","",snames)
  
  ### check file sizes
  fileSizes <- unlist(lapply(countFiles,file.size)) 
  
  ## remove files with no coverage
    failedSamples <- countFiles[which(fileSizes ==0)]
    if(length(failedSamples) > 0 ){
    countFiles <- countFiles[-which(fileSizes ==0)]
    snames <- snames[-which(fileSizes ==0)]
    write.table(failedSamples,paste0(input.yaml$output.directory,"/results/",cohort.name,".samples.failed.coverage.txt"),quote=F,row.names=F)
    message(paste0("Wrote failed samples to the file ",cohort.name,".samples.failed.coverage.txt. Proceeding with the rest of the samples"))
    manifest <- manifest %>% filter( sampleID %in% snames)
    }
    ## failed coverage
    if(nrow(manifest %>% filter(! sampleID %in% snames)) >0 ){
    
    write.table(manifest %>% filter(! sampleID %in% snames),paste0(input.yaml$output.directory,"/results/",cohort.name,".samples.failed.coverage.2.txt"),quote=F,row.names=F)
    message(paste0("Wrote failed samples to the file ",cohort.name,".samples.failed.coverage.2.txt. Proceeding with the rest of the samples"))
    manifest <- manifest %>% filter( sampleID %in% snames)
    }
    write.table(manifest ,paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,"_manifest.txt"),row.names =F,quote =F, sep ="\t")
   
    if(length(countFiles) > 0){
    firstSample <- as.data.frame(data.table::fread(countFiles[1]))
    nCols <- dim(firstSample)[2]
    
    names(firstSample)[nCols] <- snames[1]
    snames <- snames[-1]
    annotations <- firstSample[,1:(nCols-1)]
    names(annotations)[1:3] <- c("chromosome","start","end")
    annotations$name <- paste0(annotations$chromosome,":",annotations$start,"-",annotations$end)
    
    countFiles <- countFiles[-1];
    firstSample <- round(firstSample[,nCols,drop=F])
    
    counts <- purrr::map_dfc(countFiles, function(x) {round(data.table::fread(x,drop=1:(nCols-1)))})
    
    names(counts) <- snames
    countmat <- as.matrix(cbind(firstSample,counts))
    
    selected.exons <- EDM::select.exons(cohort.counts = countmat,bin.length = (annotations$end-annotations$start)/1000,n.bins.reduced = 10000)
    countmat.selected <- countmat[selected.exons,]
    
    countmat.cor <- cor(countmat.selected)
    countData <- list(countmat = countmat, 
                      annotations=annotations, 
                      selected.exons=selected.exons,
                      correlations = countmat.cor)
    
    # saveRDS(countData,file=paste0(input.yaml$output.directory,"/results/",cohort.name,".exomedepth.cohort.rds"))
    
        men <- manifest$sampleID[manifest$sex == "M" | manifest$sex == "1"]
        women <- manifest$sampleID[manifest$sex == "F" |manifest$sex == "2"]
        
        ## cleanup CountData
        countData[[2]]$chromosome <- gsub("chr","",countData[[2]]$chromosome)
        ## autosome
        countmat.auto <- countmat[which(countData[[2]]$chromosome %in% c(1:22)),]
        annotations.auto <- annotations[which(countData[[2]]$chromosome %in% c(1:22)),]
        selected.exons.auto <- EDM::select.exons(cohort.counts = countmat.auto,bin.length = (annotations.auto$end-annotations.auto$start)/1000,n.bins.reduced = 10000)
        countmat.selected.auto <- countmat[selected.exons.auto,]
        
        countmat.cor.auto <- cor(countmat.selected.auto)
        countData.auto <- list(countmat = countmat.auto,
                               annotations=annotations.auto,
                               selected.exons=selected.exons.auto,
                               correlations = countmat.cor.auto)
        
        saveRDS(countData.auto,file=paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,".exomedepth.cohort.auto.rds"))
        
        ## chrX - Male
        if(length(men) > 0 ){
        countmat.men <- countmat[which(countData[[2]]$chromosome %in% "X" ),colnames(countmat) %in% men]
        annotations.men <- annotations[which(countData[[2]]$chromosome %in% "X"),]
        selected.exons.men <- EDM::select.exons(cohort.counts = countmat.men,bin.length = (annotations.men$end-annotations.men$start)/1000,n.bins.reduced = 10000)
        countmat.selected.men <- countmat.men[selected.exons.men,]
        
        countmat.cor.men <- cor(countmat.selected.men)
        countData.men <- list(countmat = countmat.men,
                              annotations=annotations.men,
                              selected.exons=selected.exons.men,
                              correlations = countmat.cor.men)
        
        saveRDS(countData.men,file=paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,".exomedepth.cohort.men.rds"))
        }
        ## chrX - Female
        if(length(women) > 0){
        countmat.women <- countmat[which(countData[[2]]$chromosome %in% "X" ),colnames(countmat) %in% women]
        annotations.women <- annotations[which(countData[[2]]$chromosome %in% "X"),]
        selected.exons.women <- EDM::select.exons(cohort.counts = countmat.women,bin.length = (annotations.women$end-annotations.women$start)/1000,n.bins.reduced = 10000)
        countmat.selected.women <- countmat.women[selected.exons.women,]
        
        countmat.cor.women <- cor(countmat.selected.women)
        countData.women <- list(countmat = countmat.women,
                                annotations=annotations.women,
                                selected.exons=selected.exons.women,
                                correlations = countmat.cor.women)
        
        saveRDS(countData.women,file=paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,".exomedepth.cohort.women.rds"))
        }
   }
}

