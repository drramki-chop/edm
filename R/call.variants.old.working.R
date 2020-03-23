#' @importFrom yaml  read_yaml
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%"
#' @export
call.variants.old <- function(columnIndex,input){

  input.yaml <- yaml::read_yaml(input)
  
  manifest <- read.table(paste0(input.yaml$output.directory,"/results/",input.yaml$manifest), header = T, sep ="\t", stringsAsFactors = F)
  names(manifest) <- c("bam","sampleID","sex")
  
  exons <- read.table(paste0(input.yaml$output.directory,"/results/bed_file.bed"),stringsAsFactors=F)
  names(exons) <- c("chromosome","start","end")
  men <- manifest$sampleID[manifest$sex == "M" | manifest$sex == "1"]
  women <- manifest$sampleID[manifest$sex == "F" |manifest$sex == "2"]
  
 cohort.object.auto <- readRDS(paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,".exomedepth.cohort.auto.rds"))
  sample_name <- colnames(cohort.object.auto[["countmat"]])[columnIndex];

  reference_list <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]],columnIndex],
                                                reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]],-columnIndex],
                                                correlations = cohort.object.auto[["correlations"]][-columnIndex,columnIndex])
  if(is.null(reference_list$reference.choice) == F  ) {
  reference_set = apply(
    X = as.matrix(cohort.object.auto[["countmat"]][, reference_list$reference.choice]),
    MAR=1, FUN=sum)
  all_exons_auto = new('ExomeDepth', test=cohort.object.auto[["countmat"]][,columnIndex],
                       reference=reference_set,
                       formula = 'cbind(test,reference) ~ 1')
  all_exons_auto = ExomeDepth::CallCNVs(x = all_exons_auto, transition.probability=as.numeric(input.yaml$transition.probability),
                                        chromosome= cohort.object.auto[["annotations"]]$chromosome, start=cohort.object.auto[["annotations"]]$start,
                                        end=cohort.object.auto[["annotations"]]$end, name=cohort.object.auto[["annotations"]]$name)

  # all_exons_auto@CNV.calls$sample <- sample_name;
  } else {
    message(paste0("\n *********** Sample ", sample_name," failed for Autosome variant calling. ************* \n"))
  }
  
  if(sample_name %in% men){
    cohort.object.x <- readRDS(paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,".exomedepth.cohort.men.rds"))
    
    # sample_name <- colnames(cohort.object.x[["countmat"]])[columnIndex];
    
    reference_list <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], which(sample_name %in% men)],
                                                  reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]],- which(sample_name %in% men)],
                                                  correlations = cohort.object.x[["correlations"]][- which(sample_name %in% men), which(sample_name %in% men)])
    if(is.null(reference_list$reference.choice) == F  ) {
    reference_set = apply(
      X = as.matrix(cohort.object.x[["countmat"]][, reference_list$reference.choice]),
      MAR=1, FUN=sum)
    all_exons_x = new('ExomeDepth', test=cohort.object.x[["countmat"]][, which(sample_name %in% men)],
                      reference=reference_set,
                      formula = 'cbind(test,reference) ~ 1')
    all_exons_x = ExomeDepth::CallCNVs(x = all_exons_x, transition.probability=as.numeric(input.yaml$transition.probability),
                                       chromosome= cohort.object.x[["annotations"]]$chromosome, start=cohort.object.x[["annotations"]]$start,
                                       end=cohort.object.x[["annotations"]]$end, name=cohort.object.x[["annotations"]]$name)
    
   # all_exons_x@CNV.calls$sample <- sample_name;
  } else {
    message(paste0("\n *********** Sample ", sample_name," failed for chrX variant calling. ************* \n"))
  }
  } else {
    cohort.object.x <- readRDS(paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,".exomedepth.cohort.women.rds"))
    
   #  sample_name <- colnames(cohort.object.x[["countmat"]])[columnIndex];
    
    reference_list <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], which(sample_name %in% women)],
                                                  reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]],- which(sample_name %in% women)],
                                                  correlations = cohort.object.x[["correlations"]][- which(sample_name %in% women), which(sample_name %in% women)])
    
    if(is.null(reference_list$reference.choice) == F  ) {
    reference_set = apply(
      X = as.matrix(cohort.object.x[["countmat"]][, reference_list$reference.choice]),
      MAR=1, FUN=sum)
    all_exons_x = new('ExomeDepth', test=cohort.object.x[["countmat"]][, which(sample_name %in% women)],
                      reference=reference_set,
                      formula = 'cbind(test,reference) ~ 1')
    all_exons_x = ExomeDepth::CallCNVs(x = all_exons_x, transition.probability=as.numeric(input.yaml$transition.probability),
                                       chromosome= cohort.object.x[["annotations"]]$chromosome, start=cohort.object.x[["annotations"]]$start,
                                       end=cohort.object.x[["annotations"]]$end, name=cohort.object.x[["annotations"]]$name)
    
   #  all_exons_x@CNV.calls$sample <- sample_name;   
    } else {
      message(paste0("\n *********** Sample ", sample_name," failed for chrX. ************* \n"))
    }
  }
  if (dim(all_exons_x@CNV.calls)[1] > 0 & dim(all_exons_auto@CNV.calls)[1] > 0){
    all_exons_auto@CNV.calls$sample <- sample_name;
    all_exons_x@CNV.calls$sample <- sample_name;
    calls <- rbind(all_exons_auto@CNV.calls,all_exons_x@CNV.calls);
    annotated <- EDM::variant.annotations(calls)
    calls_all <- cbind(calls,annotated[,4:dim(annotated)[2]])
    all_exons <- list(auto_calls = all_exons_auto,x_calls = all_exons_x)
    write.table(calls_all,paste0(input.yaml$output.directory,"/results/individual.edm.calls/",sample_name,".edm.calls.txt"), row.names = F,quote=F, sep ="\t")
    saveRDS(all_exons,paste0(input.yaml$output.directory,"/results/individual.ed.objects/",sample_name,".ed.object.rds"))
    cor.data <- data.frame(cor = cor(all_exons$auto_calls@test, all_exons$auto_calls@reference) , sample_name = sample_name)
    write.table(cor.data,paste0(input.yaml$output.directory,"/results/individual.edm.calls/",sample_name,".edm.cor.txt"), row.names = F,quote=F,sep="\t") 
    return(all_exons)
  } 
  if (dim(all_exons_auto@CNV.calls)[1] > 0){
    all_exons_auto@CNV.calls$sample <- sample_name;
    calls <- rbind(all_exons_auto@CNV.calls);
    annotated <- EDM::variant.annotations(calls)
    calls_all <- cbind(calls,annotated[,4:dim(annotated)[2]])
    all_exons <- list(auto_calls = all_exons_auto)
    write.table(calls_all,paste0(input.yaml$output.directory,"/results/individual.edm.calls/",sample_name,".edm.calls.txt"), row.names = F,quote=F, sep ="\t")
    saveRDS(all_exons,paste0(input.yaml$output.directory,"/results/individual.ed.objects/",sample_name,".ed.object.rds"))
    cor.data <- data.frame(cor = cor(all_exons$auto_calls@test, all_exons$auto_calls@reference) , sample_name = sample_name)
    write.table(cor.data,paste0(input.yaml$output.directory,"/results/individual.edm.calls/",sample_name,".edm.cor.txt"), row.names = F,quote=F,sep ="\t")
    return(all_exons)
  } else {
    return(NULL)
  } 
  
}
