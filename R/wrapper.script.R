#' @export
#' @importFrom yaml  read_yaml
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%"
#' @importFrom clustermq Q
wrapper.script <- function(input){ 
  ## read yaml file  
  input.yaml <- yaml::read_yaml(input)
  
  ## create directories
  system(paste0("mkdir -p ",input.yaml$output.directory,"/coverage"))
  system(paste0("mkdir -p ",input.yaml$output.directory,"/results"))
  system(paste0("mkdir -p ",input.yaml$output.directory,"/logs"))
  system(paste0("mkdir -p ",input.yaml$output.directory,"/results/individual.edm.calls/"))
  system(paste0("mkdir -p ",input.yaml$output.directory,"/results/individual.ed.objects/"))
  ##check and  read  manifest
  if(is.null(input.yaml[["manifest"]])){
    message('Cannot proceed without a  manifest. Please look at the sample manifest provided in the documentation.\n')
    break;
  }else if(!is.null(input.yaml[["manifest"]]) & file.exists(input.yaml$manifest)){
    message('Reading the manifest...')
    manifest <- read.table(input.yaml$manifest, header = T, sep ="\t", stringsAsFactors = F)
    names(manifest) <- c("bam","sampleID","sex")
    
    nSamples <- NROW(manifest)
    nFemales <- NROW(manifest[manifest$sex == "M",])
    nMales <- NROW(manifest[manifest$sex == "F",])
    
    message(paste0('Found:\n * ',
                   nSamples,' individuals (',
                   nMales,' Males and ',nFemales,' females.)\n'))
    
  }else if(sum(is.na(manifest)) > 0){
    # check for empty entries in the manifest
    message('Entries in manifest cannot be empty. Please check your manifest before running.\n');
    break;
  }
  
  # Check for duplicate individual IDs
  if(sum(duplicated(sort(manifest$sampleID))) > 0){
    message('Found these duplicates in individual IDs\n', manifest$sampleID[duplicated(manifest$sampleID)] );
    break;
  }
  
  ## check bam files
  message('Checking the bam files...\n');
  if(sum(file.exists(as.character(manifest$bam))) < length(manifest$bam)){
    message('List of missing bam files:\n\n');
    print(paste(manifest$bam[!file.exists(as.character(manifest$bam))],"\n"))
    stop('Looks like some of your bam files are missing. Check your bam paths.\n\n\n')
  }
  
  ## check for index files
  message('Checking the index files for bams...\n');
  manifest$index <- NA;
  manifest$index[which(file.exists(paste0(as.character(manifest$bam),".bai")))] <- paste0(as.character(manifest$bam[which(file.exists(paste0(as.character(manifest$bam),".bai")))]),".bai")
  manifest$index[which(file.exists(paste0(gsub(pattern = "bam","bai",manifest$bam))))] <- paste0(gsub(pattern = "bam","bai",manifest$bam[which(file.exists(paste0(gsub(pattern = "bam","bai",manifest$bam))))]))
  manifest$index[which(file.exists(paste0(as.character(manifest$bam),".crai")))] <- paste0(as.character(manifest$bam[which(file.exists(paste0(as.character(manifest$bam),".crai")))]),".crai")
  manifest$index[which(file.exists(paste0(gsub(pattern = "cram","crai",manifest$bam))))] <- paste0(gsub(pattern = "cram","crai",manifest$bam[which(file.exists(paste0(gsub(pattern = "cram","crai",manifest$bam))))]))

  message(paste0('Samples with missing index files are written to results/missing.index.samples.txt \n sample_names:')) 
  message(paste0(manifest$sampleID[is.na(manifest$index)],'\n'));
  write.table( paste0('Index file not found for ',manifest[is.na(manifest$index),]),paste0(input.yaml$output.directory,"/results/missing.index.samples.txt"),row.names=F,quote=F, sep ="\t") 
  manifest <- manifest[!is.na(manifest$index),]
  # write.table(manifest[,1:3],paste0(input.yaml$output.directory,"/results/",input.yaml$manifest),row.names =F,quote =F, sep ="\t")
  write.table(manifest[,1:3],paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,"_manifest.txt"),row.names =F,quote =F, sep ="\t")
  # for (i in 1:length(manifest$bam)){
  #  if(file.exists(gsub(pattern = "bam","bai",manifest$bam[i]))){
  #    manifest$bai[i] <- gsub(pattern = "bam","bai",manifest$bam[i])
  #  } else if(file.exists(paste0(as.character(manifest$bam[i]),".bai"))){
  #    manifest$bai[i] <- paste0(as.character(manifest$bam[i]),".bai")
  #  } else{
  #    message(paste0('Samples with missing index files are written to temp/missing.index.samples.txt\n'));
  #    write.table( paste0('Index file not found for ',manifest$bam[i]),paste0(input.yaml$output.directory,"/temp/missing.index.samples.txt"),row.names=F,quote=F)
      ## Edit manifest based on index
  #    manifest <- manifest[
  #  }
  #}
 
 
  ## edit this based on cluster template
  ## clustermq template
  if(is.null(input.yaml[["hpc.scheduler"]])){
    message('Cannot proceed without a Scheduler type\n')
    break;
  }else{
    if(is.null(input.yaml[["clustermq.template"]])){ 
      options(
        clustermq.scheduler = paste0(input.yaml$hpc.scheduler),
        clustermq.template = paste0(find.package("EDM"),"/template/",input.yaml$hpc.scheduler,"_template"),
        clustermq.defaults = list(conda=paste0(input.yaml$env.name))
      )
    } else {
      options(
        clustermq.scheduler = paste0(input.yaml$hpc.scheduler),
        clustermq.template = paste0(input.yaml$clustermq.template),
        clustermq.defaults = list(conda=paste0(input.yaml$env.name))
      )
    }
  }
  
 if(is.null(input.yaml[["env.name"]])){
    message('########### Cannot proceed without a environment name. ###########\n')
    break;
  }

  if(is.null(input.yaml[["memory"]])){
    message('########### Cannot proceed without memory value. ###########\n')
    break;
  }

  
  EDM::check.exons.def(input=input)

  if( is.null(input.yaml[["precomputed.coverage"]]) | input.yaml[["precomputed.coverage"]] ==F){
    ## get Counts step
    #  source("clustermq_getCounts.R")
    nSamples = NROW(manifest)
    message(' Count step...   \n')
    clustermq::Q(get.bam.counts.mosdepth,pkgs="EDM",task_id=1:nSamples,input = input,n_jobs=nSamples, max_calls_worker = 1, template=list(job_name="compute_coverage", env.name = input.yaml$env.name ,memory = (input.yaml$memory*1024), mem = paste0(input.yaml$memory,"G")),export = list(input = input),timeout = 72000)
    
    ## move log files
    system(paste0("mv compute_coverage* ",input.yaml$output.directory,"/logs/."))
    
    ## gather counts
    message(' Gathering Counts...   \n')
    clustermq::Q(gather.mosdepth.coverage,pkgs="EDM",task_id=1,input = input,n_jobs=1, max_calls_worker = 1,template=list(job_name="gather_coverage",env.name = input.yaml$env.name ,memory = (input.yaml$memory*1024), mem = paste0(input.yaml$memory,"G")), export = list(input = input))
    
    #  source("gather_counts.R") 
    ## move log files
    system(paste0("mv gather_coverage* ",input.yaml$output.directory,"/logs/."))
  } else {
    
    if(is.null(input.yaml[["precomputed.coverage.autosomes"]]) | input.yaml[["precomputed.coverage.autosomes"]] ==F ){
      message('########### Cannot proceed without a file name. ###########\n')
      break;
    } else {
      if(file.exists(input.yaml$precomputed.coverage.autosomes)){
        precomputed_auto <- readRDS(input.yaml$precomputed.coverage.autosomes)
        if(length(colnames(precomputed_auto[["countmat"]]) %in% manifest$sampleID) != nrow(manifest)){
          message(' Not all samples from Manifest are present in precomputed.coverage.autosomes')
          message(paste0('These samples are missing: \n',manifest$sampleID[!manifest$sampleID %in% colnames(precomputed_auto[["countmat"]]) ]))
        } else {
          saveRDS(precomputed_auto, file = paste0(input.yaml$output.directory, 
                                                "/results/", input.yaml$cohort.name, ".exomedepth.cohort.auto.rds"))
        }
      }
    }
    
    if(is.null(input.yaml[["precomputed.coverage.men"]]) | input.yaml[["precomputed.coverage.men"]] ==F ){
      message('########### Cannot proceed without a file name. ###########\n')
      break;
    } else {
      if(file.exists(input.yaml$precomputed.coverage.men)) {
        precomputed_men <- readRDS(input.yaml$precomputed.coverage.men)
        if(length(colnames(precomputed_men[["countmat"]]) %in% manifest$sampleID) != nrow( manifest[manifest$sex == "M",]) ) {
          message(' Not all samples from Manifest are present in precomputed.coverage.autosomes')
          message(paste0('These samples are missing: \n',manifest[manifest$sex == "M",]$sampleID[! manifest[manifest$sex == "M",]$sampleID %in% colnames(precomputed_men[["countmat"]]) ]))
        } else {
          saveRDS(precomputed_men, file = paste0(input.yaml$output.directory, 
                                                  "/results/", input.yaml$cohort.name, ".exomedepth.cohort.men.rds"))
        }
      }
    }
    
    if(is.null(input.yaml[["precomputed.coverage.women"]]) | input.yaml[["precomputed.coverage.autosomes"]] ==F ){
      message('########### Cannot proceed without a file name. ###########\n')
      break;
    } else {
      if(file.exists(input.yaml$precomputed.coverage.women)) {
        precomputed_women <- readRDS(input.yaml$precomputed.coverage.women)
        if(length(colnames(precomputed_women[["countmat"]]) %in% manifest$sampleID) != nrow( manifest[manifest$sex == "F",]) ) {
          message(' Not all samples from Manifest are present in precomputed.coverage.autosomes')
          message(paste0('These samples are missing: \n',manifest[manifest$sex == "F",]$sampleID[! manifest[manifest$sex == "F",]$sampleID %in% colnames(precomputed_women[["countmat"]]) ]))
        } else {
          saveRDS(precomputed_women, file = paste0(input.yaml$output.directory, 
                                                 "/results/", input.yaml$cohort.name, ".exomedepth.cohort.women.rds"))
        }
      }
    }
    
    write.table(manifest, paste0(input.yaml$output.directory, 
                                 "/results/", input.yaml$cohort.name, "_manifest.txt"), 
                row.names = F, quote = F, sep = "\t")
  }
  
  
  ## call variants
  message('Calling variants...   \n')
  manifest <-  read.table(paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,"_manifest.txt"),header=T, stringsAsFactors=F,sep ="\t")
   nSamples = NROW(manifest)
   call.variants <- suppressWarnings(call.variants)
  clustermq::Q(call.variants,pkgs=list("EDM","ExomeDepth"),columnIndex=1:nSamples,n_jobs=nSamples,input=input, max_calls_worker = 1,template=list(job_name="call_variants",env.name = input.yaml$env.name ,memory = (input.yaml$memory*1024), mem = paste0(input.yaml$memory,"G")), export = list(input = input),timeout = 1440000)
  
  ## move log files
  system(paste0("mv call_variants* ",input.yaml$output.directory,"/logs/."))
  
  ## gather variants
  message(' Gathering Calls...   \n')
  clustermq::Q(gather.variant.calls,pkgs=list("EDM"),task_id=1,input = input,n_jobs=1, max_calls_worker = 1,template=list(job_name="gather_calls",env.name = input.yaml$env.name ,memory = (input.yaml$memory*1024), mem = paste0(input.yaml$memory,"G")), export = list(input = input))
  system(paste0("mv gather_calls* ",input.yaml$output.directory,"/logs/."))

  message(' Plotting correlation...   \n')
  clustermq::Q(plot.sample.max.correlation,pkgs=list("EDM"),input = input,n_jobs=1, max_calls_worker = 1,template=list(job_name="plotting_correlation",env.name = input.yaml$env.name ,memory = (input.yaml$memory*1024), mem = paste0(input.yaml$memory,"G")), export = list(input = input))

  system(paste0("mv plotting_correlation* ",input.yaml$output.directory,"/logs/."))
  system(paste0("rm ",input.yaml$output.directory,"/results/individual.edm.calls/*.edm.cor.txt"))
  
  ## summary data
  cohort.summary <- read.table(paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,".edm.summary.cohort.txt"),header=T,stringsAsFactors=F)
   
  ## summary
  message('Summary of number of calls in cohort \n')
  message('Mean of cohort: ',round(mean(cohort.summary$Freq)),'\n')
  message('Median of cohort: ',median(cohort.summary$Freq),'\n')
  message('Min calls: ',min(cohort.summary$Freq),' Max calls:',max(cohort.summary$Freq),'\n')
  message('Standard deviation: ',paste0(sd(cohort.summary$Freq),'\n'))
  return(message('Pipeline ran to complete'))
}   
