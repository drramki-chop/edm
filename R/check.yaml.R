#' @importFrom yaml  read_yaml
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%"
#' @export
check.yaml <- function(input){
  input.yaml <- yaml::read_yaml(input)
  
  yaml_fields <- c("cohort.name",
  "manifest",
  "output.directory",
  "transition.probability",
  "bed.file",
  "hpc.scheduler",
  "env.name",
  "fasta",
  "clustermq.template",
  "memory",
  "precomputed.controls",
  "precomputed.controls.autosomes",
  "precomputed.controls.men",
  "precomputed.controls.women",
  "reproducibility",
  "random.controls",
  "iterations",
  "ped",
  "precomputed.coverage",
  "precomputed.coverage.autosomes",
  "precomputed.coverage.men",
 "precomputed.coverage.women")
  
  if(is.null(input.yaml[["manifest"]])){
    message('Cannot proceed without a  manifest. Please look at the sample manifest provided in the documentation.\n')
    break;
  }
  if(is.null(input.yaml[["hpc.scheduler"]])){
    message('Cannot proceed without a Scheduler type\n')
    break;
  }
  
  if(is.null(input.yaml[["cohort.name"]])){
    input.yaml[["cohort.name"]] <- "cnv_run"
  }
  
  if(is.null(input.yaml[["output.directory"]])){
    input.yaml[["output.directory"]] <- getwd()
  }
  
  if(is.null(input.yaml[["transition.probability"]])){
    input.yaml[["transition.probability"]] <- "1e-4"
  }
  
  if(is.null(input.yaml[["env.name"]])){
    input.yaml[["env.name"]] <- "edm_env"
  }
  
  if(is.null(input.yaml[["memory"]])){
    input.yaml[["memory"]] <- 24
  }
  
  system(paste0("mkdir -p ",input.yaml$output.directory))
  system(paste0("mkdir -p ",input.yaml$output.directory,"/coverage"))
  system(paste0("mkdir -p ",input.yaml$output.directory,"/results"))
  system(paste0("mkdir -p ",input.yaml$output.directory,"/logs"))
  system(paste0("mkdir -p ",input.yaml$output.directory,"/results/individual.edm.calls/"))
  system(paste0("mkdir -p ",input.yaml$output.directory,"/results/individual.ed.objects/"))
  
  if( length(yaml_fields[!yaml_fields %in% names(input.yaml)]) > 0){
    
      foo <- yaml_fields[!yaml_fields %in% names(input.yaml)]
      
    for(i in 1:length(foo) ){
      tmp <- foo[i]
      input.yaml[[tmp]] <- NA
    }
   yaml::write_yaml(input.yaml, paste0(input.yaml$output.directory,"/results/input.yaml"))
    
  } else{
   yaml::write_yaml(input.yaml, paste0(input.yaml$output.directory,"/results/input.yaml"))
  }
    return(input.yaml)
}
