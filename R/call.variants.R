#' @importFrom yaml  read_yaml
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%"
#' @export
call.variants <- function(columnIndex,input){
  input.yaml <- yaml::read_yaml(input)
  input.columnIndex <- columnIndex
  # manifest <- read.table(paste0(input.yaml$output.directory, 
   #                             "/results/", input.yaml$manifest), header = T, sep = "\t", 
   #                      stringsAsFactors = F)
  manifest <- read.table(paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,"_manifest.txt"),header=T, sep = "\t",stringsAsFactors=F)
  names(manifest) <- c("bam", "sampleID", "sex")
  sample_name <- manifest$sampleID[columnIndex]
  exons <- read.table(paste0(input.yaml$output.directory, "/results/bed_file.bed"), 
                      stringsAsFactors = F)
  names(exons) <- c("chromosome", "start", "end")
  men <- manifest$sampleID[manifest$sex == "M" | manifest$sex == 
                             "1"]
  women <- manifest$sampleID[manifest$sex == "F" | manifest$sex == 
                               "2"]
  if (input.yaml$precomputed.controls) {
    cohort.object.auto <- readRDS(input.yaml[["precomputed.controls.autosomes"]])
    nControls <- dim(cohort.object.auto[["countmat"]])[2]
    batch <- readRDS(paste0(input.yaml$output.directory, 
                            "/results/", input.yaml$cohort.name, ".exomedepth.cohort.auto.rds"))
    cohort.object.auto$countmat <- cbind(cohort.object.auto$countmat, 
                                         batch$countmat[, columnIndex, drop = FALSE])
    columnIndex = nControls + 1
    cohort.object.auto$selected.exons <- EDM::select.exons(cohort.counts = cohort.object.auto$countmat, 
                                                           bin.length = (cohort.object.auto$annotations$end - 
                                                                           cohort.object.auto$annotations$start)/1000, n.bins.reduced = 10000)
    countmat.selected <- cohort.object.auto$countmat[cohort.object.auto$selected.exons, 
                                                     ]
    cohort.object.auto$correlations <- cor(countmat.selected)
  } else {
    cohort.object.auto <- readRDS(paste0(input.yaml$output.directory, 
                                         "/results/", input.yaml$cohort.name, ".exomedepth.cohort.auto.rds"))
  }
  reference_list <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                               columnIndex], reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                                                                 -columnIndex], correlations = cohort.object.auto[["correlations"]][-columnIndex, 
                                                                                                                                                                                                                                    columnIndex])
  if (is.null(reference_list$reference.choice) == F) {
    reference_set = apply(X = as.matrix(cohort.object.auto[["countmat"]][, 
                                                                         reference_list$reference.choice]), MAR = 1, FUN = sum)
    all_exons_auto = new("ExomeDepth", test = cohort.object.auto[["countmat"]][, 
                                                                               columnIndex], reference = reference_set, formula = "cbind(test,reference) ~ 1")
    all_exons_auto = ExomeDepth::CallCNVs(x = all_exons_auto, 
                                          transition.probability = as.numeric(input.yaml$transition.probability), 
                                          chromosome = cohort.object.auto[["annotations"]]$chromosome, 
                                          start = cohort.object.auto[["annotations"]]$start, 
                                          end = cohort.object.auto[["annotations"]]$end, name = cohort.object.auto[["annotations"]]$name)
  } else {
    message(paste0("\n *********** Sample ", sample_name, 
                   " failed for Autosome variant calling. ************* \n"))
  }
  if (sample_name %in% men) {
    if (input.yaml$precomputed.controls) {
      cohort.object.x <- readRDS(input.yaml[["precomputed.controls.men"]])
      nControls <- dim(cohort.object.x[["countmat"]])[2]
      batch <- readRDS(paste0(input.yaml$output.directory, 
                              "/results/", input.yaml$cohort.name, ".exomedepth.cohort.men.rds"))
      columnIndex <- which(colnames(batch$countmat) %in% 
                             sample_name)
      cohort.object.x$countmat <- cbind(cohort.object.x$countmat, 
                                        batch$countmat[, columnIndex, drop = FALSE])
      columnIndex = nControls + 1
      cohort.object.x$selected.exons <- EDM::select.exons(cohort.counts = cohort.object.x$countmat, 
                                                          bin.length = (cohort.object.x$annotations$end - 
                                                                          cohort.object.x$annotations$start)/1000, n.bins.reduced = 10000)
      countmat.selected <- cohort.object.x$countmat[cohort.object.x$selected.exons, 
                                                    ]
      cohort.object.x$correlations <- cor(countmat.selected)
    } else {
      cohort.object.x <- readRDS(paste0(input.yaml$output.directory, 
                                        "/results/", input.yaml$cohort.name, ".exomedepth.cohort.men.rds"))
    }
    reference_list <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                              columnIndex], reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                                                                             -columnIndex], correlations = cohort.object.x[["correlations"]][-columnIndex, 
                                                                                                                                                                                                                             columnIndex])
    if (is.null(reference_list$reference.choice) == F) {
      reference_set = apply(X = as.matrix(cohort.object.x[["countmat"]][, 
                                                                        reference_list$reference.choice]), MAR = 1, FUN = sum)
      all_exons_x = new("ExomeDepth", test = cohort.object.x[["countmat"]][, 
                                                                           columnIndex], reference = reference_set, formula = "cbind(test,reference) ~ 1")
      all_exons_x = ExomeDepth::CallCNVs(x = all_exons_x, 
                                         transition.probability = as.numeric(input.yaml$transition.probability), 
                                         chromosome = cohort.object.x[["annotations"]]$chromosome, 
                                         start = cohort.object.x[["annotations"]]$start, 
                                         end = cohort.object.x[["annotations"]]$end, name = cohort.object.x[["annotations"]]$name)
    } else {
      message(paste0("\n *********** Sample ", sample_name, 
                     " failed for chrX variant calling. ************* \n"))
    }
  } else {
    if (input.yaml$precomputed.controls) {
      cohort.object.x <- readRDS(input.yaml[["precomputed.controls.women"]])
      nControls <- dim(cohort.object.x[["countmat"]])[2]
      batch <- readRDS(paste0(input.yaml$output.directory, 
                              "/results/", input.yaml$cohort.name, ".exomedepth.cohort.women.rds"))
      columnIndex <- which(colnames(batch$countmat) %in% 
                             sample_name)
      cohort.object.x$countmat <- cbind(cohort.object.x$countmat, 
                                        batch$countmat[, columnIndex, drop = FALSE])
      columnIndex = nControls + 1
      cohort.object.x$selected.exons <- EDM::select.exons(cohort.counts = cohort.object.x$countmat, 
                                                          bin.length = (cohort.object.x$annotations$end - 
                                                                          cohort.object.x$annotations$start)/1000, n.bins.reduced = 10000)
      countmat.selected <- cohort.object.x$countmat[cohort.object.x$selected.exons, 
                                                    ]
      cohort.object.x$correlations <- cor(countmat.selected)
    } else {
      cohort.object.x <- readRDS(paste0(input.yaml$output.directory, 
                                        "/results/", input.yaml$cohort.name, ".exomedepth.cohort.women.rds"))
    }
    reference_list <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                              columnIndex], reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                                                                             -columnIndex], correlations = cohort.object.x[["correlations"]][-columnIndex, 
                                                                                                                                                                                                                             columnIndex])
    if (is.null(reference_list$reference.choice) == F) {
      reference_set = apply(X = as.matrix(cohort.object.x[["countmat"]][, 
                                                                        reference_list$reference.choice]), MAR = 1, FUN = sum)
      all_exons_x = new("ExomeDepth", test = cohort.object.x[["countmat"]][, 
                                                                           columnIndex], reference = reference_set, formula = "cbind(test,reference) ~ 1")
      all_exons_x = ExomeDepth::CallCNVs(x = all_exons_x, 
                                         transition.probability = as.numeric(input.yaml$transition.probability), 
                                         chromosome = cohort.object.x[["annotations"]]$chromosome, 
                                         start = cohort.object.x[["annotations"]]$start, 
                                         end = cohort.object.x[["annotations"]]$end, name = cohort.object.x[["annotations"]]$name)
    } else {
      message(paste0("\n *********** Sample ", sample_name, 
                     " failed for chrX. ************* \n"))
    }
  }
  if (dim(all_exons_x@CNV.calls)[1] > 0 & dim(all_exons_auto@CNV.calls)[1] > 
      0) {
    all_exons_auto@CNV.calls$sample <- sample_name
    all_exons_x@CNV.calls$sample <- sample_name
    calls <- rbind(all_exons_auto@CNV.calls, all_exons_x@CNV.calls)
    annotated <- EDM::variant.annotations(calls)
    calls_all <- cbind(calls, annotated[, 4:dim(annotated)[2]])
    all_exons <- list(auto_calls = all_exons_auto, x_calls = all_exons_x)
    write.table(calls_all, paste0(input.yaml$output.directory, 
                                  "/results/individual.edm.calls/", sample_name, ".edm.calls.txt"), 
                row.names = F, quote = F, sep = "\t")
    saveRDS(all_exons, paste0(input.yaml$output.directory, 
                              "/results/individual.ed.objects/", sample_name, ".ed.object.rds"))
    cor.data <- data.frame(cor = cor(all_exons$auto_calls@test, 
                                     all_exons$auto_calls@reference), sample_name = sample_name)
    write.table(cor.data, paste0(input.yaml$output.directory, 
                                 "/results/individual.edm.calls/", sample_name, ".edm.cor.txt"), 
                row.names = F, quote = F, sep = "\t")
    return(all_exons)
  } else {
  if (dim(all_exons_auto@CNV.calls)[1] > 0) {
    all_exons_auto@CNV.calls$sample <- sample_name
    calls <- rbind(all_exons_auto@CNV.calls)
    annotated <- EDM::variant.annotations(calls)
    calls_all <- cbind(calls, annotated[, 4:dim(annotated)[2]])
    all_exons <- list(auto_calls = all_exons_auto)
    write.table(calls_all, paste0(input.yaml$output.directory, 
                                  "/results/individual.edm.calls/", sample_name, ".edm.calls.txt"), 
                row.names = F, quote = F, sep = "\t")
    saveRDS(all_exons, paste0(input.yaml$output.directory, 
                              "/results/individual.ed.objects/", sample_name, ".ed.object.rds"))
    cor.data <- data.frame(cor = cor(all_exons$auto_calls@test, 
                                     all_exons$auto_calls@reference), sample_name = sample_name)
    write.table(cor.data, paste0(input.yaml$output.directory, 
                                 "/results/individual.edm.calls/", sample_name, ".edm.cor.txt"), 
                row.names = F, quote = F, sep = "\t")
    return(all_exons)
  } else {
    return(NULL)
  }
}
}

