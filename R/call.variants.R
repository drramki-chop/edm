#' @importFrom yaml  read_yaml
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%"
#' @export
call.variants <- function(columnIndex,input){
  input.yaml <- yaml::read_yaml(input)
  input.columnIndex <- columnIndex
  manifest <- read.table(paste0(input.yaml$output.directory, 
                                "/results/", input.yaml$cohort.name, "_manifest.txt"), 
                         header = T, sep = "\t", stringsAsFactors = F)
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
    columnIndex_n = which(colnames(cohort.object.auto[["countmat"]]) %in% 
                            sample_name)
    cohort.object.auto$countmat <- cbind(cohort.object.auto$countmat, 
                                         batch$countmat[, columnIndex_n, drop = FALSE])
    columnIndex_n = nControls + 1
    cohort.object.auto$selected.exons <- EDM::select.exons(cohort.counts = cohort.object.auto$countmat, 
                                                           bin.length = (cohort.object.auto$annotations$end - 
                                                                           cohort.object.auto$annotations$start)/1000, n.bins.reduced = 10000)
    countmat.selected <- cohort.object.auto$countmat[cohort.object.auto$selected.exons, 
                                                     ]
    cohort.object.auto$correlations <- cor(countmat.selected)
    reference_list <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                 columnIndex_n], reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                                                                     -columnIndex_n], correlations = cohort.object.auto[["correlations"]][-columnIndex_n, 
                                                                                                                                                                                                                                          columnIndex_n])
  } else {
    cohort.object.auto <- readRDS(paste0(input.yaml$output.directory, 
                                         "/results/", input.yaml$cohort.name, ".exomedepth.cohort.auto.rds"))
    if (is.null(input.yaml$ped) == F) {
      ped <- read.table(paste0(input.yaml$ped), stringsAsFactors = F, 
                        header = F, fill = T)
      columnIndex_fam <- which(colnames(cohort.object.auto[["countmat"]]) %in% 
                                 (ped %>% filter(V1 %in% (ped %>% filter(V2 %in% 
                                                                           sample_name) %>% pull(1))) %>% pull(2)))
      columnIndex_n = which(colnames(cohort.object.auto[["countmat"]]) %in% 
                              sample_name)
      reference_list <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                   columnIndex_n], reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                                                                       -columnIndex_fam], correlations = cohort.object.auto[["correlations"]][-columnIndex_fam, 
                                                                                                                                                                                                                                              columnIndex_n])
    } else {
      columnIndex_n = which(colnames(cohort.object.auto[["countmat"]]) %in% 
                              sample_name)
      reference_list <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                   columnIndex_n], reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                                                                       -columnIndex_n], correlations = cohort.object.auto[["correlations"]][-columnIndex_n, 
                                                                                                                                                                                                                                            columnIndex_n])
    }
  }
  
  if (is.null(reference_list$reference.choice) == F) {
    reference_set = apply(X = as.matrix(cohort.object.auto[["countmat"]][, 
                                                                         reference_list$reference.choice]), MAR = 1, FUN = sum)
    all_exons_auto = new("ExomeDepth", test = cohort.object.auto[["countmat"]][, 
                                                                               columnIndex_n], reference = reference_set, formula = "cbind(test,reference) ~ 1")
    all_exons_auto = ExomeDepth::CallCNVs(x = all_exons_auto, 
                                          transition.probability = as.numeric(input.yaml$transition.probability), 
                                          chromosome = cohort.object.auto[["annotations"]]$chromosome, 
                                          start = cohort.object.auto[["annotations"]]$start, 
                                          end = cohort.object.auto[["annotations"]]$end, name = cohort.object.auto[["annotations"]]$name)
  } else {
    message(paste0("\n *********** Sample ", sample_name, 
                   " failed for Autosome variant calling. ************* \n"))
  }
  references <- data.frame(matrix(nrow =3,ncol=2))
  names(references) = c("test","value")
  references$test[1] <- "autosomes_reference"
  references$test[2] <- "chrx_reference"
  references$test[3] <- "correlation"
  references$value[1] <- paste(reference_list$reference.choice,collapse=";")
  
  if (sample_name %in% men) {
    if (input.yaml$precomputed.controls) {
      cohort.object.x <- readRDS(input.yaml[["precomputed.controls.men"]])
      nControls <- dim(cohort.object.x[["countmat"]])[2]
      batch <- readRDS(paste0(input.yaml$output.directory, 
                              "/results/", input.yaml$cohort.name, ".exomedepth.cohort.men.rds"))
      columnIndex_x <- which(colnames(batch$countmat) %in% 
                               sample_name)
      cohort.object.x$countmat <- cbind(cohort.object.x$countmat, 
                                        batch$countmat[, columnIndex_x, drop = FALSE])
      columnIndex_x = nControls + 1
      cohort.object.x$selected.exons <- EDM::select.exons(cohort.counts = cohort.object.x$countmat, 
                                                          bin.length = (cohort.object.x$annotations$end - 
                                                                          cohort.object.x$annotations$start)/1000, n.bins.reduced = 10000)
      countmat.selected <- cohort.object.x$countmat[cohort.object.x$selected.exons, 
                                                    ]
      cohort.object.x$correlations <- cor(countmat.selected)
      reference_list <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                columnIndex_x], reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                                                                                 -columnIndex_x], correlations = cohort.object.x[["correlations"]][-columnIndex_x, 
                                                                                                                                                                                                                                   columnIndex_x])
    } else {
      cohort.object.x <- readRDS(paste0(input.yaml$output.directory, 
                                        "/results/", input.yaml$cohort.name, ".exomedepth.cohort.men.rds"))
      if (is.null(input.yaml$ped) == F) {
        ped <- read.table(paste0(input.yaml$ped), stringsAsFactors = F, 
                          header = F, fill = T)
        columnIndex_fam <- which(colnames(cohort.object.x[["countmat"]]) %in% 
                                   (ped %>% filter(V1 %in% (ped %>% filter(V2 %in% 
                                                                             sample_name) %>% pull(1))) %>% pull(2)))
        columnIndex_x <- which(colnames(cohort.object.x[["countmat"]]) %in% 
                                 sample_name)
        reference_list <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                  columnIndex_x], reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                                                                                   -columnIndex_fam], correlations = cohort.object.x[["correlations"]][-columnIndex_fam, 
                                                                                                                                                                                                                                       columnIndex_x])
      } else {
        columnIndex_x <- which(colnames(cohort.object.x[["countmat"]]) %in% 
                                 sample_name)
        reference_list <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                  columnIndex_x], reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                                                                                   -columnIndex_x], correlations = cohort.object.x[["correlations"]][-columnIndex_x, 
                                                                                                                                                                                                                                     columnIndex_x])
      }
    }
    if (is.null(reference_list$reference.choice) == F) {
      reference_set = apply(X = as.matrix(cohort.object.x[["countmat"]][, 
                                                                        reference_list$reference.choice]), MAR = 1, FUN = sum)
      all_exons_x = new("ExomeDepth", test = cohort.object.x[["countmat"]][, 
                                                                           columnIndex_x], reference = reference_set, formula = "cbind(test,reference) ~ 1")
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
      columnIndex_x <- which(colnames(batch$countmat) %in% 
                               sample_name)
      cohort.object.x$countmat <- cbind(cohort.object.x$countmat, 
                                        batch$countmat[, columnIndex_x, drop = FALSE])
      columnIndex_x = nControls + 1
      cohort.object.x$selected.exons <- EDM::select.exons(cohort.counts = cohort.object.x$countmat, 
                                                          bin.length = (cohort.object.x$annotations$end - 
                                                                          cohort.object.x$annotations$start)/1000, n.bins.reduced = 10000)
      countmat.selected <- cohort.object.x$countmat[cohort.object.x$selected.exons, 
                                                    ]
      cohort.object.x$correlations <- cor(countmat.selected)
      reference_list <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                columnIndex_x], reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                                                                                 -columnIndex_x], correlations = cohort.object.x[["correlations"]][-columnIndex_x, 
                                                                                                                                                                                                                                   columnIndex_x])
    } else {
      cohort.object.x <- readRDS(paste0(input.yaml$output.directory, 
                                        "/results/", input.yaml$cohort.name, ".exomedepth.cohort.women.rds"))
      if (is.null(input.yaml$ped) == F) {
        ped <- read.table(paste0(input.yaml$ped), stringsAsFactors = F, 
                          header = F, fill = T)
        columnIndex_fam <- which(colnames(cohort.object.x[["countmat"]]) %in% 
                                   (ped %>% filter(V1 %in% (ped %>% filter(V2 %in% 
                                                                             sample_name) %>% pull(1))) %>% pull(2)))
        columnIndex_x <- which(colnames(cohort.object.x[["countmat"]]) %in% 
                                 sample_name)
        reference_list <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                  columnIndex_x], reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                                                                                   -columnIndex_fam], correlations = cohort.object.x[["correlations"]][-columnIndex_fam, 
                                                                                                                                                                                                                                       columnIndex_x])
      }else {
        columnIndex_x <- which(colnames(cohort.object.x[["countmat"]]) %in% 
                                 sample_name)
        reference_list <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                  columnIndex_x], reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                                                                                   -columnIndex_x], correlations = cohort.object.x[["correlations"]][-columnIndex_x, 
                                                                                                                                                                                                                                     columnIndex_x])
      }
    }
    if (is.null(reference_list$reference.choice) == F) {
      reference_set = apply(X = as.matrix(cohort.object.x[["countmat"]][, 
                                                                        reference_list$reference.choice]), MAR = 1, FUN = sum)
      all_exons_x = new("ExomeDepth", test = cohort.object.x[["countmat"]][, 
                                                                           columnIndex_x], reference = reference_set, formula = "cbind(test,reference) ~ 1")
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
  
  references$value[2] <- paste(reference_list$reference.choice,collapse=";")
  
  if (dim(all_exons_x@CNV.calls)[1] > 0 & dim(all_exons_auto@CNV.calls)[1] > 0) {
    all_exons_auto@CNV.calls$sample <- sample_name
    all_exons_x@CNV.calls$sample <- sample_name
    all_exons <- list(auto_calls = all_exons_auto, x_calls = all_exons_x)
    calls.first <- rbind(all_exons_auto@CNV.calls, all_exons_x@CNV.calls)
    annotated <- EDM::variant.annotations(calls.first)
    calls.annotated <- cbind(calls.first, annotated[, 4:dim(annotated)[2]])
    
    if (input.yaml$reproducibility == T) {
      calls.first$reproducibility <- 0
      sample.del <- calls.first[calls.first$type %in% c("deletion"), 
                                ]
      sample.dup <- calls.first[calls.first$type %in% c("duplication"), 
                                ]
      sample.del.gr <- GenomicRanges::GRanges(sample.del$id)
      sample.dup.gr <- GenomicRanges::GRanges(sample.dup$id)
      if (is.null(input.yaml$iterations) == F) {
        iterations = as.numeric(input.yaml$iterations)
      } else {
        iterations <- 1000
      }
      
      for (iter in 1:iterations) {
        n.random.controls <- 100
        if (input.yaml$precomputed.controls) {
          rand.samples <- sort(sample(c(1:nControls), 
                                      n.random.controls))
          reference_list_auto <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                            columnIndex_n], reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                                                                                rand.samples], correlations = cohort.object.auto[["correlations"]][rand.samples, 
                                                                                                                                                                                                                                                   columnIndex_n])
          nControls <- dim(cohort.object.x[["countmat"]])[2]
          rand.samples <- sort(sample(c(1:nControls), 
                                      n.random.controls))
          reference_list_x <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                      columnIndex_x], reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                                                                                       rand.samples], correlations = cohort.object.x[["correlations"]][rand.samples, 
                                                                                                                                                                                                                                       columnIndex_x])
        } else {
          nControls <- dim(cohort.object.auto[["countmat"]])[2]
          rand.samples <- sort(sample(c(1:nControls), 
                                      n.random.controls))
          if (is.null(input.yaml$ped) == F) {
            ped <- read.table(paste0(input.yaml$ped), 
                              stringsAsFactors = F, header = F, fill = T)
            columnIndex_fam <- which(colnames(cohort.object.auto[["countmat"]]) %in% 
                                       (ped %>% filter(V1 %in% (ped %>% filter(V2 %in% 
                                                                                 sample_name) %>% pull(1))) %>% pull(2)))
            columnIndex_n = which(colnames(cohort.object.auto[["countmat"]]) %in% 
                                    sample_name)
            while (length(which(columnIndex_fam %in% 
                                rand.samples == T)) > 0) {
              rand.samples <- sort(sample(c(1:nControls), 
                                          n.random.controls))
            }
            reference_list_auto <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                              columnIndex_n], reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                                                                                  rand.samples], correlations = cohort.object.auto[["correlations"]][rand.samples, 
                                                                                                                                                                                                                                                     columnIndex_n])
            nControls <- dim(cohort.object.x[["countmat"]])[2]
            rand.samples <- sort(sample(c(1:nControls), 
                                        n.random.controls))
            columnIndex_fam <- which(colnames(cohort.object.x[["countmat"]]) %in% 
                                       (ped %>% filter(V1 %in% (ped %>% filter(V2 %in% 
                                                                                 sample_name) %>% pull(1))) %>% pull(2)))
            columnIndex_x <- which(colnames(cohort.object.x[["countmat"]]) %in% 
                                     sample_name)
            while (length(which(columnIndex_fam %in% 
                                rand.samples == T)) > 0) {
              rand.samples <- sort(sample(c(1:nControls), 
                                          n.random.controls))
            }
            reference_list_x <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                        columnIndex_x], reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                                                                                         rand.samples], correlations = cohort.object.x[["correlations"]][rand.samples, 
                                                                                                                                                                                                                                         columnIndex_x])
          } else {
            columnIndex_n = which(colnames(cohort.object.auto[["countmat"]]) %in% 
                                    sample_name)
            while (columnIndex_n %in% rand.samples == 
                   T) {
              rand.samples <- sort(sample(c(1:nControls), 
                                          n.random.controls))
            }
            reference_list_auto <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                              columnIndex_n], reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                                                                                  rand.samples], correlations = cohort.object.auto[["correlations"]][rand.samples, 
                                                                                                                                                                                                                                                     columnIndex_n])
            nControls <- dim(cohort.object.x[["countmat"]])[2]
            rand.samples <- sort(sample(c(1:nControls), 
                                        n.random.controls))
            while (columnIndex_x %in% rand.samples == 
                   T) {
              rand.samples <- sort(sample(c(1:nControls), 
                                          n.random.controls))
            }
            reference_list_x <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                        columnIndex_x], reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]], 
                                                                                                                                                                         rand.samples], correlations = cohort.object.x[["correlations"]][rand.samples, 
                                                                                                                                                                                                                                         columnIndex_x])
          }
        }
        if (is.null(reference_list_auto$reference.choice) == 
            F) {
          reference_set_auto = rowSums(cohort.object.auto[["countmat"]][, 
                                                                        reference_list_auto$reference.choice, drop= F])
          all_exons_auto = new("ExomeDepth", test = cohort.object.auto[["countmat"]][, 
                                                                                     columnIndex_n], reference = reference_set_auto, 
                               formula = "cbind(test,reference) ~ 1", subset.for.speed = 10000)
          all_exons_auto = ExomeDepth::CallCNVs(x = all_exons_auto, 
                                                transition.probability = as.numeric(input.yaml$transition.probability), 
                                                chromosome = cohort.object.auto[["annotations"]]$chromosome, 
                                                start = cohort.object.auto[["annotations"]]$start, 
                                                end = cohort.object.auto[["annotations"]]$end, 
                                                name = cohort.object.auto[["annotations"]]$name)
        } else {
          message(paste0("\n *********** Sample ", sample_name, 
                         " failed for Autosome variant calling. ************* \n"))
        }
        if (is.null(reference_list_x$reference.choice) == 
            F) {
          reference_set_x = rowSums(cohort.object.x[["countmat"]][, 
                                                                  reference_list_x$reference.choice, drop =F])
      
          all_exons_x = new("ExomeDepth", test = cohort.object.x[["countmat"]][, 
                                                                               columnIndex_x], reference = reference_set_x, 
                            formula = "cbind(test,reference) ~ 1")
          all_exons_x = ExomeDepth::CallCNVs(x = all_exons_x, 
                                             transition.probability = as.numeric(input.yaml$transition.probability), 
                                             chromosome = cohort.object.x[["annotations"]]$chromosome, 
                                             start = cohort.object.x[["annotations"]]$start, 
                                             end = cohort.object.x[["annotations"]]$end, 
                                             name = cohort.object.x[["annotations"]]$name)
        } else {
          message(paste0("\n *********** Sample ", sample_name, 
                         " failed for chr X variant calling. ************* \n"))
        }
        
        all_exons_auto@CNV.calls$sample <- sample_name
        if (dim(all_exons_x@CNV.calls)[1] > 0) {
          all_exons_x@CNV.calls$sample <- sample_name
          my.calls.final <- rbind(all_exons_auto@CNV.calls, 
                                  all_exons_x@CNV.calls)
        } else {
          my.calls.final <- all_exons_auto@CNV.calls
        }
        
        my.calls.final.del.gr <- suppressWarnings( GenomicRanges::GRanges(my.calls.final[my.calls.final$type %in% 
                                                                         c("deletion"), ]$id) )
        my.calls.final.dup.gr <- suppressWarnings( GenomicRanges::GRanges(my.calls.final[my.calls.final$type %in% 
                                                                         c("duplication"), ]$id) )
        sample.del$reproducibility <- sample.del$reproducibility + 
          suppressWarnings(  as.numeric(GenomicRanges::countOverlaps(sample.del.gr, 
                                                  my.calls.final.del.gr) > 0) )
        sample.dup$reproducibility <- sample.dup$reproducibility + 
          suppressWarnings( as.numeric(GenomicRanges::countOverlaps(sample.dup.gr, 
                                                  my.calls.final.dup.gr) > 0) )
      }
      sample.cnv.calls <- rbind(sample.del, sample.dup)
    } else {
      sample.cnv.calls <- calls.annotated
    }
    
    # annotated <- EDM::variant.annotations(sample.cnv.calls)
    # calls_all <- cbind(sample.cnv.calls, annotated[, 4:dim(annotated)[2]])
    write.table(sample.cnv.calls, paste0(input.yaml$output.directory, 
                                  "/results/individual.edm.calls/", sample_name, ".edm.calls.txt"), 
                row.names = F, sep = "\t")
    saveRDS(all_exons, paste0(input.yaml$output.directory, 
                              "/results/individual.ed.objects/", sample_name, ".ed.object.rds"))
    references$value[3] <- cor(all_exons$auto_calls@test, 
                               all_exons$auto_calls@reference)
    references$ID <- sample_name
    write.table(references, paste0(input.yaml$output.directory, 
                                 "/results/individual.edm.calls/", sample_name, ".edm.stat.txt"), 
                row.names = F, sep = "\t")
    return(all_exons)
  } else {
    if (dim(all_exons_auto@CNV.calls)[1] > 0) {
      all_exons_auto@CNV.calls$sample <- sample_name
      calls.first <- rbind(all_exons_auto@CNV.calls)
      all_exons <- list(auto_calls = all_exons_auto)
      annotated <- EDM::variant.annotations(calls.first)
      calls.annotated <- cbind(calls.first, annotated[, 4:dim(annotated)[2]])
      
      if (input.yaml$reproducibility == T) {
        calls.first$reproducibility <- 0
        sample.del <- calls.first[calls.first$type %in% 
                                    c("deletion"), ]
        sample.dup <- calls.first[calls.first$type %in% 
                                    c("duplication"), ]
        sample.del.gr <- GenomicRanges::GRanges(sample.del$id)
        sample.dup.gr <- GenomicRanges::GRanges(sample.dup$id)
        
        if (is.null(input.yaml$iterations) == F) {
          iterations = as.numeric(input.yaml$iterations)
        } else {
          iterations <- 1000
        }
        
        for (iter in 1:iterations) {
          n.random.controls <- 100
          if (input.yaml$precomputed.controls) {
            rand.samples <- sort(sample(c(1:nControls), 
                                        n.random.controls))
            reference_list_auto <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                              columnIndex_n], reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                                                                                  rand.samples], correlations = cohort.object.auto[["correlations"]][rand.samples, 
                                                                                                                                                                                                                                                     columnIndex_n])
          } else {
            nControls <- dim(cohort.object.auto[["countmat"]])[2]
            rand.samples <- sort(sample(c(1:nControls), 
                                        n.random.controls))
            if (is.null(input.yaml$ped) == F) {
              ped <- read.table(paste0(input.yaml$ped), 
                                stringsAsFactors = F, header = F, fill = T)
              columnIndex_fam <- which(colnames(cohort.object.auto[["countmat"]]) %in% 
                                         (ped %>% filter(V1 %in% (ped %>% filter(V2 %in% 
                                                                                   sample_name) %>% pull(1))) %>% pull(2)))
              columnIndex_n = which(colnames(cohort.object.auto[["countmat"]]) %in% 
                                      sample_name)
              while (length(which(columnIndex_fam %in% 
                                  rand.samples == T)) > 0) {
                rand.samples <- sort(sample(c(1:nControls), 
                                            n.random.controls))
              }
              reference_list_auto <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                columnIndex_n], reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                                                                                    rand.samples], correlations = cohort.object.auto[["correlations"]][rand.samples, 
                                                                                                                                                                                                                                                       columnIndex_n])
            } else {
              columnIndex_n = which(colnames(cohort.object.auto[["countmat"]]) %in% 
                                      sample_name)
              while (columnIndex_n %in% rand.samples == 
                     T) {
                rand.samples <- sort(sample(c(1:nControls), 
                                            n.random.controls))
              }
              reference_list_auto <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                columnIndex_n], reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]], 
                                                                                                                                                                                    rand.samples], correlations = cohort.object.auto[["correlations"]][rand.samples, 
                                                                                                                                                                                                                                                       columnIndex_n])
            }
          }
          if (is.null(reference_list_auto$reference.choice) == 
              F) {
            reference_set_auto = rowSums(cohort.object.auto[["countmat"]][, 
                                                                          reference_list_auto$reference.choice, drop= F])
            all_exons_auto = new("ExomeDepth", test = cohort.object.auto[["countmat"]][, 
                                                                                       columnIndex_n], reference = reference_set_auto, 
                                 formula = "cbind(test,reference) ~ 1", 
                                 subset.for.speed = 10000)
            all_exons_auto = ExomeDepth::CallCNVs(x = all_exons_auto, 
                                                  transition.probability = as.numeric(input.yaml$transition.probability), 
                                                  chromosome = cohort.object.auto[["annotations"]]$chromosome, 
                                                  start = cohort.object.auto[["annotations"]]$start, 
                                                  end = cohort.object.auto[["annotations"]]$end, 
                                                  name = cohort.object.auto[["annotations"]]$name)
          } else {
            message(paste0("\n *********** Sample ", 
                           sample_name, " failed for Autosome variant calling. ************* \n"))
          }
          
          all_exons_auto@CNV.calls$sample <- sample_name
          my.calls.final <- all_exons_auto@CNV.calls
          my.calls.final.del.gr <- suppressWarnings( GenomicRanges::GRanges(my.calls.final[my.calls.final$type %in% 
                                                                           c("deletion"), ]$id) )
          my.calls.final.dup.gr <- suppressWarnings(  GenomicRanges::GRanges(my.calls.final[my.calls.final$type %in% 
                                                                           c("duplication"), ]$id) )
          sample.del$reproducibility <- sample.del$reproducibility + 
            as.numeric(suppressWarnings( GenomicRanges::countOverlaps(sample.del.gr, 
                                                    my.calls.final.del.gr) ) > 0)
          sample.dup$reproducibility <- sample.dup$reproducibility + 
            as.numeric(suppressWarnings( GenomicRanges::countOverlaps(sample.dup.gr, 
                                                    my.calls.final.dup.gr) ) > 0)
        }
        sample.cnv.calls <- rbind(sample.del, sample.dup)
      } else {
        sample.cnv.calls <- calls.annotated
      }
      
      write.table(sample.cnv.calls, paste0(input.yaml$output.directory, 
                                    "/results/individual.edm.calls/", sample_name, 
                                    ".edm.calls.txt"), row.names = F, sep = "\t")
      saveRDS(all_exons, paste0(input.yaml$output.directory, 
                                "/results/individual.ed.objects/", sample_name, 
                                ".ed.object.rds"))
      references$value[3] <- cor(all_exons$auto_calls@test, 
                                 all_exons$auto_calls@reference)
      references$ID <- sample_name
      write.table(references, paste0(input.yaml$output.directory, 
                                   "/results/individual.edm.calls/", sample_name, 
                                   ".edm.stat.txt"), row.names = F, sep = "\t")
      return(all_exons)
    } else {
      return(NULL)
    }
  }
}
