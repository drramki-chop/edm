## gather coverage files mosdepth to create a cohort for use in EDM
## This is a custom script for a specific application at Children's Hospital of Philadelphia.
## This can be edited for your requirements.

suppressMessages(library(optparse))
suppressMessages(library(EDM))

options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("--manifest"), action="store",
                type='character', help="coverage computed by mosdepth"),
  make_option(c("--output-folder"),action = "store",
                type = "character", help = "output directory"),
  make_option(c("--cohort-name"), action = "store",
                type = "character", help = "Name of the cohort")
)

opt = parse_args(OptionParser(option_list=option_list))

expected.columns <- c("familyID","sampleID","paternalID","maternalID","sex","phenotype","mosdepth.coverage","batchID","sample")

manifest <- read.table(opt$manifest, header = T, sep = "\t", stringsAsFactors = F)
output.directory <- as.character(opt$`output-folder`)

### Check for permissions
system(paste0("mkdir -p ",output.directory))

cohort.name = opt$`cohort-name`

if(all.equal(sort(expected.columns),sort(names(manifest)))){
  countFiles <- as.character(manifest$mosdepth.coverage)

  ## Check file sizes
  fileSizes <- unlist(lapply(countFiles, file.size))
  failedSamples <- countFiles[which(fileSizes == 0)]
  if(length(failedSamples) > 0){
    write.table(as.data.frame(failedSamples),paste0(output.directory,"/",cohort.name,".zero.coverage.files.txt"))
    stop("Looks like some of the mosdepth output has file size zero. Please check log.")
  }

  snames <- manifest$sampleID
  firstSample <- as.data.frame(data.table::fread(countFiles[1]))
  nCols <- dim(firstSample)[2]
  names(firstSample)[nCols] <- snames[1]
  snames <- snames[-1]

  annotations <- firstSample[, 1:3]
  names(annotations)[1:3] <- c("chromosome", "start", "end")
  annotations$name <- paste0(annotations$chromosome, ":", annotations$start, "-", annotations$end)
  countFiles <- countFiles[-1]
  firstSample <- round(firstSample[, nCols, drop = F])
  counts <- purrr::map_dfc(countFiles,
                          function(x) {
                            round(data.table::fread(x, drop = 1:(nCols - 1)))
                          })
  names(counts) <- snames
  countmat <- as.matrix(cbind(firstSample, counts))
  selected.exons <- EDM::select.exons(cohort.counts = countmat,
      bin.length = (annotations$end - annotations$start)/1000,
      n.bins.reduced = 10000)
  countmat.selected <- countmat[selected.exons, ]
  countmat.cor <- cor(countmat.selected)
  countData <- list(countmat = countmat, annotations = annotations,
                    selected.exons = selected.exons, correlations = countmat.cor)

  men <- manifest$sampleID[manifest$sex %in% c("male")]
  women <- manifest$sampleID[manifest$sex %in% c("female")]

  countmat.auto <- countmat[which(countData[[2]]$chromosome %in% c(1:22)), ]
  annotations.auto <- annotations[which(countData[[2]]$chromosome %in% c(1:22)), ]
  selected.exons.auto <- EDM::select.exons(cohort.counts = countmat.auto,
                                          bin.length = (annotations.auto$end - annotations.auto$start)/1000,
                                          n.bins.reduced = 10000)
  countmat.selected.auto <- countmat[selected.exons.auto,]
  countmat.cor.auto <- cor(countmat.selected.auto)
  countData.auto <- list(countmat = countmat.auto, annotations = annotations.auto,
                        selected.exons = selected.exons.auto, correlations = countmat.cor.auto, manifest = manifest)
  saveRDS(countData.auto, file = paste0(output.directory,"/",cohort.name, ".exomedepth.cohort.auto.rds"))



  ### Sex chromosome / ChrX
  if (length(men) > 0) {
      manifest.men <- manifest[manifest$sex %in% "male",]
      countmat.men <- countmat[which(countData[[2]]$chromosome %in% "X"),
                              colnames(countmat) %in% men]
      annotations.men <- annotations[which(countData[[2]]$chromosome %in% "X"), ]
      selected.exons.men <- EDM::select.exons(cohort.counts = countmat.men,
                                              bin.length = (annotations.men$end - annotations.men$start)/1000,
                                              n.bins.reduced = 10000)
      countmat.selected.men <- countmat.men[selected.exons.men,]
      countmat.cor.men <- cor(countmat.selected.men)
      countData.men <- list(countmat = countmat.men, annotations = annotations.men,
                            selected.exons = selected.exons.men, correlations = countmat.cor.men, manifest = manifest.men)
      saveRDS(countData.men, file = paste0(output.directory,"/",cohort.name, ".exomedepth.cohort.men.rds"))
  }
  if (length(women) > 0) {
      manifest.women <- manifest[manifest$sex %in% "female",]
      countmat.women <- countmat[which(countData[[2]]$chromosome %in% "X"),
                                colnames(countmat) %in% women]
      annotations.women <- annotations[which(countData[[2]]$chromosome %in% "X"), ]
      selected.exons.women <- EDM::select.exons(cohort.counts = countmat.women,
                                                bin.length = (annotations.women$end - annotations.women$start)/1000,
                                                n.bins.reduced = 10000)
      countmat.selected.women <- countmat.women[selected.exons.women,]
      countmat.cor.women <- cor(countmat.selected.women)
      countData.women <- list(countmat = countmat.women,
                              annotations = annotations.women, selected.exons = selected.exons.women,
                              correlations = countmat.cor.women, manifest = manifest.women)
      saveRDS(countData.women, file = paste0(output.directory,"/",cohort.name, ".exomedepth.cohort.women.rds"))
  }
}else{
  stop("Looks like there is an issue with column names in the manifest. Please check documentation.")
}
