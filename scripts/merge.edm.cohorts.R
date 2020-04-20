## Merge control cohorts
##

suppressMessages(library(optparse))

option_list = list(
  make_option(c("--existing-cohort"), action="store",
                type='character', help="Existing EDM cohort"),
  make_option(c("--new-cohort"),action = "store",
                type = "character", help = "New EDM cohort"),
  make_option(c("--merged-cohort-name"), action = "store",
                type = "character", help = "Name of the merged cohort"),
  make_option(c("--output-folder"),action = "store",
                type = "character", help = "output directory"),
  make_option(c("--exclude-parents", action= "store",default=TRUE,
                type="logical", help = "exclude parental samples"))
)

opt = parse_args(OptionParser(option_list=option_list))

existing.cohort <- readRDS(opt$`existing-cohort`)
new.cohort <- readRDS(opt$`new-cohort`)
merged.cohort.name <- opt$`merged-cohort-name`
output.directory <- as.character(opt$`output-folder`)

## Check annotations and count matrices

if(all.equal(existing.cohort$annotations,new.cohort$annotations) == FALSE){
  stop("Cohort annotations do not match. Stopping here.")
}

if(NROW(existing.cohort$countmat) != NROW(new.cohort$countmat)){
  stop("Count matrices do not match between cohorts. Stopping here")
}

## check column names
if(length(unique(c(colnames(existing.cohort$countmat), colnames(new.cohort$countmat)))) <
  length(c(colnames(existing.cohort$countmat), colnames(new.cohort$countmat)))){
      stop("Looks like there are duplicate sample names")
}

if(all.equal(sort(colnames(existing.cohort$manifest)),sort(colnames(new.cohort$manifest))) == F){
  stop("Manifest format does not match between the cohorts. Stopping here.")
}

if(NCOL(existing.cohort$countmat) != NROW(existing.cohort$manifest)){
  stop(paste("Manifest does not match the count matrix in the existing cohort"))
}
if(NCOL(new.cohort$countmat) != NROW(new.cohort$manifest)){
  stop(paste("Manifest does not match the count matrix in the new cohort"))
}

countmat <- cbind(existing.cohort$countmat,new.cohort$countmat)
merged.manifest <- rbind(existing.cohort$manifest, new.cohort$manifest)
if(opt$`exclude-parents`){
  merged.manifest <- merged.manifest[merged.manifest$sample.type %in% c("parent","father","mother"),]
  countmat <- countmat[,colnames(countmat) %in% merged.manifest$sampleID]
}

annotations <- existing.cohort$annotations
selected.exons <- EDM::select.exons(cohort.counts = countmat,
                                        bin.length = (annotations$end - annotations$start)/1000,
                                        n.bins.reduced = 10000)

countmat.selected <- countmat[selected.exons,]
countmat.cor <- cor(countmat.selected)
countData <- list(countmat = countmat, annotations = annotations,
                  selected.exons = selected.exons, correlations = countmat.cor, manifest = merged.manifest)

saveRDS(countData, file = paste0(output.directory,"/",merged.cohort.name, ".exomedepth.cohort.rds"))
