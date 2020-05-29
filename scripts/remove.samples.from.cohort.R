suppressMessages(library(optparse))
suppressMessages(library(EDM))

options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("--cohort-autosomes"), action="store",
                type='character', help="autosomes"),
  make_option(c("--sex-cohort-men"),action = "store",
                type = "character", help = "Cohort of male samples"),
  make_option(c("--sex-cohort-women"), action = "store",
                type = "character", help = "Cohort of female samples"),
  make_option(c("--exclude-samples"), action = "store",
                type = "character", help = "manifest with samples to exclude"),
  make_option(c("--output-cohort-name"), action = "store",
                type = "character", help = "Name of the new cohort")
)

opt = parse_args(OptionParser(option_list=option_list))

expected.columns <- c("familyID","sampleID","paternalID","maternalID","sex","phenotype","mosdepth.coverage","batchID","sample.type")

exclude.manifest <- read.table(opt$`exclude-samples`, header = T, sep = "\t", stringsAsFactors = F)

cohort.object.auto <- readRDS(opt$`cohort-autosomes`);
cohort.object.men<- readRDS(opt$`sex-cohort-men`)
cohort.object.women <- readRDS(opt$`sex-cohort-women`)


### Exclude samples
manifest.autosomes <- cohort.object.auto$manifest;

if(NROW(manifest.autosomes[manifest.autosomes$sampleID %in% exclude.manifest$sampleID,]) > 0){
	cohort.object.auto$countmat <- cohort.object.auto$countmat[,-which(colnames(cohort.object.auto$countmat) %in% exclude.manifest$sampleID)]

	cohort.object.auto$selected.exons <- EDM::select.exons(cohort.counts = cohort.object.auto$countmat,
                                                      bin.length = (cohort.object.auto$annotations$end-cohort.object.auto$annotations$start)/1000,
                                                      n.bins.reduced = 10000)
	countmat.selected <- cohort.object.auto$countmat[cohort.object.auto$selected.exons,]
	cohort.object.auto$correlations <- cor(countmat.selected)
	cohort.object.auto$manifest <- cohort.object.auto$manifest[cohort.object.auto$manifest$sampleID %in% colnames(cohort.object.auto$countmat),]
	saveRDS(cohort.object.auto,paste0(opt$`output-cohort-name`,".auto.rds"))
} else{
	stop("I did not find any sample to exclude in the autosomal cohort. Stopping here.");
}

if(NROW(exclude.manifest[exclude.manifest$sex %in% c("male"),]) > 0){
    cohort.object.men$countmat <- cohort.object.men$countmat[,-which(colnames(cohort.object.men$countmat) %in% exclude.manifest$sampleID)]

    cohort.object.men$selected.exons <- EDM::select.exons(cohort.counts = cohort.object.men$countmat,
                                                        bin.length = (cohort.object.men$annotations$end-cohort.object.men$annotations$start)/1000,
                                                        n.bins.reduced = 10000)
    countmat.selected <- cohort.object.men$countmat[cohort.object.men$selected.exons,]
    cohort.object.men$correlations <- cor(countmat.selected)
    cohort.object.men$manifest <- cohort.object.men$manifest[cohort.object.men$manifest$sampleID %in% colnames(cohort.object.men$countmat),]
    saveRDS(cohort.object.men,paste0(opt$`output-cohort-name`,".men.rds"))
}else{
    message(paste0("There were no male samples in the manifest to exclude."))
}

if(NROW(exclude.manifest[exclude.manifest$sex %in% c("female"),]) > 0){
    cohort.object.women$countmat <- cohort.object.women$countmat[,-which(colnames(cohort.object.women$countmat) %in% exclude.manifest$sampleID)]

    cohort.object.women$selected.exons <- EDM::select.exons(cohort.counts = cohort.object.women$countmat,
                                                        bin.length = (cohort.object.women$annotations$end-cohort.object.women$annotations$start)/1000,
                                                        n.bins.reduced = 10000)
    countmat.selected <- cohort.object.women$countmat[cohort.object.women$selected.exons,]
    cohort.object.women$correlations <- cor(countmat.selected)
    cohort.object.women$manifest <- cohort.object.women$manifest[cohort.object.women$manifest$sampleID %in% colnames(cohort.object.women$countmat),]
    saveRDS(cohort.object.women,paste0(opt$`output-cohort-name`,".women.rds"))
}else{
    message(paste0("There were no male samples in the manifest to exclude."))
}

