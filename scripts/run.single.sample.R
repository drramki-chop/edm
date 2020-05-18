#!/usr/bin/Rscript

## Calls CNV for a single sample
## Uses mosdepth coverage output and uses control cohort created by EDM to call variants

## List of inputs
## sample.coverage = mosdepth output (.bed.gz)
## sample.id = sample name
## sex = M / F
## exon.definition = exon bed file
## control.autosomes.rds = control cohort created by EDM


## List of inputs from config

## controls.autosomes
## controls.sex.chromosome.male
## controls.sex.chromosome.female

## Hard-coded inputs
## Parameters such as bin size and transition.probability can be hard-coded


suppressMessages(library(optparse))
suppressMessages(library(EDM))

options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("--sample-coverage"), action="store", 
              type='character', help="coverage computed by mosdepth"),
  make_option(c("--sex"), action="store", 
              type='character', help="sample sex"),  
  make_option(c("--sample-id"),action="store",
  		        type='character',help="sample id"),
  make_option(c("--control-autosomes"), action="store",
		type='character', help = "Pre-computed control cohort (autosomes)"),
  make_option(c("--control-sex-chromosome"), action = "store",
  		type='character', help = "Pre-computed control cohort (sex chromosome)"),
  make_option(c("--exclude-family-from-controls"), action = "store", default = NULL,
		type = 'character', help = "exclude familyID from controls"),
  make_option(c("--exclude-sample-from-controls"), action = "store", default = NULL,
		type = 'character', help = "exclud sample from the control cohort")
)

opt = parse_args(OptionParser(option_list=option_list))

message(paste("\n\n Provided sex is ",opt$sex,"and the control cohort provided is",basename(opt$`control-sex-chromosome`),"\n"))

controls.autosomes = opt$`control-autosomes`
controls.sex.chromosome = opt$`control-sex-chromosome`

cohort.object.auto <- readRDS(controls.autosomes)
cohort.object.x <- readRDS(controls.sex.chromosome)

auto.manifest <- cohort.object.auto$manifest
chrX.manifest <- cohort.object.x$manifest

if(!is.null(opt$`exclude-family-from-controls`)){
	exclude.family <- opt$`exclude-family`;
    
    	auto.exclude.family <- auto.manifest[auto.manifest$familyID %in% exclude.family, ]
    	chrX.exclude.family <- chrX.manifest[chrX.manifest$familyID %in% exclude.family, ]

    if(NROW(auto.exclude.family) > 0){
	message(paste0("Excluding familyID ",exclude.family," (",NROW(auto.exclude.family)," samples)"," from the autosome control cohort\n"))
	countmat <- cohort.object.auto[["countmat"]]
	countmat <- countmat[,-which(colnames(countmat) %in% auto.exclude.family$sampleID)]
    	cohort.object.auto[["countmat"]] <- countmat
    }else{
	stop(paste0("FamilyID ",exclude.family," not found in controls (autosomes)."))
    }

    if(NROW(chrX.exclude.family) > 0){
	message(paste0("Excluding familyID ",exclude.family," from the chrX control cohort\n"))
	countmat <- cohort.object.x[["countmat"]]
	countmat <- countmat[,-which(colnames(countmat) %in% chrX.exclude.family$sampleID)]
    	cohort.object.x[["countmat"]] <- countmat	
    }else{
	stop(paste0("FamilyID ",exclude.family," not found in controls (chrX)."))
    }
}

if(!is.null(opt$`exclude-sample-from-controls`)){
	exclude.sample <- opt$`exclude-sample`;
    
	auto.exclude.sample <- auto.manifest[auto.manifest$sampleID %in% exclude.sample, ]
	chrX.exclude.sample <- chrX.manifest[chrX.manifest$sampleID %in% exclude.sample, ]

    if(NROW(auto.exclude.sample) > 0){
	message(paste0("Excluding sampleID ",exclude.sample," from the autosome control cohort\n"))
	countmat <- cohort.object.auto[["countmat"]]
	countmat <- countmat[,-which(colnames(countmat) %in% auto.exclude.sample$sampleID)]
	cohort.object.auto[["countmat"]] <- countmat
    }else{
	stop(paste0("SampleID ",exclude.sample," not found in controls (autosomes)."))
    }
    
    if(NROW(chrX.exclude.sample) > 0){
	message(paste0("Excluding sampleID ",exclude.sample," from the chrX control cohort\n"))
	countmat <- cohort.object.x[["countmat"]]
	countmat <- countmat[,-which(colnames(countmat) %in% chrX.exclude.sample$sampleID)]
        cohort.object.x[["countmat"]] <- countmat
    }else{
	stop(paste0("SampleID ",exclude.sample," not found in controls (chrX)."))
    }
}

nControls.autosomes <- dim(cohort.object.auto[["countmat"]])[2]
nControls.sex.chromosome <- dim(cohort.object.x[["countmat"]])[2]

sample.coverage <- as.data.frame(data.table::fread(opt$`sample-coverage`))

sample.auto.coverage <- sample.coverage[sample.coverage$V1 %in% c(1:22),]
sample.sex.coverage <- sample.coverage[sample.coverage$V1 %in% c("X"),]

coverage.nCols <- dim(sample.coverage)[2]; 

### Insert the test sample to the control cohort
cohort.object.auto$countmat <- cbind(cohort.object.auto$countmat,as.matrix(round(sample.auto.coverage[,coverage.nCols,drop=FALSE])))
columnIndex = nControls.autosomes + 1

### re-compute correlation
cohort.object.auto$selected.exons <- EDM::select.exons(cohort.counts = cohort.object.auto$countmat,
                                                      bin.length = (cohort.object.auto$annotations$end-cohort.object.auto$annotations$start)/1000,
                                                      n.bins.reduced = 10000)

countmat.selected <- cohort.object.auto$countmat[cohort.object.auto$selected.exons,]
cohort.object.auto$correlations <- cor(countmat.selected)

reference_list <- EDM::select.reference.panel(test.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]],columnIndex],
                                                reference.counts = cohort.object.auto[["countmat"]][cohort.object.auto[["selected.exons"]],-columnIndex],
                                                correlations = cohort.object.auto[["correlations"]][-columnIndex,columnIndex])
if(is.null(reference_list$reference.choice) == F){
	
    #reference_set = apply(
    #                    X = as.matrix(cohort.object.auto[["countmat"]][, reference_list$reference.choice]),
    #                    MAR=1, FUN=sum)
	
    reference_set = rowSums(as.matrix(cohort.object.auto[["countmat"]][, reference_list$reference.choice]))
    all_exons_auto = new('ExomeDepth', test=cohort.object.auto[["countmat"]][,columnIndex],
                         reference=reference_set,
                         formula = 'cbind(test,reference) ~ 1')
    all_exons_auto = ExomeDepth::CallCNVs(x = all_exons_auto, transition.probability=1e-4,
                                          chromosome= cohort.object.auto[["annotations"]]$chromosome, start=cohort.object.auto[["annotations"]]$start,
                                          end=cohort.object.auto[["annotations"]]$end, name=cohort.object.auto[["annotations"]]$name)
    if(dim(all_exons_auto@CNV.calls)[1] > 0){
        all_exons_auto@CNV.calls$sample <- opt$`sample-id`
        saveRDS(all_exons_auto,paste0(opt$`sample-id`,".auto.edObject.rds"));
        CNV.calls <- all_exons_auto@CNV.calls
	
	## Sample stats
	sample.stats <- data.frame(matrix(nrow=1,ncol= 5))
	names(sample.stats) <- c("sample","n.deletions","n.duplications","total.cnvs","model.phi")
	sample.stats$sample[1] <- opt$`sample-id`;
	sample.stats$n.deletions[1] <- NROW(CNV.calls[CNV.calls$type %in% "deletion",])
	sample.stats$n.duplications[1] <- NROW(CNV.calls[CNV.calls$type %in% "duplication",])
	sample.stats$total.cnvs[1] <- NROW(CNV.calls);
	sample.stats$model.phi[1] <- unique(all_exons_auto@phi);
	write.csv(sample.stats,paste0(opt$`sample-id`,".sample.stats.csv"),quote=F,row.names=F)

    } else {
	message("Warning: Did not find any CNV in autosomes.")
	flush.console();
    }
}

############### Sex chromosome #####################3

sample.sex.coverage <- sample.coverage[sample.coverage$V1 %in% c("X"),]
columnIndex = nControls.sex.chromosome + 1
cohort.object.x$countmat <- cbind(cohort.object.x$countmat,as.matrix(round(sample.sex.coverage[,coverage.nCols,drop=FALSE])))
cohort.object.x$selected.exons <- EDM::select.exons(cohort.counts = cohort.object.x$countmat,
							bin.length = (cohort.object.x$annotations$end - cohort.object.x$annotations$start)/1000,
							n.bins.reduced=10000)
countmat.selected = cohort.object.x$countmat[cohort.object.x$selected.exons,]
cohort.object.x$correlations <- cor(countmat.selected)

reference_list <- EDM::select.reference.panel(test.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]],columnIndex],
                                                reference.counts = cohort.object.x[["countmat"]][cohort.object.x[["selected.exons"]],-columnIndex],
                                                correlations = cohort.object.x[["correlations"]][-columnIndex,columnIndex])
if(is.null(reference_list$reference.choice) == F){
    reference_set = apply(
                        X = as.matrix(cohort.object.x[["countmat"]][, reference_list$reference.choice]),
                        MAR=1, FUN=sum)
    all_exons_x = new('ExomeDepth', test=cohort.object.x[["countmat"]][,columnIndex],
                         reference=reference_set,
                         formula = 'cbind(test,reference) ~ 1')
    all_exons_x = ExomeDepth::CallCNVs(x = all_exons_x, transition.probability=1e-4,
                                          chromosome= cohort.object.x[["annotations"]]$chromosome, start=cohort.object.x[["annotations"]]$start,
                                          end=cohort.object.x[["annotations"]]$end, name=cohort.object.x[["annotations"]]$name)
    if(dim(all_exons_x@CNV.calls)[1] > 0){
	all_exons_x@CNV.calls$sample <- opt$`sample-id`
    	saveRDS(all_exons_x,paste0(opt$`sample-id`,".sex.edObject.rds"));
    	CNV.calls <- rbind(CNV.calls,all_exons_x@CNV.calls)
    }else{
	    message("Did not find any CNVs in ChrX");
	    flush.console();
    }
}

if(dim(CNV.calls)[1] > 0){
    write.table(CNV.calls,paste0(opt$`sample-id`,".edm.calls.txt"),sep="\t",row.names=F,quote=F)
} else {
    message("Warning: Didn't find any CNV in this sample.")
}
