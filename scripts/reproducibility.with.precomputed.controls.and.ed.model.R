## reproducibility.w.precomputed.controls.R 

## This function takes an exomedepth model (object), control cohort, and number of required replications.
## It bootstraps random combinations of controls and uses the estimated parameters (phi) from the original model
## to call CNVs to create an estimate of reproducibility.

## Notes on this version
## ---------------------

## Uses random controls with default ExomeDepth model in each iteration. Time intensive. Set --iteartions 10 for tesing purposes.

# Rscript reproducibility.with.precomputed.controls.default.R \
# --exomedepth-object CWES-0215-P-A-DGD-16-92.ExomeDepthObject.rds \
# --control-cohort cwes.exomedepth.cohort.auto.rds --iterations 10 --n-random-controls 200

suppressMessages(library(optparse))
suppressMessages(library(EDM))

options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("--exomedepth-object"), action="store", 
              type='character', help="Single sample exomedpeth object"),
  make_option(c("--control-cohort"), action="store", 
              type='character', help="Control cohort created by EDM (.rds)"),  
  make_option(c("--iterations"),action="store",default=1000,
              type='numeric',help="Number of iterations to run"),
  make_option(c("--n-random-controls",action="store",default=100),
              type="numeric",help = "Number of controls to choose")
)

start_time <- Sys.time()

opt = parse_args(OptionParser(option_list=option_list))
cohort.object <- readRDS(opt$`control-cohort`)
exomedepth.object <- readRDS(opt$`exomedepth-object`)
iterations <- opt$iterations
n.random.controls <- opt$`n-random-controls`

countmat <- cohort.object[["countmat"]]
n.controls <- dim(countmat)[2]
sample.name <- unique(exomedepth.object@CNV.calls$sample)
countmat <- cbind(countmat,exomedepth.object@test)
sample.index = n.controls + 1
colnames(countmat)[sample.index] <- sample.name

countmat.selected <- countmat[cohort.object$selected.exons,]

### Something to keep in mind for large cohorts
correlations <- cor(countmat)


############# Reproducibility ########################
## CNV calling code re-used from EoxmeDepth

sample.cnv.calls <- exomedepth.object@CNV.calls
sample.cnv.calls$reproducibility <- 0
sample.del <- sample.cnv.calls[sample.cnv.calls$type %in% c("deletion"),]
sample.dup <- sample.cnv.calls[sample.cnv.calls$type %in% c("duplication"),]
sample.del.gr <- GenomicRanges::GRanges(sample.del$id)
sample.dup.gr <- GenomicRanges::GRanges(sample.dup$id)

for (iter in 1:iterations) {
  message(paste("iteration",iter))
  flush.console()
  rand.samples <- sort(sample(c(1:n.controls), n.random.controls))
 
  reference_list <-
    EDM::select.reference.panel(
      test.counts = exomedepth.object@test[cohort.object$selected.exons],
      reference.counts = cohort.object[["countmat"]][cohort.object[["selected.exons"]], rand.samples],
      correlations = correlations[rand.samples, sample.index]
    )


  reference_set = rowSums(countmat[, reference_list$reference.choice])
  if(NROW(countmat) > 25000){
	 subset.exons.for.speed = 10000;
  }else{
	subset.exons.for.speed = NROW(countmat)
  }

  all_exons_auto = new('ExomeDepth', test=countmat[,sample.index],
                       reference=reference_set,
                       formula = 'cbind(test,reference) ~ 1', subset.for.speed = subset.exons.for.speed)
  all_exons_auto = ExomeDepth::CallCNVs(x = all_exons_auto, transition.probability=1e-4,
                                        chromosome= cohort.object[["annotations"]]$chromosome, start=cohort.object[["annotations"]]$start,
                                        end=cohort.object[["annotations"]]$end, name=cohort.object[["annotations"]]$name)
  all_exons_auto@CNV.calls$sample <- sample.name
  

  my.calls.final <- all_exons_auto@CNV.calls
  
  #write.table(my.calls.final,paste0(sample.name,".iter.",iter,".txt"), row.names=F,quote=F,sep="\t")
  my.calls.final.del.gr <- GenomicRanges::GRanges(my.calls.final[my.calls.final$type %in% c("deletion"),]$id)
  my.calls.final.dup.gr <- GenomicRanges::GRanges(my.calls.final[my.calls.final$type %in% c("duplication"),]$id)
  sample.del$reproducibility <- sample.del$reproducibility +  as.numeric(GenomicRanges::countOverlaps(sample.del.gr,my.calls.final.del.gr) > 0)
  sample.dup$reproducibility <- sample.dup$reproducibility +  as.numeric(GenomicRanges::countOverlaps(sample.dup.gr,my.calls.final.dup.gr) > 0)
}

sample.cnv.calls <- rbind(sample.del,sample.dup)
write.table(sample.cnv.calls,paste0(sample.name,".edm.reproducibility.calls.txt"),quote=F,row.names=F,sep="\t")
total_time = Sys.time() - start_time

message(paste("Took ",total_time))
