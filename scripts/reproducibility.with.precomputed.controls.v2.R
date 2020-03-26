## reproducibility.w.precomputed.controls.R 

## This function takes an exomedepth model (object), control cohort, and number of required replications.
## It bootstraps random combinations of controls and uses the estimated parameters (phi) from the original model
## to call CNVs to create an estimate of reproducibility.

## Notes on this version
## ---------------------

## Uses an random controls

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


############# Defaults ##############################
transition.probability = 1e-4
expected.CNV.length = 50000
transitions <- matrix(
  nrow = 3,
  ncol = 3,
  c(
    1. - transition.probability,
    transition.probability / 2.,
    transition.probability / 2.,
    0.5,
    0.5,
    0.,
    0.5,
    0,
    0.5
  ),
  byrow = TRUE
)

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
  expected = rep(mean(exomedepth.object@test/reference_set,na.rm=T),length(exomedepth.object@test))
  
  loglike <- .Call(
    "get_loglike_matrix",
    phi = exomedepth.object@phi,
    expected = expected,
    total = as.integer(reference_set + exomedepth.object@test),
    observed = as.integer(exomedepth.object@test),
    mixture = 1
  )
  
  shift <- 0
  my.calls.final <- data.frame()
  for (chrom in 1:22) {
    good.pos <-
      which (exomedepth.object@annotations$chromosome == chrom)
    loc.annotations <- exomedepth.object@annotations[good.pos , ]
    loc.expected <- exomedepth.object@expected[good.pos]
    loc.test <- exomedepth.object@test[good.pos]
    loc.total <-
      exomedepth.object@reference[good.pos] + exomedepth.object@test[good.pos]
    positions <- loc.annotations$start
    end.positions <- loc.annotations$end
    
    loc.likelihood <-
      rbind(c(-Inf, 0, -Inf), loglike[good.pos, c(2, 1, 3)], c(-100, 0,-100))
    
    my.calls <-
      viterbi.hmm (
        transitions,
        loglikelihood = loc.likelihood,
        positions = as.integer(
          c(
            positions[1] - 2 * expected.CNV.length,
            positions,
            end.positions[length(end.positions)] + 2 * expected.CNV.length
          )
        ),
        #include position of new dummy exon
        expected.CNV.length = expected.CNV.length
      )
    my.calls$calls$start.p <- my.calls$calls$start.p -1  ##remove the dummy exon, which has now served its purpose
    my.calls$calls$end.p <- my.calls$calls$end.p -1  ##remove the dummy exon, which has now served its purpose
    #loc.likelihood <- loc.likelihood[ -1, c(2,1, 3) ]  ##remove the dummy exon, which has now served its purpose
    loc.likelihood <- loc.likelihood[ -c(1,nrow(loc.likelihood)), c(2,1, 3), drop = FALSE ] ##remove both of the dummy exons, which have now served their purpose
    
    ################################ Now make it look better, add relevant info
    if (nrow(my.calls$calls) > 0) {
      
      my.calls$calls$start <- loc.annotations$start[ my.calls$calls$start.p ]
      my.calls$calls$end <- loc.annotations$end[ my.calls$calls$end.p ]
      my.calls$calls$chromosome <- as.character(loc.annotations$chromosome[ my.calls$calls$start.p ])
      
      my.calls$calls$id <- paste('chr', my.calls$calls$chromosome, ':',  my.calls$calls$start, '-',  my.calls$calls$end, sep = '')
      my.calls$calls$type <- c('deletion', 'duplication')[ my.calls$calls$type ]
      
      ########## make things pretty
      my.calls$calls$BF <- NA
      my.calls$calls$reads.expected <- NA
      my.calls$calls$reads.observed <- NA
      
      
      for (ir in 1:nrow(my.calls$calls)) {
        
        if (my.calls$calls$type[ir] == 'duplication') my.calls$calls$BF[ir] <-  sum(loc.likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir],3 ] - loc.likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir],2 ])
        
        if (my.calls$calls$type[ir] == 'deletion') my.calls$calls$BF[ir] <-  sum(loc.likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir], 1 ] - loc.likelihood [ my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir],2  ])
        
        my.calls$calls$reads.expected[ ir ] <-  sum( loc.total [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ] * loc.expected [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ ir ] ])
        my.calls$calls$reads.observed[ ir ] <-  sum( loc.test [my.calls$calls$start.p[ir] : my.calls$calls$end.p[ir] ] )
      }
      
      my.calls$calls$reads.expected <- as.integer( my.calls$calls$reads.expected)
      my.calls$calls$reads.ratio <-  signif(my.calls$calls$reads.observed / my.calls$calls$reads.expected, 3)
      my.calls$calls$BF <- signif( log10(exp(1))*my.calls$calls$BF, 3)
      
      #### shift the numbering properly
      my.calls$calls$start.p <- my.calls$calls$start.p + shift
      my.calls$calls$end.p <- my.calls$calls$end.p + shift
      
      if (nrow(my.calls.final) == 0) {my.calls.final <- my.calls$calls} else {my.calls.final <- rbind.data.frame(my.calls.final, my.calls$calls)}
      #message('Number of calls for chromosome ', chrom, ' : ', nrow(my.calls$calls))
    }
    shift <- shift + length(good.pos)
  }
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
