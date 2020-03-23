#' @importFrom yaml  read_yaml
#' @export
gather.variant.calls <- function(task_id,input){
input.yaml <- yaml::read_yaml(input)

cnvFiles <- list.files(paste0(input.yaml$output.directory,"/results/individual.edm.calls/"),pattern = "*.edm.calls.txt", full.names=T)

# calls <- purrr::map_dfc(cnvFiles, function(x) {data.table::fread(x)}) 

calls <- do.call(rbind,lapply(cnvFiles, function(fn) read.table(file.path(paste0(fn)),header=T,stringsAsFactors=F , sep ="\t")))


summary.cohort <- as.data.frame(table(calls$sample))

dels <- GenomicRanges::GRanges(calls[calls$type %in% "deletion",]$id)
dups <- GenomicRanges::GRanges(calls[calls$type %in% "duplication",]$id)

del.calls <- calls[calls$type %in% "deletion",]
dups.calls <- calls[calls$type %in% "duplication",]

del.calls$dels.exact.match <- GenomicRanges::countOverlaps(dels,dels, type="equal")
del.calls$dels.any.match <- GenomicRanges::countOverlaps(dels,dels)
del.calls$dups.exact.match <- NA
del.calls$dups.any.match <- NA

dups.calls$dels.exact.match <- NA
dups.calls$dels.any.match <- NA
dups.calls$dups.exact.match <- GenomicRanges::countOverlaps(dups,dups, type="equal")
dups.calls$dups.any.match <- GenomicRanges::countOverlaps(dups,dups)

all_calls <- rbind(del.calls, dups.calls)


write.table(all_calls,paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,".edm.cohort.calls.txt"),row.names=F,sep="\t",quote=F)
write.table(summary.cohort,paste0(input.yaml$output.directory,"/results/",input.yaml$cohort.name,".edm.summary.cohort.txt"),row.names=F,sep="\t",quote=F)
}
