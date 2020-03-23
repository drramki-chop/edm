#' @export
variant.annotations <- function(calls){

get.positional.annotation <- function(gr){
  # gr <- GRanges(gr);
  
  start.gr <- GRanges(seqnames(gr),IRanges(start(gr)-5,start(gr)+5))
  start.hits <- findOverlaps(start.gr,exons.hg19.forAnn.gr)
  end.gr <- GRanges(seqnames(gr),IRanges(end(gr)-5,end(gr)+5))
  end.hits <- findOverlaps(end.gr,exons.hg19.forAnn.gr) ##########|||||||| fix the hardcoded object names ||||||
  
  start.hits <- exons.hg19.forAnn[subjectHits(start.hits),]
  end.hits <- exons.hg19.forAnn[subjectHits(end.hits),]
  
  start.hits <- start.hits[,c("gene","exonPosition","exonNumber","strand")]
  end.hits <- end.hits[,c("gene","exonPosition","exonNumber","strand")]
  names(start.hits)[2] <- "BP1.exonPosition"
  names(end.hits)[2] <- "BP2.exonPosition"
  names(start.hits)[3] <- "BP1.exonNumber"
  names(end.hits)[3] <- "BP2.exonNumber"
  
  genes <- merge(start.hits,end.hits,all.x=T,all.y=T)
  # genes <- merge(genes,genes.strand)
  genes <- merge(genes,exonCounts)
  # genes$positionalInfo <- apply(genes,1,get.positional.info)
  
  gene_return = paste0(genes$gene,collapse = "|")
  # gene_position = paste0(genes$positionalInfo, collapse = "|")
  # return(paste0(gene_return,";",gene_position))
  return(gene_return)
}

all_anno = readRDS(paste0(find.package("EDM"),"/data/pli_omim_rvis.rds"))
names(all_anno) = c("genes","rvis_percentile","Mim.Number","Phenotypes","mis_z","pLI" )



exonCounts = readRDS(paste0(find.package("EDM"),"/data/exonCounts_by_gene.rds"))
exons.hg19.forAnn = readRDS(paste0(find.package("EDM"),"/data/exons.hg19.annotations.rds"))
exons.hg19.forAnn.gr = GRanges(exons.hg19.forAnn$interval)

clinvar <- readRDS(paste0(find.package("EDM"),"/data/clinvar_110718.rds")) 
clingenHI <- readRDS(paste0(find.package("EDM"),"/data/clingenHI_110718.rds")) 
clingenTS <- readRDS(paste0(find.package("EDM"),"/data/clingenTS_110718.rds")) 

query = calls

query <- GRanges(query$id);
queryDF <- as.data.frame(query)[,1:3]

pos_ann <- lapply(seq_along(query), function(X) get.positional.annotation(query[X])) %>% do.call(rbind,.) %>% as.data.frame()
names(pos_ann) = "genes"
queryDF <- cbind(queryDF,pos_ann)
queryDF$genes <- as.character(queryDF$genes)
queryDF$clinvar <- sapply(queryDF$genes,
                          FUN = function(x){
                            foo <- unique(clinvar[clinvar$GeneSymbol %in% as.character(unlist(strsplit(x,"\\|"))),])
                            if(dim(foo)[1]  < 1){return(NA)}
                            paste(foo$DiseaseName,collapse=" | ")
                          })


# clingenHI[clingenHI$HaploinsufficiencyScore > 3,]$HaploinsufficiencyScore <- NA;
queryDF$ClingenHaploInsufficiency <- sapply(queryDF$genes,
                                            FUN=function(x){
                                              foo <- unique(clingenHI[clingenHI$geneSymbol %in% as.character(unlist(strsplit(x,"\\|"))),]);
                                              if(dim(foo)[1] < 1) { return(NA)};
                                              paste(foo$HaploinsufficiencyScore,collapse=" | ")
                                            })

# clingenTS[clingenTS$TriplosensitivityScore > 3,]$TriplosensitivityScore <- NA;
queryDF$ClingenTriplosensitivity <- sapply(queryDF$genes,
                                           FUN=function(x){
                                             foo <- unique(clingenTS[clingenTS$geneSymbol %in% as.character(unlist(strsplit(x,"\\|"))),]);
                                             if(dim(foo)[1] < 1) { return(NA)};
                                             paste(foo$TriplosensitivityScore,collapse=" | ")
                                           })
# queryDF[is.na(queryDF)] <- ".";

queryDF %>% left_join(all_anno, by = "genes") -> queryDF 

return(queryDF)                      

}


