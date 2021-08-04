#' Auxillary functions
#'
#' @description get_pos: get a random postion on a chromosome
#' @param chrom A character vector of chromsomes formatted as e.g. "chr1"
#' @rdname auxillary_functions
#'
#' @export
get_pos<-function(chrom){
  map_int(chrom,~{
    lim <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)[.x]
    return(as.integer(round(runif(1L,1L,lim),digits=0)))
  })}

#' @description mut_ref: mutate a reference base semi randomly, with a Ti/Tv of roughly 2, needs to wrapped in a map/lapply for working with vectors
#' @param base an single DNA nucleotide (ATCG) to mutate
#'
#' @rdname auxillary_functions
#' @export
mut_ref<-function(base){
  pyr<-c("C","T")
  pur<-c("A","G")
  if(base=="N"){return(sample(c(pyr,pur),1))}
  flip<-runif(1)
  if(flip>(1/3)){
    if(base %in% pyr){
      return(pyr[!pyr %in% base])
    }else{
      return(pur[!pur %in% base])
    }
  }else{
    if(base %in% pyr){
      return(sample(pur,1))
    }else{
      return(sample(pyr,1))
    }
  }
}
#' @description rev_encode: reverse the encoding of genotypes from numbers to that use in VCF i.e. 0/0, 0/1, 1/1
#' @param x a vector of numbers from 0 - 2 to encode. ANy other numbers will become no calls ("./.")
#' @rdname auxillary_functions
#' @export
rev_encode<-function(x){
  dplyr::case_when(
    x==0~"0/0",
    x==1~"0/1",
    x==2~"1/1",
    TRUE ~ "./."
  )
}
