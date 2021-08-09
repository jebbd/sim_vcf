#' Read a VCF into a data.table
#' 
#' @param path Path to the VCF file to read
#' @param ... any other options which can be passed to data.table::fread
#' @export
fread_vcf<-function(path,...){
  data.table::fread(path,skip="#CHROM",
                    blank.lines.skip = TRUE,
                    keepLeadingZeros = TRUE,
                    ...)
}

#' Check for AF allele frequency in INFO field
#' Checks 1% of the supplied data.tale/frame for an INFO entry of AF=\emph{f} 
#' @param dt A data.table/frame or tibble from a VCF (see fread_vcf)
#' @export
check_allele_freq_info<-function(dt){
  all(stringr::str_detect(dt$INFO[sample(1:nrow(dt),floor(nrow(dt)*0.01))],"\\bAF=[\\d\\.]+"))
}

#' Check for AC and AN (allele count and number) in INFO field
#' Checks 1% of the supplied data.tale/frame for an INFO entry of AC=\emph{C} and AN=\emph{N} 
#' @param dt A data.table/frame or tibble from a VCF (see fread_vcf)
#' @export
check_allele_cn_info<-function(dt){
  rows<-sample(1:nrow(dt),floor(nrow(dt)*0.01))
  all(c(
    stringr::str_detect(dt$INFO[rows],"\\bAC=[\\d\\.]+"),
    stringr::str_detect(dt$INFO[rows],"\\bAC=[\\d\\.]+")
  ))
}

#' Extract or calculate Allele Frequency from a VCF
#' If the AF field is not detected, the AC (allele count) and AN (allele number) fields are checked.
#' As a last resort the genotypes can be used
#' The data.table is modified in place?? Check
#' @param dt a data.table from fread_vcf
#' @export
extract_maf<-function(dt){
  if(check_allele_freq_info(dt)){
    dt[,"AF":=lapply(.SD,function(x){as.double((stringr::str_match(x,"\\bAF=([\\d\\.]+)"))[,2])}),.SDcols=c("INFO")]
  }
  else if(check_allele_cn_info(dt)){
    dt[,"AF":=lapply(.SD,function(x){
      as.double((stringr::str_match(x,"\\bAC=([\\d\\.]+)"))[,2])/as.double((stringr::str_match(x,"\\bAN=([\\d\\.]+)"))[,2])
      }),.SDcols=c("INFO")]
  }
  else{
    stop("No allele frequencies found or calculated")
  }
}
