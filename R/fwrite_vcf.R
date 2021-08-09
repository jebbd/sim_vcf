#' Write a vVCF file from a formatted data.table
#' The formatting should be done before invoking this call
#' bcftools should be installed for this command 
#' 
#'@param infile to input VCF use in Sim
#'@param outfile The name of the output path to write
#'@param dt The data.table with formatted data 
#'@export
fwrite_vcf<-function(infile,outfile,dt){
  hfile<-tempfile()
  system(paste0("source ~/.bash_profile;bcftools view -h ",infile,"> ",hfile),ignore.stderr = TRUE)
  header<-readLines(hfile)
  header<-purrr::discard(header,~{stringr::str_detect(.x,"FORMAT")&!stringr::str_detect(.x,"ID=GT")})
  hdr[length(hdr)+1]<-"##source=SimulatedGenotypes"
  header<-header[c(1,length(hdr),2:(length(hdr)-1))]
  fwrite(as.list(header),outfile,quote=FALSE)
  fwrite(dt,file=outfile,sep = "\t",append = TRUE,col.names = TRUE)
}
