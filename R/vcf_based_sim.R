#' Simulate genotypes from a vector of minor allele frequencies
#' All genotypes are in HWE
#' @param p vector of minor(?) allele frequencies
#' @returns A matrix of genotypes encoded for a VCF file (0/0,0/1,1/1)
#' @export
simulate_geno<-function(p,num_ppl=500){
  hw_freqs<-purrr::map(p,
                       ~{
                         t(rmultinom(1,num_ppl,
                                     c((1-.x)^2,2*.x*(1-.x),.x^2)))
                       }
  )%>%
    {do.call(rbind,.)}
  
  genotypes<-do.call(rbind,
                     purrr::map(1:nrow(hw_freqs),
                                ~{
                                  sample(
                                    c(
                                      rep("0/0",hw_freqs[.x,1]),
                                      rep("0/1",hw_freqs[.x,2]),
                                      rep("1/1",hw_freqs[.x,3])
                                    )
                                  )
                                }
                     )
  )
  
  return(genotypes)
}
#' Create a fuly populated data.table
#' 
#' @param dt The data.table
#' @param num_ppl The number of individuals to simulate
#' @export
simulate_geno_dt<-function(dt,num_ppl,maf_cut_off=0){
  qt<-copy(dt)
  extract_maf(qt)
  qt[,"FORMAT":="GT"]
  genotypes<-simulate_geno(dt$AF,num_ppl=num_ppl)
  
  colnames(genotypes)<-paste0("i_",sprintf(paste0("%0",ceiling(log10(num_ppl))+1,"d"),1:num_ppl))
  qt<-cbind(qt,as.data.table(genotypes))
  qt<-qt[AF>maf_cut_off][,AF:=NULL]
  return(qt)
}
