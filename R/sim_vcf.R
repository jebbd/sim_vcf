#' Simulate Variants in Hardy Weinberg for the human genome (hg38)
#'
#' @param num_vars The number of variants to simulate, default: 500
#' @param num_ppl The number of individuals to simulate, default: 5000
#' @param path The output file path, default: sim.vcf.gz
#'
#' @import data.table
#' @export
simulate_vcf<-function(num_vars=500L,num_ppl=5000L,path="sim.vcf.gz"){

  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

  header<-paste("##fileformat=VCFv4.2",
  paste0("##fileDate=",lubridate::today()),
          "##source=mySimpleSims",
          '##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">',
          sep="\n")

  dt<-data.table(chrom=paste0("chr",sort(round(runif(num_vars,1L,22L),digits=0))))

  header<-paste(header,
                paste(
                      paste0("##contig=<ID=",
                             names(seqlengths(genome)[1:22]),
                             ",length=",
                             GenomeInfoDb::seqlengths(genome)[1:22],
                             ",assembly=hg38,species=\"Homo sapiens\">"),
                      collapse="\n"),
                sep="\n")

  dt[,POS:=get_pos(dt$chrom)]

  dt[,maf:=map_dbl(1:nrow(dt),
                  ~{
                    maf=rnorm(1,mean=0.2,sd=0.075)
                    if(maf<0){maf<-maf*-1}
                    if(maf>0.5){maf<-1-maf}
                    return(maf)}
                  )
     ]

  dt[,`:=`(
          ID=paste0(dt$chrom,":",dt$POS),
          REF=BiocGenerics::as.vector(getSeq(genome,dt$chrom,start=dt$POS,end=dt$POS))
         )
   ]

  dt[,`:=`(
         ALT=map_chr(dt$REF,mut_ref),
         QUAL=".",FILETR="PASS",INFO=".",FORMAT="GT"
         )
     ]

  ppl<-paste0("i_",sprintf("%05d",1:num_ppl))


  hw_freqs<-purrr::map(dt$maf,
                       ~{
                          t(rmultinom(1,length(ppl),
                                  c((1-.x)^2,2*.x*(1-.x),.x^2)))
                        }
                      )%>%
    {do.call(rbind,.)}

  genotypes<-do.call(rbind,
        purrr::map(1:nrow(hw_freqs),
                   ~{
                      sample(
                        c(
                          rep(0,hw_freqs[.x,1]),
                          rep(1,hw_freqs[.x,2]),
                          rep(2,hw_freqs[.x,3])
                          )
                        )
                    }
                  )
  )

  colnames(genotypes)<-ppl

  dt<-cbind(dt,as.data.table(genotypes))
  setnames(dt,c("chrom"),c("#CHROM"))
  dt[,maf:=NULL]
  dt[,paste0(ppl):=lapply(.SD,rev_encode),.SDcols=ppl]

  fwrite(strsplit(header,"\n"),path,quote=FALSE)
  fwrite(dt,file=path,sep = "\t",append = TRUE,col.names = TRUE)
}



