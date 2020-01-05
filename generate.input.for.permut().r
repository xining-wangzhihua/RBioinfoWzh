#last edited at 20200105

library(magrittr)
library(ape)
library(seqinr)

generate.input.for.permut=function(distribution.tsv,snp.fasta){
  #>>>internal functions begin>>>
  #<<<internal functions end<<<
  #>>>control the input of distribution.tsv and snp.fasta begin>>>
  if(!is.vector(distribution.tsv)|is.list(distribution.tsv)|length(distribution.tsv)!=1){
    stop("distribution.tsv must be a string.")
  }
  if(!is.vector(snp.fasta)|is.list(snp.fasta)|length(snp.fasta)!=1){
    stop("snp.fasta must be a string.")
  }
  if(!file.exists(distribution.tsv)){stop("files specified by distribution.tsv doesn't exist.")}
  if(!file.exists(snp.fasta)){stop("files specified by snp.fasta doesn't exist.")}
  #<<<control the input of distribution.tsv and snp.fasta end<<<
  #>>>import datum in distribution.tsv and snp.fasta begin>>>
  distribution.tsv=read.table(file=distribution.tsv,fileEncoding="UTF-8-BOM",sep="\t",
                              header=TRUE,stringsAsFactors=FALSE,row.names=1) %>% as.matrix()
  snp.fasta=ape::read.FASTA(file=snp.fasta,type="DNA") %>% as.matrix() %>% as.character()
  if(ncol(distribution.tsv)!=nrow(snp.fasta)){
    stop("column number in distribution.tsv must equal number of sequences in snp.fasta")
  }
  if(any(colnames(distribution.tsv)!=rownames(snp.fasta))){
    stop("haplotype/genotype names in distribution.tsv don't match with those in snp.fasta")
  }
  if(anyNA(distribution.tsv)){stop("items in distribution.tsv can't be NA")}
  if(any(distribution.tsv=="")){stop("items in distribution.tsv can't be \"\"")}
  if(nrow(snp.fasta)==0){stop("snp.fasta must contain at least one sequence")}
  if(ncol(snp.fasta)==0){stop("sequences in snp.fasta can't be empty")}
  #<<<import datum in distribution.tsv and snp.fasta end<<<
  #>>>represent residues in sequences with 0~10 begin>>>
  for(i in 1:nrow(snp.fasta)){
    #use seqinr::s2n(), because chartr() can only recognize single character but not "10" as ten.
    snp.fasta[i,] %<>% seqinr::s2n(seq=.,levels=s2c("-atgcrymwsk"),base4=TRUE,forceToLower=TRUE) %>% as.vector()
  }
  #<<<represent residues in sequences with 0~10 end<<<
  #>>>export datum in distribution and snp.fasta begin>>>
  output_file=strftime(Sys.time(),format="input for permut, %Y%m%d_%H%M%S.txt")
  if(file.exists(output_file)){stop("can't create the file: \"",output_file,"\".\t it already exists.")}
  tulip=file(description=output_file,open="wt")
  cat("\n#meanings of the following 3 integers: genotype/haplotype number, ",
      "population number, segregating site number\n",file=tulip,append=FALSE,sep="")
  cat(ncol(distribution.tsv),"\t",nrow(distribution.tsv),"\t",ncol(snp.fasta),"\n\n#row and column of ",
      "the following matrix represent population multiply genotype/haplotype\n",file=tulip,append=TRUE,sep="")
  write.table(x=distribution.tsv,file=tulip,append=TRUE,sep="\t",
              quote=FALSE,row.names=FALSE,col.names=FALSE)
  cat("\n#row and column of the following matrix represent genotype/haplotype multiply segregating site\n",
      file=tulip,append=TRUE,sep="")
  write.table(x=snp.fasta,file=tulip,append=TRUE,sep="\t",
              quote=FALSE,row.names=FALSE,col.names=FALSE)
  close(tulip);rm(tulip)
  message("\nresult is written in: \"",output_file,"\".\n")
  #<<<export datum in gcsv and gxml end<<<
  rm(i,distribution.tsv,snp.fasta)
  invisible(output_file)#rm(output_file)
}

if(FALSE){
  #here are the test codes
  distribution.tsv="~/../OneDrive/Sep13_atrata/analysis, dnasp/atrata-its, genotype, distribution in populations.tsv"
  snp.fasta="~/../OneDrive/Sep13_atrata/analysis, dnasp/atrata-its, genotype, snp sequence.fasta"
  generate.input.for.permut(distribution.tsv,snp.fasta)
  rm(distribution.tsv,snp.fasta,generate.input.for.permut)
}
