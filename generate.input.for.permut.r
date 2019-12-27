#last edited at 20191227

generate.input.for.permut=function(genotype_xml,genotype_csv){
  #>>>r package dependency, function dependency, internal functions begin>>>
  library(magrittr)
  library(xml2)
  library(stringr)
  library(Hmisc)
  library(tibble)
  library(ape)
  library(seqinr)
  source(file="https://github.com/ywd5/r-zm/raw/master/site.polymorphism().r",local=TRUE)
  f_all_equal=function(x){unname(x) %>% unique() %>% {length(.)==1} %>% return();}
  #<<<r package dependency, function dependency, internal functions end<<<
  #>>>control the input of genotype_xml and genotype_csv begin>>>
  if(!is.vector(genotype_xml)|is.list(genotype_xml)){stop("input of genotype_xml must be a string.")}
  if(!is.vector(genotype_csv)|is.list(genotype_csv)){stop("input of genotype_csv must be a string.")}
  if(!file.exists(genotype_xml)){stop("files specified by genotype_xml doesn't exist.")}
  if(!file.exists(genotype_csv)){stop("files specified by genotype_csv doesn't exist.")}
  #<<<control the input of genotype_xml and genotype_csv end<<<
  #>>>import datum in genotype_xml begin>>>
  gxml=xml2::read_xml(x=genotype_xml) %>% xml2::as_list() %>% .[[1]]
  if(length(gxml)==0){stop("genotype_xml can't be empty.")}
  if(!f_all_equal(x=lengths(gxml))){stop("each branch of genotype_xml must has equal length.")}
  if(!f_all_equal(x=names(gxml))){stop("each branch of genotype_xml must has same name.")}
  for(i in 1:length(gxml)){gxml[[i]]=unlist(gxml[[i]]);}
  ans=gxml[[1]]
  if(length(gxml)>1)for(i in 2:length(gxml)){ans=rbind(ans,gxml[[i]]);}
  gxml=as_tibble(ans)[c("genotype","sequence")]
  rm(ans)
  gxml$genotype=str_squish(gxml$genotype)
  gxml$sequence=str_squish(gxml$sequence)
  if(anyNA(gxml$genotype)|anyNA(gxml$sequence)){stop("elements in genotype_xml can't be NA.")}
  if(any(gxml$genotype=="")|any(gxml$sequence=="")){stop("elements in genotype_xml can't be \"\".")}
  rm(genotype_xml)
  #<<<import datum in genotype_xml end<<<
  #>>>import datum in genotype_csv begin>>>
  gcsv=read.csv(genotype_csv,header=TRUE,row.names=1) %>% as.matrix()
  if(nrow(gxml)!=ncol(gcsv)){stop("genotye numbers don't match between genotype_xml and genotype_csv")}
  if(any(gxml$genotype!=colnames(gcsv))){stop("genotype names don't match between genotype_xml and genotype_csv")}
  if(anyNA(gcsv)){stop("items in genotype_csv can't be NA.")}
  if(any(gcsv=="")){stop("items in genotype_csv can't be \"\".")}
  rm(genotype_csv)
  #<<<import datum in genotype_csv end<<<
  #>>>main manipulation begin>>>
  gxml=setNames(object=strsplit(gxml$sequence,split=""),nm=gxml$genotype) %>% ape::as.DNAbin() %>% as.matrix()
  tulip=tibble(polymorphism=site.polymorphism(gxml),count=0,besides_acgt=FALSE) %>% rowid_to_column(var="index")
  tulip$count=lengths(tulip$polymorphism)
  for(i in nrow(tulip)){tulip$besides_acgt[i]=any(tulip$polymorphism[[i]] %nin% s2c("acgt"));}
  if(any(tulip$count==1 & tulip$besides_acgt)){
    ans=tulip$index[tulip$count==1 & tulip$besides_acgt]
    ans=tibble(index=ans,residue=tulip$polymorphism[ans])
    if(length(ans$residue)==0){ans$residue=character();}#this won't occur, just for testing
    ans=paste0(ans$index,"th residue:\t",ans$residue,"\n")
    warning("residues at the following sites are all same but not a,c,g,t:\n",ans)
    rm(ans)
  }
  gxml=gxml[,tulip$count!=1] %>% as.character()
  for(i in 1:nrow(gxml)){
    gxml[i,]=s2n(seq=gxml[i,],levels=s2c("-atgcrymwsk"),base4=TRUE,forceToLower=TRUE)
  }
  rm(tulip)
  #<<<main manipulation end<<<
  #>>>export datum in gcsv and gxml begin>>>
  output_file=strftime(Sys.time(),format="genotype input for permut, %Y%m%d_%H%M%S.txt")
  if(file.exists(output_file)){stop("can't create the file: \"",output_file,"\".\t it already exists.")}
  cat("\n#meanings of the following 3 integers: genotypes/haplotypes number, ",
      "population number, segregating site number\n",file=output_file,append=FALSE,sep="")
  cat(ncol(gcsv),"\t",nrow(gcsv),"\t",ncol(gxml),"\n\n#row and column of the following matrix ",
      "represent population and genotype/haplotype\n",file=output_file,append=TRUE,sep="")
  write.table(x=gcsv,file=output_file,append=TRUE,sep="\t",
              quote=FALSE,row.names=FALSE,col.names=FALSE)
  cat("\n#row and column of the following matrix represent genotype/haplotype and segregating site.\n",
      file=output_file,append=TRUE,sep="")
  write.table(x=gxml,file=output_file,append=TRUE,sep="\t",
              quote=FALSE,row.names=FALSE,col.names=FALSE)
  message("\ninput for permut is written in: \"",output_file,"\".\n")
  #<<<export datum in gcsv and gxml end<<<
  rm(gcsv,gxml,f_all_equal,site.polymorphism,i)
  invisible(output_file)#rm(output_file)
}

if(FALSE){
  #here are the test codes
  genotype_xml="~/../OneDrive/Sep13_atrata/analysis, dnasp/atrata-its, genotype from nexus.xml"
  genotype_csv="~/../OneDrive/Sep13_atrata/analysis, dnasp/atrata-its, genotype distribution in population.csv"
  generate.input.for.permut(genotype_xml,genotype_csv)
  rm(genotype_xml,genotype_csv,generate.input.for.permut)
}
