#last edited at 20191212

#a function to calculate the frequencies of segregating sites in sequences

require(package="Biostrings")
require(package="tibble")
require(package="magrittr")
require(package="stringr")
require(package="questionr")
require(package="XML")
source(file="https://github.com/ywd5/r-zm/raw/master/paste.xml.elements().r")
#meta-data, retrench, exaltation, affinage, succinct, terse,brief, prune,curtail, truncate, abridge

pick.rare.sequences=function(xss,residue_index_start_from=1){
  risf=residue_index_start_from
  rm(residue_index_start_from)
  #>>>internal functions begin>>>
  f_paste_byrow=function(df){
    #for internal usage only
    l=dim(df)
    if(!identical(length(l),2L)){stop("invalid input")}
    l=l[1]
    zm=rep("",times=l)
    for(i in 1:l){zm[i]=paste0(df[i,],collapse="");}
    return(zm)#rm(df,l,i,zm)
  }
  f_create_directory=function(prefix="residue index start from ?"){
    #for internal usage only
    if(!is.vector(prefix)|is.list(prefix)){stop("invalid input")}
    if(length(prefix)!=1){stop("length of \"prefix\" must be 1.")}
    prefix=as.character(prefix)
    if(grepl(pattern="[\a\b\f\n\r\t\v\"\\\\*?|<>:/]",x=prefix)){
      stop("\"prefix\" contain characters which are invalid for file pathes.")
    }
    if(nchar(prefix)>50){stop("\"prefix\" contain too many characters.")}
    #
    zm=Sys.time() %>% as.character() %>% strsplit(split="[-: ]") %>% .[[1]]
    if(length(zm)!=6)stop("bug")
    zm=paste0(zm[1],zm[2],zm[3],"-",zm[4],"h",zm[5],"m",zm[6],"s")
    zm=paste0(prefix,", ",zm)
    #
    if(dir.exists(paths=zm)){stop("can't create the directory. it already exist.")}
    dir.create(path=zm)
    return(zm)#rm(prefix,zm)
  }
  #<<<internal functions end<<<
  #>>>control the input of xss begin>>>
  if(is.null(xss)){stop("\"xss\" can't be NULL.")}
  if(!is(object=xss,class2="XStringSet")){stop("\"xss\" must be an XStringSet.")}
  if(length(xss)<3){stop("\"xss\" must contain at least 3 sequences.")}
  if(nchar(xss) %>% unique() %>% length()!=1){
    stop("sequences in \"xss\" must have equal length, i.e. aligned.")
  }
  if(nchar(xss[[1]])==0){stop("sequences in \"xss\" must contain at least 1 residue.")}
  #
  ids=names(xss)
  if(is.null(ids)){ids=rep("not_set",times=length(xss));}
  ids[is.na(ids)]="not_set"
  ids=str_remove_all(string=ids,pattern="[\a\b\f\n\r\t\v\"\\\\*?|<>:/]")
  ids[grep(pattern="^ *$",x=ids)]="not_set"
  ids=make.unique(names=ids,sep="_")
  names(xss)=ids
  rm(ids)
  #<<<control the input of xss end<<<
  #>>>control the input of residue_index_start_from begin>>>
  if(is.null(risf)){stop("\"residue_index_start_from\" can't be NULL.");}
  if(!is.vector(risf)|is.list(risf)){stop("\"residue_index_start_from\" must be an integer.");}
  if(length(risf)!=1){stop("length of \"residue_index_start_from\" must be 1");}
  risf=as.numeric(risf)
  if(is.na(risf)|is.infinite(risf)){stop("\"residue_index_start_from\" must be an integer.");}
  if(risf!=as.integer(risf)){
    warning("\"residue_index_start_from\" is converted to an integer via as.integer()")
  }
  risf=as.integer(risf)
  #<<<control the input of residue_index_start_from end<<<
  #>>>generate the datum structure of whole manipulation begin>>>
  l=c(0,0,0,0,5)
  names(l)=c("sequence number","unique sequence number",
             "residue number","interested residue number","output matrix row number")
  l[1]=length(xss)
  l[3]=nchar(xss[[1]])
  residue_entry=rep(FALSE,times=l[3])
  info=tibble(id=names(xss),
              sequence_hash="",
              weight=0,
              small=tibble(numbers=list(0),sites=list(0),min=0,min_site=0),
              big=tibble(numbers=list(0),sites=list(0),max=0,max_site=0))
  residue=tibble(interested=as.matrix(xss) %>% as_tibble(.name_repair="minimal"),
                 complementary_hash=tibble("",.name_repair="minimal"))
  datum=tibble(frequency=tibble(""),
               other_frequencies=tibble(""),
               complementary_frequency=tibble(""))
  #<<<generate the datum structure of whole manipulation end<<<
  #>>>calculate xss$interested_entry and xss$interested begin>>>
  for(i in 1:l[3])if(residue$interested[[i]] %>% unique() %>% length()!=1){
    residue_entry[i]=TRUE
  }
  if(!any(residue_entry)){stop("all sequences in \"xss\" are same.")}
  residue$interested=residue$interested[residue_entry]
  l[4]=ncol(residue$interested)
  #<<<calculate xss$interested_entry and xss$interested end<<<
  #>>>calculate xss$sequence_hash begin>>>
  info$sequence_hash=f_paste_byrow(df=residue$interested) %>% match(x=.,table=unique(.))
  info$sequence_hash=paste0("x",info$sequence_hash)
  #<<<calculate xss$sequence_hash end<<<
  #>>>merge number of same sequences into weight begin>>>
  msu=split(x=1:l[1],f=info$sequence_hash) %>% unname()
  spbu=split(x=info$id,f=info$sequence_hash) %>% unname()
  l[2]=length(msu)
  for(i in 1:l[2]){
    msu[[i]]=msu[[i]][1]
    spbu[[i]]=str_sort(spbu[[i]],numeric=TRUE)
  }
  msu=unlist(msu)
  xss=xss[msu]
  residue=residue[msu,]
  info=info[msu,]
  info$id=spbu
  rm(msu,spbu)
  info$weight=lengths(info$id)
  #<<<merge number of same sequences into weight end<<<
  #>>>calculate xss$complementary_hash begin>>>
  if(l[4]>1)for(i in 1:l[2])for(j in 1:l[4]){
    residue$complementary_hash[i,j]=paste0(residue$interested[i,-j],collapse="")
  }
  if(l[4]==1){residue$complementary_hash[[1]]="";}
  #
  for(i in 1:l[4]){
    ans=residue$complementary_hash[[i]]
    ans=match(x=ans,table=unique(ans))
    ans=paste0("x",ans)
    residue$complementary_hash[[i]]=ans
    rm(ans)
  }
  #<<<calculate xss$complementary_hash end<<<
  #>>>calculate the 3 frequencies begin>>>
  blank_tibble=matrix(data="",nrow=l[2],ncol=l[4]) %>% as_tibble(.name_repair="minimal")
  datum=tibble(frequency=blank_tibble,other_frequencies=blank_tibble,complementary_frequency=blank_tibble)
  rm(blank_tibble)
  for(i in 1:l[4]){
    itmo=residue$interested[[i]]
    msu=questionr::wtd.table(x=itmo,weights=info$weight)
    msu=setNames(object=as.integer(msu),nm=names(msu))
    spbu=setNames(object=rep("",times=length(msu)),nm=names(msu))
    for(j in 1:length(msu)){
      spbu[j]=paste0(as.integer(msu[-j]),names(msu[-j])) %>% paste0(collapse=",")
    }
    #
    datum$frequency[[i]]=msu[itmo] %>% unname
    datum$other_frequencies[[i]]=spbu[itmo] %>% unname
    rm(itmo,msu,spbu)
  }
  for(i in 1:l[4]){
    itmo=residue$complementary_hash[[i]]
    msu=questionr::wtd.table(x=itmo,weights=info$weight)
    datum$complementary_frequency[[i]]=msu[itmo] %>% unname
    rm(itmo,msu)
  }
  #<<<calculate the 3 frequencies end<<<
  #>>>calculate info$small and info$big begin>>>
  for(i in 1:l[2]){
    itmo=(0.02*info$weight[i]) %>% ceiling()
    msu=datum$frequency[i,] %>% unlist() %>% unname()
    spbu=datum$complementary_frequency[i,] %>% unlist() %>% unname()
    #
    info$small$sites[[i]]=(1:l[4])[msu<=3 | (msu<10&msu<itmo)]
    info$small$numbers[[i]]=msu[ info$small$sites[[i]] ]
    info$small$min[i]=min(msu)
    info$small$min_site[i]=which.min(msu)
    #
    info$big$sites[[i]]=(1:l[4])[ spbu>info$weight[i] ]
    info$big$numbers[[i]]=spbu[ info$big$sites[[i]] ]
    info$big$max[i]=max(spbu)
    info$big$max_site[i]=which.max(spbu)
    #
    rm(itmo,msu,spbu)
  }
  #<<<calculate info$small and info$big end<<<
  #>>>change sites in info$small and info$big from relative site to absolute site begin>>>
  ans=(risf:(risf+l[3]-1))[residue_entry] %>% as.character()
  for(i in 1:l[2]){
    info$small$sites[[i]]=ans[ info$small$sites[[i]] ]
    info$big$sites[[i]]=ans[ info$big$sites[[i]] ]
  }
  info$small$min_site=ans[info$small$min_site]
  info$big$max_site=ans[info$big$max_site]
  rm(ans)
  #<<<change sites in info$small and info$big from relative site to absolute site end<<<
  #>>>sort items in xss, info, residue and datum begin>>>
  ans=rep("",times=l[2])
  for(i in 1:l[2]){ans[i]=paste0(info$id[[i]],collapse=", ");}
  ans=match(x=ans,table=str_sort(ans,numeric=TRUE))
  ans=order(-info$weight,-info$small$min,info$big$max,ans)
  xss=xss[ans]
  info=info[ans,]
  residue=residue[ans,]
  datum=datum[ans,]
  rm(ans)
  #<<<sort items in xss, info, residue and datum end<<<
  #>>>generate print_vector begin>>>
  print_vector=tibble(file_name=rep("",times=l[2]),small="",big="",id="")
  #
  print_vector$file_name=paste0("sequence frequency ",info$weight,
                                ", min of is ",info$small$min,
                                ", max without is ",info$big$max)
  for(i in 1:l[2])if(info$weight[i]<=3){
    print_vector$file_name[i]=paste0(print_vector$file_name[i],"; ",paste0(info$id[[i]],collapse=", "))
  }
  ans=f_create_directory(prefix=paste0("residue index start from ",risf))
  print_vector$file_name=paste0(ans,"/",1:l[2],"; ",print_vector$file_name,".txt")
  ans=file.create(print_vector$file_name)
  if(!all(ans)){stop("can't create text files. this may due to too long file names.")}
  rm(ans)
  #
  for(i in 1:l[2]){
    print_vector$small[i]=paste0(info$small$numbers[[i]],"@",info$small$sites[[i]],collapse=" ")
    print_vector$big[i]=paste0(info$big$numbers[[i]],"@",info$big$sites[[i]],collapse=" ")
    print_vector$id[[i]]=paste0(info$id[[i]],collapse="\n")
  }
  print_vector$small=paste0(" ",print_vector$small," ")
  print_vector$big=paste0(" ",print_vector$big," ")
  print_vector$id=unlist(print_vector$id) %>% paste0("\n",.,"\n")
  #
  print_vector$small=paste.xml.elements(names="small_frequencies_of",values=as.list(print_vector$small))
  print_vector$big=paste.xml.elements(names="big_frequencies_without",values=as.list(print_vector$big))
  print_vector$id=paste.xml.elements(names="sequence_ids",values=as.list(print_vector$id))
  print_vector$small=paste0(print_vector$small,"\n")
  print_vector$big=paste0(print_vector$big,"\n")
  print_vector$id=paste0(print_vector$id,"\n")
  #<<<generate print_vector end<<<
  #>>>generate print_matrix begin>>>
  print_matrix=matrix(data="",nrow=l[5],ncol=l[4]) %>% as_tibble(.name_repair="minimal")
  print_matrix=rep(x=list(print_matrix),times=l[2])
  #
  for(i in 1:l[2]){
    print_matrix[[i]][1,]=(risf:(risf+l[3]-1))[residue_entry] %>% as.character()
    print_matrix[[i]][2,]=residue$interested[i,] %>% as.character()
    print_matrix[[i]][3,]=datum$frequency[i,] %>% as.character()
    print_matrix[[i]][4,]=datum$other_frequencies[i,] %>% as.character()
    print_matrix[[i]][5,]=datum$complementary_frequency[i,] %>% as.character()
    #each print_matrix[[i]] is a 5-row tibble, nrow(print_matrix[[i]])==l[5]
  }
  for(i in 1:l[2])for(j in 1:l[4]){
    print_matrix[[i]][[j]]=print_matrix[[i]][[j]] %>% str_pad(string=.,width=max(nchar(.))+2,side="both")
  }
  blank=matrix(data=" ",nrow=l[5],ncol=l[3]) %>% as_tibble(.name_repair="minimal")
  blank[2,]=as.character(xss[[1]]) %>% unname() %>% strsplit(split="") %>% .[[1]]
  for(i in 1:l[2]){
    ans=blank
    ans[residue_entry]=print_matrix[[i]]
    print_matrix[[i]]=f_paste_byrow(df=ans) %>% paste0("~",.,"~")
    rm(ans)
  }
  rm(blank)
  #
  #set xml tags (2 of which have attributes) around print_matrix
  msu=c("residue_index","residue","frequency_of_residue","frequencies_of_other_residues",
        "refequency_without_residue")
  spbu=list(character(),character(),c(min=""),character(),c(max=""))
  for(i in 1:l[2]){
    ans=spbu
    ans[[3]]["min"]=paste0(" ",info$small$min[i],"@",info$small$min_site[i]," ")
    ans[[5]]["max"]=paste0(" ",info$big$max[i],"@",info$big$max_site[i]," ")
    print_matrix[[i]]=paste.xml.elements(names=msu,attrs=ans,values=as.list(print_matrix[[i]]),
                                         pad_start_tag=TRUE) %>% paste0("\n")
    rm(ans)
  }
  #<<<generate print_matrix end<<<
  #>>>export print_vector and print_matrix begin>>>
  for(i in 1:l[2]){
    cat(print_vector$small[i],print_vector$big[i],print_vector$id[i],print_matrix[[i]],
        file=print_vector$file_name[i],sep="",append=FALSE)
  }
  message("\n\nthe following files are created:\n",paste0(print_vector$file_name,collase="\n"),"\n\n")
  #<<<export print_vector and print_matrix end<<<
  i=0;j=0;rm(xss,risf,residue_entry,l,i,j,f_paste_byrow,f_create_directory,info,residue,datum,print_matrix)
  invisible(print_vector$file_name)#rm(print_vector)
}

if(FALSE){
  #here are the test codes
  ans=DNAStringSet(c("acgt","acct","acgt"))
  ans=DNAStringSet(c("a","c","a"))
  tempo=pick.rare.sequences(xss=ans)
}
