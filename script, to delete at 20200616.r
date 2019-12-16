#last edited at 20191216
library(Biostrings)


file_list=list.files(path="?must specify?",pattern="?must specify?",full.names="TRUE",recursive=TRUE)


l=length(file_list)
datum=rep("",times=l)
names(datum)=gsub(pattern="^.+/(.+)\\.seq$",replacement="\\1",x=file_list)
for(i in 1:l){
  datum[[i]]=paste0(readLines(con=file_list[i]),collapse="")
}


datum=DNAStringSet(x=datum,use.names=TRUE)
writeXStringSet(x=datum,filepath="?must specify?",format="fasta",width=80)


rm(file_list,datum,l,i)
