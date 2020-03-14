#last edited at 20200314

#a function similar to base::as.Date.character(), but handle more regular expressions while trying formats
#and it only success if all strings match the regular expression and no NA produced, otherwise it will return NULL

as.Date.character2 = function(
  x, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y", "%Y/%m/%d", "%Y%m%d")
){
  i = 0; pending = TRUE; while(pending){
    i = i + 1
    if(!anyNA( as.Date(x=x, format=tryFormats[i]) )){pending = FALSE;}
    if(i == length(tryFormats)){pending = FALSE;}
  }
  x = as.Date(x=x, format=tryFormats[i])
  rm(tryFormats, pending, i)
  if(anyNA(x)){return(NULL);}else{return(x);}
}

if(FALSE){
  #here are the test codes
  ans = c("2020-11-12", "2020-12-11")
  ans = c("11/12/2020", "12/11/2020")
  ans = c("12/13/2020", "11/12/2020")
  ans = c("13/12/2020", "12/11/2020")
  ans = c("2020/11/12", "2020/12/11")
  ans = c("20200202", "202022")
  as.Date.character2(x=ans)
}
