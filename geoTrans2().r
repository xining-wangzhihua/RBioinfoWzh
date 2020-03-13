#last edited at 20200313
#this is a stable version

#a function similar to pegas::geoTrans() but more robust
#a function to transform (or convert) geographic cooridinate string to number

require(magrittr)
require(stringr)
require(tibble)

geoTrans2=function(x){
  #>>>control the input of x begin>>>
  if(!is.character(x) | length(x) == 0){stop("x must be a non-empty character vector")}
  if(anyNA(x)){stop("x can't contain NA")}
  x %<>% enc2utf8() %>% toupper() %>% str_squish()
  #<<<control the input of x end<<<
  #>>>main manipulation begin>>>
  x %<>% str_remove_all(pattern="[\n\r\t]")
  x %<>% gsub(pattern="([NSEW])", replacement="\\1 ", x=.)
  x %<>% chartr(old="°′″", new="DMS", x=.) %>% chartr(old="\'\"", new="MS", x=.)
  x %<>% str_squish()
  #
  tulip = c("^(N |S |E |W |)([[:digit:]]+)D([[:digit:]]+)M([[:digit:]]{3})S$",
            "^(N |S |E |W |)([[:digit:]]+)D([[:digit:]]+)M([[:digit:]]{3})M$",
            "^(N |S |E |W |)([[:digit:]]+)D([[:digit:]]+)M([[:digit:].]+)M$",
            "^(N |S |E |W |)([[:digit:]]+)D([[:digit:].]+)M$",
            "^(N |S |E |W |)([[:digit:].]+)D$")
  x %<>% sub(pattern=tulip[1], replacement="\\1\\2D\\3.\\4M0S", x=.)
  x %<>% sub(pattern=tulip[2], replacement="\\1\\2D\\3.\\4M0S", x=.)
  x %<>% sub(pattern=tulip[3], replacement="\\1\\2D\\3M\\4S", x=.)
  x %<>% sub(pattern=tulip[4], replacement="\\1\\2D\\3M0S", x=.)
  x %<>% sub(pattern=tulip[5], replacement="\\1\\2D0M0S", x=.)
  remove(tulip)
  #
  tulip = "^(N |S |E |W |)([[:digit:].]+)D([[:digit:].]+)M([[:digit:].]+)S$"
  tulip = strcapture(proto=data.frame(direction="", degree=0, minute=0, second=0, stringsAsFactors=FALSE),
                     pattern=tulip, x=x)#stringr::str_match() won't work for non-ACSII character
  tulip$direction %<>% str_trim(side="right")
  for(i in 1:nrow(tulip)) if(!is.na(tulip[[1]][i])){
    ans = tulip[i,]
    ans %<>% {(.$degree<(-180))|(.$degree>180)|(.$minute<0)|(.$minute>60)|(.$second<0)|(.$second>60);}
    if(ans){tulip[i,]=NA;}
    remove(ans)
  }
  x = tibble(direction=tulip$direction, degree=tulip$degree + tulip$minute / 60 + tulip$second / 3600)
  remove(tulip, i)
  #<<<main manipulation end<<<
  return(x)
}

if(FALSE){
  #here are the test codes
  x=c("N 43°27'30\"",
      "N43°27'30\"",
      "43°27'30\"",
      "43°27.5'",
      "n 43d27'30\"")
  geoTrans2(x)
}
