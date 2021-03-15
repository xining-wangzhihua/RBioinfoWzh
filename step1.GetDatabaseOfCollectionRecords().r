#this is a stable version
#last edited at 20200326

#get dependencies
if(FALSE){
  Sys.setlocale(locale="Chinese (Simplified)")
  require(magrittr); require(stringr); require(tibble); require(dplyr);
  require(readxl)
  #require(gdata)
  source(file="./geoTrans3().r", encoding="UTF-8")
  source(file="./as.Date.character2().r")
  source(file="./SplitSpeciesNames().r")
  #
  format.of.db.collection = read.table(file="./auxiliary files/format of database of collection records.txt",
                                       header=TRUE, sep="\t", row.names=NULL, stringsAsFactors=FALSE,
                                       fileEncoding="UTF-8", colClasses="character",
                                       quote="\"\'", allowEscapes=FALSE, skip=4, comment.char="#")
  if(dim(format.of.db.collection) %>% identical(., c(28L, 5L)) %>% !.){stop("auxiliary files bug");}
  #names(format.of.db.collection) %<>% enc2utf8()
  for(i in 1:length(format.of.db.collection)){format.of.db.collection[[i]] %<>% enc2utf8();}; remove(i);
  f_ans = function(datum, is.function.list){
    if(!is.function.list) for(i in 1:length(datum)){
      datum[i] %<>% paste0("\"", ., "\"") %>% str2expression() %>% eval() %>% enc2utf8()
    }
    if(is.function.list) for(i in 1:length(datum)) if(datum[[i]] %>% {. != "" & . != "f_makeshift"}){
      datum[[i]] %<>% enc2utf8() %>% str2expression() %>% eval()
    }
    return(datum)
    remove(i, datum)
  }
  format.of.db.collection$alias %<>% strsplit(split=";")
  format.of.db.collection$valid_characters %<>% f_ans(is.function.list=FALSE)
  format.of.db.collection$function_to_apply %<>% as.list() %>% f_ans(is.function.list=TRUE)
  format.of.db.collection$valid_format %<>% f_ans(is.function.list=FALSE)
  remove(f_ans)
  #
  f_import_db = function(db.file, db.name){
    db.file = read.table(file=db.file, sep="\t", row.names=NULL, col.names=c("full", "short"),
                         header=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8")
    if(ncol(db.file) !=2){stop("table specified by ", db.name, "must contain 2 columns");}
    for(i in 1:length(db.file)){
      db.file[[i]] %<>% enc2utf8() %>% stringr::str_squish()
      if(anyNA(db.file[[i]])){stop("table specified by ", db.name, " can't contain NA");}
    }; remove(i);
    if(anyDuplicated(db.file$full) != 0){
      stop("the first column of ", db.name, " is a key column, thus can't contain duplicated items")
    }
    if(intersect(x=db.file$full, y=db.file$short) %>% length() %>% {. != 0;}){
      stop("2 columns of the reference table can't contain same items")
    }
    return(db.file)
    remove(db.file, db.name)
  }
  db.scientific.names = f_import_db(db.file="./auxiliary files/db of species names.tsv", db.name="species")
  db.province.names = f_import_db(db.file="./auxiliary files/db of province names.tsv", db.name="province")
  db.city.names = f_import_db(db.file="./auxiliary files/db of city names.tsv", db.name="city")
  db.district.names = f_import_db(db.file="./auxiliary files/db of district names.tsv", db.name="district")
  remove(f_import_db)
}

#messsage() and warning()
if(TRUE){
  warning("in the part outputing x_invalid_rows, try to find a function escapes \\n\\r\\f\\t\\v\n")
  warning("try to find a function to validate people names in (original/revised)_analyser\n")
  warning("the part checking consistency between 原始鉴定 & 修订鉴定 hasn't been written\n")
  warning("try to use ape::label2table() to validate (原始/修订)种名")
  warning("use a new function to replace NA with \"\"")
  message("\nyou can use readxl::read_excel(path=, sheet=), readxl::read_xls(path=, sheet=)",
          " or readxl::read_xlsx(path=, sheet=) to read in an excel file\n")
  message("\nhow to read in an .ods file is unknown at present\n")
  message("\nwhen you use date format in excel, libreoffice, or r, etc, please notice:\n",
          "excel count date from 1899-12-31, and treat 1900-02-29 as existing (in my office suite), ",
          "so change the format of \"60\" from general to date in excel will generate 1900-02-29.\n",
          "libreoffice count date from 1899-12-30, and doesn't treat 1900-02-29 as existing\n",
          "r doesn't specify a date from which it count, and doesn't treat 1900-02-29 as existing\n",
          "if there is an error from date in integer format, try minusing or adding the integer by 1\n")
  #
  #synonym of excise: snip, remove, trim, resect, erase, eliminate, efface, expunge, expurgate
  #synonym of standardized: normalize, convert, homogenize, assimilate, uniform
  #synonym of forbid: prohibit, ban, illegal, disallow, taboo, proscribe
  #
  #in zm.ValidateDatabaseOfCollection, check if all herbarium_id and collection_id is unique
  #in zm.ValidateDatabaseOfCollection, message (unique) items in each column
  #in zm.ValidateDatabaseOfCollection, sort segements in each herbarium_id, habitat
  #in zm.ValidateDatabaseOfCollection, sort records according to several columns
}


step1.GetDatabaseOfCollectionRecords = function(x){
  ans = "x must be a data.frame. db.*.name should be a file name or NULL."; remove(ans);
  #get constants (not being allowed being set via parameters) in parent environment-------
  format.of.db.collection = get(x="format.of.db.collection", envir=parent.frame())
  db.scientific.names = get(x="db.scientific.names", envir=parent.frame())
  db.province.names = get(x="db.province.names", envir=parent.frame())
  db.city.names = get(x="db.city.names", envir=parent.frame())
  db.district.names = get(x="db.district.names", envir=parent.frame())
  #inner (internal) constants----------------------------
  destination = strftime(x=Sys.time(),format="r, zm.step1.GetDatabaseOfCollectionRecords, %Y%m%d-%H%M%S")
  dir.create(path=destination)
  message("the directory:\n", destination, "\nis created. output files will be written there\n")
  #inner (internal) expression constants--------------------------
  MoveInvalidRecords = c("message(\"the following invalid items are (generated) in column \", names(x[i]), \":\")",
                         "message(paste0(x[!is_valid,i,drop=TRUE], collapse=\"\n\"), \"\n\")",
                         "x_invalid_rows <- rbind(x_invalid_rows, x_backup[!is_valid,])",
                         "x <- x[is_valid,]",
                         "x_backup <- x_backup[is_valid,]",
                         "l <- nrow(x)",
                         "if(l == 0){stop(\"all records (rows) in the table are invalid\");}"
                         ) %>% str2expression()
  #check input: x, and get basic information----------------------
  if(!is.data.frame(x)){stop("x must be a data.frame")}
  if(nrow(x) == 0){stop("x can't be empty");}
  for(i in 1:length(x)) if(!is.character(x[[i]])){stop("each column of x must be in character mode");}
  if(is.null( names(x) )){stop("names(x) can't be NULL");}
  if(anyNA( names(x) )){stop("names(x) can't contain NA");}
  names(x) %<>% enc2utf8() %>% str_squish()
  if(any(names(x) == "")){stop("each items of names(x) can't be empty");}
  if(anyDuplicated(names(x)) != 0){stop("names(x) can't contain duplicated items");}
  #
  for(i in 1:length(x)){
    x[[i]] %<>% {replace(x=., list=is.na(.), values="") %>% enc2utf8();}
  }
  x_backup = x
  x_invalid_rows = x[integer(),]
  l = nrow(x)
  #
  tulip = format.of.db.collection %>% {setNames(object=rep(.$column_name, times=lengths(.$alias)), nm=unlist(.$alias));}
  tulip = tulip[names(x)] %>% unname()
  tulip[is.na(tulip)] = names(x)[is.na(tulip)]
  names(x) = tulip
  remove(tulip)
  #
  if(setdiff(names(x), format.of.db.collection$column_name) %>% length() %>% {. > 5}){
    stop("there are too many (> 5) columns with unrecognized names")
  }
  #check if there is any invalid characters in x---------------------
  for(i in 1:length(x)) if(names(x[i]) %in% format.of.db.collection$column_name){
    ans = which(format.of.db.collection$column_name == names(x[i]))
    ans = format.of.db.collection$valid_characters[ans]
    is_valid = grepl(pattern=ans, x=x[[i]])
    if(!all(is_valid)){eval(expr=MoveInvalidRecords);}
    remove(ans, is_valid)
  }
  for(i in 1:length(x)){x[[i]] %<>% stringr::str_squish();}
  #replace "f_makeshift" with complex functions in format.of.db.collection------------------
  if(FALSE){
    #used for test
    ans = rep(FALSE, times=nrow(format.of.db.collection))
    for(i in 1:length(ans)) if(identical(format.of.db.collection$function_to_apply[[i]], "f_makeshift")){ans[i] = TRUE;}
    format.of.db.collection$column_name[ans]#these column need a specified function
    remove(ans)
    format.of.db.collection[format.of.db.collection$column_name=="", 3:5]
  }
  #
  #a function to apply to "植株高度"
  format.of.db.collection$function_to_apply[format.of.db.collection$column_name == "植株高度"][[1]] = function(datum){
    datum %<>% tolower() %>% chartr(old=enc2utf8("-。厘米"), new="~.cm", x=.)
    datum %<>% gsub(pattern=enc2utf8("大约"),replacement="ca.",x=.)
    datum %<>% gsub(pattern=enc2utf8("约"),replacement="ca.",x=.)
    datum %<>% gsub(pattern="about",replacement="ca.",x=.)
    datum %<>% gsub(pattern="circa",replacement="ca.",x=.)
    datum %<>% gsub(pattern="c\\.",replacement="ca.",x=.)
    datum %<>% gsub(pattern="ca\\.",replacement="ca. ",x=.)
    return(str_squish(datum))
  }
  #a function to apply to "采集编号"
  format.of.db.collection$function_to_apply[format.of.db.collection$column_name == "采集编号"][[1]] = function(datum){
    datum %<>% gsub(pattern="gao", replacement="Gao", x=., ignore.case=TRUE)
    datum %<>% gsub(pattern="chen", replacement="Chen", x=., ignore.case=TRUE)
    datum %<>% gsub(pattern="zhang", replacement="Zhang", x=., ignore.case=TRUE)
    return(datum)
  }
  #a function to apply to "date_clutter"
  format.of.db.collection$function_to_apply[format.of.db.collection$column_name == "date_clutter"][[1]] = function(datum){
    string_index = grepl(pattern="[^[:digit:]]", x=datum)
    strings = datum[string_index] %>% as.Date.character2()
    if(is.null(strings)){
      stop("date_clutter(采集日期) formats error. please split the table according to date formats.")
    }else{strings = as.character(strings);}
    numbers = datum[!string_index] %>% as.integer() %>% as.Date.numeric(origin="1899-12-31") %>% as.character()
    datum[string_index] = strings
    datum[!string_index] = numbers
    rm(string_index, strings, numbers)
    return(datum)
  }
  #a function to apply to "省"
  format.of.db.collection$function_to_apply[format.of.db.collection$column_name == "省"][[1]] = function(datum){
    destination = get("destination", envir=parent.frame())
    db = get("db.province.names", envir=parent.frame()) %>% tibble::deframe()
    db.increment = datum
    datum = unname(db[datum])
    #
    db.increment = db.increment[is.na(datum)]
    if(length(db.increment) > 0) if(any( grepl(pattern=enc2utf8("省"), x=db.increment) )){
      ans = tibble(full=db.increment, short=str_remove_all(string=db.increment, pattern=enc2utf8("省"))) %>% unique()
      ans = ans[ans$full != ans$short,]
      write.table(x=ans, file=paste0(destination, "/db.province.names, increment.tsv"), sep="\t",
                  quote=TRUE, row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8")
      remove(ans)
      message("increment of db.province.names are written to:\n",
              paste0(destination, "/db.province.names, increment.tsv"),
              "\nplease check this file for erroneous conversions\n")
      db.increment %<>% str_remove_all(pattern=enc2utf8("省"))
    }
    #
    if(length(db.increment) > 0){datum[is.na(datum)] = db.increment;}
    return(datum)
    remove(db, db.increment, datum)
  }
  #a function to apply to "市"
  format.of.db.collection$function_to_apply[format.of.db.collection$column_name == "市"][[1]] = function(datum){
    destination = get("destination", envir=parent.frame())
    db = get("db.city.names", envir=parent.frame()) %>% tibble::deframe()
    db.increment = datum
    datum = unname(db[datum])
    #
    db.increment = db.increment[is.na(datum)]
    if(length(db.increment) > 0) if(any( grepl(pattern=enc2utf8("地区|市"), x=db.increment) )){
      ans = tibble(full=db.increment, short=str_remove_all(string=db.increment, pattern=enc2utf8("地区|市"))) %>% unique()
      ans = ans[ans$full != ans$short,]
      write.table(x=ans, file=paste0(destination, "/db.city.names, increment.tsv"), sep="\t",
                  quote=TRUE, row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8")
      remove(ans)
      message("increment of db.city.names are written to:\n",
              paste0(destination, "/db.city.names, increment.tsv"),
              "\nplease check this file for erroneous conversions\n")
      db.increment %<>% str_remove_all(pattern=enc2utf8("地区|市"))
    }
    #
    if(length(db.increment) > 0){datum[is.na(datum)] = db.increment;}
    return(datum)
    remove(db, db.increment, datum)
  }
  #a function to apply to "县"
  format.of.db.collection$function_to_apply[format.of.db.collection$column_name == "县"][[1]] = function(datum){
    destination = get("destination", envir=parent.frame())
    db = get("db.district.names", envir=parent.frame()) %>% tibble::deframe()
    db.increment = datum
    datum = unname(db[datum])
    #
    db.increment = db.increment[is.na(datum)]
    if(length(db.increment) > 0) if(any( grepl(pattern=enc2utf8("县"), x=db.increment) )){
      ans = tibble(full=db.increment, short=str_remove_all(string=db.increment, pattern=enc2utf8("县"))) %>% unique()
      ans = ans[ans$full != ans$short,]
      write.table(x=ans, file=paste0(destination, "/db.district.names, increment.tsv"), sep="\t",
                  quote=TRUE, row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8")
      remove(ans)
      message("increment of db.district.names are written to:\n",
              paste0(destination, "/db.district.names, increment.tsv"),
              "\nplease check this file for erroneous conversions\n")
      db.increment %<>% str_remove_all(pattern=enc2utf8("县"))
    }
    #
    if(length(db.increment) > 0){datum[is.na(datum)] = db.increment;}
    return(datum)
    remove(db, db.increment, datum)
  }
  #a function to apply to "纬度"
  format.of.db.collection$function_to_apply[format.of.db.collection$column_name == "纬度"][[1]] = function(datum){
    datum = geoTrans3(x=datum, allowed_directions=c("N", "S"))
    datum$degree = format(datum$degree, digits=16)
    datum = paste0(datum$direction, " ", datum$degree) %>% {replace(x=., list=(. == "NA NA"), values=NA);}
    return(datum)
  }
  #a function to apply to "经度"
  format.of.db.collection$function_to_apply[format.of.db.collection$column_name == "经度"][[1]] = function(datum){
    datum = geoTrans3(x=datum, allowed_directions=c("E", "W"))
    datum$degree = format(datum$degree, digits=16)
    datum = paste0(datum$direction, " ", datum$degree) %>% {replace(x=., list=(. == "NA NA"), values=NA);}
    return(datum)
  }
  #a function to apply to "(原始/修订)种名"
  f_ans = function(datum){
    destination = get("destination", envir=parent.frame())
    db = get("db.scientific.names", envir=parent.frame()) %>% tibble::deframe()
    db.increment = datum
    datum = unname(db[datum])
    #
    db.increment = db.increment[is.na(datum)]
    if(length(db.increment) > 0){
      ans = SplitSpeciesNames(x=db.increment) %>% {.[!is.na(.[,1]),,drop=FALSE];}
      #remove NAs generated by SplitSpeciesNames(), that's why don't use SplitSpeciesNames(db.increment)==db.increment before entering this condition
      ans = tibble(full=paste(ans[,1], ans[,2], sep=" "), short=ans[,1]) %>% unique()
      ans = ans[ans$full != ans$short,]
      if(nrow(ans) > 0){
        write.table(x=ans, file=paste0(destination, "/db.scientific.names, increment.tsv"), sep="\t",
                    quote=TRUE, row.names=FALSE, col.names=FALSE, fileEncoding="UTF-8")
        message("increment of db.scientific.names are written to:\n",
                paste0(destination, "/db.scientific.names, increment.tsv"),
                "\nplease check this file for erroneous conversions\n")
        db.increment_backup = db.increment
        db.increment = tibble::deframe(ans)[db.increment] %>% unname()
        db.increment[is.na(db.increment)] = db.increment_backup[is.na(db.increment)]
        remove(db.increment_backup)
      }
      remove(ans)
    }
    #
    if(length(db.increment) > 0){datum[is.na(datum)] = db.increment;}
    return(datum)
    remove(db, db.increment, datum)
  }
  format.of.db.collection$function_to_apply[format.of.db.collection$column_name == "原始种名"][[1]] = f_ans
  format.of.db.collection$function_to_apply[format.of.db.collection$column_name == "修订种名"][[1]] = f_ans
  remove(f_ans)
  #a function to apply to "boolean采集个体数"
  format.of.db.collection$function_to_apply[format.of.db.collection$column_name == "boolean采集个体数"][[1]] = function(datum){
    datum[datum %in% c("1", "yes")] = "TRUE"
    datum[datum %in% c("0", "no")] = "FALSE"
    datum[!(datum %in% c("1", "yes", "TRUE", "0", "no", "FALSE"))] = NA_character_
    return(datum)
  }
  #
  for(i in 1:nrow(format.of.db.collection)) if(!is.function(format.of.db.collection$function_to_apply[[i]])){stop("bug");}
  #apply functions to each column--------------
  for(i in 1:length(x)) if(names(x[i]) %in% format.of.db.collection$column_name){
    f_ans = format.of.db.collection %>% {.$function_to_apply[.$column_name == names(x[i])][[1]];}
    x[[i]][x[[i]] != ""] %<>% f_ans(datum=.)
    if(!is.character(x[[i]])){stop("bug: a functions to apply didn't return character vectors.");}
    x[[i]] %<>% str_squish()#here str_squish() may be redundant, but ensure results
    if(anyNA(x[[i]])){
      is_valid = !is.na(x[[i]])
      eval(expr=MoveInvalidRecords)
      remove(is_valid)
    }
    remove(f_ans)
  }
  #check if there is any invalid format in x-----------------------------------------------------
  for(i in 1:length(x)) if(names(x[i]) %in% format.of.db.collection$column_name){
    ans = format.of.db.collection$valid_format[names(x[i]) == format.of.db.collection$column_name]
    is_valid = grepl(pattern=ans, x=x[[i]]) | (x[[i]] == "")
    if(!all(is_valid)){eval(expr=MoveInvalidRecords);}
    remove(is_valid, ans)
  }
  #check the consistency between repetition columns------------
  if(all(c("date", "date_clutter") %in% names(x))){
    ans = (x$date_clutter != "") & (x$date != "")
    ans = any(x$date_clutter[ans] != x$date[ans])
    if(ans){stop("date_clutter(采集日期) doesn't match with date(采集日期_character)");}
    remove(ans)
  }
  #check the consistency between original_analysis and revised_analysis-------
  #original and revised data can't be the same
  #leave this part of codes for lack of time
  #sort segments in 标本馆编号 and 生境-----------------
  if("标本馆编号" %in% names(x)) for(i in 1:l){
    ans = x$`标本馆编号`[i] %>% strsplit(split=";") %>% .[[1]]
    ans =ans[ans != ""]#this step is useless, just ensures results
    x$`标本馆编号`[i] = str_sort(x=ans, numeric=TRUE) %>% paste0(collapse=";")
    remove(ans)
  }
  if("生境" %in% names(x)) for(i in 1:l){
    ans = x$`生境`[i] %>% strsplit(split=";") %>% .[[1]]
    ans =ans[ans != ""]#this step is useless, just ensures results
    x$`生境`[i] = sort(x=ans) %>% paste0(collapse=";")
    remove(ans)
  }
  #sort records in x according to `采集编号` and then `标本馆编号`--------------------------------
  x = arrange(.data=x, `采集编号`, `标本馆编号`)
  #message() if there are invalid records or not-------------------------------------------------
  if(nrow(x_invalid_rows) > 0){
    f_ans = function(p, r, datum){gsub(pattern=p, replacement=r, x=datum, fixed=TRUE);}
    for(i in 1:length(x_invalid_rows)){
      x_invalid_rows[[i]] %<>% f_ans(p="\\", r="\\\\", datum=.) %>% f_ans(p="\n", r="\\n") %>%
        f_ans(p="\r", r="\\r") %>% f_ans(p="\f", r="\\f") %>% f_ans(p="\t", r="\\t") %>% f_ans(p="\v", r="\\v")
    }
    remove(f_ans)
    write.table(x=x_invalid_rows, file=paste0(destination, "/invalid records.tsv"), sep="\t",
                quote=TRUE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")
    message("invalid records (rows) are written to:\n", paste0(destination, "/invalid records.tsv\n"))
  }else{message("all records (rows) are valid\n");}
  #fill x to be a full table, and put un-recognized columns into remark-------------------------------------
  if(all(c("date", "date_clutter") %in% names(x))){x = x[names(x) != "date_clutter"];}
  if(!("date" %in% names(x)) & ("date_clutter" %in% names(x))){
    warning("please check the \"date\" column, as there is inconsistency among excel, libreoffic, r, etc\n")
    names(x)[names(x) == "date_clutter"] = "date"
  }
  #
  tulip = list(rep("", times=l)) %>% rep(times=nrow(format.of.db.collection))
  tulip %<>% setNames(object=., nm=format.of.db.collection$column_name) %>% as_tibble()
  for(i in 1:length(x)) if(names(x[i]) %in% format.of.db.collection$column_name){
    tulip[[names(x[i])]] = x[[i]]
  }
  #
  if(!all(names(x) %in% format.of.db.collection$column_name)){
    x = x[!(names(x) %in% format.of.db.collection$column_name)]
    for(i in 1:length(x)){x[[i]] %<>% paste0(names(x[i]), ": ", .);}
    names(x) = paste0("remark_", 2:(length(x)+1))
    tulip = dplyr::bind_cols(tulip, x)
  }
  x = tulip
  remove(tulip)
  #output x---------------------------------------------------------------------------
  for(i in 1:length(x)) if(anyNA(x[[i]])){stop("bug");}
  for(i in 1:length(x)) if(any( grepl(pattern="[\n\r\f\t\v]", x=x[[i]]) )){stop("bug");}
  write.table(x=x, file=paste0(destination, "/db of collection.tsv"), sep="\t",
              quote=TRUE, row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")
  message("the new generated database of collection records is written in:\n",
          paste0(destination, "/db of collection.tsv\n"))
  remove(db.scientific.names, db.province.names, db.city.names, db.district.names,
         format.of.db.collection, MoveInvalidRecords, destination, l, i, x_backup, x_invalid_rows)
  invisible(x)
}

#here are the test codes
if(FALSE){
  #test if all empty cells are allowed in each column
  x = "./test/arrange again, Micranthes records before 2018.xls"
  x = "./test/20180725 look for Micrathes specimen by records.tsv"
  x = "./test/db of collection.tsv"
  x = "./r, zm.step1.GetDatabaseOfCollectionRecords/invalid records.tsv"
  x = "./r, zm.step1.GetDatabaseOfCollectionRecords/db of collection.tsv"
  #x = readxl::read_xls(path=x, sheet=excel_sheets(path=x)[1], .name_repair="minimal")
  x = read.table(file=x, header=TRUE, sep="\t", quote="\"", row.names=NULL,
                 allowEscapes=FALSE, colClasses="character",
                 check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8")
  names(x) %<>% enc2utf8()
  for(i in 1:length(x)){
    cat(i, "\t", class(x[[i]]), "\n")
    x[[i]] %<>% as.character() %>% enc2utf8() %>% {replace(x=., list=is.na(.), values="");}
  }; remove(i);
  #
  ans = x$`采集日期`
  ans[!grepl(pattern="[^[:digit:]]", x=ans)] %<>% as.integer() %>% {. - 1} %>% as.character()
  x$`采集日期` = ans
  x$`采集编号` %<>% gsub(pattern="chen", replacement="Chen") %>% gsub(pattern="gao", replacement="Gao") %>% gsub(pattern="zhang", replacement="Zhang")
  x$`纬度`
  x$`经度`
  remove(ans)
  #
  x = arrange(.data=x, `采集编号`)
  y = x
  write.table(x=x, file="./test/db of collection - new.tsv", quote=TRUE,
              sep="\t", row.names=FALSE, col.names=TRUE, fileEncoding="UTF-8")
  x = read.table(file="./test/db of collection - new.tsv", header=TRUE, sep="\t",
                 row.names=NULL, colClasses="character", allowEscapes=FALSE,
                 check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8")
  identical(x=x, y=y)
  remove(y)
  #
  zm.step1.GetDatabaseOfCollectionRecords(x=x)
  remove(x)
}
