#! /usr/bin/Rscript
# license: BSD 3-Clause License
# author: https://github.com/ywd5
# version: 20210314

# this function don't accept PHASE command-line options, it uses specifically defined ones
# adjusting smaller_for_faster[1] may has more effect than smaller_for_faster[2]

if("" == "get dependencies"){
  require(magrittr); require(tibble); require(dplyr); require(stringr); require(purrr);
  require(parallel)
}

if("" == "test codes"){
  setwd("C:/Users/wangz/OneDrive/github/RBioinfoZm/new phase/")
  input="./input.txt"; input_d="./input_d.txt"; input_k="./input_k.txt";
  smaller_for_faster=c(0.5, 20) # for 1~3
  smaller_for_faster=c(0.2, 5) # for 4
  smaller_for_faster=c(0.05, 1) # for 5
}

RunPhase <- function(input="./input.txt", input_d="./input_d.txt", input_k="./input_k.txt",
                     smaller_for_faster=c(5, 200)){
  # check if PHASE is installed---------------------------------------------------------------------
  Sys.which("PHASE")[1] %>% {if(. == ""){stop("can't find PHASE in shell path");}}
  .Platform$OS.type %>% grepl(pattern="windows", x=., ignore.case=TRUE) %>%
    {if(.){warning("the author of PHASE said (see website), PHASE in Windows isn't very stable.");}}
  # check the input---------------------------------------------------------------------------------
  if(!is.character(input)){stop("input is invalid");}
  if(length(input) != 1){stop("input is invalid");}
  if(is.na(input)){stop("input is invalid");}
  if(!is.character(input_d)){stop("input_d is invalid");}
  if(length(input_d) != 1){stop("input_d is invalid");}
  if(is.na(input_d)){stop("input_d is invalid");}
  if(!is.character(input_k)){stop("input_k is invalid");}
  if(length(input_k) != 1){stop("input_k is invalid");}
  if(is.na(input_k)){stop("input_k is invalid");}
  #
  if(!file.exists(input)){stop("the file specified by input doesn't exist");}
  #
  if(!is.numeric(smaller_for_faster)){stop("smaller_for_faster is invalid");}
  if(length(smaller_for_faster) != 2L){stop("smaller_for_faster is invalid");}
  if(anyNA(smaller_for_faster)){stop("smaller_for_faster is invalid");}
  if(any(smaller_for_faster <= 0)){stop("smaller_for_faster must be > 0");}
  # get output file names---------------------------------------------------------------------------
  output <- base::expand.grid(option=c("-MR0", "-MR1 1", "-MR3"), repeatation=as.character(1:10),
                              stringsAsFactors=FALSE) %>%
    as_tibble() %>% dplyr::arrange(option) %>%
    dplyr::mutate(directory=paste0("./output ", option, ", ", repeatation),
                  file=paste0(directory, "/output.nexus"))
  output$directory %>% dir.exists() %>% any() %>%
    {if(.){warning("Some of the output directories already exist.");}}
  output$directory %>% lapply(FUN=dir.create, showWarnings=FALSE)
  output$directory %>% dir.exists() %>% {!all(.)} %>%
    {if(.){stop("direcotory creation error");}}
  # get PHASE parameters----------------------------------------------------------------------------
  phase_parameters <- readLines(con=input) %>% stringr::str_trim(string=., side="both") %>%
    {.[. != ""]} %>% .[1:2] %>% as.integer() %>% {1 / (.[1] * .[2] * 2) / smaller_for_faster[1]} %>%
    base::signif(x=., digits=3) %>% base::format(scientific=FALSE) %>% paste0(c("-F", "-O"), .) %>%
    c("-f0", ., paste0("-d\"", input_d, "\""), paste0("-k\"", input_k, "\""),
      paste0("\"", input, "\"")) %>% paste0(collapse=" ")
  phase_parameters <- readLines(con=input) %>% stringr::str_trim(string=., side="both") %>%
    {.[. != ""]} %>% .[2] %>% as.integer() %>% {. * smaller_for_faster[2]} %>%
    base::cbind(output$option, phase_parameters, paste0("\"", output$file, "\""), .) %>%
    purrr::array_tree(array=., margin=1) %>% lapply(FUN=base::unname)
  # run PHASE---------------------------------------------------------------------------------------
  a_cluster <- parallel::detectCores(logical=TRUE) %>% {. * 0.75;} %>%
    ceiling() %>% parallel::makeCluster(spec=., outfile="")
  message("\nusing ", length(a_cluster), " cluster(s)\n")
  f_show_process <- function(xx){
    cat(xx[3], " begins running\n") # cat() inside clusterApplyLB() won't print out information
    base::system2(command="PHASE", args=xx) # system(show.output.on.console=FALSE)
    cat(xx[3], " ends running\n")
    invisible(NULL)
  }
  cat("\nPHASE begins running\n")
  parallel::clusterApplyLB(cl=a_cluster, x=phase_parameters, fun=f_show_process)
  parallel::stopCluster(cl=a_cluster)
  cat("\nPHASE ends runing\n")
  remove(a_cluster, f_show_process)
  # return result-----------------------------------------------------------------------------------
  remove(input, input_d, input_k, smaller_for_faster, phase_parameters)
  message("PHASE outputs were written to the following files:\n",
          paste0(output$file, collapse="\n"), "\n")
  message("PHASE outputs were written to the following directories:\n",
          paste0(output$directory, collapse="\n"), "\n")
  return(output$file)
}
