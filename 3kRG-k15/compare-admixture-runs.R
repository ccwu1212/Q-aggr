
#' This R file contains functions that process Admixture / fastStructure output such as:
#' 1) preparing files for CLUMPP and executing CLUMPP
#' 2) creating admixture barplot


########################### CONFIG (EDIT AS NECESSARY) ########################
# Please edit the following line to point to the path to your CLUMPP executable
CLUMPP_exe = '/soft/CLUMPP/CLUMPP_MacOSX.1.1.2/CLUMPP'

###############################
# Check if CLUMPP exists 
if( !file.exists(CLUMPP_exe)){
  cat("Please provide a path to your CLUMPP executable:\n")
  while(TRUE) { CLUMPP_exe = readline(); if(nchar( CLUMPP_exe) > 0) break}
  if( !file.exists(CLUMPP_exe)){
    stop( paste0("\nThe file\n", CLUMPP_exe, "\nwas not found.\n"))
  }
}

########################
library(fpc)
library(matrixStats)
library(igraph)

library(tools)

#########################
#' Admixture plot
graph_allQ <- function(permQ_file,  # CLUMPP-generated permutedQfile
                       Ks = "",  # list of values of K used
                       pops,  # population assignment for samples
                       reorder.by.max = FALSE,  # do reordering of samples within pop by max admix coef
                       hclust.reorder = TRUE,  # do reordering of samples by hierarchical clustering
                       spacing_between_groups = 10, # size of blanks between populations on the graph
                       min.group.size.show = 1, # minimal size of a group to appear on graph
                       admix_comp_colors=NULL,
                       admix_color_table = NULL,
                       ylabels = NULL,
                       save_Q_file = "",
                       pop_sequence = NULL,
                       ...    # additional arguments to barplot you might want to give
                       )
{

  if(is.character(permQ_file)){
    clumpQ = read.table( permQ_file)
    K = ncol(clumpQ)-5    # maximal K
    clumpQ = clumpQ[, -c(1,3,5)]    # remove dummy columns
    
    names(clumpQ)= c("Ind", "Pop", paste("C", 1:(ncol(clumpQ)-2), sep="_" ))
    pQ = clumpQ
    I = length(unique(pQ[,1]))    # number of individuals
    R = nrow(clumpQ)/I   # number of runs
    
    Ks = rep_len(Ks, R) # so you can use Ks = scalar if you have the same K 
    
    QQ = clumpQ[,-c(1,2) ]   # all Q coefs
  } else if(length(dim(permQ_file))==3){
    aQ = permQ_file  #array
    R = dim(aQ)[[3]]
    I = dim(aQ)[[1]]
    K = dim(aQ)[[2]]
    QQ = array2permQ(aQ)
  } 
  else if( length(dim(permQ_file))==2){
    QQ = permQ_file
    R = 1
    I = nrow(QQ)
    K = ncol(QQ)
   }
  
  cat("I: ", I, "  K:", K, " R:", R, "\n")
  QQ = as.data.frame(QQ)
  byRun = rep(1:R, each=I)  # column indicating the run number
  Qs = split(QQ, byRun)
  Qs = lapply(Qs, as.matrix)
  
  Qbar = do.call(cbind, Qs)
  
  #lastQ = Qs[[ length(Qs) ]] #should not change, only used for pop assignments
  pops0 = pops
  
  # Diagnostic image: 
  #pdf("diag_Qperm_test.pdf")
  #image(as.matrix(dist(t(Qbar))), useRaster=T)
  #dev.off()

  rownames(Qbar) = pops
  colnames(Qbar) = rep( paste0("Q",1:K), len=ncol(Qbar) )
  if(save_Q_file != ""){
    write.csv(Qbar, file.path( save_Q_file ))
  } else write.csv(Qbar, file.path( "Q-to_graph.csv" ))
  
  if(min.group.size.show > 1){
    pt = table(pops)
    good_pops = names(pt)[pt >= min.group.size.show]
    keep = pops %in% good_pops
    Qbar = Qbar[keep,]
    pops = pops[keep]
    Qs = lapply(Qs, function(x) x[keep, ])
    QQ = QQ[keep, ] # recycle   #Q[keep + (1:R - 1)(length(pops)), ]
  }
  
  ## Reordering samples
  if(is.null(pop_sequence)){
    pop_int = as.integer( as.factor(pops))
  } else {
    keep = pops %in% pop_sequence
    Qbar = Qbar[keep, ]
    pops = pops[keep]
    Qs = lapply(Qs, function(x) x[keep, ])
    QQ = QQ[keep, ] # recycle   #Q[keep + (1:R - 1)(length(pops)), ]
    pop_int = as.integer( factor(pops, levels = pop_sequence))
  }
 
  ro = order(pop_int)
  Qbar = Qbar[ro, ]
  pops = pops[ro]
  Qs = lapply(Qs, function(Q)Q[ro,])
  pop_int = pop_int[ro]

  
  if(hclust.reorder){
    #Qbars = split(Qbar, pops)
    perm0 = 1:nrow(Qbar)
    for(i in unique(pop_int)){
      hcQ = Qbar[ pop_int == i, ]
      hc= hclust(dist(hcQ), method="ward.D")
      rperm = hc$order
      #Qbar[ pop_int == i, ] = Qbar[ pop_int == i, ][ rperm, ]
      perm0[ pop_int == i] = perm0[ pop_int == i ][ rperm ]
    }
    
    #hcQ = cbind(Qbar, 1e9 * pop_int)
    #tree=hclust(dist(hcQ), method="ward.D")
    #rperm1 = tree$order
    rperm1 = perm0
    Qbar = Qbar[rperm1,]
    Qs = lapply(Qs, function(Q)Q[rperm1,])
    pops = pops[rperm1]
  }
  
  # Reordering by maximum
  if( reorder.by.max){
    sumQ = Reduce(`+`, Qs , 0)
    avgQ = sumQ / R
    # Compute rankings
    ords = apply(avgQ, 1, function(x) order(x, decreasing=TRUE))
    ords = t(ords)
    imax = ords[,1]
    Qord = apply(avgQ, 1, function(v) sort(v, decreasing=TRUE))
    Qord = t(Qord)
    ro = order(pops, ords[,1], ords[,2], -Qord[,1],  -Qord[,2] )
    
    # # Actual reordering of samples
    Qbar = Qbar[ro, ]
    pops = pops[ro]
    Qs = lapply(Qs, function(Q)Q[ro,])
    #lastQ = lastQ[ro,]
  }
  #else {
  #  ro = order(pops)
  #}
  lastQ = Qs[[ length(Qs) ]] 

  #write.csv( data.frame(id=rownames(lastQ), pop=pops, lastQ), 
  write.csv(as.data.frame(pops),
             paste0("ReorderedId_", format(Sys.time(), "%b_%d_%Hh%M"),".csv"), row.names=T, quote=F)
  

  
  ### Compute spaces between groups
  PopSet = unique(pops)
  sp = rle(pops)
  space_ixs = cumsum(sp$lengths) + 1
  space_ixs = space_ixs[-length(space_ixs)]
  spaces = integer( nrow(Qbar))
  spaces[space_ixs] = spacing_between_groups
  ### Prepare graphics
  oldpar = par()  
  par(lwd=0.8)
  par(mar=par()$mar + c(4.5,0,0,0))
  ### Determine colors
  if(is.null(admix_comp_colors)){
    if(is.null(admix_color_table)){
      H = as.integer(K/2)
      admix_comp_colors = c(  hsv((1:H)/H, s=0.8, v=1),
                              hsv((1:(K-H))/(K-H), s=0.9, v=0.5) ) [order(colSums(QQ))]
    } else {
      # col_nm = guess_Q_colnames( lastQ, pops0)
      
      col_nm = guess_Q_colnames3( lastQ, pops)
      admix_comp_colors = admix_color_table[ col_nm ]
      # Add colors for unidentified populations
      if(any(is.na(admix_comp_colors))){
        theNA = which(is.na(admix_comp_colors))
        numNA = length(theNA)
        candidate_colors = hsv( (1:numNA)/numNA, s=0.4, v=0.55)
        na_names = names(admix_comp_colors)[ theNA ] 
        admix_comp_colors[ theNA ] = candidate_colors
        cat("Please choose colors for:\n", paste(na_names), "\n")
        # Remember current device
        dev_c = dev.cur()
        # and go to a new one
        dev.new()
        barplot(rep(1,numNA), col=candidate_colors, names.arg=na_names)
        # getcolors()
        dev.set(dev_c)
        
     
        
        }
    }
   
  }

  set.seed(4321)   
  x_at = barplot( t(Qbar), 
                  #horiz=FALSE, 
                  col = rep( admix_comp_colors, len=ncol(Qbar)),
                  space = spaces,
                  xaxt = "n",
                  axes = FALSE   , border=NA  #,ylab="K"
  )
  
  if(!is.null(ylabels)){
    axis(2, at=-0.5+(1:R), labels=ylabels, tick=F, las=2, cex.axis=0.8)  
  } else {
    # If we use different Ks, plot them on y-axis
    if(exists("Ks")) 
      if(length(unique(Ks)) > 1){ 
      axis(2, at=-0.5+(1:R), labels=paste( Ks), tick=F, las=3)    
    }  
  }

  ## Add population labels
  pops_at = sapply(1:length(PopSet), function(i){
    ii = which(PopSet[[i]] == pops)
    median(x_at[ii])
  })
  axis(1, at=pops_at, labels=PopSet, tick=FALSE, las=2, cex.axis=0.7)
  invisible(Qbar)
}


############################################################
#' Load Q files into a list
#' Input : character vector of file paths
#' Output: a list of Q matrices ordered by number of columns
allQs <- function(files){
  if( length(files) < 1) stop("No files given.\n")
  Qs = list()
  for(i in seq(along=files)){
    #cat("Importing ", files[[i]], "..\n" )
    Qs[[i]] = read.table(files[[i]] )
  }
  Ks = sapply(Qs, ncol)
  o = order(Ks)
  Qs = Qs[ o ]
  attr(Qs, "Ks") = Ks[ o ]
  attr(Qs, "files") = files[ o ]
  Qs
}


##############################################
#' Make a single "wide" matrix out of a list of Q matrices by cbind
list_to_mat <- function(Qs){
  Ks = sapply(Qs, ncol)
  N_ind=  nrow(Qs[[1]])
  QQ = unlist(Qs)
  QQ = matrix(QQ, nrow= N_ind)
  stopifnot(all(QQ[,1] == Qs[[1]][,1] ))
  colnames(QQ) = unlist( lapply(Ks, function(i) paste(i, 1:i, sep="_")))
  list(Q=QQ, K=Ks)
}

####################################################
#' Make a "long" matrix out of a list of Qs by adding zero columns and stacking up
#' @param Qs List of matrices with same number of rows but possibly different ncol
fill_columns <- function(Qs){
  K = sapply(Qs, ncol)
  I = nrow(Qs[[1]])
  maxK = max(K)
  for(i in seq(along=Qs)){
    Qs[[i]] = as.matrix(Qs[[i]])
    if(K[[i]] < maxK ){
      Qs[[i]] = cbind( Qs[[i]], matrix(0, nrow=I, ncol = maxK - K[[i]]) )
    }
    colnames(Qs)<- NULL
  }
  Qs
}

######################################################

#' Make a CLUMPP indfile
#' Qs: list of Q matrices
#' indfile_name: path to output CLUMPP indfile
#' sample_ids: sample ids to put in the indfile
make_CLUMPP_indfile <- function(Qs, 
                                indfile.name="indfile.txt"
                                #,sample_ids=NULL
                                ){
  # We expect a list, but if it's only one matrix, we will convert it to a list
  if( !is.null(dim(Qs)) ) Qs  = list(Qs)
  R = length(Qs)
  cat("Creating CLUMPP indfile ", indfile.name, "\n")
  cat("\tK=", ncol(Qs[[length(Qs) ]]), "  R=", R,"\n")
  I = nrow(Qs[[1]])
  # 
  sample_ids=1:I
  #if(is.null(sample_ids)){
  #  sample_ids = 1:I
  #} else {
  #  stop("create indfile: bad usage\n")
  #  if(class(sample_ids) != "character" ){
      sample_ids = as.character(sample_ids)
  #  }
  #}
  
  ## Rbind all Q matrices
  # If Qs have different K values, need to fill columns by zeroes first
  Qs = fill_columns(Qs)
  allQ = do.call(rbind, Qs)
  
  
  inds = rep(sample_ids, R)
  C1 =       1:(R * I)
  C3 = rep( "(x)", R*I)
  C4 = rep(1,      R*I)
  C5 = rep(":",    R*I)
  # DEBUG: if the lengths are not equal, recover() to see what's going on
  if( length(sample_ids) != I ){
    cat("There is some inconsistency in the supplied data. (Possible reason: not all samples are assigned a population).\n")
    recover()
  }
  df = cbind(C1, inds, C3,C4,C5, allQ)
  write.table(df, file=indfile.name, quote=F, col.names=F, row.names=F)
  invisible(df)
}

######################################################

#' Make a CLUMPP paramfile
#' @param paramfile: path to output paramfile
#' @param repeats: a CLUMPP parameter
make_CLUMPP_paramfile <- function(Qs, paramfile, repeats=2000){
  R = length(Qs)   # number of runs
  C = nrow(Qs[[1]])   # number of individuals
  K = max( sapply(Qs, ncol)) # number of clusters = max number of clusters among the runs
  
  cat("DATATYPE 0 # individual
      INDFILE indfile.txt
      OUTFILE outfile.txt
      MISCFILE miscfile.txt
      C ", C, "
      K ", K, "
      R ", R, "
      
      M 3  # 1: FullSearch 2:Greedy 3:LargeKGreedy
      W 0
      S 1
      
      GREEDY_OPTION 2 # 1:all orders of runs 2: some number 3: prespecified orders
      REPEATS ", repeats, "
      
      OVERRIDE_WARNINGS 1
      
      PRINT_RANDOM_INPUTORDER 0
      PRINT_PERMUTED_DATA 1  # 0,1,2   1-all Q matrices in one file, 2-separate files for distruct
      #EVERY_PERMFILE permQ_
PERMUTED_DATAFILE permutedQ.txt
ORDER_BY_RUN 0
PRINT_EVERY_PERM 0
      ", file=paramfile )
  
}


######################################################
loadPermutedQfile <- function(permQ_file){
  clumpQ = read.table( permQ_file )
  clumpQ = clumpQ[, -c(1,3,5)]    # remove dummy columns
  names(clumpQ)= c("Ind", "Pop", paste("C", 1:(ncol(clumpQ)-2), sep="_" ))
  clumpQ
}

######################################################
permuted2list <- function(permQ_file, Ks=NULL){
  clumpQ = loadPermutedQfile(permQ_file)
  pQ = clumpQ
  I = length(unique(clumpQ[,1]))    # number of individuals
  R = nrow(clumpQ)/I   # number of runs
  
  QQ = clumpQ[,-c(1,2) ]   # all Q coefs
  byRun = rep(1:R, each=I)  # column indicating the run number
  Qs = split(QQ, byRun)
  Qs = lapply(Qs, as.matrix)
  Qs
}

######################################################
permuted2array <- function(clumpQ, Ks){
  pQ = clumpQ
  I = length(unique(clumpQ[,1]))    # number of individuals
  R = nrow(clumpQ)/I   # number of runs
  
  QQ = clumpQ[,-c(1,2) ]   # all Q coefs
  QQ = as.matrix(QQ)
  K = ncol(QQ)
  byRun = rep(1:R, each=I)  # column indicating the run number
  Qs = split(QQ, byRun)
  Qs = lapply(Qs, as.matrix)
  aQ = array( NA_real_, dim = c(I, K,  R))      #length(QQ)/(I*K)) )
  for(i in 1:R){
    aQ[,,i] = Qs[[i]]
  }  
  aQ
}

######################################################
array2permQ <- function(aQ){
  t( matrix( aperm(aQ,c(2,1,3)), nr=dim(aQ)[[2]]))
}

######################################################
writeArrayQ <- function(aQ, file){
  f = open(file, "w+")
  QQ = array2permQ(aQ)
  write.table(QQ, file, quote=F, row.names=F, col.names=F )
}

######################################################
#' Distance between two Q matrices  - sum of absolute differences
distanceBetweenRuns <- function(Q1, Q2){
  # Q1, Q2 assumed after CLUMPP reordering
  d = abs(Q1 - Q2)
  d = as.matrix(d)
  sum(d)
}

######################################################
#' Number of outlier samples for each run
#' @param aQ An array
#' @return Vector 
numOutliersPerRun <- function(aQ){
  I  = dim(aQ)[[1]]
  K = dim(aQ)[[2]]
  R = dim(aQ)[[3]]
  if(R<2) return( integer(R) )
  N_out = integer(R)
  # Loop over samples
  # In each sample find the "central" cluster and
  for(i in 1:I){
    sample = aQ[i,, ,drop=TRUE]
    # Find median value for each component across runs
    meds = apply(sample, 1, median)
    # so meds is a vector of what we think as approximately true Q values for the sample
    
    # Define the error 
    dmed = abs(sample - meds)
    err = colMeans(dmed) # sum sq median dev error for each run
    
    #Threshold:
    h = median(err) + 6 * IQR(err)   # too conservative ?
    #h = min(0.1, h)
    
    if(i %% 5000 == 0){
      par(ask=TRUE) 
      plot(err, main=i)
      abline(h = h)			
    }
    aberrant = err > h
    # Count the number of aberrant samples per run
    N_out[ which(aberrant) ] = N_out[ which( aberrant) ] + 1
    if( length(N_out) > R ) browser()
  }
  
  N_out
}


######################################################
# 
getAverageQ <- function(aQ, outlierFraction = 0.02){
  # browser()
  R =  dim(aQ)[[3]]
  I = dim(aQ)[[ 1 ]]
  N_out = numOutliersPerRun(aQ)
  if(R > 2){
    good_runs = N_out < outlierFraction * I
    while( length(good_runs)==0){
      outlierFraction = outlierFraction + 0.05
      cat("*** getAverageQ: There are no good runs. Parameters are being adjusted: ")
      cat(" outlierFraction = ", outlierFraction, "\n")
      good_runs = N_out < outlierFraction * I
    }    
  } else good_runs = !logical(R) # all are good when R<=2

  good_Q = aQ[, , good_runs , drop=FALSE]
  av = apply(good_Q, c(1,2), mean)
  sd =  apply(good_Q, c(1,2), sd)
  ans = list(meanQ=av, sdQ=sd)
  attr(ans, "runs")=which(good_runs)
  attr(ans, "N_out") = N_out
  ans
}


######################################################
readSampleIDfile = function(sample_id_file, expectedNumSamples){
  sidf = read.table(sample_id_file, stringsAsFactors=FALSE, colClasses="character")
  I = expectedNumSamples
  if(nrow(sidf) != I){
    if(nrow(sidf)== I+1){
      sidf = sidf[-1, ] # skip first row
    } else {
      stop("Sample ID file has ", nrow(sidf), "rows, but Q matrices have ", I , "rows\n")
    }
  }
  if(length(sidf)==1){
    sample_id = sidf[[1]]
  } else {
    sidf = sidf[,1:2]
    ii = apply(sidf,2, function(v) length(unique(v)) )
    sample_id = sidf[[ which.max( rev(ii)) ]]
  }
  as.character(sample_id)
}


########################################################
#' Comparison across K
#' 
#compareDifferentK <- 
compare_runs_different_K <-
  function(input_dir, 
           Q_pattern="[.]+.*Q$", 
           output_dir=input_dir, 
           sample_id_file=NULL, 
           pop_file = NULL,
                recalculate.permQ = FALSE,
                 MAKE_GRAPH=TRUE,
                GRAPH_FILE = "admixture-plot.pdf",
           admix_color_table = NULL,
           ...
){
  
  Qfiles = dir(input_dir, recursive = TRUE, pattern = Q_pattern, full.names=TRUE)
  cat("Found", length(Qfiles), "files:\n" )
  print(Qfiles)
  
  Qs = allQs(Qfiles)
  if(length(Qs)==0){
    stop("Could not load Q matrices.\n")
  }
  I = nrow(Qs[[1]])
  
  ### Sample ids and pop info
  sample_id = 1:I
  if(!is.null(sample_id_file)){
    sample_id = readSampleIDfile(sample_id_file, expectedNumSamples = I)
    fam = read.table(sample_id_file, stringsAsFactors=F, colClasses="character")
  }
  
  if(!is.null(pop_file)){
    sample_pop_info = read.table(pop_file, stringsAsFactors=FALSE)
    sample_pop_info[, 1] = gsub(" ", "_", sample_pop_info[, 1] )
    #pop_assn = sample_pop_info[,3]
    #names(pop_assn) = sample_pop_info[,1]  
    fampop = merge(fam[,1:2], sample_pop_info, by=1, sort=FALSE )
    if(nrow(fampop)>0){
      sample_id = fampop[,2]
      pop_assn = fampop[,3]
      names(pop_assn) = sample_id
    } else {
      stop("Could not find any corresponidng entries between popfile and fam file! Please check yor popfile.\n")
    }
  } else {
    sample_pop_info = data.frame(id=sample_id, pop=rep("",length(sample_id)),  stringsAsFactors=F)
    pop_assn = sample_pop_info[,2]
    names(pop_assn) = sample_pop_info[,1]  
  }

  samples_with_pop = which(sample_id %in% names(pop_assn))
  
  if( length(samples_with_pop) < length(sample_id )){
    warning(paste( length(sample_id)-length(samples_with_pop), "samples have no population assignment. " ))
  }
  #pops = rep("", length(sample_id))
  #pops[samples_with_pop] =   pop_assn[ sample_id[samples_with_pop]  ]

  pops = pop_assn[ fam[[1]] ]
  pops[ is.na(pops) ] = ""
  
  #### 
  CLUMPP_dir = file.path(output_dir, "clumpp")
  dir.create(CLUMPP_dir, rec=TRUE, showW = FALSE)
  
  make_CLUMPP_indfile(Qs, indfile.name = file.path(CLUMPP_dir, "indfile.txt")
                      #, sample_ids = NULL
                      )
  make_CLUMPP_paramfile(Qs, paramfile = file.path(CLUMPP_dir, "paramfile"), repeats = 3000)
  
  permQfile = file.path(CLUMPP_dir, "permutedQ.txt")
  
  # Creating CLUMPP command file and running CLUMPP
  if(!exists("CLUMPP_exe")){
    warning("The variable CLUMPP_exe should be set in the config.R file. Skipping running CLUMPP.")
    cat("Please run CLUMPP in the directory ", CLUMPP_dir, "\n")
  } else {
    clumppCommandFile = file.path(output_dir, "runClummpForEachK.sh")
    #sink(file=clumppCommandFile)
    file.create(clumppCommandFile)
    
    if( file.exists(permQfile) && ! recalculate.permQ ){
      cat("Skipping CLUMPP...\n")
    } else {
      cat("Running CLUMPP in the directory ", CLUMPP_dir,  "\n")
      command = paste("cd ", CLUMPP_dir, "; ",  CLUMPP_exe, " 2> clumpp.log 1> clumpp.out " )
      cat(command,"\n",  file = clumppCommandFile, append = TRUE)
      tryCatch( system(command),
                error = function(e){
                  msg = paste0("Command\n", command, "\nfailed to run\n")
                  cat(msg)
                } ) 
    }    
  } # end running CLUMPP
  
  QQ = list_to_mat(Qs)

  # QQ is a cbind of all Q matrices. dim = 3023 * 145 
  
  ylabels = gsub("average.mode-", "", file_path_sans_ext(basename(attr(Qs, "files"))))
  ylabels = gsub(".*[.]", "", ylabels)
  library(stringr)
  K1 = str_match(file_path_sans_ext(basename(attr(Qs, "files"))), "K([[:digit:]]+)" )[,2]
  K2 = str_match(file_path_sans_ext(basename(attr(Qs, "files"))), "([[:digit:]]+)[.]Qmode" )[,2]
  K1[ is.na(K1)] = K2[is.na(K1)]
  K_int = sapply(Qs, ncol)
  ylabels = paste0("K", K_int)
  ylabels = make.unique(ylabels)
  cat("Note K labels:\n", paste(ylabels), "\n")
  qqTable = cbind.data.frame( id = fam[[1]], subpop_v3 = pops, QQ$Q, stringsAsFactors=F)
  write.csv(qqTable, file.path(output_dir, "Qtable.csv"), row.names=F)
  write.table(qqTable, file.path(output_dir, "Qtable.txt"), row.names=F)
  
  tryCatch({
    GRAPH_DIR  = file.path(output_dir, "Graphs")
    dir.create(GRAPH_DIR, rec=TRUE, showWarn=FALSE)
    graph_file = file.path(GRAPH_DIR, GRAPH_FILE)
    pdf(file = graph_file, width = 12, height = 8)	
    graph_allQ( permQfile, Ks = QQ$K, pops = pops, 
               reorder.by.max = F, 
               hclust.reorder = TRUE,
               #spacing_between_groups=10, 
               ylabels = ylabels,
               min.group.size.show=1,
               admix_color_table = admix_color_table,
               save_Q_file = file.path(GRAPH_DIR, "Q-to_graph.csv"),
               ...)
    dev.off()
    cat("Created file ", graph_file, "\n")
    return(graph_file)
  }, error = function(e){
    graphics.off()
    e
  } , 
  finally = {}
  )
  
}

#########################################################

# Computes Q modes for different runs with same value of K:
# 1) For each K , collect Q matrices from different runs that have K 
# 2) prepare and run CLUMPP
# 3) find modes for Q and remove aberrant runs
# 4) graph all runs and "good" runs

#' Compute Q modes for different runs with same value of K:
#'
#' @param sample_id_file FAM file that was used when running ADMIXTURE
#' @param pop_file File with population information (formatted FID IID POP )
#' 

  compare_runs_same_K <-
  function(    input_dir, # will look for Q files in all subdirs of this directory
               Q_pattern = "Q$", # Q file pattern    # "[.]+.*Q$", 
               output_dir=input_dir, 
               sample_id_file=NULL,
               pop_file = NULL, 
               recalculate.permQ = FALSE,
               reorder.runs = "none",
               clumpp.repeats=1500,
               run_clust_method="igraph"
  ){
    Qfiles = dir(input_dir, recursive = TRUE, pattern = Q_pattern, full.names=TRUE)
    #Qfiles_dirs = dirname(Qfiles)
    
    # Loading Q files
    Qs = allQs(Qfiles)
    if(length(Qs)==0){
      stop("Could not load Q matrices.\n")
    }
    
    I = nrow(Qs[[1]])
    #Loading sample IDs
    sample_id = 1:I
    if(!is.null(sample_id_file)){
      sample_id = readSampleIDfile(sample_id_file, expectedNumSamples = I)
      fam = read.table(sample_id_file, stringsAsFactors=F, colClasses="character")
    }
    
    if(!is.null(pop_file)){
      sample_pop_info = read.table(pop_file, stringsAsFactors=FALSE)
      if(all(sample_pop_info[,1]==sample_pop_info[,2])){
        sample_pop_info[[2]]=NULL
      }
      sample_pop_info[, 1] = gsub(" ", "_", sample_pop_info[, 1] )
      #pop_assn = sample_pop_info[,3]
      #names(pop_assn) = sample_pop_info[,1]  
      fampop = merge(fam[,1:2], sample_pop_info, by=c(1:(ncol(sample_pop_info)-1)), sort=FALSE )
      if(nrow(fampop)>0){
        sample_id = fampop[,2]
        pop_assn = fampop[,3]
        names(pop_assn) = sample_id
      } else {
        
        stop("Could not find any corresponidng entries between popfile and fam file! Please check yor popfile.\n")
      }
    } else {  # no pop file
      sample_pop_info = data.frame(id=sample_id, pop=rep("",length(sample_id)),  stringsAsFactors=F)
      pop_assn = sample_pop_info[,2]
      names(pop_assn) = sample_pop_info[,1]  
    }
    if(F){
      samples_with_pop = which(sample_id %in% names(pop_assn))
      
      if( length(samples_with_pop) < length(sample_id )){
        warning(paste( length(sample_id)-length(samples_with_pop), "samples have no population assignment. " ))
      }
      pops = rep("", length(sample_id))
      pops[samples_with_pop] =   pop_assn[ sample_id[samples_with_pop]  ]
      pops2= pops      
    }

    # try this:
    pops = pop_assn[ fam[[1]] ]
    pops[ is.na(pops) ] = ""
    #if( any(pops2 != pops) )browser()
    
    
    
    # Creating CLUMPP directories
    per_K_dir_name <- function(output_dir, K) file.path(output_dir, paste0("clumpp/clumpp-K", K))
    Klist = attr(Qs, "Ks")
    uK = unique(Klist)
    
    for(K in uK){
      per_K_dir =  per_K_dir_name(output_dir, K)
      dir.create( per_K_dir, recursive = TRUE, showWarn=FALSE)
      cat(" Created directory ", per_K_dir, "\n")
      
      Qs_perK = Qs[ Klist %in% K ]
      Qs_perK_files = attr(Qs, "files")[Klist %in% K  ]
      Qs_dirs = dirname(Qs_perK_files)
      adm_out_files = file.path( Qs_dirs, paste0( K, ".out" ))
      
      indfile.name=file.path(per_K_dir, "indfile.txt")
      if(!file.exists(indfile.name)){
        make_CLUMPP_indfile(Qs_perK, 
                            indfile.name = indfile.name
                            #  ,sample_ids = 1:length(sample_id) 
        )        
      } else {
        cat("Indfile ", indfile.name, " already exists, skipping.\n")
      }

      make_CLUMPP_paramfile(Qs_perK, 
                            paramfile = file.path(per_K_dir, "paramfile"),
                            repeats=clumpp.repeats)
      
    }
    
    # Creating CLUMPP command file and running CLUMPP
    if(!exists("CLUMPP_exe")){
      warning("The variable CLUMPP_exe should be set in the config.R file. Skipping running CLUMPP.")
    } else {
      clumppCommandFile = file.path(output_dir, "runClummpForEachK.sh")
      #sink(file=clumppCommandFile)
      file.create(clumppCommandFile)
      for(K in uK){
        cat("K=", K, "\n")
        per_K_dir =  per_K_dir_name(output_dir, K)
        # First check if the permQfile exists:
        #permQfile = dir(per_K_dir, patt="permutedQ.txt", full.names=TRUE)
        permQfile = file.path(per_K_dir, "permutedQ.txt")
        if( file.exists(permQfile) && ! recalculate.permQ ){
          cat("Skipping CLUMPP...\n")
          next
        }
        
        cat("Running CLUMPP for K=", K, "in the directory ", per_K_dir,  "\n")
        command = paste("cd ", per_K_dir, "; ",  CLUMPP_exe, " 2> clumpp.log 1> clumpp.out " )
        cat(command,"\n",  file = clumppCommandFile, append = TRUE)
        tryCatch( system(command),
                  error = function(e){
                    msg = paste0("Command\n", command, "\nfailed to run\n")
                    cat(msg)
                  } )
      }
      #system( paste0( "bash ") )
    }
    

    
    GRAPH_DIR = file.path(output_dir, "Graphs")
    dir.create(GRAPH_DIR, showW = FALSE, rec= TRUE)
    
    # Computing averages / modes
    AVG_Q_DIR = file.path(output_dir, "Best_Q")
    dir.create(AVG_Q_DIR, showW = FALSE, rec=TRUE)
    
    cat("Computing Q modes \n")
    for(K in uK){
      cat("###  K=", K, " ### \n")
      per_K_dir =  per_K_dir_name(output_dir, K)
      
      permQfile = dir(per_K_dir, patt="permutedQ.txt", full.names=TRUE)
      if(length(permQfile)==0){
        cat("Warning: No permQfile for K=", K, " - skipping.\n")
        next
      }
      permQ = loadPermutedQfile( permQfile)
      aQ = permuted2array(permQ)
      
      ##
      Qs_perK = Qs[ Klist %in% K ]
      Qs_perK_files = attr(Qs, "files")[Klist %in% K  ]
      Qs_dirs = dirname(Qs_perK_files)
      adm_out_files = file.path( Qs_dirs, paste0( K, ".out" ))
      # Extract CV error info
      
      WHICH_K_DIR = file.path(output_dir, "Which_K")
      if(reorder.runs != "none"){
        WHICH_K_DIR = file.path(output_dir, paste0("Runs-", reorder.runs))
        dir.create(WHICH_K_DIR, showW = FALSE, rec= TRUE)
      }
      CV_file = file.path( WHICH_K_DIR, paste0("cv-error-K",K,".txt") )
      loglike_file = file.path( WHICH_K_DIR, paste0("loglike-K",K,".txt") )
      

      
      if(reorder.runs == "cv"){
        # Extract info from .out files
        cv_info = system(paste0(
          "find ", input_dir , " -name '", K, ".out' | xargs grep CV | sed 's/:CV//g' | cut -f1,4 -d' ' | sort -n -k1  > ", CV_file ), 
          intern = F)  
        tryCatch({     
          
          if( file.exists(CV_file)){
            cv_info = read.table(CV_file, stringsAsFactors=FALSE)
            cv_dirs = dirname(cv_info[,1])
            rownames(cv_info) = cv_dirs
            # What if some Q matrices lack .out files ?
            with_cv = (Qs_dirs  %in% cv_dirs)
            #cv_info = cv_info[ rownames(cv_info) %in% Qs_dirs, ]
            cv_info = cv_info[ Qs_dirs[with_cv], ]
            o_cv = order(cv_info[,2])
            new_order = c( which(!with_cv),  which(with_cv)[o_cv]  )   
          } else {
            stop("Could not find CV error file")
          }
          aQ = aQ[,,new_order]
        }, error=function(e){
          warning("Could not reorder by CV error")
          e
        })
      }
      
      if(reorder.runs == "loglike"){
        
        loglike_info = system(paste0(
          "find ", input_dir , " -name '", K, ".out' | xargs grep '^Loglikelihood' | sed s/:Loglikelihood://g  > ",   loglike_file
        ),  intern = F)
        
        tryCatch({
          
          if( file.exists(loglike_file)){
            LL_info = read.table(loglike_file, stringsAsFactors=FALSE)
            LL_dirs = dirname(LL_info[,1])   
            rownames(LL_info) = LL_dirs
            # What if some Q matrices lack .out files ?
            with_LL = (Qs_dirs %in% LL_dirs)
            LL_info = LL_info[ Qs_dirs[with_LL], ]
            o_LL = order(LL_info[,2])
            new_order_LL = c( which(!with_LL),  which(with_LL)[o_LL])
          } else {
            stop("Could not find loglike file.")
          }
          
          aQ = aQ[,,new_order_LL]
        }, error=function(e){
          warn("Could not reorder by loglikelohood.")
          e
        })
      }
      
      
      if(reorder.runs == "hclust"){
        HeightCutoffFraction = 0.29
        Runs = apply(aQ, 3, as.numeric)
        # DRuns = dist(t(Runs))
        #hc = hclust(DRuns, method="ward")
        
        
        #heightCutoff=  HeightCutoffFraction * max(hc$height)
        #heightCutoff = mean(heightCutoff, quantile(hc$height, 0.8) )
        #cls = cutree(hc, h= heightCutoff)
        
        rcls = clusterRuns(aQ, runClustMethod=run_clust_method)
        
        result = rcls[[ run_clust_method ]]
        cls = result$cls
        
        ord =  result$order
        hc = rcls$hc
        used_runs = which( result$used_runs )
        
        
        nR = length(Qs_perK_files)
        ## if(length(cls)!= nR) browser() # not checking anymore, because of outlier runs
        run_info = data.frame(n=1:nR,  file=Qs_perK_files, cls=integer(nR))
        run_info[ used_runs, "cls"] = cls
        if(!is.null(ord)){
          run_info[ used_runs, ]  = run_info[ used_runs[ ord ], ]  
        }
        write.table(run_info, 
                    file.path(WHICH_K_DIR, paste0("run_info-",run_clust_method,"-K", K,".txt")), 
                    quote=F, row.names=F, sep="\t")
        
        tryCatch({
          pdf(file= file.path(WHICH_K_DIR, paste0("runs-dendro-K", K,".pdf")) )
          plot(hc, sub='')
          dev.off()          
        }, error = function(e) {cat("Error, cannot plot the dendrogram of runs.\n")} 
        , finally = dev.off
        )

        
        #tryCatch({
        #  pdf(file= file.path(WHICH_K_DIR, paste0("runs-pam-K", K,".pdf")) )
        #  plot(rcls$pamRes, which=1)
        #  plot(rcls$pamRes, which=2)
        #  dev.off()
        #}, error = function(e){ dev.off() } )

        aQ = aQ[ ,, used_runs]
        if(!is.null(ord)){
          aQ = aQ[ , ,ord ]
        }
      }
     
      Qmodes_dir = file.path(output_dir, paste0("Q-modes_", run_clust_method))
      Qmodes_K_dir = file.path(Qmodes_dir, paste0("K", K))
      Qmodes_K_dir_small = file.path(Qmodes_K_dir, "small")
      
      dir.create(Qmodes_K_dir, rec=T, showW=F)
      dir.create(Qmodes_K_dir_small, rec=T, showW=F)
      
                               
      jpeg(file=file.path(Qmodes_dir, "", paste0("all_runs.K",K,".png")), 
           width=800, height=800)
      graph_allQ(aQ, Ks=K, pops=pops )
      dev.off()
      
     
      
      for(m in unique(cls)){
        memb = which(m == cls)
        cat(".. Mode ", m, " (", length(memb) ,"runs )... \n")
        aQm = aQ[ ,,memb , drop=F] # get runs of this cluster
        avgQ = getAverageQ(aQm)
        if( dim(aQm)[[3]] > 2){
          good_aQm = aQm[ , ,attr(avgQ, "runs"), drop=F]        
        } else good_aQm = aQm
       
        #browser()
        MGRAPH_DIR = file.path(Qmodes_K_dir, paste0("run-cluster-",m, "_", length(memb)))
        if(length(memb)<2){
          MGRAPH_DIR = file.path(Qmodes_K_dir, "small", paste0("run-cluster-",m, "_", length(memb)))
        }
        dir.create(MGRAPH_DIR, rec=T, showW=F)
        
        pdf(file=file.path(MGRAPH_DIR, paste0("goodruns-m",m, ".K",K,".pdf")))
        graph_allQ(good_aQm, Ks=K, pops=pops )
        dev.off()
        
        pdf(file=file.path(MGRAPH_DIR, paste0("all_runs-m",m, ".K",K,".pdf")))
        graph_allQ(aQm, Ks=K, pops=pops )
        dev.off()
        
        #pdf(file=file.path(MGRAPH_DIR, paste0("avgGoodRuns-m",m, ".K",K,".pdf")))
        #graph_allQ(avgQ, Ks=K, pops=pops )
        #dev.off()
        
        # Write avg Q matrix
        # avg_Q_file1 = file.path(per_K_dir, paste0( "average.mode-",m, ".K", K, ".Q"))
        avg_Q_file2 = file.path( MGRAPH_DIR, paste0( "average.mode-",m, ".K", K, ".Qmode"))
        
        #write.table(avgQ$meanQ, avg_Q_file1 , row.names=F, col.names=F, quote=F)
        write.table(avgQ$meanQ, avg_Q_file2 , row.names=F, col.names=F, quote=F)
        
        # Save the array as well
        
        saveRDS( aQm, file=file.path( MGRAPH_DIR, paste0( "allruns-",m, ".K", K, ".rds") ) )
        
      }

    }
    invisible(Qs)
  }


#################################
### Run clustering  #####
#' 
clusterRuns = function(aQ, runClustMethod="pam", 
                       HeightCutoffFraction=0.25,
                       removeOutlierRuns=TRUE){
  mQ = apply(aQ, 3, as.numeric)
  R = dim(aQ)[[3]]
  K = dim(aQ)[[2]]
  max_clust = as.integer(R^0.5 + 1)
  
  dm = dist( t(mQ), method="manhattan")
 
  # Outlier removal
  if(removeOutlierRuns){
    sdm = apply(as.matrix(dm),2,sort)
    out1 = sdm[2,] > quantile(sdm[2,],0.9) 
    # distance to second closest neighbour 
    # snd = secondNeibDist = colMeans( sdm[2:3,] )
    # outlierRuns = (snd > max( median(snd) + 6/1.35*snd, quantile(snd, 0.8) ) )# -beyod the "normal" IQR but not more than 20% of observations
    acc_runs = ! out1      # accepted runs
    if(sum(acc_runs) < 0.25*R){
      cat("Outliers not removed.\n ")
      acc_runs = !logical(R)    # all true
    }
  }

  dm0 = dm
  dm = as.dist( as.matrix(dm)[ acc_runs, acc_runs ] )
  
  # 1. Hier clustering
  hc = hclust(dm, method="single") 
  heightCutoff=  HeightCutoffFraction * max(hc$height)
  #heightCutoff = mean(heightCutoff, quantile(hc$height, 0.8) )
  cls = cutree(hc, h = heightCutoff)

  if( length(unique(cls))> max_clust){
    cls = cutree(hc, k=max_clust)
  }
  
  hier = list(
    order=hc$order, 
    #hc=hc,
    cls = cls,
    aQ = aQ[,,acc_runs, drop=F],
    used_runs = acc_runs
  )
  
  ############# IGRAPH
  S = 60 # samples to misclassify
  clusterSeparDistance = 2*S + 0.1*(2*S)*K
  
  admD = matrix(NA, nc=R, nr=R)
  admD[lower.tri(admD)] = apply( combn(R,2),2, function(ij){
    i = ij[1]; j=ij[2]
    dr = abs(aQ[,,i] - aQ[,,j])
    admD[j,i] = admD[i,j] = sum( rowMaxs(dr) > 0.5 )
  })
  #diag(admD) = 0
  cSD = 50
  g = graph.adjacency(adjmatrix=as.matrix(admD) < cSD ,diag=FALSE, mode="lower")
  cc = clusters(g)
  igraph = list(cls = cc$membership, 
                used_runs = !logical(R), 
                graph=g)
  
  ### PAM
  if(F){
    tryCatch( {
      pamk.best = pamk(dm)
      pamRes = pam( dm , k=pamk.best$nc , keep.diss=T)   
      PAM = list(pam = pamRes) } , error=function(e){ cat("Error in pam , skipping pam.\n")}  )
  }
  
#   if(F){
#     pamResult=list()
#     for(k in 2:(R-1)){
#       pamResult[[k]] = pam(dm, k )  # medioids, clustering
#     }
#   }

  ans=list()
  ans$hier = hier
  ans$hc = hc
  #if(exists("pamRes")) ans$pamRes= pamRes
  ans$igraph = igraph

  ans
}


#### END #####



