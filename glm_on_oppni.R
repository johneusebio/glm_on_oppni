# TODO LIST ----------------------------------------------------------------------------------------
# TODO parse.split_info.ons3col: allow to write the 3 columns to a diff text file for each cond
# 

# main wrapper

glm_on_oppni <- function(oppni_path, output_dir, pipe="FIX", sNorm=T, contrast=NULL) {
  oppni.data  <- import_oppni(oppni_path)
  oppni.input <- oppni.data$input_data
  oppni.pipe  <- oppni.data$pipe_config
  
  output_dir.fsf <- file.path(output_dir, "feat_results")
  output_dir.ons <- file.path(output_dir, "onsets")
  
  if (is.windows()) {
    oppni.input[] <- lapply(oppni.input, FUN = function(X) {gsub(x=X, pattern="/global/home/hpc3586/", replacement="Z://", fixed=T)})
  }
  
  # loop through each row
  for (op.row in nrow(oppni.input)) {
    tmp.row     <- oppni.input[op.row,]
    tmp.name    <- basename(tmp.row$OUT)
    tmp.pattern <- paste0(tmp.name, "_*")
    
    oppni.task <- oppni_task.subj(oppni_row    = oppni.input[op.row,], 
                                  onsets.dir   = output_dir.ons, 
                                  file_pattern = tmp.pattern)
    oppni_fsf.subj(oppni_task = oppni.task, oppni_pipe_config = oppni.pipe, 
                   output_dir = output_dir.fsf, contrast = contrast, oppni_row = tmp.row, oppni_pipe = pipe)
  }
}

# importing oppni data -----------------------------------------------------------------------------

import_oppni <- function(oppni_path) {
  
  # check if directories exist
  if (!dir.exists(oppni_path)) {
    stop("Specified OPPNI output directory does not exist.")
  }
  if (!dir.exists(file.path(oppni_path,"optimization_results"))) {
    stop("optimization_results directory does not exist.")
  }
  
  #check if key files exist
  if (!file.exists(file.path(oppni_path,"input_file.txt"))) {
    stop("Unable to find OPPNI input_file.txt")
  }
  
  # collect the oppni data
  oppni_input_table <- parse.oppni_input(file.path(oppni_path, "input_file.txt"))
  
  # get pipeline info
  oppni_pipe_config <- pickle_reader(pickle_reader.path = pickle_reader.path, 
                                     pkl.data = file.path(oppni_path, "pronto-options.pkl"))
  oppni_pipe_config <- oppni_pipe_config[[2]]
  
  # compile into list
  oppni_data             <- list()
  oppni_data$input_data  <- oppni_input_table
  oppni_data$pipe_config <- oppni_pipe_config
  
  return(oppni_data)
}

parse.oppni_input <- function(oppni_input_path) {
  require("plyr")
  
  con        <- file(oppni_input_path, open="r")
  line.count <- 0
  while(T) {
    line <- readLines(con, n=1)
    
    if (length(line)==0) {break}
    
    line.count <- line.count+1
    # print(line)
    line.0 <- parse.oppni_input.extract_line(line)
    line.0 <- parse.oppni_input.as_table(line.0, sort=T)
    
    parse.oppni_input.check_names(oppni_input_row  = line.0, 
                                  oppni_input_line = line, 
                                  line_num         = line.count)
    
    if (line.count==1) {
      oppni_input_table <- line.0
    } else {
      oppni_input_table <- plyr::rbind.fill(oppni_input_table, line.0)
    }
  }
  
  close(con)
  return(oppni_input_table)
}

parse.oppni_input.extract_line <- function(line) {
  line.split <- unlist(strsplit(line, split = " "))
  return(line.split)
}

parse.oppni_input.as_table <- function(line, sort=T) {
  line.tab <- lapply(X = line, FUN = parse.oppni_input.split_entry, split = "=")
  line.tab <- do.call("cbind", line.tab)
  if (sort) {
    line.tab <- line.tab[, order(names(line.tab))]
  }
  return(line.tab)
}

parse.oppni_input.split_entry <- function(line.entry, split) {
  split_entry <- unlist(strsplit(line.entry, split=split))
  entry_col   <- data.frame(split_entry[2], stringsAsFactors = F)
  
  colnames(entry_col) <- split_entry[1]
  return(entry_col)
}

parse.oppni_input.check_names <- function(oppni_input_row, oppni_input_line=NULL, line_num = NULL) {
  acceptable_names <- c("IN", "OUT", "DROP", "TASK", "PHYSIO", "STRUCT", "CUSTOMREG")
  
  if (any(!names(oppni_input_row) %in% acceptable_names)) {
    invalid_names <- names(oppni_input_row)[!names(oppni_input_row) %in% acceptable_names]
    stop("Invalid column name(s): ", invalid_names, "\n  Check row ", line_num, ": \n\n", oppni_input_line)
  }
}

# extract task information associated with the provided oppni input.txt file -----------------------

oppni_task.subj <- function(oppni_row, onsets.dir=NULL, file_pattern=NULL) {
  # extract task info
  oppni.task <- parse.split_info(split_info=oppni_row$TASK, DROP=oppni_row$DROP, 
                                 onsets.dir=onsets.dir, file_pattern=file_pattern)
  
  return(oppni.task)
}

# get task info from split_info file ---------------------------------------------------------------

parse.split_info <- function(split_info, DROP = NULL, onsets.dir=NULL, file_pattern=NULL) {
  
  # find the number of scans to drop
  if (DROP!=0) {
    DROP.bkup <- DROP
    DROP      <- regmatches(DROP, gregexpr("(?<=\\[).+?(?=\\])", DROP, perl=T))[[1]]
    DROP      <- as.numeric(unlist(strsplit(x = DROP, split = ",")))
    DROP      <- DROP[1]
  } else {
    DROP <- 0
  }
  
  if (!is.numeric(DROP)) {stop("Invalid DROP provided: ", DROP.bkup)}
  
  split_info.list      <- list()
  split_info.list$meta <- list()
  split_info.list$cond <- list()
  
  con        <- file(split_info, open="r")
  line.count <- 0
  cc <- 0
  
  while(T) {
    line <- readLines(con, n=1)
    if (length(line)==0) {break}
    line.count <- line.count + 1
    
    # check for meta-data
    line.meta_data.check <- lapply(X=c("UNIT","TR_MSEC","TYPE"), 
                                   FUN=function(X,line) {startsWith(x=line, prefix=X)}, line=line)
    line.meta_data.check <- any(unlist(line.meta_data.check))
    
    # check for a task condition
    line.condition.check <- lapply(X=c("NAME","ONSETS","DURATION"), 
                                   FUN=function(X,line) {startsWith(x=line, prefix=X)}, line=line)
    line.condition.check <- unlist(line.condition.check)
    
    if (line=="") {
      next
    } else if (line.meta_data.check) {
      tmp.meta                              <- parse.split_info.line(line)
      split_info.list$meta[names(tmp.meta)] <- tmp.meta[[names(tmp.meta)]]
    } else if (any(line.condition.check)) {
      
      tmp.cond <- parse.split_info.line(line)
      
      if (names(tmp.cond)=="NAME") {
        cc <- cc + 1
        
        split_info.list$cond[[cc]]          <- list()
        split_info.list$cond[[cc]]$NAME     <- tmp.cond[[names(tmp.cond)]]
      } else if (names(tmp.cond)=="ONSETS") {
        split_info.list$cond[[cc]]$ONSETS   <- tmp.cond[[names(tmp.cond)]]
        split_info.list$cond[[cc]]$ONSETS   <- parse.split_info.drop(
          drop_tr = DROP, 
          onsets  = split_info.list$cond[[cc]]$ONSETS, 
          tr_msec = split_info.list$meta$TR_MSEC)
        
      } else if (names(tmp.cond)=="DURATION") {
        split_info.list$cond[[cc]]$DURATION <- tmp.cond[[names(tmp.cond)]]
      }
    } else  {
      stop("Invalid text on line #", line.count, ": ", line, "\n\nSee file: '", split_info, "'")
    }
    
  }
  
  close(con)
  
  split_info.list <- parse.split_info.time2tr(split_info.list)
  split_info.list <- parse.split_info.ons2sec(split_info.list)
  split_info.list <- parse.split_info.ons3col(split_info.list)
  split_info.list <- parse.split_info.mk_onsFile(split_info.list,
                                                 onsets.dir = onsets.dir,
                                                 file_pattern = file_pattern)
  
  return(split_info.list)
}

parse.split_info.line <- function(split_info_line) {
  split_info_line <- unlist(strsplit(split_info_line, split = "="))
  
  label <- split_info_line[1]
  
  data  <- split_info_line[2]
  data  <- regmatches(data, gregexpr("(?<=\\[).+?(?=\\])", data, perl=T))[[1]]
  data  <- unlist(strsplit(data, split = ","))
  
  # try to convert to numeric. No change if there's a warning
  tryCatch(
    {
      data = as.numeric(data)
    },
    warning = function(data){return(data)},
    error   = function(data){return(data)}
  )
  
  out          <- list()
  out[[label]] <- data
  
  return(out)
}

parse.split_info.drop <- function(onsets, drop_tr, tr_msec) {
  check.varList <- c("onsets", "drop_tr", "tr_msec")
  check.numeric <- !unlist(lapply(X = list(onsets, drop_tr, tr_msec), FUN = is.numeric))
  if (any(check.numeric)) {
    stop("This variable(s) must be numeric: ", paste(check.varList[check.numeric], collapse = ", "))
  }
  
  new.ons <- onsets - (tr_msec * drop_tr)
  
  return(new.ons)
}

parse.split_info.time2tr <- function(split_info.list) {
  require("pracma")
  
  if (!is.list(split_info.list)) {stop("split_info.list must be a list.")}
  
  for (cc in 1:length(split_info.list$cond)) {
    
    if (pracma::strcmpi(split_info.list$meta$UNIT, "TR")) {
      split_info.list$cond[[cc]]$ONSETS_TR <- split_info.list$cond[[cc]]$ONSETS_TR
    } else if (pracma::strcmpi(split_info.list$meta$UNIT, "SEC")) {
      split_info.list$cond[[cc]]$ONSETS_TR <- time2tr(onsets   = split_info.list$cond[[cc]]$ONSETS * 1000, 
                                                      tr       = split_info.list$meta$TR_MSEC, 
                                                      rounding = "up")
    } else if (pracma::strcmpi(split_info.list$meta$UNIT, "MSEC")) {
      split_info.list$cond[[cc]]$ONSETS_TR <- time2tr(onsets   = split_info.list$cond[[cc]]$ONSETS, 
                                                      tr       = split_info.list$meta$TR_MSEC, 
                                                      rounding = "up")
    }
    
  }
  
  return(split_info.list)
}

parse.split_info.ons2sec <- function(split_info.list) {
  
  for (cc in 1:length(split_info.list$cond)) {
    tmp.ons.sec <- ons2sec(onsets     = split_info.list$cond[[cc]]$ONSETS, 
                           orig_units = split_info.list$meta$UNIT, 
                           tr_msec    = split_info.list$meta$TR_MSEC)
    
    split_info.list$cond[[cc]]$ONSETS_SEC <- tmp.ons.sec
  }
  
  return(split_info.list)
}

parse.split_info.ons3col <- function(split_info.list, output.file = NULL) {
  
  for (cc in 1:length(split_info.list$cond)) {
    tmp.3col <- mk.ons.3col(onsets   = split_info.list$cond[[cc]]$ONSETS_SEC, 
                            duration = split_info.list$cond[[cc]]$DURATION, 
                            input = rep(x = 1, length(split_info.list$cond[[cc]]$ONSETS_SEC)))
    
    split_info.list$cond[[cc]]$ONSETS_3COL <- tmp.3col
    
    if (!is.null(output.file)) {
      write.table(x=tmp.3col, file=output.file, append=F, sep=" ", row.names=F, col.names=F)
      
      split_info.list$cond[[cc]]$ONSETS_3COL.file <- output.file
    }
  }
  
  return(split_info.list)
}

parse.split_info.mk_onsFile <- function(split_info.list, onsets.dir, file_pattern) {
  if (!is.null(onsets.dir)) {
    if (!dir.exists(onsets.dir)) {dir.create(onsets.dir, showWarnings = T, recursive = T)}
    
    if (is.null(file_pattern)) {
      stop("No file_pattern specified.")
    } else if (!is.character(file_pattern)) {
      stop("file_pattern must be a character.")
    } else if (length(file_pattern)>1) {
      stop("You may only specify one file pattern.")
    } else if (!grepl(pattern="*", x = file_pattern, fixed = T)) {
      stop("The file pattern must contain a wildcard '*'")
    }
  }
  
  for (cc in 1:length(split_info.list$cond)) {
    file.name <- gsub(x=file_pattern, pattern="*", replacement=split_info.list$cond[[cc]]$NAME, fixed=T)
    file.name <- paste0(file.name, ".txt")
    
    abs.path <- file.path(onsets.dir, file.name)
    
    write.table(x = split_info.list$cond[[cc]]$ONSETS_3COL, 
                file = abs.path, append = F, sep = "\t", row.names = F, col.names = F, eol = "\n")
    
    split_info.list$cond[[cc]]$ONSETS_3COL_FILENAME <- abs.path
  }
  
  return(split_info.list)
}

# the main attraction - this will make the .fsf file for the GLM -----------------------------------

oppni_fsf.subj <- function(oppni_row, oppni_task, oppni_pipe, oppni_pipe_config, output_dir, nNorm, file_pattern, contrast = NULL) {
  
  # make basic single-condition contrasts if no contrast is provided
  if (is.null(contrast)) {
    contrast <- oppni_fsf.mk_con(oppni_task)
  } else if (!is.list(contrast) & all(unlist(lapply(X = contrast, FUN = is.numeric)))) {
    contrast <- list(contrast)
  } else if (!is.list(contrast)) {
    stop("contrast must be specified as a numeric vector or a list of numeric vectors.")
  }
  
  file_pattern <- basename(oppni_row$OUT)
  
  # check if the contrasks are valid 
  oppni_fsf.check_con(oppni_task=oppni_task, contrast=contrast)
  oppni_fsf.write_txt(oppni_task        = oppni_task, 
                      oppni_pipe_config = oppni_pipe_config, 
                      contrast          = contrast, 
                      output_dir        = output_dir, oppni_row = oppni_row, 
                      oppni_pipe        = oppni_pipe, 
                      oppni_sNorm = sNorm, file_pattern = file_pattern)
  
} ################################################################################################## FINISH THIS

oppni_fsf.write_txt <- function(oppni_row, oppni_task, oppni_pipe, oppni_sNorm, oppni_pipe_config, contrast, output_dir, file_pattern) {
  if (!oppni_pipe %in% c("IND","FIX","CON")) {
    stop("Invalid oppni_pipe. Select IND, FIX, or COND.")
  }
  
  fsf_dir <- file.path(output_dir, "fsf")
  if (!dir.exists(fsf_dir)) {
    dir.create(fsf_dir,showWarnings = F, recursive = T)
  }
  
  feat_dir <- file.path(output_dir, "feat_results_lvl1") 
  if (!dir.exists(feat_dir)) {
    dir.create(feat_dir, showWarnings = F, recursive = T)
  }
  
  fsf.path    <- file.path(fsf_dir, gsub(x=file_pattern, pattern="*", replacement="feat", fixed=T))
  feat.output <- file.path(feat_dir, gsub(x=file_pattern, pattern="*", replacement="lvl1"))
  
  nifti_file <- oppni_fsf.get_outputNifti(oppni_row, pipe = oppni_pipe)
  
  print(nifti_file)
  
  n_timepts  <- get.nifti.hdr(file = nifti_file, args = "-nv")
  n_voxels   <- get.nifti.hdr(file = nifti_file, args = "-n4")
  n_voxels   <- prod(as.numeric(unlist(strsplit(n_voxels, split = " "))), na.rm = T)
  
  # basic setup
  cat('set fmri(version) 6.00', file = fsf.path, append = F)
  cat('set fmri(level) 1', file=fsf.path, append=T)
  cat('set fmri(analysis) 7', file=fsf.path, append=T)
  cat('set fmri(relative_yn) 0', file=fsf.path, append=T)
  cat('set fmri(help_yn) 1', file=fsf.path, append=T)
  cat('set fmri(featwatcher_yn) 0', file=fsf.path, append=T)
  cat('set fmri(sscleanup_yn) 0', file=fsf.path, append=T)
  
  # define analysis_specific lines
  cat(paste0('set fmri(otuputdir) "', feat.output, '"'), file=fsf.path, append=T)
  cat('set fmri(tr) ', oppni_task$meta$TR_MSEC/1000, file=fsf.path, append=T)
  cat('set fmri(npts) ', ntimepts, append = T)
  cat('set fmri(ndelete) 0', file=fsf.path, append=T)
  cat('set fmri(multiple) 1', file=fsf.path, append=T)
  cat('set fmri(inputtype) 2', file=fsf.path, append=T)
  cat('set fmri(filtering_yn) 1', file=fsf.path, append=T)
  cat('set fmri(brain_thresh) 10', file=fsf.path, append=T)
  cat('set fmri(critical_z) 5.3', file=fsf.path, append=T)
  cat('set fmri(noise) 0.000000', file=fsf.path, append=T)
  cat('set fmri(noisear) 0.000000', file=fsf.path, append=T)
  cat('set fmri(mc) 0', file=fsf.path, append=T)
  cat('set fmri(sh_yn) 0', file=fsf.path, append=T)
  cat('set fmri(regunwarp_yn) 0', file=fsf.path, append=T)
  cat('set fmri(gdc) ""', file=fsf.path, append=T)
  cat('set fmri(dwell) 0.7', file=fsf.path, append=T)
  cat('set fmri(te) 35', file=fsf.path, append=T)
  cat('set fmri(signallossthresh) 10', file=fsf.path, append=T)
  cat('set fmri(unwarp_dir) y-', file=fsf.path, append=T)
  cat('set fmri(st) 0', file=fsf.path, append=T)
  cat('set fmri(st_file) ""', file=fsf.path, append=T)
  cat('set fmri(bet_yn) 0', file=fsf.path, append=T)
  cat('set fmri(smooth) 0.0', file=fsf.path, append=T)
  cat('set fmri(norm_yn) 0',file=fsf.path, append=T)
  cat('set fmri(perfsub_yn) 0', file=fsf.path, append=T)
  cat('set fmri(temphp_yn) 0', file=fsf.path, append=T)
  cat('set fmri(templp_yn) 0', file=fsf.path, append=T)
  cat('set fmri(melodic_yn) 0', file=fsf.path, append=T)
  cat('set fmri(stats_yn) 1', file=fsf.path, append=T)
  cat('set fmri(prewhiten_yn) 1', file=fsf.path, append=T)
  cat('set fmri(motionevs) 0', file=fsf.path, append=T)
  cat('set fmri(motionevsbeta) ""', file=fsf.path, append=T)
  cat('set fmri(scriptevsbeta) ""', file=fsf.path, append=T)
  cat('set fmri(robust_yn) 0', file=fsf.path, append=T)
  cat('set fmri(mixed_yn) 2', file=fsf.path, append=T)
  cat('set fmri(randomisePermutations) 5000', file=fsf.path, append=T)
  
  # parameters specific to this analysis
  cat('set fmri(evs_orig)', length(oppni_task$cond)  , file=fsf.path, append=T)
  cat('set fmri(evs_orig)', length(oppni_task$cond)*2, file=fsf.path, append=T)
  cat('set fmri(ncon_orig)', length(contrast), file=fsf.path, append=T)
  cat('set fmri(ncon_real)', length(contrast), file=fsf.path, append=T)
  cat('set fmri(nftests_orig) 0', file=fsf.path, append=T)
  cat('set fmri(nftests_real) 0', file=fsf.path, append=T)
  cat('set fmri(constcol) 0', file=fsf.path, append=T)
  cat('set fmri(poststats_yn) 1', file=fsf.path, append=T)
  cat('set fmri(threshmask) ""', file=fsf.path, append=T)
  cat('set fmri(thresh) 3', file=fsf.path, append=T)
  cat('set fmri(prob_thresh) 0.05', file=fsf.path, append=T)
  cat('set fmri(z_thresh) 2.33', file=fsf.path, append=T)
  cat('set fmri(zdisplay) 0', file=fsf.path, append=T)
  cat('set fmri(zmin) 2', file=fsf.path, append=T)
  cat('set fmri(zmax) 8', file=fsf.path, append=T)
  cat('set fmri(rendertype) 1', file=fsf.path, append=T)
  cat('set fmri(bgimage) 1', file=fsf.path, append=T)
  cat('set fmri(tsplot_yn) 1', file=fsf.path, append=T)
  cat('set fmri(reginitial_highres_yn) 0', file=fsf.path, append=T)
  cat('set fmri(reginitial_highres_search) 90', file=fsf.path, append=T)
  cat('set fmri(reginitial_highres_dof) 3', file=fsf.path, append=T)
  cat('set fmri(reghighres_yn) 0', file=fsf.path, append=T)
  cat('set fmri(reghighres_search) 90', file=fsf.path, append=T)
  cat('set fmri(reghighres_dof) BBR', file=fsf.path, append=T)
  cat('set fmri(regstandard_yn) 1', file=fsf.path, append=T)
  cat('set fmri(alternateReference_yn) 0', fsf.path, append = T)
  cat(paste0('set fmri(regstandard) "', oppni_pipe_config$reference,'"'), fsf.path, append = T)
  cat('set fmri(regstandard_search) 90', file=fsf.path, append=T)
  cat('set fmri(regstandard_dof) 3', file=fsf.path, append=T)
  cat('set fmri(regstandard_nonlinear_yn) 0', file=fsf.path, append=T)
  cat('set fmri(regstandard_nonlinear_warpres) 0', file=fsf.path, append=T)
  cat('set fmri(paradigm_hp) 100', file=fsf.path, append=T)
  cat('set fmri(totalVoxels) ', n_voxels, file=fsf.path, append=T)
  cat('set fmri(ncopeinputs) 0', file=fsf.path, append=T)
  cat(paste0('set feat_files(1) "', nifti_file,'"'), file=fsf.path, append=T)
  cat('set fmri(confoundevs) 0', file=fsf.path, append=T)
  
  # add each condition as a seperate EV
  for (cc in 1:length(oppni_task$cond)) {
    cond.name <- oppni_task$cond[[cc]]$NAME
    cond.file <- oppni_task$cond[[cc]]$ONSETS_3COL.file
    
    cat(paste0('set fmri(evtitle',cc,') "', cond.name, '"'), file=fsf.path, append=T)
    cat(paste0('set fmri(shape'  ,cc,') 3'), file=fsf.path, append=T)
    cat(paste0('set fmri(convolve'  ,cc,') 3'), file=fsf.path, append=T) # double-gamma HRF
    cat(paste0('set fmri(convolve_phase'  ,cc,') 0'), file=fsf.path, append=T) 
    cat(paste0('set fmri(tempfilt_yn'  ,cc,') 0'), file=fsf.path, append=T) 
    cat(paste0('set fmri(deriv_yn'  ,cc,') 1'), file=fsf.path, append=T) # temporal derivative 
    cat(paste0('set fmri(custom'  ,cc,') "', cond.file, '"'), file=fsf.path, append=T)
    
    # Orthogonalise EV cc wrt EV ii
    for (ii in 0:length(oppni_task$cond)) {
      cat('set fmri(ortho',cc,'.',ii,') 0', file=fsf.path, append=T)
    }
  }
  
  # contrast and F-test mode
  cat('set fmri(con_mode_old) orig', file=fsf.path, append=T)
  cat('set fmri(con_mode) orig', file=fsf.path, append=T)
  
  # add contrasts to the .fsf file
  for (cc in 1:length(contrast)) {
    con.title <- oppni_fsf.name_con(con_vector = contrast[[cc]], oppni_task = oppni_task)
    
    # Display images for contrast_real cc
    cat(paste0('set fmri(conpic_real.',cc,') 1'), file=fsf.path, append=T)
    
    # Title for contrast_real cc
    cat(paste0('set fmri(conname_real.',cc,') "',con.title,'"'), file=fsf.path, append=T)
    
    rl.con.ind <- 0
    # the real contrasts
    for (cci in 1:length(oppni_task$cond)) {
      for (cci.i in 1:2) {
        rl.con.ind <- rl.con.ind + 1
        rl.con.val <- ifelse(cci.i==1, yes = contrast[[cc]][cci], no = 0)
        
        # Real contrast_real vector cc element rl.con.ind
        cat(paste0('set fmri(con_real',cc,'.',rl.con.ind,') ',rl.con.val), file=fsf.path, append=T)
      }
    }
  }
  
  # add contrast_orig to .fsf file
  for (cc in 1:length(contrast)) {
    con.title <- oppni_fsf.name_con(con_vector = contrast[[cc]], oppni_task = oppni_task)
    
    cat(paste0('set fmri(conpic_orig.',cc,') 1'), file=fsf.path, append=T)
    cat(paste0('set fmri(conname_orig.',cc,') "',con.title,'"'), file=fsf.path, append=T)
    
    for (cc.i in 1:length(contrast[[cc]])) {
      cat(paste0('set fmri(con_orig',cc,'.',cc.i,') ',contrast[[cc]][cc.i]), file=fsf.path, append=T)
    }
  }
  
  # contrast masking
  cat('set fmri(conmask_zerothresh_yn) 0', file=fsf.path, append=T)
  for (cc in 1:length(contrast)) {
    for (dd in 1:length(contrast)) {
      if (cc==dd) {next}
      cat(paste0('set fmri(conmask',cc,'_',dd,') 0'), file=fsf.path, append=T)
    }
  }
  cat('set fmri(conmask1_1) 0', file=fsf.path, append=T) # do any contrast masking at all?
  
  # Now options that don't appear in the GUI
  cat('set fmri(alternative_mask) ""', file=fsf.path, append=T)
  cat('set fmri(init_initial_highres) ""', file=fsf.path, append=T)
  cat('set fmri(init_highres) ""', file=fsf.path, append=T)
  cat('set fmri(init_standard) ""', file=fsf.path, append=T)
  cat('set fmri(overwrite_yn) 1', file=fsf.path, append=T)
}

oppni_fsf.get_outputNifti <- function(oppni_row, pipe, sNorm=T) {
  # check if specified pipe is valid
  if (!pipe %in% c("IND", "FIX", "CON")) {stop("Invalid pipe option. Select IND, FIX, or COND.")}
  
  # read in the data
  oppni.subj <- basename(oppni_row[,"OUT"])
  oppni.out  <- dirname(oppni_row[,"OUT"])
  
  oppni.out_dir   <- file.path(oppni.out, "optimization_results", "processed")
  oppni.out_nifti <- grep(x = list.files(oppni.out_dir), pattern = oppni.subj, value = T)
  oppni.out_nifti <- grep(x = oppni.out_nifti, pattern = pipe, value = T)
  oppni.out_nifti <- file.path(oppni.out_dir, oppni.out_nifti)
  
  if (sNorm) {
    oppni.out_nifti <- grep(x = oppni.out_nifti, pattern = "sNorm", value = T)
  } else {
    oppni.out_nifti <- oppni.out_nifti[!grepl(x = oppni.out_nifti, pattern = "sNorm")]
  }
  
  if (length(oppni.out_nifti)>1) {
    stop("More than one matching file found: ", paste(oppni.out_nifti, collapse = ", "))
  }
  
  return(oppni.out_nifti)
}

oppni_fsf.mk_con <- function(oppni_task) {
  n.cond   <- length(oppni_task$cond)
  contrast <- list()
  for (cc in 1:n.cond) {
    tmp.con        <- rep(0, n.cond)
    tmp.con[cc]    <- 1
    contrast[[cc]] <- tmp.con
  }
  
  return(contrast)
}

oppni_fsf.name_con <- function(con_vector, oppni_task) {
  if (!is.numeric(con_vector)) {
    stop("con_vector must be a numeric vector")
  } else if (!is.list(oppni_task)) {
    stop("oppni_task must be a list")
  } else if (!length(con_vector) != length(oppni_task)) {
    stop("con_vector must have the same number of entries as oppni_task$cond")
  }
  
  nm.ls <- NULL
  for (cc in 1:length(con_vector)) {
    if (con_vector[cc]==0) {next}
    
    tmp.factor <- abs(con_vector[cc])
    tmp.factor <- ifelse(test = tmp.factor==1, yes = "", no = tmp.factor)
    tmp.valence<- ifelse(test = con_vector[cc] < 0, yes = "-", no = "+")
    tmp.nm     <- oppni_task$cond[[cc]]$NAME
    tmp.nm     <- ifelse(test = tmp.factor > 1, yes = paste0("(",tmp.nm,")"), no = tmp.nm)
    
    tmp.nm <- paste0(tmp.valence, tmp.factor, tmp.nm)
    
    nm.ls  <- c(nm.ls,tmp.nm)
  }
  
  contrast.name <- paste0(nm.ls, collapse = "")
}

oppni_fsf.check_con <- function(oppni_task, contrast) {
  # check to make sure the contrasts are all valid
  con.len <- unique(unlist(lapply(X = contrast, FUN = length)))
  if (length(con.len) != 1) {
    stop("All contrasts must have the same length")
  }
  
  con.class <- !unlist(lapply(contrast, is.numeric))
  if (any(con.class)) {
    stop("All contrasts must be numeric vectors.")
  }
  
  con.all_zero <- unlist(lapply(contrast, FUN = function(X) {return(all(X==0))}))
  if (any(con.all_zero)) {
    tmp.contrast     <- contrast
    con.all_zero.ind <- which(con.all_zero)
    for (zero in con.all_zero.ind) {
      tmp.contrast[[zero]] <- NULL
    }
    contrast <- tmp.contrast
    warning("The following contrast(s) was removed for containing all zero values: ", 
            paste(con.all_zero.ind, sep = ", "))
  }
}

# NIFTI file operations ----------------------------------------------------------------------------

get.nifti.hdr <- function(file, args = NULL) {
  if (!is.character(file)) {
    stop("Filepath to nifti file, 'file', must be a character.")
  } else if (!file.exists(file)) {
    stop("The specified file does not exist: ", file)
  }
  
  
  if (!is.null(args) & !is.character(args)) {
    stop("Specified 3dinfo arguments, 'args', must be a character.")
  }
  
  hdr.info <- system(command = paste("3dinfo", args, file, sep = " "), intern = T)
  
  return(hdr.info)
}

# convert between units ----------------------------------------------------------------------------

time2tr <- function(onsets, tr_dur, rounding = "up") {
  # Onsets and TR must be in the same units
  
  if (!rounding %in% c("up", "down")) {
    stop("Rounding must be set to either 'up' or 'down'")
  }
  if (!is.numeric(onsets) | !is.numeric(tr_dur)) {
    stop("Both onsets and tr must be numeric.")
  }
  if (length(tr_dur) != 1) {
    stop("tr must be a single number.")
  }
  
  tr_ons <- onsets / tr_dur
  
  if (rounding == "up") {
    tr_ons <- ceiling(tr_ons)
  } else if (rounding == "down") {
    tr_ons <- floor(tr_ons)
  }
  
  return(tr_ons)
}

ons2sec <- function(onsets, orig_units = "msec", tr_msec = NULL) {
  require("pracma")
  
  if (pracma::strcmpi(orig_units, "SEC")) {
    ons.sec <- onsets
  } else if (pracma::strcmpi(orig_units, "MSEC")) {
    ons.sec <- onsets / 1000
  } else if (pracma::strcmpi(orig_units, "TR")) {
    if (is.null(tr_msec)) {stop("tr_msec must be specified")}
    else if (!is.numeric(tr_msec)) {stop("tr_msec must be a numeric value")}
    
    if (length(tr_msec) > 1) {
      warning("More than one tr_msec value provided. Will only use the first value.")
      tr_msec <- tr_msec[1]
    }
    
    ons.sec <- onsets * tr_msec / 1000
  } else {
    stop("orig_units must be one of the following (case-insensitive): SEC, MSEC, TR")
  }
  
  return(ons.sec)
}

mk.ons.3col <- function(onsets, duration, input) {
  #' Onsets must be in seconds.
  
  # check to make sure the input is good
  if (length(onsets) != length(duration) | length(onsets) != length(input)) {
    stop("All three variables (onsets, duration, input) must be the same length.")
  }
  
  var.names   <- c("onsets", "duration", "input")
  var.numeric <- unlist(lapply(X = list(onsets, duration, input), FUN = is.numeric))
  
  if (!any(var.numeric)) {
    stop("These variables must be numeric: ", paste(var.names[!var.numeric], collapse = ", "))
  }
  
  ons.3col <- data.frame(onsets = onsets, duration = duration, input = input)
  
  return(ons.3col)
}

# pickle-reader ------------------------------------------------------------------------------------

pickle_reader <- function(pickle_reader.path, pkl.data) {
  require("reticulate")
  
  reticulate::source_python(pickle_reader.path)
  pickle_data <- read_pickle_file(pkl.data)
  
  return(pickle_data)
}

# system functions ---------------------------------------------------------------------------------

convert.winPath <- function(winPath) {
  file.parts    <- unlist(strsplit(winPath, ":"))
  file.parts[1] <- tolower(file.parts[1])
  
  linPath <- file.path("/mnt", file.parts[1], file.parts[2], fsep = "/")
  while (as.logical(grepl(x = linPath, pattern = "//"))) {
    linPath <- gsub(pattern = "//", replacement = "/", x = linPath)
  }
  
  return(linPath)
}

is.windows <- function() {
  require("pracma")
  
  os.name    <- Sys.info()["sysname"]
  win.status <- pracma::strcmpi(s1 = os.name, "Windows")
  
  return(win.status)
}

################
# TESTING ZONE # -----------------------------------------------------------------------------------
################

# test variables -----------------------------------------------------------------------------------

if (is.windows()) {
  # windows
  oppni_path         <- "Z:\\SART_data\\oppni\\output_dir\\GO\\young\\processing_GO_sart_young_replicate_erGNB"
  pickle_reader.path <- "Z:\\JE_packages\\glm_on_oppni\\pickle_reader.py"
  template_fsf.path  <- "Z:\\JE_packages\\glm_on_oppni\\template.fsf"
} else {
  # cac
  oppni_path         <- "/global/home/hpc3586/SART_data/oppni/output_dir/GO/young/processing_GO_sart_young_replicate_erGNB"
  pickle_reader.path <- "/global/home/hpc3586/JE_packages/glm_on_oppni/pickle_reader.py"
  template_fsf.path  <- "/global/home/hpc3586/JE_packages/glm_on_oppni/template.fsf"
}

oppni_input_path <- file.path(oppni_path, "input_file.txt")
oppni_pipe_path  <- file.path(oppni_path, "pronto-options.pkl")

output_dir <- "C:/Users/enter/Desktop/test_feat"

