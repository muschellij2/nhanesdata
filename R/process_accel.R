#' Process NHANES 2003-2004 and 2005-2006 accelerometry data
#'
#' @description
#' This function processes the raw NHANES 2003-2004 and 2005-2006 accelerometry data
#' as provided by the NHANES study on the CDC website.
#' Note that due to the large file size of the unzipped .xpt files, this function uses a non-trivial
#' amount of RAM (~12 GB peak). To avoid crashing your computer when running, please ensure you have
#' enough RAM available.
#'
#' @param write logical argument indicating whether a .rda file should be created for each wave of processed data.
#'   Defaults to FALSE.
#'
#' @param local logical argument indicating whether the zippped raw .xpt accelerometry files are stored locally.
#' If FALSE, will download the data into a temporary file from the CDC website and process the data. If TRUE,
#' the zipped data will be sourced locally. Defaults to FALSE.
#'
#' @param localpath character string indicating where the locally zipped raw .xpt files are.
#'
#' @param deleteraw logical argument indicating whether to delete the unzipped .xpt files after reading them into R.
#'
#' @examples
#' \dontrun{
#' library("nhanesaccel")
#' process_accel()
#' process_flags()
#' }
#'
#' @importFrom haven read_xpt
#'
#' @importFrom utils write.csv unzip
#'
#' @export
process_accel <- function(write=FALSE, local=FALSE, localpath=NULL, deleteraw=TRUE){
        waves_accel <- paste0("PAXRAW_", c("C","D"))
        names_accel <- c("0304","0506")
        urls <- c("https://wwwn.cdc.gov/Nchs/Nhanes/2003-2004/PAXRAW_C.ZIP",
                  "https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/PAXRAW_D.ZIP")
        for(i in seq_along(waves_accel)){
                datapath <- c()
                if(!local){
                        datapath <- system.file("extdat/act/",package="nhanesdata")
                        temp <- tempfile()
                        download.file(urls[i], temp)
                        sim.data <- read_xpt(unzip(temp,
                                                   tolower(paste0(waves_accel[i],".xpt"))))
                        unlink(temp)
                }
                if(local & !is.null(localpath)){
                        datapath <- localpath
                        sim.data <- read_xpt(unzip(paste0(datapath, waves_accel[i],".ZIP"),
                                                   tolower(paste0(waves_accel[i],".xpt")),
                                                   exdir=datapath)
                                             )
                        file.remove(paste0(datapath, waves_accel[i],".xpt"))
                }
                if(local & is.null(localpath)){
                        datapath <- system.file("extdat/act/",package="nhanesdata")
                        sim.data <- read_xpt(unzip(paste0(datapath,waves_accel[i],".ZIP"),
                                                   tolower(paste0(waves_accel[i],".xpt")),
                                                   exdir=datapath)
                                             )
                        file.remove(paste0(datapath, waves_accel[i],".xpt"))
                }

                uid <- as.character(unique(sim.data$SEQN))

                ## read_xpt reads all variables in as numeric
                ## check for truncation on reading in
                stopifnot(all(nchar(uid)==5))

                ## create empty data frame with full 7 days of data
                ## for each subject (10,080 rows/subject)
                n    <- length(uid)
                seqn <- rep(uid,each=10080)
                paxn <- rep(c(1:10080),n)
                full.list <- data.frame(SEQN=seqn,PAXN=paxn)
                rm(list=c("n","seqn","paxn"))

                ## merge data sets to create data will NAs for missing days/times
                inx <- match(paste0(full.list$SEQN, "_", full.list$PAXN),
                             paste0(sim.data$SEQN, "_", sim.data$PAXN))
                full.na <- cbind(full.list, sim.data[inx,-c(1,5)])
                rm(list=c("full.list","inx"))


                ## create id and day of the week variables to fill in
                ## note: this assumes PAXCAL/PAXSTAT do not change from with subjects
                ## this was verified outside of this function
                u_inx  <- which(!duplicated(sim.data$SEQN))
                u_data <- sim.data[u_inx,c('SEQN','PAXCAL','PAXSTAT','PAXDAY'), drop=FALSE]

                cal  <- rep(u_data$PAXCAL,each=7)
                stat <- rep(u_data$PAXSTAT,each=7)

                rm(list=c("u_inx"))

                weekday <- rep(NA, length(uid)*7)
                inx_cur <- 1
                for(k in seq_along(uid)){
                        d <- u_data$PAXDAY[k]
                        if (d==1) {x<-c(1:7)}
                        else if (d==2) {x<-c(2:7,1)}
                        else if (d==3) {x<-c(3:7,1:2)}
                        else if (d==4) {x<-c(4:7,1:3)}
                        else if (d==5) {x<-c(5:7,1:4)}
                        else if (d==6) {x<-c(6:7,1:5)}
                        else if (d==7) {x<-c(7,1:6)}
                        weekday[inx_cur:(inx_cur+6)] <- x
                        inx_cur <- inx_cur + 7
                }
                rm(list=c("k","x","inx_cur","sim.data","u_data"))

                id2       <- rep(uid,each=7)
                idweekday <- data.frame(SEQN=id2,PAXCAL=cal,PAXSTAT=stat,WEEKDAY=weekday,wave=i, stringsAsFactors = FALSE)
                rm(list=c("weekday","id2","uid"))

                col.name <- paste0("MIN",1:1440)

                pax      <- full.na$PAXINTEN
                pax.wide <- data.frame(matrix(pax,ncol=1440,byrow=T))
                colnames(pax.wide)<-col.name

                out.name <- paste0("PAXINTEN", "_", LETTERS[i+2])
                assign(out.name, data.frame(idweekday,pax.wide,stringsAsFactors = FALSE), envir=parent.frame())


                if(write){
                        eval(parse(text=paste0('save(', out.name,',file="', datapath, out.name, '.rda")')))
                        message(paste0("Wave ", i, " Saved as: ", datapath, out.name))
                }

                rm(list=c("pax","pax.wide","col.name","out.name","idweekday"))

                message(paste("Wave", i, "Processed"))
        }
}





#' Process wear/non-wear flags for NHANES 2003-2004 and 2005-2006 accelerometry data
#'
#' @description
#' This function creates wear/non-wear flag matrices for processed NHANES 2003-2004 and 2005-2006 accelerometry data
#'
#'
#'
#' @param write logical argument indicating whether a .rda file of wear/non-wear flags
#' should be created for each wave of processed data. Defaults to FALSE.
#'
#' @param local logical argument indicating whether the processed .rda accelerometry files are stored locally.
#' If FALSE, will load them from the package data.
#'
#' @param localpath character string indicating where the locally processed .rda files are if you don't want to use those processed and included
#' in the package.
#'
#' @param window size of the moving window used to assess non-wear in minutes. Defaults to 90 minutes.
#' See \code{\link{accel.weartime}} for more details.
#'
#' @param tol maximum number of minutes with counts greater than 0 within the window allowed before a particular
#' minute is considered "wear". That is, if for a given minute, the window around that minute has tol + 1
#' activity counts greater than 0, this miute is considered "wear".
#' See \code{\link{accel.weartime}} for more details.
#'
#' @param tol.upper maximium activity count for any minute within the window
#' Defaults to 99. That is, for a given minute, if the window contains any minutes with
#' activity counts greater than tol.upper, this minute is considered "wear".
#' See \code{\link{accel.weartime}} for more details.
#'
#' @param ... aditional arguments to be passed to \code{\link{accel.weartime}}.
#'
#'
#' @examples
#' \dontrun{
#' process_accel(write=FALSE)
#' }
#'
#' @references
#'
#' @importFrom accelerometry accel.weartime
#'
#' @importFrom utils write.csv
#'
#' @export
process_flags <- function(write=FALSE, local=FALSE,localpath=NULL,
                          window=90, tol=2, tol.upper=99, ...){
        waves_accel <- paste0("PAXINTEN_", c("C","D"))
        out.name    <- paste0("Flags_", LETTERS[3:4])
        for(i in seq_along(waves_accel)){
                loaded <- waves_accel[i] %in% ls(globalenv())
                if(!loaded){
                        if(!local){
                                data(list=waves_accel[i], envir = environment(), package="nhanesdata")
                        }
                        if(local){
                                load(paste0(localpath, waves_accel[i]))
                        }
                        eval(parse(text=paste0("full_data = ",waves_accel[i])))
                }
                if(loaded){
                        message(paste( waves_accel[i],"found in the Global Environment. Using this object to create WNW flags."))
                        eval(parse(text=paste0("full_data = globalenv()[[\"", waves_accel[i], "\"]]")))
                }

                activity_data <- as.matrix(full_data[,paste0("MIN",1:1440)])
                activity_data[is.na(activity_data)] = 0    # replace NAs with zeros

                WMX = matrix(NA,nrow = nrow(activity_data), ncol = ncol(activity_data))

                t1 = Sys.time()
                pb <- txtProgressBar(min = 1, max = nrow(activity_data), style = 3)
                for (j in 1 : nrow(activity_data)){

                        activity_data_j = activity_data[j,]
                        wearMark = accel.weartime(activity_data_j,window = window,
                                                  tol = tol, tol.upper = tol.upper)
                        WMX[j,] = wearMark
                        setTxtProgressBar(pb, j)

                }

                t2 = Sys.time()
                print(paste('total time:', as.character(round(t2 - t1,2))))

                out = data.frame(full_data[,-which(colnames(full_data)%in%paste0("MIN",1:1440)),drop=FALSE], WMX,
                                 stringsAsFactors = FALSE)
                names(out) = names(full_data)

                out[is.na(full_data)] = NA ## put NAs back where they belong

                assign(out.name[i], out, envir=parent.frame())
                if(write){
                        eval(parse(text=paste0('save(', out.name[i],',file="', out.name[i],
                                               '.rda", envir=parent.frame())')))
                }



        }
}








#' Process mortality data for NHANES 2003-2004 and 2005-2006 waves
#'
#' @description
#' This function creates a clean mortality dataset which can be combined with data from the
#' NHANES 2003-2004/2005-2006 waves.
#'
#'
#'
#' @param write logical argument indicating whether a .rda file of wear/non-wear flags
#' should be created for each wave of processed data. Defaults to FALSE.
#'
#'
#'
#' @examples
#' \dontrun{
#'
#'
#' }
#'
#' @references
#'
#' @importFrom accelerometry accel.weartime
#'
#' @importFrom utils write.csv
#'
#' @export
process_mort <- function(write=FALSE){
        waves_mort <- c("NHANES_2003_2004_MORT_2011_PUBLIC.dat",
                        "NHANES_2005_2006_MORT_2011_PUBLIC.dat")

        for(i in seq_along(waves_mort)){
                path <- system.file("extdat/mort", waves_mort[i],
                                    package = "nhanesdata")
                raw.data = readLines(path)

                N = length(raw.data)

                seqn = NULL
                eligstat = NULL
                mortstat = NULL
                causeavl = NULL
                ucod_leading = NULL
                diabetes = NULL
                hyperten = NULL

                permth_int = NULL
                permth_exm = NULL
                mortsrce_ndi = NULL
                mortsrce_cms = NULL
                mortsrce_ssa = NULL
                mortsrce_dc = NULL
                mortsrce_dcl= NULL


                for (j in 1:N){

                        seqn = c(seqn,substr(raw.data[j],1,5))
                        eligstat = c(eligstat,as.numeric(substr(raw.data[j],15,15)))
                        mortstat = c(mortstat,as.numeric(substr(raw.data[j],16,16)))
                        causeavl = c(causeavl,as.numeric(substr(raw.data[j],17,17)))
                        ucod_leading = c(ucod_leading,substr(raw.data[j],18,20))
                        diabetes = c(diabetes,as.numeric(substr(raw.data[j],21,21)))
                        hyperten = c(hyperten,as.numeric(substr(raw.data[j],22,22)))

                        permth_int = c(permth_int,as.numeric(substr(raw.data[j],44,46)))
                        permth_exm = c(permth_exm,as.numeric(substr(raw.data[j],47,49)))

                        mortsrce_ndi = c(mortsrce_ndi,substr(raw.data[j],50,50))
                        mortsrce_cms = c(mortsrce_cms,substr(raw.data[j],51,51))
                        mortsrce_ssa = c(mortsrce_ssa,substr(raw.data[j],52,52))
                        mortsrce_dc = c(mortsrce_dc,substr(raw.data[j],53,53))
                        mortsrce_dcl = c(mortsrce_dcl,substr(raw.data[j],54,54))

                }


                out.name <- paste0("Mortality_",LETTERS[i+2])
                out <- data.frame(seqn, eligstat,
                                  mortstat, causeavl,
                                  ucod_leading, diabetes,
                                  hyperten,
                                  permth_exm, permth_int,
                                  mortsrce_ndi, mortsrce_cms,
                                  mortsrce_ssa, mortsrce_dc,
                                  mortsrce_dcl,
                                  stringsAsFactors=FALSE)

                colnames(out) = c('SEQN', 'eligstat',
                                  'mortstat', 'causeavl',
                                  'ucod_leading', 'diabetes',
                                  'hyperten',
                                  'permth_exm', 'permth_int',
                                  'mortsrce_ndi', 'mortsrce_cms',
                                  'mortsrce_ssa', 'mortsrce_dc',
                                  'mortsrce_dcl')

                assign(out.name, out, envir=parent.frame())

                if(write) {
                        eval(parse(text=paste0('save(', out.name,',file="', out.name,
                                               '.rda", envir=parent.frame())')))
                }

        }



}








#' Process covariate data for NHANES 2003-2004 and 2005-2006 waves
#'
#' @description
#' This function merges covariate data
#' from several NHANES data files.
#' The resulting data matrix includes survey weights as well as a number of demographic and
#' comorbidity responses.
#'
#'
#' @param cohorts numeric vector with entires corresponding to the first year of a given
#' NHANES wave of interest. Will only accept cohorts 2003-2004 and beyond.
#'
#' @param varnames character vector indicating which column names are to be searched for.
#' Will check all files in dataPath
#'
#' @param dataPath file path where covariate data are saved. Covariate data must be in .XPT format,
#' and should be in their own folder. For example, PAXRAW_C.XPT should not be located in the folder with
#' your covariate files. This will not cause an error, but the code will take much longer to run.
#'
#' @param write logical argument indicating whether a .rda file of covariate data. Defaults to FALSE.
#'
#'
#'
#' @examples
#' \dontrun{
#'
#'
#' }
#'
#' @references
#'
#'
#' @export
process_covar <- function(cohorts=c(2003,2005),
                          varnames = c('WTMEC2YR','WTINT2YR','SDMVPSU','SDMVSTRA',
                                       'RIDAGEMN','RIDAGEEX','RIDRETH1','RIAGENDR',
                                       'BMXWT','BMXHT','BMXBMI','DMDEDUC2',
                                       'ALQ101', 'ALQ110','ALQ120Q','ALQ120U','ALQ130',
                                       'SMQ020','SMD030','SMQ040',
                                       'SMD070',
                                       'MCQ220','MCQ160F','MCQ160B','MCQ160C','PFQ061B','PFQ061C', 'DIQ010'),
                          dataPath=NULL,
                          write=FALSE){
        years    <- seq(2003, 2023, by=2)
        if(!all(cohorts %in% years)) stop("One or more cohorts invalid")
        stopifnot(length(cohorts) >= 1)
        stopifnot(is.vector(cohorts))

        cohorts <- sort(cohorts)
        varnames <- c('SEQN',varnames)
        ## find all files of .xpt structure that correspond to the year specified in the data
        if(is.null(dataPath)){
                dataPath <- system.file("extdat/covar/",package="nhanesdata")
        }
        files_full   <- list.files(dataPath)

        for(i in seq_along(cohorts)){
                cohort <- cohorts[i]

                pathExt <- paste('_', LETTERS[3:26], '.XPT', sep='')[which(years == cohort)]
                files   <- files_full[substr(files_full, (nchar(files_full) - 5), nchar(files_full)) == pathExt]

                if(length(files) == 0) next

                ## find which files contain the variables requested by the user
                covarMats <- lapply(files, function(x){
                        mat <- read_xpt(paste0(dataPath,x))
                        mat <- mat[,colnames(mat)%in%varnames,drop=FALSE]
                        if(!is.null(dim(mat))) mat
                        else NULL
                })
                ##
                covarMats <- covarMats[!vapply(covarMats, is.null,logical(1))]
                matchedNames <- lapply(covarMats, colnames)
                numMatched   <- length(unlist(matchedNames))

                if(numMatched == 0) stop('Error: No Variable Names Recognized for this Year/Variable Combination')
                if(numMatched > 0)  message(
                        paste("For", cohort, "cohort,",
                              (numMatched - length(matchedNames)),
                              'Covariates Found of', (length(varnames)-1),'specified.',
                              'Missing covariates:',
                              paste(setdiff(varnames, unlist(matchedNames)),collapse=", ") ))

                ## Merge the covariate data
                ids       <- sort(unique(unlist(lapply(covarMats, function(x) x[,'SEQN']))))
                totalCols <- sum(vapply(covarMats, ncol, numeric(1))) - length(covarMats) + 1
                CovarMat           <- matrix(NA,ncol=totalCols,nrow=length(ids))
                colnames(CovarMat) <- c('SEQN', unlist(sapply(covarMats,function(x) colnames(x)[-1])))
                CovarMat[,'SEQN']  <- ids
                CovarMat <- data.frame(CovarMat)
                invisible(lapply(covarMats, function(x) CovarMat[,colnames(x)[-1]] <<- as.matrix(x[match(CovarMat$SEQN,x$SEQN),-1]) ) )



                CovarMat$SEQN <- as.character(CovarMat$SEQN)

                out.name <- paste0("Covariate_",substr(pathExt,2,2))

                assign(out.name, CovarMat, envir=parent.frame())

                if(write) {
                        eval(parse(text=paste0('save(', out.name,',file="', out.name,
                                               '.rda", envir=parent.frame())')))
                }

        }


}









#' Remove days with too few/too much weartime and NHANES data quality flags.
#'
#' @description
#' This function subsets accelerometry data by wear/non-wear time criteria, as well as
#' NHANES data quality flags.
#'
#'
#'
#' @examples
#' \dontrun{
#' process_accel(write=FALSE)
#' }
#'
#' @references
#'
#'
#' @export
exclude_accel <- function(act, flags, threshold = 600, rm_PAXSTAT = TRUE, rm_PAXCAL = TRUE,
                          return_act = FALSE){

        stopifnot(all(is.finite(act$PAXSTAT),is.finite(act$PAXCAL)))
        flag_nonwear <- rowSums(flags[, paste('MIN', 1:1440, sep='')], na.rm = TRUE) < threshold

        ## remove nonwear days and days flagged by NHANES
        cond <- c("flag_nonwear", "act$PAXSTAT!=1","act$PAXCAL!=1")[c(TRUE, rm_PAXSTAT, rm_PAXCAL)]
        cond <- paste(cond, collapse="|")

        if(return_act) return(eval(parse(text=paste("return(act[!(",cond, "),])"))))

        eval(parse(text=paste("return(which(!(",cond,")))")))
}




#' Reweight NHANES accelerometry data
#'
#' @description
#' This function re-weights accelerometry data for NHANES 2003-2004,2005-2006 waves
#'
#'
#'
#' @examples
#' \dontrun{
#' process_accel(write=FALSE)
#' }
#'
#' @references
#'
#'
#' @export
reweight_accel <- function(data){
        data$valid <- 1
        data.wave1 <- subset(data, wave == 1)
        data.wave2 <- subset(data, wave == 2)

        data.wave1.adjusted <- nhanes.accel.reweight(acceldata=data.wave1,
                                                     wave=1, seqn.column=1,
                                                     include.column=which(colnames(data.wave1) == 'valid'))
        # Rename previous 2-year weight variable so new one does not overwrite it

        data.wave2.adjusted <- nhanes.accel.reweight(acceldata=data.wave2,
                                                     wave=2, seqn.column=1,
                                                     include.column=which(colnames(data.wave1) == 'valid'))

        all.data <- rbind(data.wave1.adjusted, data.wave2.adjusted)
        if(nrow(data.wave1) == 0 | nrow(data.wave2) == 0) {
                all.data$wtmec4yr_adj <- all.data$wtmec2yr_adj
        } else {
                all.data$wtmec4yr_adj <- all.data$wtmec2yr_adj/2
        }

        all.data$NormWts <- (all.data$wtmec4yr_adj/sum(all.data$wtmec4yr_adj))*nrow(all.data)

        all.data
}

