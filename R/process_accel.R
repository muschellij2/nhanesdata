#' Process NHANES 2003-2004 and 2005-2006 accelerometry data
#'
#' @description
#' This function processes the raw NHANES 2003-2004 and 2005-2006 accelerometry data
#' as provided by the NHANES study on the CDC website.
#' Note that due to the large file size of the unzipped .xpt files, this function uses a non-trivial
#' amount of RAM (~12 GB peak). To avoid crashing your computer when running, please ensure you have
#' enough RAM available.
#'
#' @param write logical argument indicating whether a csv file should be created for each wave of processed data.
#'   Defaults to FALSE.
#'
#' @examples
#' \dontrun{
#'
#' process_accel(write=FALSE)
#'
#' }
#'
#' @references
#'
#' @ importFrom haven read_xpt
#'
#' @ importFrom utils write.csv
#'
#' @export
process_accel <- function(write=FALSE){
        waves_accel <- paste0("PAXRAW_", c("C","D"))
        names_accel <- c("0304","0506")
        for(i in seq_along(waves_accel)){
                sim.data <- read_xpt(unzip(paste0(waves_accel[i],".zip"),
                                           tolower(paste0(waves_accel[i],".xpt"))))
                uid      <- unique(sim.data$SEQN)

                print('loaded data')
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
                # full.na  <- left_join(full.list, sim.data, by=c("SEQN","PAXN"))
                # rm(list=c("full.list"))

                inx <- match(paste0(full.list$SEQN, "_", full.list$PAXN),
                             paste0(sim.data$SEQN, "_", sim.data$PAXN))
                full.na <- cbind(full.list, sim.data[inx,-c(1,5)])
                rm(list=c("full.list","inx"))

                print('made it past merging')

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

                print('made it to dates')

                id2       <- rep(uid,each=7)
                idweekday <- data.frame(SEQN=id2,PAXCAL=cal,PAXSTAT=stat,WEEKDAY=weekday)
                rm(list=c("weekday","id2","uid"))

                col.name <- paste0("MIN",1:1440)

                pax      <- full.na$PAXINTEN
                pax.wide <- as.data.frame(matrix(pax,ncol=1440,byrow=T))
                colnames(pax.wide)<-col.name

                out.name <- paste0("PAXINTEN", "_", LETTERS[i+2])
                assign(out.name, cbind(idweekday,pax.wide), envir=parent.frame())


                if(write){
                        eval(parse(text=paste0('save(', out.name,',file="', out.name, '.Rda")')))
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
#' @param write logical argument indicating whether a csv file of wear/non-wear flags
#' should be created for each wave of processed data. Defaults to FALSE.
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
#'
#' process_accel(write=FALSE)
#'
#' }
#'
#' @references
#'
#' @ importFrom accelerometry accel.weartime
#'
#' @ importFrom utils write.csv
#'
#' @export
process_flags <- function(write=FALSE, window, tol, tol.upper, ...){

}








#' Process mortality data for NHANES 2003-2004 and 2005-2006 waves
#'
#' @description
#' This function creates a clean mortality dataset which can be combined with data from the
#' NHANES 2003-2004/2005-2006 waves.
#'
#'
#'
#' @param write logical argument indicating whether a csv file of wear/non-wear flags
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
#' @ importFrom accelerometry accel.weartime
#'
#' @ importFrom utils write.csv
#'
#' @export
process_mort <- function(write=FALSE){

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
#'
#' @param write logical argument indicating whether a csv file of wear/non-wear flags
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
#' @ importFrom utils write.csv
#'
#' @export
process_covar <- function(write=FALSE){

}











