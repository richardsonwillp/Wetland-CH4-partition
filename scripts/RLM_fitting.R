# Written to independently do some of the processing steps from the Iwata program, namely:

# 1) Load wavelet coefficient files, filter by frequency 
# 2) Fit RLM between Wm (dependent var) and Wx (independent var)
# 3) Save slopes of lines of best fit, calculate RMSDs for Wm - Wx and save

# in order to fit dynamic ebullition threshold without needing to have done an entire partitioning run previously

#### 2022-02-01 UPDATE: RLM fitting between Wc and Wq added to assess c-q similarity for quality assurance

## CHOOSE METHOD FOR CALCULATING fn 
fn_type <- 'include_zc' # alternative would be 'original'

args = commandArgs(trailingOnly = TRUE)
# get working directory
base_path <- getwd()

site_yr <- Sys.getenv('site_yr')

## fn bounds
### think about loading numeric arguments using as.double(args[x])
fn_LB_str <- Sys.getenv('LB'); fn_UB_str <- Sys.getenv('UB')
fn_LB <- as.double(fn_LB_str) #[-] lower frequency bound
fn_UB <- as.double(fn_UB_str) #[-] upper frequency bound

# user parameters to change between different site-years of data (configured for WAY 3 SUMMER 2017 currently)
zm <- as.double(Sys.getenv('Zm')) #[m] measurement height
dt <- as.double(Sys.getenv('DT'))#[s] temporal resolution of EC measurements
N_wvlt <- 2**(floor(log2((30*60)/dt))) # number of wavelet coefficients in each half-hour period

# filepaths
wvlt_file_loc <- paste(base_path, '/wvlet/', site_yr, sep='')
zc_loc <- Sys.getenv('zc_file')
output_loc <- paste(base_path, '/ref_data/Output/', site_yr, '/', site_yr,'_RLMfitRMSDs_fnLB_', fn_LB_str,'_fnUB_', fn_UB_str, '.csv', sep='')

# libraries
library(MASS)
library(Metrics)

# load canopy height data
zc_df <- read.delim(zc_loc, header=FALSE, sep='\t', col.names=c('date', 'canopy_height'), stringsAsFactors=FALSE)
# get list of wvlt coef filepaths
wvlt_filepaths = list.files(wvlt_file_loc, all.files=FALSE, full.names=TRUE, recursive=TRUE)

# main routine of script
main <- function(filepaths){
  
  n_files <- length(filepaths)
  
  # Make variables to hold results
  N_coeff_out <- rep(NA, n_files); tstamp_out <- rep(NA, n_files)
  q_mods <- list(); mq_slope_out <- rep(NA, n_files); q_RMSD_out <- rep(NA, n_files); q_MAD_out <- rep(NA, n_files); q_R2_out <- rep(NA, n_files)
  c_mods <- list(); mc_slope_out <- rep(NA, n_files); c_RMSD_out <- rep(NA, n_files); c_MAD_out <- rep(NA, n_files); c_R2_out <- rep(NA, n_files)
  T_mods <- list(); mT_slope_out <- rep(NA, n_files); T_RMSD_out <- rep(NA, n_files); T_MAD_out <- rep(NA, n_files); T_R2_out <- rep(NA, n_files)
  cq_mods <- list(); cq_slope_out <- rep(NA, n_files); cq_RMSD_out <- rep(NA, n_files); cq_MAD_out <- rep(NA, n_files); cq_R2_out <- rep(NA, n_files)
  
  # iterate over all half hour periods
  for(i in 1:n_files){
    
    # get info for this period
    file <- filepaths[i]
    ## date and time
    day <- substr(strsplit(strsplit(file, '/')[[1]][8], '-')[[1]][2], start=1, stop=8)
    date_str <- gsub("(\\d{4})(\\d{2})(\\d{2})$","\\1\\2\\3", day)
    time <- strsplit(strsplit(strsplit(file, '/')[[1]][8], '_')[[1]][2], '[.]')[[1]][1]
    time_str <- gsub("(\\d{2})(\\d{2})$","\\1:\\2", time)
    datetime_str <- paste(date_str, time_str, sep=' ')
    
    ## zc and zd
    zc <- zc_df$canopy_height[zc_df$date == date_str]
    zd <- 0.67*zc
    
    # column names for wavelet coefficient files
    cols <- c('Index_wvlt', 'u', 'w', 'T', 'q', 'c', 'm')
    
    # read wavelet file
    wvlt_df <- read.table(file, col.names=cols, skip=1)

    # calculate mean wind speed
    u_mean <- mean(wvlt_df$u)

    # calculate normalized frequency, add it to frame
    wvlt_df$j = (N_wvlt/wvlt_df$Index_wvlt)*dt
    if (fn_type == 'original'){
      wvlt_df$fn <- (zm)/(wvlt_df$j*u_mean)
    }else{
      wvlt_df$fn <- (zm - zd)/(wvlt_df$j*u_mean)
    }
    
    # filter wvlt_df based on fn_LB, fn_UB
    wvlt_df_filt <- wvlt_df[fn_LB <= wvlt_df$fn & wvlt_df$fn <= fn_UB, ]
    
    ## get number of coeffs used to fit model
    N_coeff_mod <- nrow(wvlt_df_filt)   
    
    # make sure that file has adequate observations for model fitting
    if (nrow(wvlt_df_filt) > 3){  
      
      # fit models, calculate all stats of interest and save
      ## m vs. q
      q_mod <- rlm(m ~ 0 + q, data=wvlt_df_filt, psi=psi.bisquare, method='MM', maxit=100)
      ## m vs. c
      c_mod <- rlm(m ~ 0 + c, data=wvlt_df_filt, psi=psi.bisquare, method='MM', maxit=100)    
      ## m vs. T
      T_mod <- rlm(m ~ 0 + T, data=wvlt_df_filt, psi=psi.bisquare, method='MM', maxit=100)
      ## c vs. q
      cq_mod <- rlm(c ~ 0 + q, data=wvlt_df_filt, psi=psi.bisquare, method='MM', maxit=100)
      
      ## get slopes
      slope_mq <- q_mod$coefficients; slope_mc <- c_mod$coefficients; slope_mT <- T_mod$coefficients; slope_cq <- cq_mod$coefficients
      
      ## Calculate RMSDs And R2
      ### q_mod
      RMSD_q <- rmse(q_mod$model$m, q_mod$fitted.values)
      MAD_q <- mad(q_mod$residuals)
      R2_q <- 1 - (sum(q_mod$residuals^2))/(sum((q_mod$model$m - mean(q_mod$model$m))^2))
      ### c_mod
      RMSD_c <- rmse(c_mod$model$m, c_mod$fitted.values)
      MAD_c <- mad(c_mod$residuals)
      R2_c <- 1 - (sum(c_mod$residuals^2))/(sum((c_mod$model$m - mean(c_mod$model$m))^2))
      ### T_mod
      RMSD_T <- rmse(T_mod$model$m, T_mod$fitted.values)
      MAD_T <- mad(T_mod$residuals)
      R2_T <- 1 - (sum(T_mod$residuals^2))/(sum((T_mod$model$m - mean(T_mod$model$m))^2))
      ### cq_mod
      RMSD_cq <- rmse(cq_mod$model$c, cq_mod$fitted.values)
      MAD_cq <- mad(cq_mod$residuals)
      R2_cq <- 1 - (sum(cq_mod$residuals^2))/(sum((cq_mod$model$c - mean(cq_mod$model$c))^2))
      
      ## save models and stats 
      tstamp_out[i] <- datetime_str; N_coeff_out[i] <- N_coeff_mod
      q_mods[[i]] <- q_mod; c_mods[[i]] <- c_mod; T_mods[[i]] <- T_mod; cq_mods[[i]] <- cq_mod
      mq_slope_out[i] <- slope_mq; q_RMSD_out[i] <- RMSD_q; q_MAD_out[i] <- MAD_q; q_R2_out[i] <- R2_q
      mc_slope_out[i] <- slope_mc; c_RMSD_out[i] <- RMSD_c; c_MAD_out[i] <- MAD_c; c_R2_out[i] <- R2_c
      mT_slope_out[i] <- slope_mT; T_RMSD_out[i] <- RMSD_T; T_MAD_out[i] <- MAD_T; T_R2_out[i] <- R2_T
      cq_slope_out[i] <- slope_cq; cq_RMSD_out[i] <- RMSD_cq; cq_MAD_out[i] <- MAD_cq; cq_R2_out[i] <- R2_cq
      
      print(paste(datetime_str, 'done', sep=' '))
    }else{
      tstamp_out[i] <- datetime_str; N_coeff_out[i] <- NA
      q_mods[[i]] <- NA; c_mods[[i]] <- NA; T_mods[[i]] <- NA; cq_mods[[i]] <- NA
      mq_slope_out[i] <- NA; q_RMSD_out[i] <- NA; q_MAD_out[i] <- NA; q_R2_out[i] <- NA
      mc_slope_out[i] <- NA; c_RMSD_out[i] <- NA; c_MAD_out[i] <- NA; c_R2_out[i] <- NA
      mT_slope_out[i] <- NA; T_RMSD_out[i] <- NA; T_MAD_out[i] <- NA; T_R2_out[i] <- NA
      cq_slope_out[i] <- NA; cq_RMSD_out[i] <- NA; cq_MAD_out[i] <- NA; cq_R2_out[i] <- NA
      print(paste(datetime_str, 'not enough data to fit model', sep=' '))
    }
  }
  output <- data.frame(Timestamp=tstamp_out, N_coeff=N_coeff_out,
                       slope_mq=mq_slope_out, RMSD_mq=q_RMSD_out, MAD_mq=q_MAD_out, R2_mq=q_R2_out,
                       slope_mc=mc_slope_out, RMSD_mc=c_RMSD_out, MAD_mc=c_MAD_out, R2_mc=c_R2_out, 
                       slope_mT=mT_slope_out, RMSD_mT=T_RMSD_out, MAD_mT=T_MAD_out, R2_mT=T_R2_out,
                       slope_cq=cq_slope_out, RMSD_cq=cq_RMSD_out, MAD_cq=cq_MAD_out, R2_cq=cq_R2_out,
                       row.names = 'Timestamp')
  return(output)
}

# execute main routine
df_output <- main(wvlt_filepaths)

# save output
write.csv(df_output, output_loc, row.names=TRUE)
