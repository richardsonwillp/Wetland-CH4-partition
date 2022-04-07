# this script determines the minimum reference scalar flux value at which scalar similarity is 
# adequate for partitioning
  ## fractional ebullition increases erroneously under low reference scalar fluxes; fit a
  ## broken line regression between fractional ebullition and each reference scalar flux
  ## using the R package 'segmented' (Muggeo 2003, 2008)

# libraries
library(segmented)

# get working directory and other environment variables
base_path <- getwd()
site_yr <- Sys.getenv('site_yr')
run_ID <- Sys.getenv('run_ID')

# load frame with all stats
master_file <- paste(base_path, '/ref_data/Output/', site_yr, '/', run_ID, '/', site_yr, '_Interm_RefDf_allHH_PostProcStats.csv', sep='')
df_master <- read.csv(master_file, stringsAsFactors=FALSE)
## Timestamp column is being labeled as 'X' by in post-processing stats, change to Timestamp
names(df_master)[names(df_master) == 'X'] <- 'Timestamp'

# load qc filtered frames; use their index values to keep desired rows from all stats frame
qfilt_file <- paste(base_path, '/ref_data/Output/', site_yr, '/', site_yr, '_Interm_RefDf_qPart_filtHH_NLfilt.csv', sep='')
qPart_filt <- read.csv(qfilt_file, stringsAsFactors=FALSE)
df_qPart <- subset(df_master, Timestamp %in% qPart_filt$Timestamp)

cfilt_file <- paste(base_path, '/ref_data/Output/', site_yr, '/', site_yr, '_Interm_RefDf_cPart_filtHH_NLfilt.csv', sep='')
cPart_filt <- read.csv(cfilt_file, stringsAsFactors=FALSE)
df_cPart <- subset(df_master, Timestamp %in% cPart_filt$Timestamp)

# filter for low CH4 fluxes (0.01 umol m-2 s-1) prior to doing analysis; scalar similarity not necessarily maintained
df_qPart <- df_qPart[df_qPart$CH4_tot >= 0.01,]
df_cPart <- df_cPart[df_cPart$CH4_tot >= 0.01,]
# filter for negative fractional ebullition
df_qPart <- df_qPart[df_qPart$frac_eb_q >= 0.0,]
df_cPart <- df_cPart[df_cPart$frac_eb_c >= 0.0,]

##########################################################################################################################
#############                         Fit broken line regressions, save parameters                           #############
##########################################################################################################################
# FracEbq - LE
FracEbq_LE_glm <- glm(frac_eb_q ~ LE, data=df_qPart)
FracEbq_LE_seg <- segmented(FracEbq_LE_glm, seg.Z = ~LE, psi=25)
summary(FracEbq_LE_seg)
# FracEbc - CO2 flux magnitude
FracEbc_Fco2_glm <- glm(frac_eb_c ~ co2_flux_mag, data=df_cPart)
FracEbc_Fco2_seg <- segmented(FracEbc_Fco2_glm, seg.Z = ~co2_flux_mag, psi=10)
summary(FracEbc_Fco2_seg)

output = data.frame(FracEbq_LE_brkpt = FracEbq_LE_seg$psi[2], FracEbq_LE_brkpt_SE = FracEbq_LE_seg$psi[3], 
                    FracEbq_LE_LHS_slope = FracEbq_LE_seg$coefficients[2], FracEbq_LE_LHS_int = FracEbq_LE_seg$coefficients[1],
                    FracEbq_LE_RHS_slope = slope(FracEbq_LE_seg)[[1]][2],
                    FracEbc_Fco2_brkpt = FracEbc_Fco2_seg$psi[2], FracEbc_Fco2_brkpt_SE = FracEbc_Fco2_seg$psi[3], 
                    FracEbc_Fco2_LHS_slope = FracEbc_Fco2_seg$coefficients[2], FracEbc_Fco2_LHS_int = FracEbc_Fco2_seg$coefficients[1],
                    FracEbc_Fco2_RHS_slope = slope(FracEbc_Fco2_seg)[[1]][2],
                    row.names = NULL)

out_path <- paste(base_path, '/ref_data/Output/', site_yr, '/', run_ID, '/', site_yr, '_FracEb_RefScalFlux_brkpts.csv', sep='')
write.csv(output, out_path, row.names=FALSE)

