## IMPORTANT: Please set code_root variable properly. 
## code_root should be set to the directory where the repository README file is located. 
## For more information, please read the repository README file
code_root="C:\\Users\\Jonathan\\OneDrive - Harvard University\\School\\Harvard\\BST249\\Project\\BST249 Project Code\\SAPHIRE-master\\"

setwd(paste0(code_root, "scripts_main"))
library(vioplot)
library("corrplot")
library(readr)
library(cairoDevice)

##
source(paste0(code_root, "R/fun_SEIRpred.R"))
source(paste0(code_root, "R/fun_SEIRsimu.R"))
source(paste0(code_root, "R/fun_SEIRfitting_BST249.R"))
source(paste0(code_root, "R/fun_BTSetup_BST249.R"))
source(paste0(code_root, "R/init_cond.R"))
source(paste0(code_root, "R/fun_R0estimate.R"))
source(paste0(code_root, "R/fun_SEIRplot.R"))
source(paste0(code_root, "R/fun_Findzero.R"))
##

init_sets_list=get_init_sets_list(r0 = 0.23)
SEIRfitting(init_sets_list, randomize_startValue = T,
            run_id = "main_analysis", output_ret = T, skip_MCMC=F)


# Simulated data
sim1 <- read.table("C:\\Users\\Jonathan\\Dropbox\\BST249 Group Project\\Simulated_data\\data_Replicates_SAPHIRE\\main_est_SAPHIRE.txt",
                   header = TRUE)
init_sets_list$daily_new_case <- sim1$Onset_expect
init_sets_list$daily_new_case_all <- c(sim1$Onset_expect, rep(0,8))
