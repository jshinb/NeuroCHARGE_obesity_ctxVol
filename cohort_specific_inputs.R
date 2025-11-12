# COHORT INFORMATION ------------------------------------------------------------------- #
# Analyst/PI information for contact and study information: put double or single quotation
# marks
#
analysis_date =     # YYYY-MM-DD
analyst = 
analyst_email = 
PI =  #use c('PI1 name','PI2 name',...) if more than 1 
PI_email = #use c('PI1 email','PI2 email',...) if more than 1 
cohort_full_name = #'Saguenay Youth Study' 
cohort_short_name = #'SYS_parents' 
ethnicities = # or c('White', 'African American') 
ancestries = # (AFR, AMR (American-Hispanic?), EAS, EUR, SAS)
FreeSurfer_version = #e.g., '5.3.0'
ICV_based_on =  #e.g., FreeSurfer/ customized pipeline (ref:PMID)
ICV_pipeline_reference = #'PMIDXXX'

# DATA INFORMATION --------------------------------------------------------------------- #
# working directory(wd), file names for the brain and non-brain datasets, and missing value code
wd=#directory where you downloaded the R scripts
input_neuro=#path to the data file with brain image variables "../data/neuro_SPS.txt"
input_non.neuro=# path to the data file with non-brain variables'../data/non.neuro_hypothetical_SPS.txt'
missing_value_code='NA' #how missing values are coded
#----------------------------------------------------------------------------------------#
