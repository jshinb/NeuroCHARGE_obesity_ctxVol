# COHORT INFORMATION ------------------------------------------------------------------- #
# Analyst/PI information for contact and study information: put double or single quotation
# marks
#
analysis_date = #*[YYYY-MM-DD]
analyst = 
analyst_email = 
PI = #*[use c('PI1 name','PI2 name',...) if more than 1]
PI_email = #*[use c('PI1 email','PI2 email',...) if more than 1]
cohort_full_name = #*[e.g., 'Saguenay Youth Study']
cohort_short_name = #*[e.g., 'SYS_parents']
ethnicities = #*["White" or use c('White', 'African American') in the presence of multiple ethnicities]
ancestries = #*["EUR" or c('AFR', 'AMR', 'EAS', 'EUR', 'SAS',...) in the presence of multiple ancestries]
FreeSurfer_version = #*[e.g., '5.3.0']
ICV_based_on =  #*[e.g., FreeSurfer/ customized pipeline (ref:PMID)]
ICV_pipeline_reference = #*['PMIDXXX']

# DATA INFORMATION --------------------------------------------------------------------- #
# working directory(wd), file names for the brain and non-brain datasets, and missing value code
wd = #*[directory where the R scripts and data files are stored]
input_neuro = #*[name of data file with brain image variables: e.g., "neuro_SPS.txt"]
input_non.neuro = #*[name of the data file with non-brain variables: e.g., 'non.neuro_SPS.txt']
missing_value_code ='NA' #*[how missing values are coded]
#----------------------------------------------------------------------------------------#

