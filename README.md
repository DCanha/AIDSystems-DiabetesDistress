# Automated Insulin Delivery Systems effect on Diabetes Distress

Automated insulin delivery systems (AID) provide an alternative to multiple daily injections, open-loop insulin systems, and other traditional forms of treatment requiring patient calculations. The constant need for monitoring becomes burdensome for the patient and contributes to adverse psychological outcomes. These include diabetes distress (DD), which encompasses the psychological impact of the disease stemming from feelings of fear, being overwhelmed, and strained social relationships, among others. DD is prevalent among people with diabetes (PwD). While AIDs improve glycemic control and are expected to alleviate diabetes burden, their impact on DD is unclear. Therefore, this repository tries to uncover whether AIDs can effectively reduce DD in PwD, providing the required R syntax to reproduce the meta-analysis results presented in the scientific work that arose from this research question.

## Repository Structure

In this repository you will find the ```code``` folder that contains all the syntax needed to reproduce the meta-analysis results:
 - ```ma-before-after.R```: Contains the R code for conducting meta-analysis on Before-After Studies (BAS) examining the impact of AID on DD. Three separate analyses are conducted for each population: adult, pediatric, and caregiver.
 - ```ma-rct.R```: Contains the R code for analyzing Randomized Controlled Trials (RCTs) to assess the effect of AID on DD. Three separate analyses are conducted for each population: adult, pediatric, and caregiver.

## How to use

To replicate the analysis or adapt it to similar research, follow the comments and structure within each R script. Ensure all dependencies are installed. The initial datasets for each analysis are created as data frames with feature values taken directly from the extracted papers. The meta-analysis is performed on standardized mean differences (SMDs) and corresponding standard errors (SEs). 

## Dependencies 

 - R (4.2.2)
 - R packages: *dplyr*, *esc*, *meta*, *metafor*, *estmeansd*
