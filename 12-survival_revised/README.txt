In this folder there are results from survival analysis in particular:
  
- script_survival.R performs the Kaplan-Meier analysis plotting the survival curves and the log-rank test p-value.

- cox_fdr.R performs the cox regression model including some explanatory variables (tumor stages, age, gender) except from high and low group. This analysis was ran for our candidates genes and fdr was calculated for each p-values. The output was saved into a csv file (cox_table_LUAD_final_table.csv, cox_table_LUSC_final_table.csv).

- functions_survival.R includes all the functions used in the previous scripts.
- the tables with metadata need to be downloaded by GDCPortal or with TCGAbiolinks functions
LUAD and LUSC folders contain clinical data and Kaplan-Meier survival plots.
