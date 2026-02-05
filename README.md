This github repository contains the necessary indices and outputs from the ML classifier for bladder cancer recurrence and corresponding Rmarkdowns. 

# Table of Contents
- 4 separate .md files generated from their corresponding Rmarkdown (.Rmd) files. (01-04 should be run sequentially)
- plots/ contains all the graphs exported from these notebooks
- model/
    - train-test split indices
    - features used in the ML classifier \
    *Caveat*: The raw models were not exported as an rds but it should be reproducible regardless
- data_cleaned/
    - clinical and expression data after preprocessing in 01-preprocessing.Rmd
- figures/
    - figures put together for the report
- UROMOL_TaLG.teachingcohort.rds
    - UROMOL training dataset
- knowles_matched_TaLG_final.rds
    - Knowles validation dataset