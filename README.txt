The Rmarkdown file main.Rmd is used to generate the main.PDF file containing all the results that are used in the report.
The raw data found in "Data" is used to generate the processed files in "Data/Processed". After MOFA is run, the resulting models are saved in "Models".
If the following files are already detected:

- "Data/Processed/GSEA_results.RData"				Runtime: ~15min
- "Data/Processed/processed_data.RData"				Runtime: ~10min
- "Data/Processed/results_MOFA_bootstrap.RData"		Runtime: ~4H
- "Data/Models/MOFA_data.RData"						Runtime: ~1H

the processing and MOFA steps are skipped and the analysis is performed on the imported data. Otherwise, the required files are generated on the spot.