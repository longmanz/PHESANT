## Full PHEASANT pipeline for UK-Biobank Nealelab GWAS.

This page describes the steps we applied in going from the raw data provided by UK-biobank, to that used to run GWAS.

* Create a new VM on google cloud, install R and git. 
* Clone the PHESANT repository.
* Run `reengineering_phenofile_pharma_exomes.r`
  * Remove cancer variables - these have a up to 31 visits and should be considered separately.
  * Restrict to the first visit (we made the assumption that the first visit contained the least missing data).
  * Prepend each column name with `x` and replace all instances of `-` and `.` with `_`, for compatibility with PHESANT.
* Restrict to the subset of samples that are output by the genetic data QC pipeline.
* The resulant phenotype file, restricted to the individuals for whom we have clean genetic data is then passed to PHESANT. A summary of how PHESANT parses the raw phenotype data is shown [here](https://github.com/astheeggeggs/PHESANT/blob/master/pipeline.pdf).
* PHESANT inputs
  * `variablelistfile` [outcome_info_final_round2.tsv](https://github.com/astheeggeggs/PHESANT/tree/master/variable-info/outcome_info_final_round2.tsv).
    * As described in the [README for PHESANT](https://github.com/astheeggeggs/PHESANT/), see the Variable information file heading.
      * The most important column is `EXCLUDED`. If there is any text in a cell in this column it is excluded.
      * `CAT_MULT_INDICATOR_FIELDS` is also important, it describes who to include as controls for this phenotype according to the rules explained in the [README for PHESANT](https://github.com/astheeggeggs/PHESANT/).
  * `datacodingfile` [data-coding-ordinal-info.txt](https://github.com/astheeggeggs/PHESANT/tree/master/variable-info/data-coding-ordinal-info.txt). See the [README for PHESANT](https://github.com/astheeggeggs/PHESANT/) for details of its format.
  * Our default settings for filters.
    * `--catmultcutoff` 50. The cutoff for exclusion when creating dichotomous variables for CAT-MULTIPLE.
    * `--catordnacutoff` 500. The cutoff for exclusion for number of non-NAs in ordered categorical variables.
    * `--catunordnacutoff` 5000. The cutoff for exclusion for number of non-NAs in unordered categorical variables.
    * `--contnacutoff` 5000. The cutoff for exclusion for number of non-NAs in continuous variables.
    * `--binnacutoff` 5000. The cutoff for exclusion for number of non-NAs in binary-variables.
    * `--bintruecutoff` 100. The cutoff for exclusion for numbers of members of a category in binary-variables.
    * `--mincategorysize` 10. The minimum number of samples in a category for categorical single, integer, and continous variables.
    * `--maxunorderedcategories` 1000. The maximum number of categories in an unordered categorical variable.
    * `--propforcontinuous` 0.2. The cutoff for proportion of samples with the same value for the variable to not be considered continuous.
* We run PHESANT in chunks on the cloud. 
* We then run the `summarise_phenotypes.r` functions to obtain a summary file of the phenotypes that made it into the final file.
