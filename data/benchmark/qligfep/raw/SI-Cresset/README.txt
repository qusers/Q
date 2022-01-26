The subfolders in wang_2015_datasets contain csv files with ddG values for all perturbations of the different approaches by
 - "cresset_valid" as supplied by Kuhn et al. (this study)
 - "schrodinger" as supplied in the supporting information of Wang et al., J. Am. Chem. Soc. 2015, 137, 7, 2695-2703; https://doi.org/10.1021/ja512751q
 - "amber" as supplied in the supporting information of Song et al., J. Chem. Inf. Model. 2019, 59, 7, 3128-3135; https://doi.org/10.1021/acs.jcim.9b00105
 - "gaff", "cgenff" and "consensus" as supplied online by Gapsys et al., Chem. Sci., 2020, Advance Article; https://doi.org/10.1039/C9SC03754C on
    https://github.com/deGrootLab/pmx/tree/master/protLig_benchmark/ddg_data, downloaded on 23 January 2020
Other files in the subfolders contain the experimental dG values as specified in the sdf files from the supporting information by Wang et al. ("experimental_dG.csv"). 

The file "final_results_all.xlsx" contains all dG values of the different approaches. Results for Kuhn et al. and Gapsys. et al. were calculated from the ddG values provided using Flare's internal free energy analysis tool. dG values for Wang et al. and Song et al. were taken from the supporting information of the corresponding papers (cf. above).

The Jupyter Notebook "recalculate_ddG.ipynb" contains the results for R, MUE, RMSE and Kendall's tau based on the input files supplied in the subfolders.
Note that the error estimates were obtained by bootstrapping and may deviate slightly upon rerunning the notebook.

The Jupyter Notebook "Statistics.ipynb" presents per-link statistics for each data set (as opposed to the per-compound statistics presented elsewhere in this paper)