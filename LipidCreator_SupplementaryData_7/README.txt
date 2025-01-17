Requirements:
 - R v3.5 or higher
 - if you run on linux, install
   - libcurl4-openssl-dev
   - libxml2-dev
   - xml2
   - libssl-dev
 - tidyverse (package for R)
 - svglite (package for R)
 - ggbeeswarm (package for R)
 - ggpmisc (package for R)

Run the script:
 1) Open a console or the R prompt or Rstudio
 2) Change into this directory
 3) Run the script with:
  - on Windows: Rscript.exe plasma-evaluation.R
  - on Linux or OS X: Rscript plasma-evaluation.R

 4) The script creates five plots in total within this directory:
    two for Figure 6a and b and three for Supplementary Information Figure 9 (pdf and svg and png version).

The script uses two files: 

 1) "Bowden-MEDM-all-species.csv", which is a list of Median of Means values and standard uncertainty values
    from the Bowden et al. paper's supplementary tables, containing results of the NIST SRM 1950 ring trial.
    The MEDM values were measured by at least five labs.
 2) "LipidCreatorVal_QC-filtered_SAMPLES-NIST_uM_04022020v1", this is the table of measurements with absolute lipid quantities also including the NIST SRM 1950 sample quantities. The file "LipidCreatorVal_BakerPanel_minArea1000_CV20_nmol-L_V1.csv", which was submitted to MetaboLights within the dataset MTBLS1375 only contains the human sample quantities.
