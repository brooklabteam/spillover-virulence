# spillover-virulence

This repository compiles code for [Brook et al. 2022](https://doi.org/10.1101/2021.10.06.463372), "Reservoir host immunology and life history shape virulence evolution in zoonotic viruses", currently in review. You should be able to reproduce all analyses by cloning this repo to your computer and editing the "homewd" address to your own working directory at the top of each R-script. 

There are nine main folders:

- duration-infection-refs: This folder contains raw .pdf files and a formatted Bibliography (in both Word and .pdf form) summarizing all citations used to procure the 'duration of infection' data with which we supplemented the "stringent_data_guth2021.csv" file (in "source") to convert reported human case fatality rates from zoonoses to measures of infection-induced mortality ('virulence').

- figure 1 : This folder contains the raw powerpoint file used to generate Figure 1 of the main text.

- figure-2: This folder contains all scripts and generates all data to produce the Figure 2 of the main text and Figure S2 in the Supplementary information.  Script "build_data_Fig2.R" codes the equations for $r^ * $, $\beta ^ * $, $\alpha ^ * $, prevalence, and total absolute infected mortality  and generates the two ".Rdata" files in the same folder and script "Fig2_S2.R" then uses those datasets to plot the figures which, respectively, are saved to the "main-figs" and "supp-figs" folders. 

- figure-3: This folder contains all scripts and generates all data to produce the Figure 3 of the main text and Figure S3 in the Supplementary information. See [this detailed brief](/figure-3/within-host-scaling.md) of how we used data from the literature to determine the parameter values used for our within-host model. The script "Fig3.R" loads two outside datasets, "Healy.csv" from [Healy et. al. 2014](https://doi.org/10.1098/rspb.2014.0298)  and "Pantheria.csv" from [Jones et al. 2009](https://doi.org/10.1890/08-1494.1) to produce Figure 3A. The script then loads neutrophil data from the [Species360 database](https://zims.species360.org/) to formulate panel 3B, and publicly available phylogenetic data from [TimeTree](http://www.timetree.org/) to generate panel 3C. After formulating panels A,B, and C, the script "Fig3.R" then loads data produced in the sub-script "brief_script.R" to plot virulence estimates by mammalian order from both the nested model and the literature in panel 3D.

- figure-s1: This folder contains one script only ('FigS1-make-PIP.R') which is self-contained and produces the supplementary figure S1, a pairwise invasibility plot (save to the supp-figs folder).

- main-figs:  This folder contains all of the main figures for the paper, produced in figure-2 and figure-3 folders above, as well as the composie image of our conceptual figure.

- phylo-tree: This folder contains scripts and data from timeTree.org that produce the phylogenetic tree seen in Figure 3 and the measures of phylogenetic distance used throughout the manuscript.

- source: This folder contains source code from [Mollentze and Streicker 2020](https://doi.org/10.1073/pnas.1919176117) for computing and visualizing partial effects from GAMs ('mollentze-streicker-2020-functions.R'). Additionally, it contains a script from [Guth et al. 2022](https://doi.org/10.1073/pnas.211362811) ('guth-2021-cfr-virus-exclude.R') that produces GAM-derived estimates of average case fatality rates from the literature across zoonoses derived from different mammalian orders, then converts these to virulence estimates using duration of infection data which we collected in part with this analysis. The duration of infection data is listed in the dataset called in the source script and included in the folder ("stringent_data_guth2021.csv"). A compiled folder of references cited for those estimates is available within this repository. Finally, the 'guth-2021-cfr-virus-exclude.R' produces two ".Rdata" files which are called by scripts in the figure-4 folder above: "gam.dat.Guth.et.al.2021.Rdata" which includes all zoonoses and "SI.gam.dat.no.rabies.Guth.et.al.2021.Rdata" which excludes rabies lyssavirus.

- supp-figs:  This folder contains all of the supplementary figures for the paper, produced in figure-2, figure-3, figure-4, and figure-s1 folders above.