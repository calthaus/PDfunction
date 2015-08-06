<<<<<<< HEAD
PDfunction
==============
Contains data and scripts to reproduce the figures of the paper "Time-kill curve analysis and pharmacodynamic functions for in vitro evaluation of antimicrobials against Neisseria gonorrhoeae" (by Sunniva Förster, Magnus Unemo, Lucy J. Hathaway, Nicola Low, Christian L. Althaus)
==============
Software Versions Used

R version 3.2.0 (2015-04-16)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 8 x64 (build 9200)
locale:
[1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.1252

=======
<html>
<h6>Contains data and scripts to reproduce the figures of the paper </h6><br>
<h3>"Time-kill curve analysis and pharmacodynamic functions for in vitro evaluation of antimicrobials against <i>Neisseria gonorrhoeae"</i> </h3><br>
<h6>Sunniva Förster, Magnus Unemo, Lucy J. Hathaway, Nicola Low, Christian L. Althaus</h6><br>
<body><br>
<b>Software Versions Used</b><br>
R version 3.2.0 (2015-04-16)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 8 x64 (build 9200)
<br>
locale:
[1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.1252
<br>
>>>>>>> 0dae818d3fb6451804dabd721603ee6639460e5f
attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     
<br>
other attached packages:
 [1] censReg_0.5-20      maxLik_1.2-4        miscTools_0.6-16    magicaxis_1.9.4     sm_2.2-5.4          plotrix_3.5-12      metafor_1.9-7       Matrix_1.2-0       
 [9] XLConnect_0.2-11    XLConnectJars_0.2-9 drc_2.5-12          MASS_7.3-40        
<br>
loaded via a namespace (and not attached):
[1] Rcpp_0.12.0      nloptr_1.0.4     plyr_1.8.3       tools_3.2.0      lme4_1.1-8       nlme_3.1-120     lattice_0.20-31  mgcv_1.8-6       glmmML_1.0      
[10] parallel_3.2.0   mvtnorm_1.0-3    SparseM_1.6      rJava_0.9-7      gtools_3.5.0     grid_3.2.0       nnet_7.3-9       survival_2.38-1  multcomp_1.4-1  
[19] TH.data_1.0-6    minqa_1.2.4      car_2.0-25       scales_0.2.5     codetools_0.2-11 splines_3.2.0    pbkrtest_0.4-2   colorspace_1.2-6 quantreg_5.11   
[28] sandwich_2.3-3   munsell_0.4.2    zoo_1.7-12
<<<<<<< HEAD
==============
Figure 1
=======


<b>Figure 1</b><br>
>>>>>>> 0dae818d3fb6451804dabd721603ee6639460e5f
Example to illustrate the parameters of the pharmacodynamic function
Example to illustrate the parameters of the pharmacodynamic function
<br>
<b>Figure 2</b><br>
script Figure2.R fits linear regression to data from WHO G, WHO K, WHO L, WHO M, WHO N treated with 11 doubling concentrations of ciprofloxacin and untreated control.
<br>
<b>Figure 3</b><br>
script Figure3.R fits linear regression to one representative experiment of DOGK18 treated with gentamicin, spectinomycin, azithromycin, benzylpenicillin, ceftriaxone, cefixime, chloramphenicole, tetracycline
<br>
<b>Figure 4</b><br>
A: script Figure4A.R fits linear regression to DOGK18 treated with cefixime. Regressions are shown as dotted lines.<br> 
B: script Figure4B.R fits pharmacodynamic function to DOGK18 treated with cefixime. Dotted line shows fit without outliers removed.<br>
C: script Figure4C.R takes model parameters of plots from figure 2 from file singleplot_statistics_panel.xls and plots them as solid lines in different colors for comparison.<br>
D: script Figure4D.R takes model parameters of plots from figure 3 from file summary_statistics_9AB (arithmetic mean of parameters from 2 replicates) and plots them as solid lines in different colors for comparison.<br>
ABCD: Latex script to make panel from A,B,C,D
<br>
<b>Table 1</b><br>
script Pharmacodynamics_DOGK18 fits pharmacodynamic function to two independent experiments of DOGK18 treated with gentamicin, spectinomycin, azithromycin, benzylpenicillin, ceftriaxone, cefixime, chloramphenicole, tetracycline. Table 1 shows the arithmetic mean of the parameters and resistance determinants as published elsewhere.
script Table1 means calculates the means of the two independent experiments.
<br>
<b>Figure S1</b><br> 
script growthcurves_3curves.R fits gompertz model to Growth_GW_short.txt (data from decline phase removed, complete raw data in Growth_GW.txt)
<br>
<b>Table S1</b><br>
Pharmacodynamics_wt.R fits pharmacodynamic function to two independent experiments of DOGK18 treated with gentamicin, spectinomycin, azithromycin, benzylpenicillin, ceftriaxone, cefixime, chloramphenicole, tetracycline and summarizes parameter estimates and standard errors.
<br>
<b>Table S2</b><br>
Pharmacodynamics_panel.R fits pharmacodynamic function to to data from WHO G, WHO K, WHO L, WHO M, WHO N treated with 11 doubling concentrations of ciprofloxacin and untreated control and summarizes parameter estimates and standard errors.
</body>
