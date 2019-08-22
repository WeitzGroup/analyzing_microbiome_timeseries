module_introduction

The MATLAB and Octave tutorials are identical in their analysis and figures. There are minor differences in the syntax when computing the partial autocorrelation (MATLAB and Octave each have their own toolbox or package). In addition, Octave requires packages to be loaded explicitly.

The non-tutorial scripts (e.g. figure_autoregression.m and randomwalk_ensemble.m) are not currently Octave compatible.



MATLAB instructions
1. Navigate to the directory “module_introduction”.
2. Run the live script “tutorial_autoregression.mlx”. You can evaluate one cell at a time by clicking in the cell and using Matlab’s “Run Section” command (e.g. command+enter on Macs). Alternatively, you can run the entire script by using Matlab’s “Run” command.
3. The “helper_functions” sub-directory must be on Matlab’s path for the tutorial to run. This is taken care of automatically (line 2).

Dependencies:	
Econometrics Toolbox for the parcorr() function (line 46). You can skip this cell (lines 45-48), and the rest of the tutorial will still run.



OCTAVE instructions
1. Navigate to the directory “module_introduction”.
2. Run the script “tutorial_autoregression_octave.m”.  The entire script runs at once.
3. The “helper_functions” sub-directory must be on Octave’s path. This is taken care of automatically (line 15).

Dependencies:
“nan” and “tsa” packages for the pacf() function (line 77). If you want to skip this part, comment out lines 13-14 and lines 76-87 before running the script. 
See https://octave.sourceforge.io for instructions on installing Octave packages.



Files in this directory:
tutorial_autoregression.mlx - main tutorial script for MATLAB
tutorial_autoregression_octave.m - main tutorial script for Octave
figure_autoregression.m - generates the manuscript figure
randomwalk_ensemble/randomwalk_ensemble.m - script for generating distribution of correlation values
randomwalk_ensemble/randomwalk_ensemble_originals.mat - ensemble results for original timeseries (already generated)
randomwalk_ensemble/randomwalk_ensemble_residuals.mat - ensemble results for residual timeseries (already generated)
helper_functions/label_ax.m - adds panel labels to manuscript figure
helper_functions/plot_correlation.m - plots correlation matrix as a heatmap
helper_functions/plot_correlation_example.m - plots two timeseries and the correlation between them
helper_functions/plot_distribution.m - plots distribution of correlation values (generated via randomwalk_ensemble.m)
helper_functions/redbluecmap.mat - colormap for plot_correlation.m