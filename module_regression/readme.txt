module_regression

There are two tutorials in this module: one on autoregression and one on regression.

The MATLAB and Octave tutorials are identical in their workflow and concepts presented.* There are minor differences in syntax. In addition, MATLAB and Octave use different algorithms to estimate AR coefficients; thus, those results and downstream analysis (e.g. the estimates of residuals, regression results) will be slightly different.

*Octave currently does not support LASSO regularization (although many other open-source software do, including R). The section on LASSO regularization is skipped in the Octave tutorial "tutorial_regression_octave.m"

The non-tutorial scripts (e.g. the analysis scripts in the "code/" directory) are not currently Octave compatible.



MATLAB tutorials
1. Navigate to the directory “module_introduction/tutorials”.
2. Run the live scripts “tutorial_autoregression.mlx” and "tutorial_regression.mlx". You can evaluate one cell at a time by clicking in the cell and using Matlab’s “Run Section” command (e.g. command+enter on Macs). Alternatively, you can run the entire script by using Matlab’s “Run” command.
3. The “helper_functions” sub-directory must be on Matlab’s path for the tutorial to run. This is taken care of automatically.

Dependencies:	
Econometrics Toolbox for parcorr(), arima(), and estimate()
Statistics and Machine Learning Toolbox for lasso()



OCTAVE tutorials
1. Navigate to the directory “module_introduction/tutorials”.
2. Run the script “tutorial_autoregression_octave.m” and "tutorial_regression_octave.m". The entire script runs at once.
3. The “helper_functions” sub-directory must be on Octave’s path. This is taken care of automatically.

Dependencies:
“nan” and “tsa” packages for pacf(), acorf(), and durlev()
See https://octave.sourceforge.io for instructions on installing Octave packages.



ALOHA 1.0 timeseries analysis (MATLAB only)
1. Navigate to the directory "module_regression/".
2. Add the following directories to the MATLAB path: code/ and data_raw/. Alternatively, you can add the entire module_regression/ directory to the MATLAB path.
3. Run the script "MASTER.m". Note that this will take a while to run. The entire analysis is run from this script via smaller scripts (e.g. "importdata.m").
4. Results are visualized via figure scripts (e.g. "figure_timeseries.m"). The manuscript figures are generated via "generate_manuscript_figures.m".

Dependencies:
Econometrics Toolbox for parcorr(), arima(), and estimate()
Statistics and Machine Learning Toolbox for lasso()



Subdirectories:
tutorials - contains both MATLAB and Octave compatible tutorials
code - main code for analysis of aloha 1.0 timeseries data
data_raw - original aloha 1.0 timeseries data
data_mat - cleaned data and analysis results
figures - all figure outputs from the main code including manuscript figures