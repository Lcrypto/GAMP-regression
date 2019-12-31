This code replicates the results in Korobilis (2019) ``High-Dimensional Macroeconomic Forecasting Using Message Passing Algorithms''.

The folder FORECASTING replicates the empirical exercise in the main paper. The code FORECASTING_UR.m estimates and forecasts using the specification in equation (37) and can be used to replicate Tables 1 & 2. The code FORECASTING.m estimates and forecasts using the specification in equation (38) and can be used to replicate Table 3. The code gives the MSFEs and logPLs, and you simply need to take the averages over the forecast evaluation period. Note that in the paper I quote these metrics relative to AR(1), hence, onece the code is finished use the following commands to print MSFEs and log-APLs:
>> mean(MSFE)./mean(MSFE(:,1))
and
>> mean(logPL)- mean(logPL(:,1))
In order to obtain Figures D.3 and D.4 in the Online Appendix you simply need to run the command
>> cumsum(MSFE)
to obtain the cumulative sum of the squared forecast errors.


The folder APPENDIX_C replicates the Monte Carlo simulations in Section C of the Online Appendix. The subfolder MONTE_CARLO_TVP does the first exercise that examines the different forms of time-variation. The subfolder MONTE_CARLO_REG does the second exercise which examines the performance of GAMP in shirnking exogenous predictors in a static regression. MONTE_CARLO_AR is the third exercise that examines the ability of GAMP as an estimator in an AR(4) likelihood. In all three cases, the main document is named MONTE_CARLO.m. Once the file is run and .mat files are saved, I also provide functions that do the boxplots I plot in this Appendix.


The folder APPENDIX_D replicates the two exercises in Online Appendix D.1. Simply run the code VOLATILITY_ESTIMATES.m in order to replicate Figures D.1 and D.2.