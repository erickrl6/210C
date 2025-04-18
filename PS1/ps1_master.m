% run_all.m

disp('Plotting True IRF...');
ps1_2_true_irf;

disp('Plotting misspecified IRFs with output lags and contemporaneous monetary shock...');
ps1_3_misspec_irf;

disp('Plotting misspecified IRF with 6 lags of monetary shocks ...');
ps1_4_eps6_irf;

disp('Plotting misspecified IRF with 6 lags of monetary base...');
ps1_4_monetary_irf;

disp('Plotting Jorda Local Projections...');
ps1_5_Jorda_irf;

disp('Plotting Jorda Local Projections with one output lag...');
ps1_5_Jorda_cont_irf;

disp('All scripts completed.');
