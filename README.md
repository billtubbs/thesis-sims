# thesis-sims

MATLAB code and data to reproduce the simulation results in my Laval University masters thesis in electrical and computer engineering:

 - [Multiple-model observers for detecting ore feed disturbances in grinding operations](http://hdl.handle.net/20.500.11794/118186), Masters thesis, 2023.

<img src="images/rod_obs_sim_resp_plot1.png" width="60%">

For the code to reproduce the results of the following conference paper, refer to this repository instead [https://github.com/billtubbs/ifac-2022-mmkf](https://github.com/billtubbs/ifac-2022-mmkf).

- William Tubbs, André Desbiens, Jocelyn Bouchard, An Observer to Detect Infrequently-Occurring Disturbances in Grinding Operations, IFAC-PapersOnLine, Volume 55, Issue 21, 2022, Pages 13-18, ISSN 2405-8963, https://doi.org/10.1016/j.ifacol.2022.09.236.

The code in this repository has been tested with MATLAB versions 2019b, 2020b, and 2021b.  It may or may not work with other versions!


## 1. Generating RODD disturbances (section 3.1 of thesis report)

The files for these simulations are currently in this repository:
 - https://github.com/billtubbs/thesis-report/


## 2. Observer simulations with linear systems (section 3.2 of thesis report)

The files for these simulations are in the [`linear-sims`](linear-sims) sub-directory.  Navigate to this directory and then follow the instructions below.


### System models

The two systems considered are defined in the following files
 - [linear-sims/sys_rodin_step.m](linear-sims/sys_rodin_step.m) - SISO linear system described in section 3.2.1 and used in the simulations in section 3.2.2
 - [linear-sims/sys_rodin_step_2x2sym2.m](linear-sims/sys_rodin_step_2x2sym2.m) - 2x2 linear system used in section 3.2.3


### Process observers

 - [process-observers/obs_rodin_step.m](process-observers/obs_rodin_step.m) - observers with un-optimized parameters used for the SISO simulations at the beginning of section 3.2.2
 - [linear-sims/obs_rodin_step_opt.m](linear-sims/obs_rodin_step_opt.m) - observers with optimized parameters for the SISO simulations at the end of section 3.2.2
 - [linear-sims/obs_rodin_step_2x2_opt.m](linear-sims/obs_rodin_step_2x2_opt.m) - observers with with optimized parameters for the 2x2 linear system simulations in section 3.2.3


### Observer parameter tuning experiments with SISO linear system (section 3.2.2)

The observer parameters were chosen by running multiple simulations with different combinations of parameter values.  Run the following scripts in the sequence shown to generate the latex for the summary tables of observer parameter search results.  Note: these simulations can take a long time to run.

For the Kalman filter (KF3):
 - [linear-sims/gen_sim_specs_sim1_3KF_Q.m](linear-sims/gen_sim_specs_sim1_3KF_Q.m)
 - [linear-sims/run_obs_sims.m](linear-sims/run_obs_sims.m) with the line `sim_name = "rod_obs_sim1_3KF_seed"` uncommented
 - [linear-sims/rod_obs_sim1_plot_KF3_rmse.m](linear-sims/rod_obs_sim1_plot_KF3_rmse.m)

For the MKF_SF observer (1995 version):
 - [linear-sims/gen_sim_specs_sim1_MKF_SF95_popt.m](linear-sims/gen_sim_specs_sim1_MKF_SF95_popt.m)
 - [linear-sims/run_obs_sims.m](linear-sims/run_obs_sims.m) with the line `sim_name = "rod_obs_sim1_MKF_SF95_popt"` uncommented
 - [linear-sims/rod_obs_sim1_MKF_SF95_popt_table.m](linear-sims/rod_obs_sim1_MKF_SF95_popt_table.m)

For the MKF_SF observer (1998 version):
 - [linear-sims/gen_sim_specs_sim1_MKF_SF98_popt.m](linear-sims/gen_sim_specs_sim1_MKF_SF98_popt.m)
 - [linear-sims/run_obs_sims.m](linear-sims/run_obs_sims.m) with the line `sim_name = "rod_obs_sim1_MKF_SF_popt"` uncommented
 - [linear-sims/rod_obs_sim1_MKF_SF98_popt_table.m](linear-sims/rod_obs_sim1_MKF_SF98_popt_table.m)

For the MKF_SP observer:
 - [linear-sims/gen_sim_specs_sim1_MKF_SP_popt.m](linear-sims/gen_sim_specs_sim1_MKF_SP_popt.m)
 - [linear-sims/run_obs_sims.m](linear-sims/run_obs_sims.m) with the line `sim_name = "rod_obs_sim1_MKF_SP_popt"` uncommented
 - [linear-sims/rod_obs_sim1_MKF_SP_popt_table.m](linear-sims/rod_obs_sim1_MKF_SP_popt_table.m)


### Observer evaluation simulations with SISO linear system (section 3.2.2)

Run the following scripts (these take a while):
 - [linear-sims/gen_sim_specs_sim1_all_seed.m](linear-sims/gen_sim_specs_sim1_all_seed.m)
 - [linear-sims/run_obs_sims.m](linear-sims/run_obs_sims.m) with the line `sim_name = "rod_obs_sim1_all_seed"` uncommented
 - [linear-sims/rod_obs_sim1_MKF_SP_popt_table.m](linear-sims/rod_obs_sim1_MKF_SP_popt_table.m)


### Observer parameter tuning experiments with 2x2 linear system (section 3.2.3)

For the Kalman filter (KF3):
 - [linear-sims/gen_sim_specs_sim2_3KF_Q.m](linear-sims/gen_sim_specs_sim2_3KF_Q.m)
 - [linear-sims/run_obs_sims.m](linear-sims/run_obs_sims.m) with the line `sim_name = "rod_obs_sim2_3KF_seed"` uncommented
 - [linear-sims/rod_obs_sim2_plot_KF3_rmse.m](linear-sims/rod_obs_sim2_plot_KF3_rmse.m)

For the MKF_SF observer (1995 version):
 - [linear-sims/gen_sim_specs_sim2_MKF_SF95_popt.m](linear-sims/gen_sim_specs_sim2_MKF_SF95_popt.m)
 - [linear-sims/run_obs_sims.m](linear-sims/run_obs_sims.m) with the line `sim_name = "rod_obs_sim2_MKF_SF95_popt"` uncommented
 - [linear-sims/rod_obs_sim2_MKF_SF95_popt_table.m](linear-sims/rod_obs_sim2_MKF_SF95_popt_table.m)

For the MKF_SF observer (1998 version):
 - [linear-sims/gen_sim_specs_sim2_MKF_SF98_popt.m](linear-sims/gen_sim_specs_sim2_MKF_SF98_popt.m)
 - [linear-sims/run_obs_sims.m](linear-sims/run_obs_sims.m) with the line `sim_name = "rod_obs_sim2_MKF_SF_popt"` uncommented
 - [linear-sims/rod_obs_sim2_MKF_SF98_popt_table.m](linear-sims/rod_obs_sim2_MKF_SF98_popt_table.m)

For the MKF_SP observer:
 - [linear-sims/gen_sim_specs_sim2_MKF_SP_popt.m](linear-sims/gen_sim_specs_sim2_MKF_SP_popt.m)
 - [linear-sims/run_obs_sims.m](linear-sims/run_obs_sims.m) with the line `sim_name = "rod_obs_sim2_MKF_SP_popt"` uncommented
 - [linear-sims/rod_obs_sim2_MKF_SP_popt_table.m](linear-sims/rod_obs_sim2_MKF_SP_popt_table.m)


### Observer evaluation simulations with 2x2 linear system (section 3.2.3)

Run the following scripts (these take a while):
 - [linear-sims/gen_sim_specs_sim2_all_seed.m](linear-sims/gen_sim_specs_sim2_all_seed.m)
 - [linear-sims/run_obs_sims.m](linear-sims/run_obs_sims.m) with the line `sim_name = "rod_obs_sim2_all_seed"` uncommented
 - [linear-sims/rod_obs_sim2_MKF_SP_popt_table.m](linear-sims/rod_obs_sim2_MKF_SP_popt_table.m)


### Plots

To make the various plots of the simulation results in section 3.2.2, run the following script, with the line `sim_name = "rod_obs_sim1_all_seed";` uncommented.

Also uncomment one of the lines in the section `choose observers to include in plots` to determine which observer results to include int he plots.  For example, to make the box plot in Figure 3.20, use the following:'
```lang-matlab
obs_sel_labels = {'KF3', 'MKF_SF95', 'MKF_SF1', 'MKF_SP1', 'SKF'};
```

To produce the plot 'Effect of random variables on the RMSE results=' in Fig. A.2 in the Appendix, run the following in sequence:
 - [linear-sims/gen_sim_specs_sim1_3KF_seed.m](linear-sims/gen_sim_specs_sim1_3KF_seed.m)
 - [linear-sims/run_obs_sims.m](linear-sims/run_obs_sims.m) with the line `sim_name = "rod_obs_sim1_3KF_seed";` uncommented
 - [linear-sims/rod_obs_sim_crmse_plot.m](linear-sims/rod_obs_sim_crmse_plot.m)


## 3. Observer simulations with grinding process simulator (section 3.3 of thesis report)

The files for these simulations are in the [`grind-sims`](grind-sims) sub-directory.  Navigate to this directory and then follow the instructions below.


### Data from grinding simulation model

The grinding simulation model used in these simulations is not included.  Instead, input-output data from simulations with this model is saved in the following directory:

The [data](data) subdirectory contains time-series data sets from 15 simulations of a grinding process model (the model itself is not available here).  Data sets 1 to 5 contain short simulations of 300 time steps.  Simulations 6 to 15 contain longer simulations of 2460 time steps.  The files contain data for 7 process variables although only two were used in this work (BASE_ORE_MIX and SAG_OF_P80_M).

| #  | Filename                                    | Use          |
| -- | ------------------------------------------- | ------------ |
| 1  | sim_OL_rc_est_mix_factor_300_1_ident.csv    | Process model estimation (Fig. 4 in paper)  |
| 2  | sim_OL_rc_est_mix_factor_300_2_ident.csv    | Process model validation (model selection)    |
| 3  | sim_OL_rc_est_mix_factor_300_3_ident.csv    | Initial observer test (Fig. 5 in paper)    |
| 4  | sim_OL_rc_est_mix_factor_300_4_ident.csv    | Not used     |
| 5  | sim_OL_rc_est_mix_factor_300_5_ident.csv    | Observer parameter optimization     |
| 6 ... 15  | sim_OL_rc_est_mix_factor_2460_6_ident.csv ... sim_OL_rc_est_mix_factor_2460_15_ident.csv  | Observer evaluation (RMSE results in Table 2 and Fig. 6 in paper)    |


### Process observers

 - [grind-sims/rod_obs_P2DcTd4.m](grind-sims/rod_obs_P2DcTd4.m) - observers with the identified linear system model and optimized parameters used for the grinding simulations in section 3.2.3.
 - [grind-sims/rod_obs_P2Dcd1_T.m](grind-sims/rod_obs_P2Dcd1_T.m) - these observers have a system model identified from simulation data without measurement nosie. This model provides better predictions of the true simulation outputs than the model identified from noisy measurements. However, results using these observers were not included in the thesis results.


### Instructions to reproduce the results

Open the script [rod_obs_sim.m](rod_obs_sim.m) and specify the input data sequences to include in the simulations in lines 55-57.  For example, specify the first 5 as follows:

```lang-matlab
i_in_seqs = [1, 2, 3, 4, 5];
```

To change which observers are included in the simulations, edit line 82:

```lang-matlab
observers = {KF1, KF2, KF3, MKF_SF95, MKF_SF1, MKF_SP1, SKF};
```

The observers are defined in separate script files.

Running [rod_obs_sim.m](rod_obs_sim.m) with the above settings should produce the following output:

```lang-none
Starting observer simulations with input seq. #1 ...
Observer simulation results saved to file: rod_obs_sim_1_1.csv
MKF simulation results saved to file: rod_obs_sim_1_1_MKF_SF95.csv
MKF simulation results saved to file: rod_obs_sim_1_1_MKF_SF1.csv
MKF simulation results saved to file: rod_obs_sim_1_1_MKF_SP1.csv
                                  KF1       KF2        KF3      MKF_SF95    MKF_SF1    MKF_SP1      SKF
                                _______    ______    _______    ________    _______    _______    _______

    RMSE                         3.5629    2.9746     2.1599     2.0206       2.014     2.1431     1.2235
    RMSE in transitions           3.937    2.7732     2.6771     2.6254      2.5853      2.814     1.5633
    RMSE in steady-state         3.1534    3.1596     1.4876     1.1525      1.2169     1.1546    0.75418
    Variance in steady-state    0.67006    8.7855     1.5117      1.369      1.3116     1.1122    0.36806
    RMSD in steady-state        0.26817    3.0764    0.67384    0.59802     0.59391    0.57502    0.29097

Existing results loaded from file: rod_obs_sim_1_summary.csv
Summary results saved to file: rod_obs_sim_1_summary.csv
Step responses identified: 1
Existing step responses loaded from file: rod_obs_sim_1_resps.csv
Step responses saved to file: rod_obs_sim_1_resps.csv

Starting observer simulations with input seq. #2 ...
Observer simulation results saved to file: rod_obs_sim_1_2.csv
MKF simulation results saved to file: rod_obs_sim_1_2_MKF_SF95.csv
MKF simulation results saved to file: rod_obs_sim_1_2_MKF_SF1.csv
MKF simulation results saved to file: rod_obs_sim_1_2_MKF_SP1.csv
                                  KF1       KF2        KF3      MKF_SF95    MKF_SF1    MKF_SP1      SKF
                                _______    ______    _______    ________    _______    _______    ________

    RMSE                         3.5217    3.0487     2.3228     2.4687      2.4607     2.5192      1.7649
    RMSE in transitions          3.9721    3.3371     2.9114     3.1844      3.1581     3.2259       2.216
    RMSE in steady-state         2.8002    2.6088     1.1098    0.79101     0.87124    0.92697     0.82963
    Variance in steady-state    0.18654     6.595     0.6185    0.15753     0.24418     0.2708    0.089818
    RMSD in steady-state        0.34913    2.5483    0.59618    0.42227     0.46404    0.44366     0.37122

Existing results loaded from file: rod_obs_sim_1_summary.csv
Summary results saved to file: rod_obs_sim_1_summary.csv
Step responses identified: 1
Existing step responses loaded from file: rod_obs_sim_1_resps.csv
Step responses saved to file: rod_obs_sim_1_resps.csv

Starting observer simulations with input seq. #3 ...
Observer simulation results saved to file: rod_obs_sim_1_3.csv
MKF simulation results saved to file: rod_obs_sim_1_3_MKF_SF95.csv
MKF simulation results saved to file: rod_obs_sim_1_3_MKF_SF1.csv
MKF simulation results saved to file: rod_obs_sim_1_3_MKF_SP1.csv
                                  KF1       KF2        KF3      MKF_SF95    MKF_SF1    MKF_SP1      SKF
                                _______    ______    _______    ________    _______    _______    _______

    RMSE                         3.5484    3.0647     2.1525      2.187      2.2579     2.1788     1.5972
    RMSE in transitions          4.2352     3.178     3.1659     3.3339      3.3265     3.3047     2.3832
    RMSE in steady-state         3.0623    2.9946     1.1804    0.98834      1.2294      1.018     0.8188
    Variance in steady-state     1.3538     10.76      1.795    0.89845      1.5499    0.99461    0.26255
    RMSD in steady-state        0.21132    2.8264    0.60988    0.50948     0.66169    0.58495    0.24328

Existing results loaded from file: rod_obs_sim_1_summary.csv
Summary results saved to file: rod_obs_sim_1_summary.csv
Step responses identified: 2
Existing step responses loaded from file: rod_obs_sim_1_resps.csv
Step responses saved to file: rod_obs_sim_1_resps.csv

Starting observer simulations with input seq. #4 ...
Observer simulation results saved to file: rod_obs_sim_1_4.csv
MKF simulation results saved to file: rod_obs_sim_1_4_MKF_SF95.csv
MKF simulation results saved to file: rod_obs_sim_1_4_MKF_SF1.csv
MKF simulation results saved to file: rod_obs_sim_1_4_MKF_SP1.csv
                                  KF1       KF2        KF3      MKF_SF95    MKF_SF1    MKF_SP1      SKF
                                _______    ______    _______    ________    _______    _______    _______

    RMSE                         3.4658    3.3323     2.0001     1.5196      1.7333     1.6168     1.1457
    RMSE in transitions          4.5037    3.3486     2.9724     2.3109      2.6663     2.5423     1.7677
    RMSE in steady-state         2.8932    3.3251     1.3684    0.98738      1.0944    0.95898    0.71758
    Variance in steady-state    0.81292    11.633     1.5673    0.99212      1.1773    0.72147     0.1773
    RMSD in steady-state        0.24586    2.7882    0.61642    0.59055     0.60218    0.51876    0.27521

Existing results loaded from file: rod_obs_sim_1_summary.csv
Summary results saved to file: rod_obs_sim_1_summary.csv
Step responses identified: 2
Existing step responses loaded from file: rod_obs_sim_1_resps.csv
Step responses saved to file: rod_obs_sim_1_resps.csv

Starting observer simulations with input seq. #5 ...
Observer simulation results saved to file: rod_obs_sim_1_5.csv
MKF simulation results saved to file: rod_obs_sim_1_5_MKF_SF95.csv
MKF simulation results saved to file: rod_obs_sim_1_5_MKF_SF1.csv
MKF simulation results saved to file: rod_obs_sim_1_5_MKF_SP1.csv
                                  KF1       KF2        KF3      MKF_SF95    MKF_SF1    MKF_SP1      SKF
                                _______    ______    _______    ________    _______    _______    ________

    RMSE                         3.0909    3.0379     2.0435     2.3128      2.1737     2.414       1.6779
    RMSE in transitions          3.7779    3.2699     2.5635     3.0196      2.8178    2.9945       2.2029
    RMSE in steady-state         1.9567    2.7302     1.1266    0.88113     0.90784    1.4206      0.58628
    Variance in steady-state    0.76163    9.5025     1.3601    0.69641     0.67505    1.1633     0.070481
    RMSD in steady-state        0.35642    2.5906    0.62558    0.55192     0.52433    0.6603      0.36777

Existing results loaded from file: rod_obs_sim_1_summary.csv
Summary results saved to file: rod_obs_sim_1_summary.csv
Step responses identified: 0
Existing step responses loaded from file: rod_obs_sim_1_resps.csv
Step responses saved to file: rod_obs_sim_1_resps.csv
run rod_obs_sim_plots.m to produce plots.
run rod_obs_calc_metrics.m to calculate evaluation metrics.
run rod_obs_step_plots.m to produce step response summary plot.
```

### Simulation results

After running [rod_obs_sim.m](rod_obs_sim.m), the results of the simulations are saved as CSV files to the [results](results) subdirectory. For example:

- rod_obs_sim_1_1_MMKF.csv
- rod_obs_sim_1_1.csv
- rod_obs_sim_1_2_MMKF.csv
- rod_obs_sim_1_2.csv
- rod_obs_sim_1_3_MMKF.csv
- rod_obs_sim_1_3.csv
- rod_obs_sim_1_4_MMKF.csv
- rod_obs_sim_1_4.csv
- rod_obs_sim_1_5_MMKF.csv
- rod_obs_sim_1_5.csv
- rod_obs_sim_1_resps.csv
- rod_obs_sim_1_summary.csv

Explanation of output results files:
- The files 'rod_obs_sim_1_1.csv', 'rod_obs_sim_1_2.csv', ... etc. contain the state estimates, output estimates and output estimation errors of each observer for the duration of each simulation.
- The files 'rod_obs_sim_1_1_MMKF.csv', 'rod_obs_sim_1_2_MMKF.csv', ... etc. contain additional data on the multi-model observer (MMKF), such as the state estimates and conditional probabilities of each of the observer's Kalman filters.
- The file 'rod_obs_sim_1_resps.csv' contains the data used to create the plot of observer responses to shocks (Fig. 6 in the paper)
- The file 'rod_obs_sim_1_summary.csv' will contain a records of all the simulation parameters, model parameters, observer parameters, and overall RMSE metrics for each simualtion. This file is not over-written by 'rod_obs_sim.m'. Every time a new simulation is run, a new row is added to 'rod_obs_sim_1_summary.csv'. ***Before re-running simulations, remove the previous results from this file or erase it completely, otherwise, some results may be duplicated.***


### Plots

To produce the plots shown in the report, run simulations for all the datasets (1 to 15) and then run the scripts [rod_obs_sim_plots.m](rod_obs_sim_plots.m) and [rod_obs_step_plots.m](rod_obs_step_plots.m).

After running these scripts, images of the plot figures will be saved in the [plots](plots) folder in pdf format:

- rod_obs_sim_1_ioplot.pdf
- rod_obs_sim_3_est.pdf
- rod_obs_sim_resp_plot1.pdf
- rod_obs_sim_resp_plot2.pdf


### Evaluation metrics

To calculate the evaluation metrics in Table 3.8 of the thesis report, run the script [rod_obs_calc_metrics.m](rod_obs_calc_metrics.m).  This should produce the following output:

```lang-none
Simulation results loaded from file: rod_obs_sim_1_summary.csv
Results for the following simulations found:
     6     7     8     9    10    11    12    13    14    15

Observer performance metrics
                                  KF1       KF2        KF3      MKF_SF95    MKF_SF1    MKF_SP1      SKF
                                _______    ______    _______    ________    _______    _______    _______

    RMSE                         3.2491    3.2575     1.8059     1.6982      1.6558     1.6719     1.2441
    RMSE in transitions          4.4357    3.2636     2.6902     2.9474      2.8302     2.9277      2.041
    RMSE in steady-state         2.8276    3.2547     1.4636     1.1186      1.1225     1.0807    0.89732
    Variance in steady-state     1.8572    10.064     1.6117    0.66163     0.69619    0.57594     0.2463
    RMSD in steady-state        0.13665    2.9053    0.61527    0.43946     0.45231    0.43034    0.16967

Latex table code:
\hline
% See script rod_obs_calc_metrics.m
% 16-Jun-2023 17:26:43 results with system P2DcTd4, sigma_M = 5, tau_ss = 1.2
RMSE($\hat{\mathbf{Y}},\mathbf{Y}$) overall & 1.81 & 1.70 & 1.66 & 1.67 & 1.24 \\
RMSE($\hat{\mathbf{Y}},\mathbf{Y}$) transient & 2.69 & 2.95 & 2.83 & 2.93 & 2.04 \\
RMSE($\hat{\mathbf{Y}},\mathbf{Y}$) steady-state & 1.46 & 1.12 & 1.12 & 1.08 & 0.90 \\
Var($\hat{\mathbf{Y}}$) steady-state & 1.61 & 0.66 & 0.70 & 0.58 & 0.25 \\
RMSD($\hat{\mathbf{Y}},\mathbf{Y}$) steady-state &0.62 & 0.44 & 0.45 & 0.43 & 0.17 \\
\hline
```


### Sensitivity analysis

To run the simulations for the sensitivity analysis in section 3.3.4, run the following script:

```lang-matlab
rod_obs_sim_sens_P2DcTd4
```

There is also a script to run the sensitivity analysis with observers that have the improved system model, but these results were not included in the thesis as they were not significantly different:

```lang-matlab
rod_obs_sim_sens_P2Dcd1_T
```

To produce the heat-map plots in figures 3.33 and 3.34, run the following Python notebook:
 - [grind-sims/Heatmap-plots-of-sensitivity-results.ipynb](grind-sims/Heatmap-plots-of-sensitivity-results.ipynb)


## Unit tests

A set of test scripts are included in each of the main sub-directories, to verify that the main sub-routines are working correctly.  To run the tests run the following command from each sub-directory of the main repository.

```lang-matlab
runtests
```
