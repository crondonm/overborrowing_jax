%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information
%
% In this code: Main File. 
%               1. Replicate full set of results
%               2. This code was fully executed in the High-Performace
%               Clúster at the University of Notre Dame. 
%               3. Memory RAM requiered: 1TB
%               4. Full execution takes about 8hrs in the latest hardware
%               available on March 2025.
%               
% Authors:  Juan Herreño, jherrenolopera@ucsd.edu
%               Carlos Rondón Moreno, crondon@bcentral.cl
%
% Date: March 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping

clearvars
clear global
close all

%% Create parameter matrix and compute transition matrix...

run Step0_Parameters.m
run Step1_transition_matrix.m

%% Solve the model ...

run Step2_dcntrlzd_eqm_imp_info.m
run Step3_planner_imp_info.m
run Step4_dcntrlzd_eqm_pfct_info.m
run Step5_planner_prfct_info.m

%% Simulate welfare costs...

run Step6_Welfare_costs.m

%% Reproduce Figures 

run '../Figures/Figure1.m'
run '../Figures/Figure2.m'
run '../Figures/Figure3a_3b.m'
run '../Figures/Figure4_5.m'
run '../Figures/Figure6.m'
run '../Figures/Figure7_8.m'
run '../Figures/Figure9.m'
run '../Figures/Figure10.m'
run '../Figures/Figure11a_11b.m'
run '../Figures/Figure12.m'
run '../Figures/Figure13.m'

% Appendix Figures

run '../Figures/Appendix/FigureA1.m'
run '../Figures/Appendix/FigureA2.m'
run '../Figures/Appendix/FigureA3.m'
run '../Figures/Appendix/FigureA4_A5.m'

%% Reproduce Tables

% Table 1: 
run '../Tables/Table1.m'
% Table 2: 
run '../Tables/Table2.m'
% Table 3:
run '../Tables/Table3.m'
% Table 4:
run '../Tables/Table4.m'
% Table A1: Appendix
run '../Tables/TableA1.m'





