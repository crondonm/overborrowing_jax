# Overborrowing and Systemic Externalities in the Business Cycle Under Imperfect Information

## Authors:

Juan David Herreño and Carlos Rondón-Moreno. 

Comments and suggestions are welcome at: crondon[at]bcentral[dot]cl.

## Estimation and Calibration of the Model

This folder contains all the necessary files to replicate the estimation of the parameters in the endowment processes. In particular, it replicates table 1 in the paper.

### Instructions

Run the file `Table1_Calibration.m` in Matlab.

If you want to replicate the full optimization algorithm, set the variable `replicate` to `true`.

Matlab's global optimization toolbox is required to run the code. Full execution time is around 3 hours in a Macbook Pro M1 Max with 32 gigas of RAM.

### Files

| File Name               | Description                                                                 |
|-------------------------|-----------------------------------------------------------------------------|
| `Table1_Calibration.m`  | Main script to replicate Table 1 in the paper. Includes the calibration and optimization algorithms. |
| `LL_klm.m`        | Computes Log-likelihood function based on the Kalman Filter.          |
| `LL_klm_opt.m`    | Evaluates Log-likelihood function at the optimal point. Needed to compute the standard errors of the estimated parameters.      |
| `Hessiancsd.m`  | Cao (2008): Computes numerical hessian.          |
| `simannb.m`        | Goffe (1999) and modified by Soest (1999) and Su (2000): Executes Simulated Annealing .                  |
| `Data.xlsx`        | Database used for the paper. See file for sources |
| `Calibration.mat`  | Contains results as presented in the paper.                  |

### Table 2

Table 2 summarizes the calibrated parameters along with those selected based on Bianchi (2011).

## Notes

Ensure that all files are in the same directory before running the scripts.