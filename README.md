# rstarBrookings2017

Replication files for
[*Safety, Liquidity, and the Natural Rate of Interest*](https://www.brookings.edu/bpea-articles/safety-liquidity-and-the-natural-rate-of-interest/)
by Marco del Negro, Domenico Giannone, Marc Giannoni, and Andrea Tambalotti,
*Brookings Papers on Economic Activity*, Spring 2017: 235-294.


## Required software

- Julia v0.5.0 or above
- DSGE.jl v0.3.1
- MATLAB 16a

**Download instructions**

1. Download Julia from `https://julialang.org/downloads/`.
2. Open the Julia REPL and type:

   a. `Pkg.add("DSGE")` to install DSGE.jl

   b. `Pkg.pin("DSGE", v"0.3.1")` to use DSGE.jl v0.3.1

   c. If, after running this replication code, you would like to use the most
      current version of DSGE.jl, type `Pkg.free("DSGE")` to un-pin the version


## Installing this repository

Git users are welcome to fork this repository and clone it for local
use. Non-Git users will probably find it easiest to download the zip
file by clicking on the green `Clone or download` button on the right
hand side of this screen, and then clicking "Download ZIP".


## Directory structure

- `dsge/`: Julia replication code for DSGE model

  - `spec1010_20_2016Q3_1221.jl`: main script that generates all results

  - `input_data`:
    - `data/`: transformed input data
	- `user/`: precomputed mode and hessian files

  - `output_data/m1010/ss20/`:
    - `estimate/`:
	  - `raw/`: raw estimation outputs: parameter draws from Metropolis-Hastings
	  - `tables/`: LaTeX tables of parameter moments

    - `forecast/`:
	  - `raw/`: full distribution, raw results (in model units) from
                all post-estimation products
      - `work/`: means and bands from raw results, transformed into
                 final units, in binary format
      - `tables/`: LaTeX tables of post-estimation results (shock
                   decomopositions, parameter histories)

- `tvar/`: results for TVAR model
  - `MainModelX.m`: main scripts that generate all the results
  - `DataCompleteLatest.xls`: input data
  - `FiguresModelX/`: output figures for each model specification
  - `Routines/`: functions called to run estimation and produce figures
  - `output_data/`: output files for all TVAR models

- `plot/`: MATLAB code for plotting all results
  - `makeRstarPlots.m`: Main driver script
  - `helperFunctions/`: MATLAB functions called by `makeRstarPlots`
  - `Figures/`: Output figure directory
  - `Tables/`: Input table directory


## How to run the DSGE code

The script `spec1010_20_2016Q3_1221.jl` generates all results for the
DSGE model. As provided, the script will create a model object with
the appropriate settings to estimate the model, compute smoothed
histories of pseudoobservables, and produce shock decompositions
("pseudoobservables" is the term we use to describe linear
transformations of states that we are interested in. In the Julia REPL
with the DSGE package loaded, type `?pseudo_measurement` for more
information.

Three boolean variables at the top of the script indicate the
operations the script will perform next. If `run_estimation = true`,
the code will load the pre-computed mode and hessian files and run
Metropolis-Hastings (be aware that this runs sequentially and could
take more than 20 hours). To compute a mode and hessian from scratch,
set `reoptimize` and `calculate_hessian` to `true` under "Settings for
estimation" in the script. This will take much longer.

If `run_modal_forecast = true`, the code will load the provided
parameter mode and compute smoothed pseudoobservable histories and
shock decompositions. This runs quite quickly.

If `run_full_forecast = true`, the code will load all parameter draws
from the estimation and compute smoothed pseudoobservables and shock
decomositions for each draw. It will then compute means and bands
across all draws. Since we compute results for 20,000 parameter draws,
we parallelize this computation across 50 workers. Blocks of 500 draws
each are sent to the workers, which compute and record results
draw-by-draw. If you are not working with a cluster, set `nworkers =
1` in `spec1010_20_2016Q3_1221.jl` (under "Parallelization"). If you
have a cluster, set `addprocsfcn` appropriately for your machine (see
ClusterManagers.jl).


## How to run the TVAR code

The `MainModelX.m` scripts generate results for TVAR model (each script named
according to the model number `X` specified in the paper). As provided, each
script estimates the TVAR model for the given specification, outputs a `.mat`
file of estimation results, produces figures for the distribution of trends, and
prints the change in trends for the variables specified in the model.

If `RunEstimation = 1`, the code will run the estimation for the specified
model. If `RunEstimation = 0`, the code will load the estimation results from a
previous estimation. The first time the code is run, `RunEstimation` should be
set to `1` to run the estimation and produce the necessary results.

We ran the estimation scripts in MATLAB R2016a.


## How to run plots

See the README.md file in the `plot` directory. All results are plotted using
MATLAB 16a.