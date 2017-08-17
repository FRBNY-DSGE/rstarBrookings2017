# rstarBrookings2017

Replication files for *Safety, Liquidity and the Natural Rate of
Interest* by Marco del Negro, Domenico Giannone, Marc Giannoni, and
Andrea Tambalotti, presented at Brookings in March 2017

## Notes

The "forecast step" (which produces smoothed histories of states and
shock decompositions) for the full distribution of DSGE parameter
draws was run in parallel on 50 workers. Each worker produced results
for a block of 500 parameter draws. See `How to run the DSGE code`
below for details on how to adjust for your machine.

All results are plotted using MATLAB 16a.

## Required software
- Julia v0.5.0 or above
- DSGE.jl v0.2.0 or above
- MATLAB 16a

**Download instructions**

1. Download Julia from `https://julialang.org/downloads/`.
2. Open the Julia REPL and type `Pkg.add("DSGE")`

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

## How to run plots

See the README.md file in the `plot` directory.
