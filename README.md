# rstarBrookings2017

Replication files for
[*Safety, Liquidity, and the Natural Rate of Interest*](https://www.brookings.edu/bpea-articles/safety-liquidity-and-the-natural-rate-of-interest/)
by Marco del Negro, Domenico Giannone, Marc Giannoni, and Andrea Tambalotti,
*Brookings Papers on Economic Activity*, Spring 2017: 235-294.


## Updated r* estimates

The [VAR Excel
file](https://github.com/FRBNY-DSGE/rstarBrookings2017/blob/master/update/TVAR_Comparison.xlsx)
contains the original as well as updated estimates of the trends in the real
return for safe and liquid assets, the convenience yield, and its safety and
liquidity components computed using the VAR model. The [VAR Vintages
file](https://github.com/FRBNY-DSGE/rstarBrookings2017/blob/master/update/vintages/TVAR_Rstar_Vintages.xlsx) contains quarterly vintages of the estimates. These are the estimates shown
in Figures 1 (black lines), 4, and 5 of the paper.  The [DSGE Excel
file](https://github.com/FRBNY-DSGE/rstarBrookings2017/blob/master/update/DSGE_Rstar_Vintages.xls)
contains original as well as updated estimates of the 30-year, 10-year, and
5-year forward r* computed using the DSGE model. These are the estimates shown
in Figures 1 (blue lines), and 12 of the paper.  The figure below shows updated
estimates of r* and the liquidity/safety component.

<img src="https://raw.githubusercontent.com/FRBNY-DSGE/rstarBrookings2017/master/update/Rstar_Figure.png" width="600">


## Required software

- Julia v0.6.0 or above
- [DSGE.jl](https://github.com/FRBNY-DSGE/DSGE.jl) v0.4.1
- MATLAB 16a

**Download instructions**

1. Download Julia from `https://julialang.org/downloads/`.
2. Open the Julia REPL and type:

   a. `Pkg.add("DSGE")` to install DSGE.jl

   b. `Pkg.pin("DSGE", v"0.4.1")` to use DSGE.jl v0.4.1

   c. If, after running this replication code, you would like to use the most
      current version of DSGE.jl, type `Pkg.free("DSGE")` to un-pin the version


## Installing this repository

Git users are welcome to fork this repository or clone it for local
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

Disclaimer
------
Copyright Federal Reserve Bank of New York. You may reproduce, use, modify, make derivative works of, and distribute and this code in whole or in part so long as you keep this notice in the documentation associated with any distributed works. Neither the name of the Federal Reserve Bank of New York (FRBNY) nor the names of any of the authors may be used to endorse or promote works derived from this code without prior written permission. Portions of the code attributed to third parties are subject to applicable third party licenses and rights. By your use of this code you accept this license and any applicable third party license.

THIS CODE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT ANY WARRANTIES OR CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY WARRANTIES OR CONDITIONS OF TITLE, NON-INFRINGEMENT, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, EXCEPT TO THE EXTENT THAT THESE DISCLAIMERS ARE HELD TO BE LEGALLY INVALID. FRBNY IS NOT, UNDER ANY CIRCUMSTANCES, LIABLE TO YOU FOR DAMAGES OF ANY KIND ARISING OUT OF OR IN CONNECTION WITH USE OF OR INABILITY TO USE THE CODE, INCLUDING, BUT NOT LIMITED TO DIRECT, INDIRECT, INCIDENTAL, CONSEQUENTIAL, PUNITIVE, SPECIAL OR EXEMPLARY DAMAGES, WHETHER BASED ON BREACH OF CONTRACT, BREACH OF WARRANTY, TORT OR OTHER LEGAL OR EQUITABLE THEORY, EVEN IF FRBNY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES OR LOSS AND REGARDLESS OF WHETHER SUCH DAMAGES OR LOSS IS FORESEEABLE.
