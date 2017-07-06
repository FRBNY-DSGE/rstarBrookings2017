using DSGE, ClusterManagers, HDF5

##########################################################################################
## SETUP
##########################################################################################

# What do you want to do?
run_estimation     = false
run_modal_forecast = false
run_full_forecast  = false

# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss18")

# Settings for data, paths, etc.
dataroot = joinpath(dirname(@__FILE__()), "input_data")
saveroot = dirname(@__FILE__())
m <= Setting(:dataroot, dataroot, "Input data directory path")
m <= Setting(:saveroot, saveroot, "Output data directory path")
m <= Setting(:data_vintage, "161223")
# m <= Setting(:use_population_forecast, true)
m <= Setting(:use_population_forecast, false)

# Settings for estimation
# set to false => will load pre-computed mode and hessian before MCMC
m <= Setting(:reoptimize, false)
m <= Setting(:calculate_hessian, false)

# Settings for forecast dates
m <= Setting(:date_forecast_start,  quartertodate("2016q4"))
m <= Setting(:date_conditional_end, quartertodate("2016q4"))
m <= Setting(:shockdec_startdate,   Nullable(date_mainsample_start(m)))

# Parallelization
m <= Setting(:forecast_block_size,  500)
nworkers = 30
addprocsfcn = addprocs_sge # choose to work with your scheduler; see ClusterManagers.jl

##########################################################################################
## RUN
##########################################################################################

# Run estimation
if run_estimation

    # Use calculated ss18 mode
    if !reoptimize(m)
        mode_file = joinpath(dataroot, "user", "paramsmode_vint=161223.h5")
        specify_mode!(m, mode_file)
    end

    # Use calculated ss18 hessian
    if !calculate_hessian(m)
        hessian_file = joinpath(dataroot, "user", "hessian_vint=161223.h5")
        specify_hessian(m, hessian_file)
    end

    estimate(m)

    # Print tables of estimated parameter moments
    groupings = DSGE.parameter_groupings(m)
    moment_tables(m, groupings = groupings)
end

# Forecast step: produces smoothed histories and shock decompositions
if run_modal_forecast || run_full_forecast

    # what do we want to produce?
    output_vars = [:histpseudo, :forecastpseudo, :shockdecpseudo]

    # conditional type
    cond_type = :none

    # Forecast label: all forecast output filenames will contain this string
    forecast_string = ""

    # Modal forecast
    if run_modal_forecast
        # run modal forecasts and save all draws
        forecast_one(m, :mode, cond_type, output_vars; verbose = :high)

        # compute means and bands
        means_bands_all(m, :mode, cond_type, output_vars)
    end

    # Full-distribution forecast
    if run_full_forecast
        my_procs = addprocsfcn(nworkers)
        @everywhere using DSGE

        forecast_one(m, :full, cond_type, output_vars; verbose = :high, forecast_string = forecast_string)
        rstar_bands = [0.68, 0.95]
        means_bands_all(m, :full, cond_type, output_vars; verbose = :high, density_bands = rstar_bands,
                        forecast_string = forecast_string)
        rmprocs(my_procs)

        meansbands_matrix_all(m, :full, cond_type, output_vars; forecast_string = forecast_string)

        # print history means and bands tables to csv
        write_meansbands_tables_all(m, :full, cond_type, [:histpseudo], forecast_string = forecast_string)

        # print shockdec means and bands tables to csv
        if any(x->contains(string(x), "shockdec"), output_vars)
            shockdec_vars = [:ExpectedAvg20YearRealNaturalRate, :ExpectedAvg20YearRealRate, :OutputGap,
                             :Forward5YearRealNaturalRate, :Forward10YearRealNaturalRate,
                             :Forward20YearRealNaturalRate, :Forward30YearRealNaturalRate,
                             :RealNaturalRate, :ExAnteRealRate]

            write_meansbands_tables_all(m, :full, cond_type, [:shockdecpseudo, :trendpseudo, :dettrendpseudo],
                                        vars = shockdec_vars,
                                        forecast_string = forecast_string)

        end
    end
end

nothing
