using DSGE, ClusterManagers, HDF5, OrderedCollections, DataFrames, ModelConstructors
using Plots, Nullables, Dates
gr() # Or specify whichever plotting backend you prefer

##########################################################################################
## SETUP
##########################################################################################

# What do you want to do?
run_estimation     = false
run_forecast       = true
make_tables        = false
plot_irfs          = false
plot_shockdecs     = true

# Number of workers for parallel?
n_workers = 20

# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss20")

# Settings for data, paths, etc.
dataroot = joinpath(dirname(@__FILE__()), "input_data")
saveroot = dirname(@__FILE__())
m <= Setting(:dataroot, dataroot, "Input data directory path")
m <= Setting(:saveroot, saveroot, "Output data directory path")
m <= Setting(:data_vintage, "190829")
m <= Setting(:use_population_forecast, false)

# Settings for estimation
# set to false => will load pre-computed mode and hessian before MCMC
m <= Setting(:reoptimize, true)
m <= Setting(:calculate_hessian, true)

# Settings for forecast dates
m <= Setting(:date_forecast_start,  quartertodate("2019-Q3"))
m <= Setting(:date_conditional_end, quartertodate("2019-Q3"))
m <= Setting(:shockdec_startdate,   Nullable(date_mainsample_start(m)))

# Parallelization
m <= Setting(:forecast_block_size,  500)
#nworkers = 30
#addprocsfcn = addprocs_sge # choose to work with your scheduler; see ClusterManagers.jl

##########################################################################################
## RUN
##########################################################################################

# Run estimation
if run_estimation

    if reoptimize(m)
        # Start from ss18 mode
        mode_file = rawpath(m, "estimate", "paramsmode.h5")
        mode_file = replace(mode_file, "ss20", "ss18")
        update!(m, h5read(mode_file, "params"))
    else
        # Use calculated ss18 mode
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
if run_forecast

    # set overrides
    overrides = forecast_input_file_overrides(m)
    overrides[:mode] = joinpath(saveroot, "output_data/m1010/ss20/estimate/raw/paramsmode_vint=161223.h5")
    overrides[:full] = joinpath(saveroot, "output_data/m1010/ss20/estimate/raw/mhsave_vint=161223.h5")

    # what do we want to produce?
    output_vars = [:histobs, :histpseudo, :forecastobs, :forecastpseudo] #=, :shockdecobs, :shockdecpseudo,
                   :irfobs, :irfpseudo] =#
    push!(output_vars, [:dettrendobs, :dettrendpseudo, :trendobs, :trendpseudo,
                       :shockdecpseudo, :shockdecobs]...)

    # conditional type
    cond_type = :none

    # Forecast label: all forecast output filenames will contain this string
    forecast_string = ""

    # Modal forecast
    # run modal forecasts and save all draws
    forecast_one(m, :mode, cond_type, output_vars; verbose = :high)

    # compute means and bands
    compute_meansbands(m, :mode, cond_type, output_vars)

    # Full-distribution forecast
    my_procs = addprocs_frbny(n_workers)
    @everywhere using OrderedCollections
    @everywhere using DSGE
    @everywhere using DataFrames

    forecast_one(m, :full, cond_type, output_vars; verbose = :high, forecast_string = forecast_string)
    rstar_bands = [0.68, 0.95] # Do not change these rstar_bands, since the makeRstarPlots.m file assumes
                               # the band ordering to be hard-coded
    compute_meansbands(m, :full, cond_type, output_vars; verbose = :high, density_bands = rstar_bands,
                       forecast_string = forecast_string)
    rmprocs(my_procs)

    meansbands_to_matrix(m, :full, cond_type, output_vars; forecast_string = forecast_string)
end

if make_tables

    # input and conditional types
    input_types = [:mode, :full]
    cond_type = :none

    # Forecast label: all forecast output filenames will contain this string
    forecast_string = ""

    vars = [:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
            :RealNaturalRate, :Forward5YearRealNaturalRate,
            :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
            :Forward30YearRealNaturalRate, :ExpectedAvg20YearRealNaturalRate]

    for input_type in input_types
        # print history means and bands tables to csv
        write_meansbands_tables_all(m, input_type, cond_type, [:histpseudo, :forecastpseudo], forecast_string = forecast_string,
                                    vars = vars)

        # print shockdec means and bands tables to csv
       #= write_meansbands_tables_all(m, input_type, cond_type, [:shockdecpseudo, :trendpseudo, :dettrendpseudo],
                                    vars = vars, forecast_string = forecast_string)=#
    end

end

if plot_irfs

    #############################################
    # Change these settings as you see fit
    #############################################
    # It is recommended to only plot irfs for ≦ 6 variables so that the subplot aligns to
    # a single page properly
    vars =  [:RealNaturalRate, :Forward5YearRealNaturalRate,:Forward10YearRealNaturalRate,
             :Forward30YearRealNaturalRate, :y_t, :OutputGap]

    shocks = [:b_liqp_sh, :b_liqtil_sh, :b_safep_sh, :b_safetil_sh,
              :z_sh, :zp_sh]

    input_type = :mode # :mode or :full depending on the forecast that you ran

    all_plots = Vector{Plots.Plot}(length(shocks))
    bands = input_type == :full ? ["68.0%", "95.0%"] : Vector{String}()
    save_plots = true # Whether or not you want to save your plots
    plotroot = figurespath(m, "forecast")
    #############################################

    # These are the shocks for whom we want to flip the impulse responses
    # (i.e. give the response to a +1 SD shock instead of the -1 SD shock
    # calculated by `impulse_responses`
    flipped_shocks = [:σ_ω_sh, :z_sh, :λ_f_sh]

    for (i, shock) in enumerate(shocks)
        # Create individual plots
        plots = plot_impulse_response(m, shock, vars, :pseudo, input_type, :none;
                                      flip_sign = shock in flipped_shocks,
                                      bands_pcts = bands,
                                      plotroot = "")

        # Plot all variables together
        nvars = length(vars)
        nrows = convert(Int, ceil(nvars/2))
        varplots = if isodd(nvars)
            empty = plot(grid = false, foreground_color_axis = :white)
            vcat(collect(values(plots)), [empty])
        else
            collect(values(plots))
        end
        p = plot(varplots..., layout = grid(nrows, 2), size = (650, 775), titlefontsize = 10)
        all_plots[i] = p

        if save_plots
            # Save
            fn = joinpath(plotroot, "IRF_" * DSGE.detexify(string(shock)) * "_vint=" * data_vintage(m) * ".pdf")
            DSGE.save_plot(p, fn)
        end
    end

end

if plot_shockdecs

    #############################################
    # Change these settings as you see fit
    #############################################
    start_date = DSGE.quarter_number_to_date(1994)
    end_date = DSGE.quarter_number_to_date(2030)

    input_type = :mode

    # vars = [:RealNaturalRate, :Forward5YearRealNaturalRate,
    #         :Forward10YearRealNaturalRate, :Forward30YearRealNaturalRate]
    vars = [:Forward30YearRealNaturalRate]
    shockdec_class = :pseudo
    save_plots = true # Whether or not you want to save your plots
    #############################################

    # function paper_shock_groupings(m::Model1010)
    #     fin      = ShockGroup("FF", [:γ_sh, :μ_e_sh, :σ_ω_sh, :b_liqtil_sh, :b_liqp_sh, :b_safetil_sh, :b_safep_sh], RGB(0.29, 0.0, 0.51)) # indigo
    #     # zpe      = ShockGroup("zp", [:zp_sh], RGB(0.0, 0.3, 0.0))
    #     # tfp      = ShockGroup("z", [:z_sh], RGB(1.0, 0.55, 0.0)) # darkorange
    #     prod      = ShockGroup("Productivity", [:z_sh, :zp_sh], RGB(1.0, 0.55, 0.0)) # darkorange
    #     return [fin, prod]
    # end

    function paper_shock_groupings(m::Model1010)
        conv      = ShockGroup("Conv Yield", [:b_liqtil_sh, :b_liqp_sh, :b_safetil_sh, :b_safep_sh], RGB(0.29, 0.0, 0.51)) # indigo
        risk = ShockGroup("Risk", [:σ_ω_sh], RGB(.161, 0.486, 0.275))
        otherff = ShockGroup("Other FF", [:γ_sh, :μ_e_sh], RGB(0.0, 0.3, 0.0))
        # zpe      = ShockGroup("zp", [:zp_sh], RGB(0.0, 0.3, 0.0))
        # tfp      = ShockGroup("z", [:z_sh], RGB(1.0, 0.55, 0.0)) # darkorange
        prod      = ShockGroup("Productivity", [:z_sh, :zp_sh], RGB(1.0, 0.55, 0.0)) # darkorange
        return [conv, risk, otherff, prod]
    end


    # By default set to empty string if save_plots is false
    plotroot = save_plots ? figurespath(m, "forecast") : ""

    plots = []
    for var in vars
        p = plot_shock_decomposition(m, var, shockdec_class, input_type, :none,
                                     start_date = start_date,
                                     end_date = end_date,
                                     legend = :left,
                                     groups = paper_shock_groupings(m),
                                     forecast_label = "",
                                     plotroot = plotroot)
        push!(plots, p)
    end

    # By default set to empty string if save_plots is false
    plotroot = save_plots ? figurespath(m, "forecast") : ""

    plots = []
    for var in vars
        p = plot_shock_decomposition(m, var, shockdec_class, :full, :none,
                                     start_date = start_date,
                                     end_date = end_date,
                                     legend = :left,
                                     groups = paper_shock_groupings(m),
                                     forecast_label = "",
                                     plotroot = plotroot)
        push!(plots, p)
    end

end

nothing
