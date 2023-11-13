using DSGE, ClusterManagers, HDF5, JLD2, JLD, OrderedCollections, DataFrames, ModelConstructors
using Plots, Nullables
gr() # Or specify whichever plotting backend you prefer

##########################################################################################
## SETUP
##########################################################################################

# What do you want to do?
run_estimation     = false
run_forecast       = true
make_tables        = true
plot_irfs          = false
plot_shockdecs     = false

# Number of workers for parallel?
n_workers = 48
# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss20", custom_settings = [Setting(:n_mon_anticipated_shocks, 6)])

# Settings for data, paths, etc.
dataroot = joinpath(dirname(@__FILE__()), "input_data")
saveroot = dirname(@__FILE__())
m <= Setting(:dataroot, dataroot, "Input data directory path")
m <= Setting(:saveroot, saveroot, "Output data directory path")
m <= Setting(:data_vintage, "230830")
m <= Setting(:cond_vintage, "230830")
m <= Setting(:cond_id, 02)
m <= Setting(:cond_full_names, [:obs_gdp, :obs_corepce, :obs_AAAspread, :obs_BBBspread,
                                :obs_nominalrate, :obs_longrate,
                                :obs_nominalrate1, :obs_nominalrate2, :obs_nominalrate3,
                                :obs_nominalrate4, :obs_nominalrate5, :obs_nominalrate6])
m <= Setting(:cond_semi_names, [:obs_AAAspread,:obs_BBBspread,
                                :obs_nominalrate, :obs_longrate,
                                :obs_nominalrate1, :obs_nominalrate2, :obs_nominalrate3,
                                :obs_nominalrate4, :obs_nominalrate5, :obs_nominalrate6])
#m <= Setting(:cond_semi_names, [:obs_AAAspread,:obs_BBBspread])

m <= Setting(:use_population_forecast, false)

forecast_string = ""

# Settings for estimation
# set to false => will load pre-computed mode and hessian before MCMC
m <= Setting(:reoptimize, true)
m <= Setting(:calculate_hessian, true)

# Settings for forecast dates
m <= Setting(:date_forecast_start,  quartertodate("2023-Q3"))
m <= Setting(:date_conditional_end, quartertodate("2023-Q3"))
m <= Setting(:shockdec_startdate,   Nullable(date_mainsample_start(m)))

# Parallelization
m <= Setting(:forecast_block_size,  500)


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
    output_vars = [:histobs, :histpseudo, :forecastobs, :forecastpseudo]#=, :shockdecobs, :shockdecpseudo,
                  # :irfobs, :irfpseudo] =#
                  #[:shockdecobs, :shockdecpseudo] #:trendobs, :trendpseudo, :dettrendobs, :dettrendpseudo]

    # conditional type
    cond_type = :full

    df = load_data(m, cond_type = cond_type, check_empty_columns = false, try_disk = false)
   # @assert false


    # Modal forecast
    # run modal forecasts and save all draws
    forecast_one(m, :mode, cond_type, output_vars; verbose = :high, forecast_string = forecast_string, df = df)

    # compute means and bands
    compute_meansbands(m, :mode, cond_type, output_vars, forecast_string = forecast_string, df = df)

    # Full-distribution forecast
    #ENV["frbnyjuliamemory"] = "2G"
    my_procs = addprocs_frbny(n_workers)
    @everywhere using OrderedCollections
    @everywhere using DSGE
    @everywhere using DataFrames

    ENV["frbnyjuliamemory"] = "2G"
    forecast_one(m, :full, cond_type, output_vars; verbose = :high, forecast_string = forecast_string, df = df)
    rstar_bands = [0.68, 0.95] # Do not change these rstar_bands, since the makeRstarPlots.m file assumes
                               # the band ordering to be hard-coded
    compute_meansbands(m, :full, cond_type, output_vars; verbose = :high, density_bands = rstar_bands,
                       forecast_string = forecast_string, df = df)
    rmprocs(my_procs)

    meansbands_to_matrix(m, :full, cond_type, output_vars; forecast_string = forecast_string)
end

if make_tables
    cond_type = :full

    # input and conditional types
    input_types = [:mode, :full]
    if !isdefined(Main, :cond_type)
        cond_type = :full
    end

    vars = [:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
            :RealNaturalRate, :Forward5YearRealNaturalRate,
            :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
            :Forward30YearRealNaturalRate, :ExpectedAvg20YearRealNaturalRate]

    #vars = [:obs_AAAspread, :obs_BBBspread]

    for input_type in input_types
        # print history means and bands tables to csv
        write_meansbands_tables_all(m, input_type, cond_type, [:histpseudo, :forecastpseudo], forecast_string = forecast_string,
                                    vars = vars)

        # print shockdec means and bands tables to csv
       # write_meansbands_tables_all(m, input_type, cond_type, [:shockdecpseudo, :trendpseudo, :dettrendpseudo], #typically pseudo
                         #          vars = vars, forecast_string = forecast_string)
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
    end_date = DSGE.quarter_number_to_date(2023.25)

    input_type = :mode

    vars = [:RealNaturalRate, :Forward5YearRealNaturalRate,
            :Forward10YearRealNaturalRate, :Forward30YearRealNaturalRate]
    shockdec_class = :pseudo
    save_plots = true # Whether or not you want to save your plots
    #############################################

    function paper_shock_groupings(m::Model1010)
        conv  = DSGE.ShockGroup("Conv Yield", [:b_liqtil_sh, :b_liqp_sh, :b_safetil_sh, :b_safep_sh],
                          RGB(.667, 0.29, 0.224))
        risk  = DSGE.ShockGroup("Risk", [:σ_ω_sh],  RGB(0.161, 0.486, 0.275))
        perm_prod  = DSGE.ShockGroup("Zp", [:zp_sh],  RGB(0.18, 0.259, 0.447))
        temp_prod = DSGE.ShockGroup("Z", [:z_sh], :blue)
        other = DSGE.ShockGroup("Other", [:g_sh, :λ_f_sh, :λ_w_sh, :γ_sh, :μ_e_sh,
                                          :π_star_sh, :μ_sh, :lr_sh, :tfp_sh, :gdpdef_sh, :corepce_sh,
                                          :gdp_sh, :gdi_sh, :dettrend], :gray40)

        return [conv, risk, perm_prod, temp_prod, other]
    end

    # By default set to empty string if save_plots is false
    plotroot = save_plots ? "test_figures/" : ""
        #figurespath(m, "forecast") : ""

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
end

nothing
