using DSGE, ClusterManagers, HDF5

##########################################################################################
## SETUP
##########################################################################################

# What do you want to do?
run_estimation     = false
run_modal_forecast = false
run_full_forecast  = false
remake_tables      = true
plot_irfs          = false
plot_shockdecs     = false

# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss20")

# Settings for data, paths, etc.
dataroot = joinpath(dirname(@__FILE__()), "input_data")
saveroot = dirname(@__FILE__())
m <= Setting(:dataroot, dataroot, "Input data directory path")
m <= Setting(:saveroot, saveroot, "Output data directory path")
m <= Setting(:data_vintage, "180509")
m <= Setting(:use_population_forecast, false)

# Settings for estimation
# set to false => will load pre-computed mode and hessian before MCMC
m <= Setting(:reoptimize, true)
m <= Setting(:calculate_hessian, true)

# Settings for forecast dates
m <= Setting(:date_forecast_start,  quartertodate("2018-Q2"))
m <= Setting(:date_conditional_end, quartertodate("2018-Q2"))
m <= Setting(:shockdec_startdate,   Nullable(date_mainsample_start(m)))

# Parallelization
m <= Setting(:forecast_block_size,  500)
nworkers = 20
addprocsfcn = addprocs_sge # choose to work with your scheduler; see ClusterManagers.jl

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
if run_modal_forecast || run_full_forecast || remake_tables

    # set overrides
    overrides = forecast_input_file_overrides(m)
    overrides[:mode] = joinpath(saveroot, "output_data/m1010/ss20/estimate/raw/paramsmode_vint=161223.h5")
    overrides[:full] = joinpath(saveroot, "output_data/m1010/ss20/estimate/raw/mhsave_vint=161223.h5")

    # what do we want to produce?
    output_vars = [:histobs, :histpseudo, :shockdecobs, :shockdecpseudo]

    # conditional type
    cond_type = :none

    # Forecast label: all forecast output filenames will contain this string
    forecast_string = ""

    # Modal forecast
    if run_modal_forecast
        # run modal forecasts and save all draws
        forecast_one(m, :mode, cond_type, output_vars; verbose = :high)

        # compute means and bands
        compute_meansbands(m, :mode, cond_type, output_vars)
    end

    # Full-distribution forecast
    if run_full_forecast
        my_procs = addprocsfcn(nworkers)
        @everywhere using DSGE

        # forecast_one(m, :full, cond_type, output_vars; verbose = :high, forecast_string = forecast_string)
        rstar_bands = [0.68, 0.90, 0.95]
        compute_meansbands(m, :full, cond_type, output_vars; verbose = :high, density_bands = rstar_bands,
                           forecast_string = forecast_string)
        rmprocs(my_procs)

        meansbands_to_matrix(m, :full, cond_type, output_vars; forecast_string = forecast_string)
    end

    # print history means and bands tables to csv
    table_vars = [:ExpectedAvg20YearRealNaturalRate, :RealNaturalRate,
                  :Forward5YearRealNaturalRate, :Forward10YearRealNaturalRate,
                  :Forward20YearRealNaturalRate, :Forward30YearRealNaturalRate]
    write_meansbands_tables_all(m, :full, cond_type, [:histpseudo], forecast_string = forecast_string,
                                vars = table_vars)

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

if remake_tables
    # print history means and bands tables to csv
    table_vars = [:Forward5YearRealNaturalRate, :Forward10YearRealNaturalRate,
                  :Forward30YearRealNaturalRate]
    write_meansbands_tables_all(m, :full, cond_type, [:histpseudo], forecast_string = forecast_string,
                                vars = table_vars)

    # print shockdec means and bands tables to csv
    write_meansbands_tables_all(m, :full, cond_type, [:shockdecpseudo, :trendpseudo, :dettrendpseudo],
                                    vars = table_vars,
                                    forecast_string = forecast_string)

end

if plot_irfs
    # ant_obs = [Symbol("obs_nominalrate" * string(i)) for i = 1:n_anticipated_shocks(m)]
    # vars = setdiff(collect(keys(m.observables)), ant_obs)
    vars =  [:RealNaturalRate, :Forward5YearRealNaturalRate,:Forward10YearRealNaturalRate,
                 :Forward30YearRealNaturalRate, :y_t, :OutputGap]
    irf_class = :pseudo

    # shocks = [:b_sh, :σ_ω_sh, :μ_sh, :z_sh, :λ_f_sh, :rm_sh]
    shocks = [:b_liqp_sh, :b_liqtil_sh, :b_safep_sh, :b_safetil_sh,
              :z_sh, :zp_sh]

    # These are the shocks for whom we want to flip the impulse responses
    # (i.e. give the response to a +1 SD shock instead of the -1 SD shock
    # calculated by `impulse_responses`
    flipped_shocks = shocks

    plotroots = "plots"
    bands = Vector{String}()

    for shock in shocks
        # Create individual plots
        plots = plot_impulse_response(m, shock, vars, :pseudo, :full, :none;
                                      flip_sign = shock in flipped_shocks,
                                      bands_pcts = bands)

        # Plot all variables together
        nvars = length(vars)
        nrows = convert(Int, ceil(nvars/2))
        varplots = if isodd(nvars)
            empty = plot(grid = false, foreground_color_axis = :white)
            vcat(collect(values(plots)), [empty])
        else
            collect(values(plots))
        end
        p = plot(varplots..., layout = grid(nrows, 2), size = (650, 775))

        # Save
        fn = joinpath(plotroot, "Appendix_IRF_" * DSGE.detexify(string(shock)) * ".pdf")
        DSGE.save_plot(p, fn)
    end

end

if plot_shockdecs
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

    start_date = DSGE.quarter_number_to_date(1994)
    end_date = DSGE.quarter_number_to_date(2017.5)

    # vars = [:RealNaturalRate, :Forward5YearRealNaturalRate,
    #         :Forward10YearRealNaturalRate, :Forward30YearRealNaturalRate]
    # shockdec_class = :pseudo
    vars = [:obs_tfp]
    shockdec_class = :obs

    plots = []
    for var in vars
        p = plot_shock_decomposition(m, var, shockdec_class, :mode, :none,
                                     #start_date = start_date,
                                     end_date = end_date,
                                     fileformat = :png, legend = :left,
                                     groups = paper_shock_groupings(m),
                                     forecast_label = "")
        push!(plots, p)
    end
end

nothing
