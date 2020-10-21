using DSGE #, ModelConstructors
using ClusterManagers, HDF5, Plots, ModelConstructors, Nullables
GR.inline("pdf")
include("/data/dsge_data_dir/dsgejl/shlok/proc/includeall.jl")

@everywhere using OrderedCollections

# What do you want to do?
run_meansbands_combo    = false
make_standard_plots     = false
make_standard_packet    = false
make_rstar_tables       = true
make_snapshot_plots     = false
make_snapshot_tables    = false
auto_add_procs          = true # set to false if you add workers before launching Julia-repl
n_workers               = 80

# Data vintages
data_vint  = "200716"
cond_vint  = "200716"
fcast_date = DSGE.quartertodate("2020-Q2")

# Initialize model objects
forecast_string1 = "snapshot_scenario=1"
forecast_string2 = "snapshot_scenario=2"
forecast_string3 = "snapshot_scenario=3"

m10 = Model1002("ss10")
usual_settings!(m10, data_vint, cdvt = cond_vint, fcast_date = fcast_date)
overrides = forecast_input_file_overrides(m10)
overrides[:mode] = "/data/dsge_data_dir/dsgejl/output_data/m1002/ss10/estimate/raw/paramsmode_vint=190829.h5"
if make_snapshot_plots
    m1 = temporary_shutdown_scenario(forecast_string1, m10; set_parameters = false,
                                     fcast_date = fcast_date, cond_vintage = cond_vint, data_vintage = data_vint)
    m2 = shutdown_with_business_cycle_dynamics_scenario(forecast_string2, m10; set_parameters = false,
                                                        fcast_date = fcast_date, cond_vintage = cond_vint, data_vintage = data_vint)
    m3 = persistent_demand_shortfall_scenario(forecast_string3, m10; set_parameters = false,
                                              fcast_date = fcast_date, cond_vintage = cond_vint, data_vintage = data_vint)
else
    m1 = temporary_shutdown_scenario(forecast_string1, m10,
                                     fcast_date = fcast_date, cond_vintage = cond_vint, data_vintage = data_vint)
    m2 = shutdown_with_business_cycle_dynamics_scenario(forecast_string2, m10,
                                                        fcast_date = fcast_date, cond_vintage = cond_vint, data_vintage = data_vint)
    m3 = persistent_demand_shortfall_scenario(forecast_string3, m10,
                                              fcast_date = fcast_date, cond_vintage = cond_vint, data_vintage = data_vint)
end
weights = [0.65; 0.1; 0.25]
combo_forecast_string = "snapshot_weights=" * join(map(x->string(x), weights), "_")

for model in [m1, m2, m3]
    model <= Setting(:saveroot, "/data/dsge_data_dir/dsgejl/")
    model <= Setting(:date_conditional_end, fcast_date) # Set conditional end to first forecast start date
    model <= Setting(:cond_full_names, [:obs_gdp, :obs_corepce, :obs_spread, # Have to add anticipated nominal rates to conditional data
                                        :obs_nominalrate, :obs_longrate,
                                        :obs_nominalrate1, :obs_nominalrate2, :obs_nominalrate3,
                                        :obs_nominalrate4, :obs_nominalrate5, :obs_nominalrate6])
    model <= Setting(:cond_semi_names, [:obs_spread,
                                        :obs_nominalrate, :obs_longrate,
                                        :obs_nominalrate1, :obs_nominalrate2, :obs_nominalrate3,
                                        :obs_nominalrate4, :obs_nominalrate5, :obs_nominalrate6])
end

# This should always be copied from the last snapshot (so, above!)
m_old = Model1002("ss60")
usual_settings!(m_old, "200514", cdvt = "200514", fcast_date = quartertodate("2020-Q2"))
m_old <= Setting(:n_regimes, 4, true, "reg", "")
m_old <= Setting(:snapshot, true, true, "snapshot", "")
m_old <= Setting(:drawparam, true, true, "drawparam", "")
m_old <= Setting(:drawprior, true, true, "drawprior", "")
overrides_old = forecast_input_file_overrides(m_old)
overrides_old[:mode] = "/data/dsge_data_dir/dsgejl/output_data/m1002/ss10/estimate/raw/paramsmode_vint=190829.h5"
overrides_old[:full] = "/data/dsge_data_dir/dsgejl/output_data/m1002/ss10/estimate/raw/mhsave_vint=190829.h5"
weights_prev = [0.6, 0.25, 0.15]
forecast_string_prev = "mayscen22=125_weights=" * join(map(x -> string(x), weights_prev), "_")

# Forecast
if run_meansbands_combo
    obs_output_vars = [:histobs, :forecastobs,
                       :hist4qobs, :forecast4qobs,
                       :bddforecast4qobs, :bddforecastobs]
    pseudo_output_vars = [:histpseudo, :forecastpseudo,
                          :hist4qpseudo, :forecast4qpseudo,
                          :histutpseudo, :forecastutpseudo]

    # Run modal meansbands first
    df = load_data(m1, cond_type = :full; try_disk = false)

    # Full-distribution MeansBands
    if auto_add_procs
        my_procs = addprocs_frbny(n_workers)
        ENV["juliafrbnymemory"] = "2G"
        @everywhere using DSGE, OrderedCollections
    end

    for cond_type in [:semi, :full]
        for wts in [[.65, .1, .25], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
            combo_forecast_string = "snapshot_weights=" * join(map(x->string(x), wts), "_")
            compute_meansbands([m1, m2, m3],
                               [:full, :full, :full],
                               [cond_type, cond_type, cond_type], obs_output_vars,
                               forecast_strings = [forecast_string1, forecast_string2, forecast_string3],
                               combo_forecast_string = combo_forecast_string,
                               weights = wts, variable_names = collect(keys(m3.observables)),
                               density_bands = [.5, .6, .68, .7, .8, .9], df = df,
                               verbose = :high)
            compute_meansbands([m1, m2, m3],
                               [:full, :full, :full],
                               [cond_type, cond_type, cond_type], pseudo_output_vars,
                               forecast_strings = [forecast_string1, forecast_string2, forecast_string3],
                               combo_forecast_string = combo_forecast_string,
                               weights = wts, density_bands = [.5, .6, .68, .7, .8, .9], df = df,
                               verbose = :high)
        end
    end

    combo_forecast_string = "snapshot_weights=" * join(map(x->string(x), weights), "_")

    rmprocs(my_procs)
end


if make_standard_plots
    gr()
    output_vars = [:forecastobs, :forecastpseudo]
    for var in output_vars
        make_forecast_plots(m1, :full, :full, var;
                            forecast_string = combo_forecast_string)
    end
end

if make_standard_packet
    write_standard_packet(m1, :full, :full; forecast_string = combo_forecast_string,
                          weights = weights,
                          outdir = "/data/dsge_data_dir/SystemwideDSGE/Current/Snapshot/SupportingFiles",
                          purpose = "Pre-Blackbook Meeting")
end

if make_snapshot_plots
    plotly()
    combo_forecast_string = "snapshot_weights=" * join(map(x->string(x), weights), "_")
    snapshot_plots(m1, plotroot = "Figures",
                   forecast_string = combo_forecast_string,
                   weights = weights,
                   cond_type = :full,
                   input_type = :full, bdd_and_unbdd = true,
                   zero_shocks = false)

    for wts in [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
        combo_forecast_string = "snapshot_weights=" * join(map(x->string(x), wts), "_")
        snapshot_plots(m1, plotroot = "Figures",
                       forecast_string = combo_forecast_string,
                       weights = wts,
                       cond_type = :full,
                       input_type = :full, bdd_and_unbdd = true,
                       zero_shocks = false)
    end
end

# Snapshot tables
if make_snapshot_tables
    combo_forecast_string = "snapshot_weights=" * join(map(x->string(x), weights), "_")
    snapshot_tables(m1, m_old, "Commentary/",
                    forecast_string1 = combo_forecast_string,
                    forecast_string2 = forecast_string_prev,
                    weights = weights,
                    input_type_new = :full, input_type_old = :full,
                    cond_old = :full, cond_new = :full, bdd_and_unbdd = true, zero_shocks = false)

    for wts in [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
        combo_forecast_string = "snapshot_weights=" * join(map(x->string(x), wts), "_")
        snapshot_tables(m1, m_old, "Commentary/",
                        forecast_string1 = combo_forecast_string,
                        forecast_string2 = forecast_string_prev,
                        weights = wts,
                        input_type_new = :full, input_type_old = :full,
                        cond_old = :full, cond_new = :full, bdd_and_unbdd = true, zero_shocks = false)
    end

end

if make_rstar_tables
    input_types = [:full]
    cond_type = :semi

    # Forecast label: all forecast output filenames will contain this string
    rstar_forecast_string = combo_forecast_string

    vars = [:NaturalRate]

    for input_type in input_types
        # print history means and bands tables to csv
        write_meansbands_tables_all(m1, input_type, cond_type, [:histpseudo, :forecastpseudo],
                                    forecast_string = rstar_forecast_string, vars = vars)
    end
end

nothing
