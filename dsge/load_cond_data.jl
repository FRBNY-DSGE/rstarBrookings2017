# This script loads in conditional daily data used for spreads in the Brookings DSGE model (Model 1010 ss20).
# This script assumes that an incomplete conditional data set already exists (see line 19).

using FredData, Dates, DataFrames, CSV, Statistics
using DSGE: lastdayofquarter, missing2nan

# Change these settings to update data file
cond_id  = 02


<<<<<<< Updated upstream
vintage      = "2023-07-13"
#cond_vintage = "2023-07-13"                  # YYYYMMDD format, should be the date w/latest available vintage
quarter  = "2023-04-01"                      # start date of current quarter, so 2020-04-01 == 2020-Q2
=======
vintage      = "2023-10-23"
cond_vintage = "2023-10-23"                  # YYYYMMDD format, should be the date w/latest available vintage
quarter  = "2023-10-01"                      # start date of current quarter, so 2020-04-01 == 2020-Q2
>>>>>>> Stashed changes

use_mean = true                              # Use the mean, otherwise use the most recent observation

# Other settings, usually do not need to change
fred_names = ["DBAA", "DAAA", "DGS10", "DGS20", "DGS30"]
save_names = [:BAA,   :AAA,   :GS10,   :GS20,   :GS30]   # usually the name for the series at quarterly frequency
units      = ["lin",  "lin",  "lin",   "lin",   "lin"]
file_loc   = "input_data/cond"
cond_idno = lpad(string(cond_id), 2, string(0)) # print as 2 digits
filename   = "cond_cdid=" * cond_idno * "_cdvt=$(vintage[vcat(3:4, 6:7, 9:10)]).csv"
quarter_over = Date(vintage) > Date(quarter) + Month(3) # Has the quarter ended?

# Fetch the data!
f = Fred()
data = DataFrame()
for (name, unit) in zip(fred_names, units)
    fred_series = get_data(f, name; vintage_dates = vintage, frequency = "q", units = unit)
    if isempty(data)
        global data[!, :date] = fred_series.data[!, :date]
        global data[!, Symbol(name)] = fred_series.data[!, :value]
    else
        global data = outerjoin(data, fred_series.data[!, [:date, :value]], on = :date)
        colnames = propertynames(data)
        val_col  = findfirst(colnames .== :value)
        colnames[val_col] = Symbol(name)
        rename!(data, colnames)
    end
end
sort!(data, [:date])

# When current quarter has not ended yet
daily_data = DataFrame()
if !quarter_over
    for (name, unit) in zip(fred_names, units)
        fred_series = get_data(f, name; observation_start = quarter,
                               vintage_dates = vintage, frequency = "d", units = unit)
        if isempty(daily_data)
            global daily_data[!, :date] = fred_series.data[!, :date]
            global daily_data[!, Symbol(name)] = fred_series.data[!, :value]
        else
            global daily_data = outerjoin(daily_data, fred_series.data[!, [:date, :value]], on = :date)
            colnames = propertynames(daily_data)
            val_col  = findfirst(colnames .== :value)
            colnames[val_col] = Symbol(name)
            rename!(daily_data, colnames)
        end
    end

    # Clean data
    sort!(daily_data, [:date])
    entry_row = findfirst(data[!, :date] .== Date(quarter))
    for name in propertynames(data)
        if name != :date
            skip_val = .!(isnan.(daily_data[!, name]))
            skip_val = map(x -> ismissing(x) ? true : x, skip_val)
            if use_mean
                global data[entry_row, name] = Float64(mean(skipmissing(daily_data[skip_val, name])))
            else
                global data[entry_row, name] = Float64(skipmissing(daily_data[skip_val, name])[end])
            end
        end
    end
end

# Fix the dates in data
date_vec = map(x -> lastdayofquarter(x), data[!, :date])
data[!, :date] = date_vec

# Read in the conditional data and update it
cond_data = CSV.read(joinpath(file_loc, filename), DataFrame)
# for name in save_names         # only uncomment these lines if you want to check that this code is correct.
#     cond_data[!, name] .= NaN  # Make sure you have another copy of the original CSV file to compare against.
# end
tmp = innerjoin(cond_data, data, on = :date) # to get the right dates
for (name, save_name) in zip(propertynames(data)[.!(propertynames(data) .== :date)], save_names)
    if save_name in propertynames(cond_data)
        miss_or_nan = isnan.(missing2nan(convert(Vector{Union{Missing, Float64}}, cond_data[!, save_name])))
        cond_data[miss_or_nan, save_name] = tmp[miss_or_nan, name]
    else
        cond_data[!, save_name] .= tmp[!, name]
    end

end
CSV.write(joinpath(file_loc, filename), cond_data)
