using Dates

function makeTime(sYear, eYear, time_res)

    if time_res == "year"
        time_vector = Dates.DateTime(sYear):Dates.Year(1):Dates.DateTime(eYear); # space through yearly timestep
    elseif time_res == "month"
        time_vector = Dates.DateTime(sYear,1,1):Dates.Month(1):Dates.DateTime(eYear,1,1); # create hte monthlyly time vector
    end

    # convert to days after unix time
    time_vector = collect(time_vector);
    time_vector = [Dates.datetime2unix(i) for i in time_vector]; # convert to unix time in sec
    time_vector = time_vector/(24*60*60); # days
    return time_vector
end




