function getEms(params, tspan)

    t_length = length(tspan);
    ems = Array{Float64}(undef, t_length, 10);
    x = ones(t_length,1);

    # Units of Tg/yr
    ems[:,1] = 412.5*x;
    ems[:,2] = 137.5*x;

    # CO
    ems[:,3] = 900*x;
    ems[:,4] = 65*x;

        # OH
        ems[:,5] = 3150*x;
        ems[:,6] = 3150*x;

        # MCF
        
        ems[:,7] = 0.2*x;
    ems[:,8] = 0.2*x;

    # N@O
   
    ems[:,9] = 9.3*x;
    ems[:,10] = 3.5*x;
    return ems
end

function convertEms(ems, params, tspan)

    # define function to convert from Tg/yr to ppb/day or ppt/day
    function Tg2ppb(ems_in, conversion_factor)
        return 2/conversion_factor * ems_in / params["YrToDay"];
    end


    conv_factors = [params["mm_ch4"], params["mm_ch4"], params["mm_co"], params["mm_co"], params["mm_oh"], params["mm_oh"], params["mm_mcf"], params["mm_mcf"], params["mm_n2o"], params["mm_n2o"]];
    ems_out = Array{Float64}(undef, size(ems));

    # call Tg2ppb for all species
    for i = 1:length(conv_factors)
        ems_out[:,i] = Tg2ppb(ems[:,i], conv_factors[i]);
    end
    return ems_out
end

    
                    

