include("getEms.jl")

function MakeBox(params, tspan, flag)

    DaysToSec = 60*60*24;
    
    # Define index for tropospheres
    CH4_NH, CH4_SH, CO_NH, CO_SH, OH_NH, OH_SH, MCF_NH, MCF_SH, N2O_NH, N2O_SH = [i for i in 1:10];

    # Stratospheric boxes
    CH4_NHS, CH4_SHS, CO_NHS, CO_SHS, OH_NHS, OH_SHS, MCF_NHS, MCF_SHS, N2O_NHS, N2O_SHS = [i for i in (N2O_SH +1): (N2O_SH+10)];

    # define reaction constants
    k_ch4, k_co, k_mcf = params["k_ch4"], params["k_co"], params["k_mcf"]
    kx_nh, kx_sh = 0.595*DaysToSec, 0.827*DaysToSec;

    # stratospheric constants
    k_ch4_strat_nh, k_ch4_strat_sh = params["k_ch4_strat_nh"], params["k_ch4_strat_sh"];
    k_co_strat, k_mcf_strat = params["k_co_strat"], params["k_mcf_strat"];
    k_n2o_strat_nh, k_n2o_strat_sh = params["k_n2o_strat_nh"], params["k_n2o_strat_sh"];
    k_oh_strat = params["k_oh_strat"];
    if flag["use_strat"]
#        kx_nh = 0.62 *DaysToSec;
        #        kx_sh = 0.845 * DaysToSec;
        kx_nh, kx_sh = 0.595*DaysToSec, 0.827*DaysToSec;
        end # if statement with kx construction in stratosphere 
    

    τᵢ = params["τᵢ"];
    τᵢ_strat = params["τᵢ_strat"]; # North-South exchange rate in the stratosphere 
    
    function interpSources(s, t)
        """
function to interpolate sources in time for differential equation slutions
inputs:
s: array of parameters being optimized
t: time as float64
output:
s: the interpolated grid of parameters
"""

        

        num_time, num_species = size(s);

        species_grid = [i for i =1:num_species];
        grid = (tspan, species_grid);
        intp = interpolate(grid, s, Gridded(Linear()));
        s_out = intp(t, species_grid);
        return s_out
    end #function interpSources


        

    function dcdt(dc, c, p, t)
        s = interpSources(p, t);
        τₛₜ = s[11];
        # Conserve mass between strat & trop, so:
        # tau_ST = tau_TS/((p_surface-p_tropopause)/(p_tropopause-p_stratopause))
        # where m_fac = (p_surface-p_tropopause)/(p_tropopause-p_stratopause)
        #             = 5.7047
        τₜₛ = τₛₜ / 5.7047;




        dc[CH4_NH] = s[CH4_NH] - k_ch4 * c[CH4_NH] * c[OH_NH] + (c[CH4_SH] - c[CH4_NH])/τᵢ;
        dc[CH4_SH] = s[CH4_SH] - k_ch4*c[CH4_SH]*c[OH_SH] + (c[CH4_NH] - c[CH4_SH])/τᵢ;

        # CO reactions
        dc[CO_NH] = s[CO_NH] + k_ch4*c[CH4_NH]*c[OH_NH] - k_co*c[CO_NH]*c[OH_NH] + (c[CO_SH] - c[CO_NH])/τᵢ;
        dc[CO_SH] = s[CO_SH] + k_ch4*c[CH4_SH]*c[OH_SH] - k_co*c[CO_SH]*c[OH_SH] + (c[CO_NH] - c[CO_SH])/τᵢ;

        #OH reactions
        dc[OH_NH] = s[OH_NH] - k_ch4*c[CH4_NH]*c[OH_NH] - k_co*c[CO_NH]*c[OH_NH] - kx_nh*c[OH_NH];
        dc[OH_SH] = s[OH_SH] - k_ch4*c[CH4_SH]*c[OH_SH] - k_co*c[CO_SH]*c[OH_SH] - kx_sh*c[OH_SH];

        # MCF in Troposphere 
        dc[MCF_NH] = s[MCF_NH] - k_mcf * c[MCF_NH] * c[OH_NH] + (c[MCF_SH] - c[MCF_NH]) / τᵢ;
        dc[MCF_SH] = s[MCF_SH] - k_mcf * c[MCF_SH] * c[OH_SH] + (c[MCF_NH] - c[MCF_SH]) / τᵢ;

        dc[N2O_NH] = s[N2O_NH] + (c[N2O_SH] - c[N2O_NH]) / τᵢ;
        dc[N2O_SH] = s[N2O_SH] + (c[N2O_NH] - c[N2O_SH]) / τᵢ;

        

        ### Stratospheric boxes
        if flag["use_strat"]
            dc[CH4_NH] += (c[CH4_NHS] - c[CH4_NH]) / τₛₜ;
            dc[CH4_SH] += (c[CH4_SHS] - c[CH4_SH]) / τₛₜ;
            dc[CO_NH] += (c[CO_NHS] - c[CO_NH]) / τₛₜ;
            dc[CO_SH] += (c[CO_SHS] - c[CO_SH]) / τₛₜ;
            dc[OH_NH] += (c[OH_NHS] - c[OH_NH]) / τₛₜ;
            dc[OH_SH] += (c[CO_SHS] - c[OH_SH]) / τₛₜ;
            dc[MCF_NH] += (c[MCF_NHS] - c[MCF_NH]) / τₛₜ;
            dc[MCF_SH] += (c[MCF_SHS] - c[MCF_SH]) / τₛₜ;
            dc[N2O_NH] += (c[N2O_NHS] - c[N2O_NH]) / τₛₜ;
            dc[N2O_SH] += (c[N2O_SHS] - c[N2O_SH]) / τₛₜ;

            # CH4 in Stratosphere
            dc[CH4_NHS] = -k_ch4_strat_nh * c[CH4_NHS] + (c[CH4_NH] - c[CH4_NHS]) / τₜₛ + (c[CH4_SHS] - c[CH4_NHS]) / τᵢ_strat;
            dc[CH4_SHS] = -k_ch4_strat_sh * c[CH4_SHS] + (c[CH4_SH] - c[CH4_SHS]) / τₜₛ + (c[CH4_NHS] - c[CH4_SHS]) / τᵢ_strat;

            # CO in stratosphere
            dc[CO_NHS] = -k_co_strat * c[CO_NHS] + (c[CO_NH] - c[CO_NHS]) / τₜₛ + (c[CO_SHS] - c[CO_NHS]) / τᵢ_strat;
            dc[CO_SHS] = -k_co_strat * c[CO_SHS] + (c[CO_SH] - c[CO_SHS]) / τₜₛ + (c[CO_NHS] - c[CO_SHS]) / τᵢ_strat;

            # OH in the stratopshere 
            dc[OH_NHS] = -k_oh_strat * c[OH_NHS] + (c[OH_SHS] - c[OH_NHS]) / τᵢ_strat; #+ (c[OH_NH] - c[OH_NHS]) / τₛ;
            dc[OH_SHS] = -k_oh_strat * c[OH_SHS] + (c[OH_NHS] - c[OH_SHS]) / τᵢ_strat; #+ (c[OH_SH] - c[OH_SHS]) / τₛ;

            # MCF in the Stratosphere
            dc[MCF_NHS] = -k_mcf_strat * c[MCF_NHS] + (c[MCF_NH] - c[MCF_NHS]) / τₜₛ + (c[MCF_SHS] - c[MCF_NHS]) / τᵢ_strat;
            dc[MCF_SHS] = -k_mcf_strat * c[MCF_SHS] + (c[MCF_SH] - c[MCF_SHS]) / τₜₛ + (c[MCF_NHS] - c[MCF_NHS]) / τᵢ_strat;

            # N2O in the Stratosphere
            dc[N2O_NHS] = -k_n2o_strat_nh * c[N2O_NHS] + (c[N2O_NH] - c[N2O_NHS]) / τₜₛ + (c[N2O_SHS] - c[N2O_NHS]) / τᵢ_strat;
            dc[N2O_SHS] = -k_n2o_strat_sh * c[N2O_SHS] + (c[N2O_SH] - c[N2O_SHS]) / τₜₛ + (c[N2O_NHS] - c[N2O_SHS]) / τᵢ_strat;
        end # stratospheric box construction
    end # function dcdt
end # function MakeBox

function BoxModelWrapper(IC, ems, params, tspan, flags)

    # convert ems from Tg/yr to ppb/day
    ems = convertEms(ems, params);
    ts = (tspan[1], tspan[end]);
    model = MakeBox(params, tspan, flags);




    # run the model
    problem = ODEProblem(model, IC, ts, ems);

    out = solve(problem, alg_hints = "stiff",saveat=tspan, dtmax=params["YrToDay"]/24);
    return out, problem
end

function StructToArray(struct_in, n_air)

    num_rows = length(struct_in);
    num_cols = length(struct_in[1][:]);
    array_out = Array{Float64}(undef, num_rows, num_cols);

    
    for i = 1:num_rows
        array_out[i,:] = struct_in[i][:];
    end #for loop

    # convert [OH] from ppb to molec/cm³
    OH_inds = [5, 6, 15, 16];
    array_out[: , OH_inds] = array_out[: , OH_inds] * n_air / 1e9;
    return array_out;
end #function StructToArray

