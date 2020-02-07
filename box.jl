
using ParameterizedFunctions, Printf
include("getEms.jl")

function makeBox(params, tspan)
 
    # 1. c(1) = OH
    # 2. C(3) = CH4
    # 3. c(5) = CO
    #c(7) = MCF
    #c(9) = N2O
    #C(11) = SF6

    # defining indecies 
    CH4_NH, CH4_SH, CO_NH, CO_SH, OH_NH, OH_SH, MCF_NH, MCF_SH, N2O_NH, N2O_SH = [i for i in 1:10];
    #    CH4_NHS, CH4_SHS, CO_NHS, CO_SHS, OH_NHS, OH_SHS, MCF_NHS, MCF_SHS, N2O_NHS, N2O_SHS = [i for i = 11:20];
    CH4_s, CO_s, OH_S, MCF_S, N2O_s = [i for i = 11:15];

# CH4 + OH
R1(c,CH4, OH) = params["k_ch4"] * c[CH4] * c[OH];

    # CO + Oh
    R2(c,CO, OH) = params["k_co"] * c[CO]*c[OH];


# MCF + OH 
    R3(c,MCF, OH) = params["k_mcf"] * c[MCF] * c[OH];

    tau_i = 1*365.25;
    tau_s = 3*365.25;
    kX_NH = 0.59* 60*60*24;
    kX_SH = 0.83*60*60*24;

    inter_trop(c, con_current, con_other, tau_i) = (c[con_other] - c[con_current])/tau_i;
    strat_trop(c, con_trop, con_strat, tau_s) = (c[con_strat] - c[con_trop])/tau_s;


    # indivitual kinetic equations:
    dCH4_NH(s,c, t) = s[CH4_NH] - R1(c, CH4_NH, OH_NH) + inter_trop(c, CH4_NH, CH4_SH, tau_i);
    dCH4_SH(s,c, t) = s[CH4_SH] - R1(c, CH4_SH, OH_SH) + inter_trop(c, CH4_SH, CH4_NH, tau_i);

    dOH_NH(s,c,t) = s[OH_NH] - R1(c, CH4_NH, OH_NH) - R2(c, CO_NH, OH_NH) - kX_NH*c[OH_NH];
    dOH_SH(s,c,t) = s[OH_SH] - R1(c, CH4_SH, OH_SH) - R2(c, CO_SH, OH_SH) - kX_SH*c[OH_SH];


    dCO_NH(s,c, t) = s[CO_NH] + R1(c, CH4_NH, OH_NH) - R2(c, CO_NH, OH_NH) + inter_trop(c, CO_NH, CO_SH, tau_i);
    dCO_SH(s,c, t) = s[CO_SH] + R1(c, CH4_SH, OH_SH) - R2(c, CO_SH, OH_SH) + inter_trop(c, CO_SH, CO_NH, tau_i);
    
    dMCF_NH(s,c, t) = s[MCF_NH] - R3(c, MCF_NH, OH_NH) + inter_trop(c, MCF_NH, MCF_SH, tau_i);
    dMCF_SH(s,c, t) = s[MCF_SH] - R3(c, MCF_SH, OH_SH) + inter_trop(c, MCF_SH, MCF_NH, tau_i);

    dN2O_NH(s,c, t) = s[N2O_NH] + inter_trop(c, N2O_NH, N2O_SH, tau_i);
    dN2O_SH(s,c, t) = s[N2O_SH] + inter_trop(c, N2O_SH, N2O_NH, tau_i);
    
    # add stratospheric box




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
#        grid = (species_grid, tspan);
        intp = interpolate(grid, s, Gridded(Linear()));
        s_out = intp(t, species_grid);
#        @printf("CH4 = %f CO = %f OH = %f \n", s_out[2], s_out[4], s_out[6]);
        return s_out
    end #function interpSources
        
        
        

    function dc_dt(dc, c, s, t)



        s = interpSources(s, t); # get the gridded source terms at time t
        dc[1:10] = [dCH4_NH(s,c,t), dCH4_SH(s,c,t), dCO_NH(s,c,t), dCO_SH(s,c,t), dOH_NH(s,c,t), dOH_SH(s,c,t), dMCF_NH(s,c,t), dMCF_SH(s,c,t), dN2O_NH(s,c,t), dN2O_SH(s,c,t)];

        end

        return dc_dt

end


params = getParams();
    k_ch4, k_co, k_mcf = params["k_ch4"], params["k_co"], params["k_mcf"];



    box = @ode_def begin
        dnh_ch4 = s_nh_ch4 - k_ch4 * nh_ch4 * nh_oh;
        dsh_ch4 = s_sh_ch4 - k_ch4 * sh_ch4 * sh_oh;
        dnh_co = s_nh_co + k_ch4 * nh_ch4 * nh_oh - k_co * nh_co * nh_oh;
        dsh_co = s_sh_co + k_ch4 * sh_ch4 * sh_oh - k_co * sh_co * sh_oh;
        dnh_oh = s_nh_oh - k_ch4 * nh_ch4 * nh_oh - k_co * nh_co * nh_oh - kx_nh * nh_oh;
        dsh_oh = s_sh_oh - k_ch4 * sh_ch4 * sh_oh - k_co * sh_co * sh_oh - kx_sh * sh_oh;
    end s_nh_ch4 s_nh_co s_nh_oh kx_nh




function BoxModelWrapper(IC, ems, params, tspan)

    # convert ems from Tg/yr to ppb/day
    ems = convertEms(ems, params);
    ts = (tspan[1], tspan[end]);
    model = makeBox(params, tspan);




    # run the model
    problem = ODEProblem(model, IC, ts, ems);

    out = solve(problem, alg_hints = "stiff",saveat=tspan, dtmax=params["YrToDay"]/100);
    return out, problem
end
 
function StructToArray(struct_in)

    num_rows = length(struct_in);
    num_cols = length(struct_in[1][:]);
    array_out = Array{Float64}(undef, num_rows, num_cols);
    for i = 1:num_rows
        array_out[i,:] = struct_in[i][:];
    end #for loop
    return array_out;
end #function StructToArray

