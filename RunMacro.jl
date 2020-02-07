# Run the box model

include("getParams.jl")
include("box.jl")
include("getIC.jl")
include("makeTime.jl")
include("getEms.jl")
include("makeObs.jl")
using DifferentialEquations
using Dates
using Interpolations
using ParameterizedFunctions
using MAT




sYear = 1980; # start year
eYear = 2017; # end year
time_res = "year"; # time step resolution
tspan = makeTime(sYear, eYear, time_res); # create time vector




# create the params struct for easier input
params = getParams();
IC = getIC(params); # Get the initial conditions 
ems = getEms(params, tspan);
ems_conv = convertEms(ems, params, tspan);
obs = readMatObs("obs.mat");

params = getParams();
ems_conv = [ems_conv[:, 1] ems_conv[:,3] ems_conv[:,5]];
# define the box model macro
    k_ch4, k_co, k_mcf = params["k_ch4"], params["k_co"], params["k_mcf"];

    box = @ode_def begin
        dnh_ch4 = s_nh_ch4 - k_ch4 * nh_ch4 * nh_oh;
        dnh_co = s_nh_co + k_ch4 * nh_ch4 * nh_oh - k_co * nh_co * nh_oh;
        dnh_oh = s_nh_oh - k_ch4 * nh_ch4 * nh_oh - k_co * nh_co * nh_oh - kx_nh * nh_oh;
    end s_nh_ch4 s_nh_co s_nh_oh kx_nh

ts = (tspan[1], tspan[end]);
IC = [IC[1], IC[3], IC[5]];
    problem = ODEProblem(box, IC, ts, ems_conv);
    out = solve(problem, alg_hints = "stiff",saveat=tspan, dtmax=params["YrToDay"]/12);
