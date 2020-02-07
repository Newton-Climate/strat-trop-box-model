# Run the box model
using Printf, Plots
include("getParams.jl")
include("NewBox.jl")
include("getIC.jl")
include("makeTime.jl")
include("getEms.jl")
include("makeObs.jl")
include("PertLifetimes.jl")
#include("invert.jl")
using DifferentialEquations
using Dates
using Interpolations
using ParameterizedFunctions
using MAT
#using Turing, DiffEqBayes, Optim
using DiffEqParamEstim




sYear = 1980; # start year
eYear = 2100; # end year
time_res = "month"; # time step resolution
tspan = makeTime(sYear, eYear, time_res); # create time vector
flags = Dict{String, Bool}();
flags["use_strat"] = true;




# create the params struct for easier input
parameters = getParams();
IC = getIC(parameters); # Get the initial conditions 
ems = getEms(parameters, tspan);
ems2 = convertEms(ems, parameters);

#Run box model to steady-state
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
IC = sol.u[end];

pert_ind = 12*50; # when to apply perturbation?
ems[pert_ind,1] += 10*12;
ems[pert_ind,2] += 10*12;

# Run box model with stratosphere 
strat_sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
strat_con = StructToArray(strat_sol.u, parameters["n_air"]);

# Run box model without stratosphere
flags["use_strat"] = false;

#Run box model to steady-state
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
IC = sol.u[end];

trop_sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
trop_con = StructToArray(trop_sol.u, parameters["n_air"]);

strat_nh , strat_sh = strat_con[end , 5], strat_con[end , 6];
trop_nh , trop_sh = trop_con[end , 5], trop_con[end , 6];

println("Tropospheric-only runs:")
@printf("NH OH is %1.2e1 \n", trop_nh)
@printf("SH OH is %1.2e1 \n", trop_sh)
lifetime, ind = CalcPertLifetime(trop_con, pert_ind);
print(lifetime/12, "\n");

println("Strat and Trop run:")
@printf("NH OH is %1.2e1 \n", strat_nh)
@printf("SH OH is %1.2e1 \n", strat_sh)
lifetime, ind = CalcPertLifetime(strat_con, pert_ind);
print(lifetime/12, "\n");



time = [i for i = 1:length(tspan)];
trop_ch4 = 0.5*(trop_con[:,1] + trop_con[:,2]);
strat_ch4 = 0.5*(strat_con[:,1] + strat_con[:,2]);
con = [trop_ch4 strat_ch4];

p = plot(tspan, trop_ch4, label = "Without Stratosphere", lw=3, title = "20-Tg Perturbation Test");
plot!(tspan, strat_ch4, label="With Stratosphere", lw=3)
ylabel!("ppb")
xlabel!("time")
plot!(xlims = (pert_ind - 5*12, time[end]), xticks=([i for i = pert_ind - 5*12:5*12:length(tspan)], string.([ i for i = -5:5:60])))
savefig("figures/perturbation_test_strat.pdf")

