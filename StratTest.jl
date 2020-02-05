# OBRun the box model
using Printf, Plots
include("getParams.jl")
include("PlotCon.jl")
include("NewBox.jl")
include("getIC.jl")
include("makeTime.jl")
include("getEms.jl")
include("makeObs.jl")
include("lifetimes.jl")
#include("invert.jl")
using DifferentialEquations
using Dates
using Interpolations
using ParameterizedFunctions
using MAT
#using Turing, DiffEqBayes, Optim
using DiffEqParamEstim




sYear = 1980; # start year
eYear = 2030; # end year
time_res = "year"; # time step resolution
tspan = makeTime(sYear, eYear, time_res); # create time vector
flags = Dict{String, Bool}();
flags["use_strat"] = true;




# create the params struct for easier input
parameters = getParams();
IC = getIC(parameters); # Get the initial conditions 
ems = getEms(parameters, tspan);

#obs = readMatObs("obs.mat", ems);
#priors = defPriors(ems2);
#Run box model
#sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
#IC = sol.u[end];
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
fixed_con = StructToArray(sol.u);

ems[20:end,11] = 1.3 * ems[20:end,11];
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
pert_con = StructToArray(sol.u);

#obs = obs';
#priors = priors';

#bayesian_result = turing_inference(model,Tsit5(),tspan,obs,priors,num_samples=10)
#monte_prob = MonteCarloProblem(model)

# obj = build_loss_objective(monte_prob,SRIW1(),L2Loss(tspan,obs,differ_weight=1.0,data_weight=0.5),maxiters=1000,
#                                  verbose=false,verbose_opt=false,verbose_steps=1,num_monte=50)
#result = Optim.optimize(obj, ems2, Optim.BFGS())


n = parameters["n_air"];
fixed_ch4 = (fixed_con[:,1] + fixed_con[:,2])/2;
pert_ch4 = (pert_con[:,1] + pert_con[:,2])/2;

# convert [OH] to ppb
pert_con[:,5] = pert_con[:,5] * n/1e9;
fixed_con[:,5] = fixed_con[:,5] * n/1e9;
pert_con[:,6] = pert_con[:,6] * n/1e9;
fixed_con[:,6] = fixed_con[:,6] * n/1e9;
#plotCon(fixed_con, tspan);

#y=[fixed_ch4, pert_ch4];
#time = [i for i =1:length(tspan)]
#plot(time, y, label=["fixed lifetime", "perturbed_lifetime"]);
#plot!(xlabel = "years", ylabel="ppb")
#plot!
#plot!(xticks=([1:5:tspan[end];], string.0:5:38)))
#savefig("strat-trop_exchange.pdf")
