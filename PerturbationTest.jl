# Run the box model
using Printf, Plots
include("getParams.jl")
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

#obs = readMatObs("obs.mat", ems);
#priors = defPriors(ems2);
#Run box model
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
IC = sol.u[end];

pert_ind = 12*50;
ems[pert_ind,1] += 10*12;
ems[pert_ind,2] += 10*12;
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);

#obs = obs';
#priors = priors';

#bayesian_result = turing_inference(model,Tsit5(),tspan,obs,priors,num_samples=10)
#monte_prob = MonteCarloProblem(model)

# obj = build_loss_objective(monte_prob,SRIW1(),L2Loss(tspan,obs,differ_weight=1.0,data_weight=0.5),maxiters=1000,
#                                  verbose=false,verbose_opt=false,verbose_steps=1,num_monte=50)
#result = Optim.optimize(obj, ems2, Optim.BFGS())

nh = sol.u[end][5];
n = parameters["n_air"];
sh = sol.u[end][6];
@printf("NH OH is %1.2e1f \n", nh*n/1e9)
@printf("SH OH is %1.2e1f \n", sh*n/1e9)

con = StructToArray(sol.u);
lifetime, ind = CalcPertLifetime(con, pert_ind);
print(lifetime/12, "\n");
time = [i for i = 1:length(tspan)];

#p = plot(tspan, con[:,1:2], label=["Northern Hemisphere", "Southern Hemisphere"], lw=2, lw=3);
#ylabel!("ppb")
#xlabel!("time")
#plot!(xlims = (pert_ind - 5*12, time[end]), xticks=([i for i = pert_ind - 5*12:5*12:length(tspan)], string.([ i for i = -5:5:60])))
#savefig("perturbation_test_strat.pdf")

