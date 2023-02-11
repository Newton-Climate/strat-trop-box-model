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


#obs = readMatObs("obs.mat", ems);
#priors = defPriors(ems2);
#Run box model
println("Running model with Stratosphere")
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
IC = sol.u[end];

pert_ind = 12*50;
ems[pert_ind,1] += 10*12;
ems[pert_ind,2] += 10*12;
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);

n = parameters["n_air"];
con_with_strat = StructToArray(sol.u, parameters["n_air"]);
nh = con_with_strat[end,5]
sh = con_with_strat[end,6]
@show nh
@show sh

lifetime = CalcPertLifetime(con_with_strat, pert_ind)
println(lifetime)




### No Stratosphere
println("Only Troposphere run")
flags["use_strat"] = false;
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
IC = sol.u[end];

pert_ind = 12*50;
ems[pert_ind,1] += 10*12;
ems[pert_ind,2] += 10*12;
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);

n = parameters["n_air"];
con_only_trop = StructToArray(sol.u, parameters["n_air"]);
nh = con_only_trop[end,5]
sh = con_only_trop[end,6]
@show nh
@show sh

lifetime = CalcPertLifetime(con_only_trop, pert_ind)
println(lifetime)

time = [i for i = 1:length(tspan)] ./12 .- pert_ind/12;
start_ind, end_ind = pert_ind-12, pert_ind+20*12
time = [start_ind:end_ind]
con_strat = con_with_strat[start_ind:end_ind,:]
con_trop = con_only_trop[start_ind:end_ind,:]


p1 = plot(time, con_strat[:,1], label="with stratosphere", color="red")
plot!(time, con_trop[:,1], label="only troposphere", color="green")
title!("Northern Troposphere")

p2 = plot(time, con_strat[:,2], label=:false, color="red")
plot!(time, con_trop[:,2], label=:false, color="green")
title!("Southern Troposphere")

p3 = plot(time, con_strat[:,11], label=:false, color="red")
plot!(time, con_trop[:,11], label=:false, color="green")
title!("Northern Stratosphere")

p4 = plot(time, con_strat[:,12], label=:false, color="red")
plot!(time, con_trop[:,12], label=:false, color="green")
title!("Southern Stratosphere")

plot(p1, p2, p3, p4, layout=(2,2))
plot!(ylabel="ppb", xlabel="years since perturbation")
plot!(fontfamily="serif-roman", legendfont=font("Computer Modern", 7))
savefig("test_lifetimes.pdf")

# globally averaged rresults 
trop_only_trop = (con_trop[:,1] + con_trop[:,2])/2
trop_only_strat = (con_strat[:,1] + con_strat[:,2])/2


# the stratospheric results 
trop_only_strat = (con_trop[:,11] + con_trop[:,12])/2
trop_only_strat = (con_strat[:,11] + con_strat[:,12])/2
