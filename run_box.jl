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
eYear = 2200; # end year
time_res = "year"; # time step resolution
tspan = makeTime(sYear, eYear, time_res); # create time vector
flags = Dict{String, Bool}();
flags["use_strat"] = false;




# create the params struct for easier input
parameters = getParams();
IC = getIC(parameters); # Get the initial conditions 
ems = getEms(parameters, tspan);


#obs = readMatObs("obs.mat", ems);
#priors = defPriors(ems2);
#Run box model
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
eq = sol.u[end];

con = StructToArray(sol.u, parameters["n_air"]);
nh = con[end,5]

sh = con[end,6]
@show nh
@show sh
