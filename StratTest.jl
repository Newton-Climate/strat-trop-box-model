# OBRun the box model
using Printf, Plots
include("getParams.jl")
include("PlotCon.jl")
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

# Define parameters
sYear = 1980; # start year
eYear = 2030; # end year
time_res = "year"; # time step resolution
tspan = makeTime(sYear, eYear, time_res); # create time vector
flags = Dict{String, Bool}();
flags["use_strat"] = true;

# create the params struct for easier input
parameters = getParams();
IC = getIC(parameters); # Get the initial conditions 
ems = getEms(parameters, tspan); # Construct the state vector (emissions and exchange time)

#Run box model with fixed 7-year STE lifetime
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
fixed_con = StructToArray(sol.u, parameters["n_air"]);

# perturb the STE lifetime τₜₛ
ems[20:end,11] = 1.3 * ems[20:end,11];
sol, model = BoxModelWrapper(IC, ems, parameters, tspan, flags);
pert_con = StructToArray(sol.u, parameters["n_air"]); # array with perturbed model results 

# Global average of CH4
fixed_ch4 = (fixed_con[:,1] + fixed_con[:,2])/2;
pert_ch4 = (pert_con[:,1] + pert_con[:,2])/2;

# Plot all boxes and species to figures/
plotCon(fixed_con, tspan);
