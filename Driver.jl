# Run the box model

include("getParams.jl")
include("box.jl")
include("getIC.jl")
include("makeTime.jl")
include("getEms.jl")
using DifferentialEquations
using Dates
using Interpolations
using ParameterizedFunctions




sYear = 1980; # start year
eYear = 2017; # end year
time_res = "year"; # time step resolution
tspan = makeTime(sYear, eYear, time_res); # create time vector




# create the params struct for easier input
params = getParams();
IC = getIC(params); # Get the initial conditions 
ems = getEms(params, tspan);
ems2 = convertEms(ems, params, tspan);


#Run box model
sol, model = BoxModelWrapper(IC, ems, params, tspan);
#print(sol.u)
