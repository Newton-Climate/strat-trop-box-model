using Interpolations
using ParameterizedFunctions, DifferentialEquations
include("getParams.jl")
include("getEms.jl")
include("makeTime.jl")
include("getIC.jl")




sYear = 1980; # start year
eYear = 2017; # end year
time_res = "year"; # time step resolution
tspan = makeTime(sYear, eYear, time_res); # create time vector

params = getParams();
#time = [ i for i in range(1,1e6)];
ems = getEms(tspan, params);
ems_conv = convertEms(ems, params);
ems = ems[:,1:6];
IC = getIC(params);
IC = IC[1:6];

k_ch4, k_co, k_mcf = params["k_ch4"], params["k_co"], params["k_mcf"];
kx_nh, kx_sh = 0.99 * 3600*24, 1.23*3600*24; # arbitrary OH reaction rates 
τ_NS = 365.25 # 1 year in days 


box = @ode_def begin
        num_time, num_species = size(s);
        s = convertEms(s, params);
        species_grid = [i for i =1:num_species];
        grid = (tspan, species_grid);
        intp = interpolate(grid, s, Gridded(Linear()));
        s = intp(t, species_grid)


        dnh_ch4 = s_nh_ch4 - k_ch4 * nh_ch4 * nh_oh - (nh_ch4 - sh_ch4)/τ_NS;
        dsh_ch4 = s_sh_ch4 - k_ch4 * sh_ch4 * sh_oh - (sh_ch4 - nh_ch4)τ_NS;
        dnh_co = s_nh_co + k_ch4 * nh_ch4 * nh_oh - k_co * nh_co * nh_oh - (nh_co - sh_co)/τ_NS;
        dsh_co = s_sh_co + k_ch4 * sh_ch4 * sh_oh - k_co * sh_co * sh_oh - (sh_co - nh_co)/τ_NS;
        dnh_oh = s_nh_oh - k_ch4 * nh_ch4 * nh_oh - k_co * nh_co * nh_oh - kx_nh * nh_oh;
        dsh_oh = s_sh_oh - k_ch4 * sh_ch4 * sh_oh - k_co * sh_co * sh_oh - kx_sh * sh_oh;
    end s_nh_ch4 s_sh_ch4 s_nh_co s_sh_co s_nh_oh s_sh_oh
ts = (tspan[1], tspan[end]);
problem = ODEProblem(box, IC, ts, ems);
out = solve(problem, alg_hints = "stiff",saveat=tspan, dtmax=params["yrToDay"]/12);
