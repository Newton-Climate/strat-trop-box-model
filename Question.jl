using DifferentialEquations, DiffEqBayes, ParameterizedFunctions
using Distributions, Interpolations

end_time = 25.0;
time = collect(1.0:end_time);
num_steps = length(time);
tspan = (1.0 ,end_time);

# define some constants 
const k1, k2, k3 = 7.89, 464.8, 3600*24;

function interpSources(s, t)
    num_time, num_species = size(s)
    species_grid = [i for i =1:num_species];
    grid = (time, species_grid);
    intp = interpolate(grid, s, Gridded(Linear()));
    s_out = intp(t, species_grid);
    return s_out
end #function interpSources

    
function BoxModel(dc ,c ,p ,t)
    s = interpSources(p ,t) # interpolate the source terms (parameters) in time

    dc[1] = s[1] - k1 * c[1] * c[3];
    dc[2] = s[2] + k1 * c[1] * c[3] - k2 * c[2] * c[3];
    dc[3] = s[3] - k1 * c[1] * c[3] - k2 * c[2] * c[3] - k3 * c[3];
end # function BoxModel

# define initial conditions
c_0 = [1700.0 100.0 3.71e-5];

# defining the source-terms (parameters) where row i is time i and col j is species j
sources = [1.05 * ones(num_steps, 1) 1.06 * ones(num_steps, 1) 7.74 * ones(num_steps, 1)];

# define and solve the forward problem
problem = ODEProblem( BoxModel ,c_0 ,tspan ,sources);
sol = solve(problem, alg_hints = "stiff",saveat=time);

# inversion
priors = [ Normal( sources[1,1] ,0.1 * sources[1,1] ) ,Normal(sources[ 1 ,2], 0.1 * sources[ 1,2]) ,Normal( sources[1 ,3] ,0.1 * sources[1,3])];
obs = sol.u;
result_turing = turing_inference( problem ,Tsit5(),tspan,obs,priors;num_samples=1000)

