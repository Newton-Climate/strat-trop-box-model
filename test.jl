using ParameterizedFunctions, DiffEqBayes, DifferentialEquations, Stan, Distributions
f1 = @ode_def LotkaVolterraTest1 begin
 dx = a*x - x*y
 dy = -3*y + x*y
end a

p = [1.5]
u0 = [1.0,1.0]
tspan = (0.0,10.0)
prob1 = ODEProblem(f1,u0,tspan,p)

σ = 0.01                         # noise, fixed for now
t = collect(range(1,stop=10, length=10))   # observation times
sol = solve(prob1,Tsit5())

randomized = VectorOfArray([(sol(t[i]) + σ * randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)

bayesian_result_hmc = dynamichmc_inference(prob1, data, [Normal(1.5, 1)], t, [bridge(ℝ, ℝ⁺, )])
