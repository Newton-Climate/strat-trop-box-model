using Distributions, Turing, Random
include("getEms.jl")

function defPriors(ems_prior, params)
    """ms
Defines the priors of the parameters to be optimized
Arguements:
ems_prior: the array of parameters to be optimized
output:
array of distributions corresponding to the parameters 
"""

    (rows, cols) = size(ems_prior);
    priors = Array{Any}(undef, rows, cols);
    Random.seed!(123) # Setting the seed
    CH4_NH, CH4_SH, CO_NH, CO_SH, OH_NH, OH_SH, MCF_NH, MCF_SH, N2O_NH, N2O_SH, SF6_NH, SF6_SH = [i for i in 1:12];

    # define the σ on each prior
    

    σ_CH4 = 20; # Tg/yr
    σ_OH = 315; # Tg/yr
    σ_CO = 300; # Tg/yr
    σ_N2O = 2; # Tg/yr
    #    σ_mcf_nh  = max([1.5*ones(size(ems_prior[:,7])),.2*ems_prior[:,7]],[],2).^2;
    σ_MCF_NH = 10;
    σ_MCF_SH  =   0.5;
    σ_τ = 2.5; # yr

    

    σ = [σ_CH4, σ_CH4,  σ_CO, σ_CO, σ_OH, σ_OH, σ_MCF_NH, σ_MCF_SH, σ_N2O, σ_N2O, σ_τ]';
#    σ = convertEms(σ, params);


    # define the μ
    μ = ems_prior[1,:];
    σ = 0.1 * μ;


    for i=1:rows
        for j = 1:cols
            priors[i,j] = Normal(μ[j], σ[j]);
        end #for loop
    end # for loop
    return priors
end #function DefPriors


    
@model function infer_ems(ode_model, priors, obs)
    println("running MCMC chain")
    p = priors
    prob = remake(ode_model, p=p)
    #predicted = solve(prob) #,Tsit5(),saveat=0.1)
    predicted = solve(problem, alg_hints = "stiff") #,saveat=tspan, dtmax=params["YrToDay"]/100);
end

#hmc_model = infer_ems(model, priors, obs)
#sample(hmc_model, NUTS(.65),1000)
#chain = mapreduce(c -> sample(hmc_model, NUTS(.65),1000), chainscat, 1:3)
