function makeBox(params, tspan)


    # 1. c(1) = OH
# 2. C(3) = CH4
# 3. c(5) = CO
#c(7) = MCF
#c(9) = N2O
#C(11) = SF6

# defining indecies 
    CH4_NH, CH4_SH, CO_NH, CO_SH, OH_NH, OH_SH, MCF_NH, MCF_SH, N2O_NH, N2O_SH, SF6_NH, SF6_SH = [i for i in 1:12];

# CH4 + OH
R1(c,CH4, OH) = params["k_ch4"] * c[CH4] * c[OH];

    # CO + Oh
    R2(c,CO, OH) = params["k_co"] * c[CO]*c[OH];


# MCF + OH 
    R3(c,MCF, OH) = params["k_mcf"] * c[MCF] * c[OH];

    tau_i = 1*365.25;
    tau_s = 3*365.25;
    kX_NH = 0.99* 60*60*24;
    kX_SH = 1.23*60*60*24;

    inter_trop(c, con_current, con_other, tau_i) = (c[con_other] - c[con_current])/tau_i;

    strat_trop(c, con_trop, con_strat, tau_s) = (c[con_strat] - c[con_trop])/tau_s;


    # indivitual kinetic equations:

    dCH4_NH(s,c, t) = s[CH4_NH] - R1(c, CH4_NH, OH_NH) + inter_trop(c, CH4_NH, CH4_SH, tau_i);
    dCH4_SH(s,c, t) = s[CH4_SH] - R1(c, CH4_SH, OH_SH) + inter_trop(c, CH4_SH, CH4_NH, tau_i);

    dOH_NH(s,c,t) = s[OH_NH] - R1(c, CH4_NH, OH_NH) - R2(c, CO_NH, OH_NH) - kX_NH*c[OH_NH];
    dOH_SH(s,c,t) = s[OH_SH] - R1(c, CH4_SH, OH_SH) - R2(c, CO_SH, OH_SH) - kX_SH*c[OH_SH];


    dCO_NH(s,c, t) = s[CO_NH] + R1(c, CH4_NH, OH_NH) - R2(c, CO_NH, OH_NH) + inter_trop(c, CO_NH, CO_SH, tau_i);
    dCO_SH(s,c, t) = s[CO_SH] + R1(c, CH4_SH, OH_SH) - R2(c, CO_SH, OH_SH) + inter_trop(c, CO_SH, CO_NH, tau_i);
    dMCF_NH(s,c, t) = s[MCF_NH] - R3(c, MCF_NH, OH_NH) + inter_trop(c, MCF_NH, MCF_SH, tau_i);
    dMCF_SH(s,c, t) = s[MCF_SH] - R3(c, MCF_SH, OH_SH) + inter_trop(c, MCF_SH, MCF_NH, tau_i);

    dN2O_NH(s,c, t) = s[N2O_NH] + inter_trop(c, N2O_NH, N2O_SH, tau_i);
    dN2O_SH(s,c, t) = s[N2O_SH] + inter_trop(c, N2O_SH, N2O_NH, tau_i);

    dSF6_NH(s,c, t) = s[SF6_NH] + inter_trop(c, SF6_NH, SF6_SH, tau_i);
    dSF6_SH(s,c, t) = s[SF6_SH] + inter_trop(c, SF6_SH, SF6_NH, tau_i);

    function dc_dt(dc, c, s, t)

        num_time, num_species = size(s);
        species_grid = [i for i =1:num_species];
        grid = (tspan, species_grid);
        intp = interpolate(grid, s, Gridded(Linear()));
        s = intp(t, species_grid)


        dc[1:end] = [dCH4_NH(s,c,t), dCH4_SH(s,c,t), dCO_NH(s,c,t), dCO_SH(s,c,t), dOH_NH(s,c,t), dOH_SH(s,c,t), dMCF_NH(s,c,t), dMCF_SH(s,c,t), dN2O_NH(s,c,t), dN2O_SH(s,c,t)];
        end

        return dc_dt

end



function BoxModelWrapper(IC, ems, params, tspan)

    # convert ems from Tg/yr to ppx/day
    ems_conv = convertEms(ems, params, tspan);
    ts = (tspan[1], tspan[end]);
    model = makeBox(params, tspan);



    # run the model
    problem = ODEProblem(model, IC, ts, ems_conv);

    out = solve(problem, alg_hints = "stiff",saveat=tspan, dtmax=params["YrToDay"]/12);
    return out, model
end


