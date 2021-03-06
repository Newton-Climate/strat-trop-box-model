function getParams()
    #1. Define the parameters
    # 2. Remember to use the unitful package to establsh units
    # 3. try using unitful to convert to ppb/day

    #    @unit molec "molec" Molecules 5u"yd" false

    τᵢ = 1; # North-South troposphere exchagne time
    τᵢ_strat = 3.3;             # Interhemispheric exchange rate in the stratosphere (Fabian et al., 1968)


            # Unit conversions
    DaysToS = 60 * 60 * 24;         # Days to Seconds
    YrToDay = 365.25;               # Years to Days

        
    # define rate constants
    k_ch4 = 3.395e-15; # molec/cm^3/s
    k_co = 2e-13;         # reaction rate (cm3/molec/s) from Prather (1993)
    k_mcf = 6.05e-15; # molec/cm^3/s
    k_co_strat     = 12;            # (yr^-1) Assuming lifetime of 1 month in the stratosphere
    k_oh_strat     = 12;            # (yr^-1) Assuming lifetime of 1 month in the stratosphere
    k_co_strat   = k_co_strat   / YrToDay;
    k_oh_strat   = k_oh_strat   / YrToDay;

    k_ch4_other    = 1/151;         # (yr^-1) Soil uptake and tropospheric chlorine (6% of total loss; Kirschke et al., 2013)
    k_co_other     = 0;             # (yr^-1) Other CO losses (currently neglecting)

        
    
# Masses
m       = 5.15e21;              # Rough: Total mass of atmosphere in g
n_air   = 2.69e19;              # Rough: number density of air (molec/cm3)
m_air   = 28.8;                 # Average molar mass of the atmosphere in g/mol (mostly N2 and O2)
m_ch4   = 16;                   # Molar mass of CH4
m_mcf   = 133.4;                # Molar mass of methyl-chloroform
m_n2o   = 44.013;               # Molar mass of nitrous oxide
m_c2h6  = 30.07;                # Molar mass of ethane
m_oh    = 17.008;               # Molar mass of hydroxyl
m_co    = 28.01;                # Molar mass of carbon monoxide
mConv   = m/m_air/1e12*1e-9;	# Factor to convert mixing ratios to Tg
RxNconv = n_air / 1e9;          # Convert ppb to molec/cm3
mm_ch4  = m_ch4  * mConv;       # Convert CH4 mixing ratios (ppb) to Tg
mm_mcf  = m_mcf  * mConv;       # Convert MCF mixing ratios (ppt) to Gg
mm_n2o  = m_n2o  * mConv;       # Convert N2O mixing ratios (ppb) to Tg
mm_c2h6 = m_c2h6 * mConv;       # Convert C2H6 mixing ratios (ppt) to Gg
mm_oh   = m_oh  * mConv;        # Convert OH number densities (ppb) to Tg
mm_co   = m_co  * mConv;        # Convert CO mixing ratios (ppb) to Tg

t_n2oNH        = 163;           # N2O lifetime (yr) in the NH (Brown et al., 2013; Table 4)
t_n2oSH        = 109;           # N2O lifetime (yr) in the SH (Brown et al., 2013; Table 4)
t_ch4_strat_nh = 188;           # CH4 lifetime (yr) in the NH (Brown et al., 2013; Table 4)
t_ch4_strat_sh = 200;           # CH4 lifetime (yr) in the SH (Brown et al., 2013; Table 4)
t_mcf_strat    = 31;            # 15# of loss is in stratosphere, using MCF lifetime of 5.5 yr (ratio of losses is proportional to inverse ratio of lifetimes Makide & Rowland, 1981)
gmOH = 1.0e6;                   # Global mean OH (molec/cm^3)
RxNconv      = RxNconv * DaysToS;
        
### Make a structure with the parameters
    # Unit conversions
    params = Dict{String, Float64}();
params["mm_ch4"]  = mm_ch4;
params["mm_mcf"]  = mm_mcf;
params["mm_n2o"]  = mm_n2o;
params["mm_c2h6"] = mm_c2h6;
params["mm_oh"]   = mm_oh;
params["mm_co"]   = mm_co;
params["YrToDay"] = YrToDay;
params["n_air"]   = n_air;
    params["gmOH"]    = gmOH;
    
# Rate constants/lifetimes
params["k_ch4"]        = RxNconv * k_ch4;   
params["k_mcf"]          = RxNconv * k_mcf;
params["k_n2o_strat_nh"]       = 1/(t_n2oNH * YrToDay);
params["k_n2o_strat_sh"]       = 1/(t_n2oSH * YrToDay);
params["k_ch4_strat_nh"] = 1/(t_ch4_strat_nh * YrToDay);
params["k_ch4_strat_sh"] = 1/(t_ch4_strat_sh * YrToDay);
params["k_mcf_strat"]    = 1/(t_mcf_strat * YrToDay);



params["k_co"]           = RxNconv * k_co;
params["k_co_strat"]     = k_co_strat;
params["k_co_other"]     = k_co_other;
params["k_oh_strat"]     = k_oh_strat;
params["RxNconv"]        = RxNconv;
    params["DaysToS"]        = DaysToS;
    params["τᵢ"] = τᵢ * YrToDay; # 365.25 days
    params["τᵢ_strat"] = τᵢ_strat * YrToDay;


        return params;
end #function
