function getIC(params)




    # Troposphere values 
    nh_ch4 = 1575; # ppb
    sh_ch4 = 1520; # ppb

    nh_mcf = 78; # ppt
    sh_mcf = 68; # ppt

    nh_n2o = 302; # ppb
    sh_n2o = 301; # ppb

    nh_oh    = (params["gmOH"]/ params["n_air"])*1e9;              # ppb
    sh_oh    = (params["gmOH"]/params["n_air"])*1e9;              # ppb

    nh_co = 95; # ppb
    sh_co = 67; # ppb

    ### Stratospheric conditions
    nh_ch4_S = 1570;                 # ppb
    sh_ch4_S = 1515; # ppb
    nh_oh_S    = (1.0e6/params["n_air"])*1e9;  # ppb
    sh_oh_S    = (1.0e6/params["n_air"])*1e9;           # ppb
    nh_co_S    = 70;                          # ppb
    sh_co_S    = 50;                          # ppb

    nh_mcf_S   = 30;                          # ppt
    sh_mcf_S   = 20;                          # ppt
    nh_n2o_S   = 302.0;                       # ppb
    sh_n2o_S   = 301.0;                       # ppb


    IC = [nh_ch4, sh_ch4, nh_co, sh_co, nh_oh, sh_oh, nh_mcf, sh_mcf, nh_n2o, sh_n2o, nh_ch4_S, sh_ch4_S, nh_co_S, sh_co_S, nh_oh_S, sh_oh_S, nh_mcf_S, sh_mcf_S, nh_n2o_S, sh_n2o_S];
    return IC
end


