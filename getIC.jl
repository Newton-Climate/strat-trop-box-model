function getIC(params)




    # Troposphere values 
    nh_ch4 = 1575; # ppb
    sh_ch4 = 1525; # ppb

    nh_mcf = 78; # ppt
    sh_mcf = 68; # ppt

    nh_n2o = 300; # ppb
    sh_n2o = 300; # ppb

    nh_oh    = (params["gmOH"]/ params["n_air"])*1e9;              # ppb
    sh_oh    = (params["gmOH"]/params["n_air"])*1e9;              # ppb

    nh_co = 95; # ppb
    sh_co = 67; # ppb

    IC = [nh_ch4, sh_ch4, nh_co, sh_co, nh_oh, sh_oh, nh_mcf, sh_mcf, nh_n2o, sh_n2o];
    return IC
end


