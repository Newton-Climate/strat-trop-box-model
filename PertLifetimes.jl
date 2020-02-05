function CalcPertLifetime(con_array, pert_ind)

    global_ch4 = (con_array[:,1] + con_array[:,2])/2;


    ch4_baseline = global_ch4[1] *ones(size(global_ch4));
    delta_ch4 = global_ch4 - ch4_baseline;
    max1, max_ind = findmax(delta_ch4);
    m1, ind1 = findmin(abs.(max1/exp(1)*ones(size(ch4_baseline)) - delta_ch4));

    lifetime = ind1 - pert_ind;
    return lifetime, ind1;
end #function CalcPertLifetime

    
                               
