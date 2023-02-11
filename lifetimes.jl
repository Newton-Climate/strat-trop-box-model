function CalcPertLifetime(con_array, pert_ind)

    global_ch4 = (con_array[:,1] + con_array[:,2])/2;
    max_con, max_ind = findmax(global_ch4)
    e = exp(1);
    ch4_baseline = global_ch4[max_ind-12] *ones(size(global_ch4));
    delta_ch4 = global_ch4 - ch4_baseline;
    max1 = maximum(delta_ch4)
    m1, ind1 = findmin(abs.(max1/exp(1) .- delta_ch4));
    return (ind1-max_ind)/12;
end

                               
