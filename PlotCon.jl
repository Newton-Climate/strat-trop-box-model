using Plots, LaTeXStrings


function plotCon(con_array, tspan)
    """
Plots the concentrations of each box
Inputs:
con_arraych
: An array of the concentrations
Outputs:
A plot object
"""
    species_label = [L"CH_4", "CO", "OH", "MCF", L"N_2O"];
    NH_trop = [i for i = 1:2:9];
    SH_trop = [i for i = 2:2:10];
    NH_strat = [i for i = 11:2:19];
    SH_strat = [i for i = 12:2:20];

    # create an array of indexes for the for-loop
    # order is species_indexes[box_number, species]

    box_names = ["NH Troposphere", "SH Troposphere", "NH Stratosphere", "SH Stratosphere"];
    num_boxes = length(box_names);
    num_species = length(NH_trop);
    species_indexes_array = [NH_trop SH_trop NH_strat SH_strat];
    
    time = [i for i = 1:length(tspan)];
    
    y_units = ["ppb", "ppb", L"molecules/cm^3", "ppt", "ppb"];

    # iterate through each species 
    # 
    for box = 1:num_boxes
        species_indexes = species_indexes_array[:, box]
        P = Array{Any}(undef, num_species, 1);
        iteration = 1;
        for species in species_indexes
            P[iteration] = plot(time, con_array[:, species], label= species_label[iteration], ylabel = y_units[iteration], xlabel = "years", lw = 3);
            iteration += 1;
        end # species for-loop

        fig = plot(P[1], P[2], P[3], P[4], P[5], layout = (5,1), title = box_names[box]);
        filename = string("figures/", box_names[box], ".pdf");
        println(filename)
#        display(fig)
        savefig(filename)
    end # box for-loop
end # function plotCon

        

    
