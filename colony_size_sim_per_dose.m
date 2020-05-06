function [raw_dose_sim, fused_dose_sim] = colony_size_sim_per_dose(ips, dd)
% This is the simulation function performing on the level of doses, given
% the dose index dd and the input parameters ips.


D = ips.Ds(dd);
%disp("Dose "+num2str(D));


%Here we choose between the different types of underlying simulation
%procedures. At the moment only "initial" is a viable option
if ips.sim_type == "nonspatial"
    fun = @(ips,w,D) colony_size_sim_nonspatial(ips,w,dd);
elseif ips.sim_type == "initial"
    fun = @(ips,w,D) colony_size_sim_per_well_arrayed(ips,w,dd);
    %fun = @(ips,w,D) colony_size_sim_per_well(ips,w,dd);
elseif ips.sim_type == "gradual"
    fun = @(ips,w,D) colony_size_sim_per_well2_arrayed(ips,w, dd);
    %fun = @(ips,w,D) colony_size_sim_per_well2(ips,w,dd);
else
    disp("Please choose a valid simulation style");
    return
end


first = true;
%iterate over all wells, for each well run the function as defined by
%ips.sim_type
for w = 1:ips.N_w(dd)
    %disp("well "+num2str(w)+" of "+num2str(ips.N_w)+" wells");
    [raw_well_sim, fused_well_sim] = fun(ips, w, dd);
    if first
        raw_dose_sim =raw_well_sim;
        raw_dose_sim2 = raw_dose_sim;
        fused_dose_sim =fused_well_sim;
        fused_dose_sim2 = fused_dose_sim;
        first = false;
    else
        raw_dose_sim2 = cat(1,raw_dose_sim2, raw_well_sim);
        for f = fields(raw_well_sim)'
            fn = f{:};
            raw_dose_sim.(fn) = cat(1,raw_dose_sim.(fn), raw_well_sim.(fn));
        end
        fused_dose_sim2 = cat(1,fused_dose_sim2, fused_well_sim);
        for f = fields(fused_well_sim)'
            fn = f{:};
            fused_dose_sim.(fn) = cat(1,fused_dose_sim.(fn), fused_well_sim.(fn));
        end
    end
end
end