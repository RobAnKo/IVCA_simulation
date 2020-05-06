function [raw_sim_output, fused_sim_output, tracked_sim_output] = colony_size_sim_main(ips)

%input: 

%experimental setup parameters
%N -                number of objects per well
%N_w -              number of wells per dose (can be a scalar or a vector 
%                   of scalars, one for each dose)
%times -            vector of time points (arbitrary units, but need to be 
%                   coherent with other values, usually hours)
%Ds -               vector of doses


%"biological" parameters
%PE -               plating efficiency: which fraction of colonies stays 
%                   alive and growing without irradiation (growth rate != 0)
%alpha -            an optional parameter inspired by the LQM: linear 
%                   exponential dose-dependent decrease in survival 
%                   probability 
%beta -             an optional parameter inspired by the LQM: quadratic 
%                   exponential dose-dependent decrease in survival 
%                   probability 

%doub_time -        one way to set growth rates in absence of irradiation: 
%                   using doubling time (inverse of growth rate) given in 
%                   hours [h]
%gamma_dt -         how does doubling time change with dose? [hGy⁻¹]
%eta0 -             another way: using the growth rate. Both .doub_time and
%                   .eta0 indicate the mean value in absence of irradiation
%                   given in inverse hours [h⁻¹]
%gamma_eta -        how does the mean growth rate change per Gy dose? 
%                   [h⁻¹Gy⁻¹]


%std_type -         an option indicating whether the growth rate standard 
%                   deviation .gr_sd is to be interpreted as a relative 
%                   value ("rel") or as an absolute value ("abs")
%gr_sd -            dose-independent standard deviation of growth rate. 
%                   Dimensionless in [0,1] if .std_type is "rel", given in 
%                   [h⁻¹] if .std_type is "abs"
%gr_sd_f -          how does the growth rate standard deviation change with
%                   dose [Gy⁻¹] or [h⁻¹Gy⁻¹] depending on .std_type

%mit_del -          mitotic delay: how much is onset of growth delayed per 
%                   dose [h⁻¹Gy⁻¹]
%in_size_type -     an option to change the source of initial sizes. 
%                   If "norm" is given, the values in .in_size_range are 
%                   inpterpreted as mean and std of a normal distribution 
%                   from which to draw sizes. If "distr" is given, sizes 
%                   are drawn from the .in_size_distr, which can either be 
%                   a distribution object or a vector of sizes (e.g. from 
%                   experimental data) from which to draw.
%in_size_range -    mean and std of a potential distribution from which to 
%                   draw initial sizes. Only has an effect if .in_size_type
%                   is "norm"
%in_size_distr -    a distribution or a vector of size values from which 
%                   to draw. Only has an effect if .in_size_type is "distr"
%in_size_min -      a lower threshold on initial sizes. If sizes below this
%                   threshold are drawn as initial sizes, they are redrawn 
%                   until they are larger than .in_size_min
%shrink -           a boolean flag which indicates whether dead colonies 
%                   should shrink (with a hardcoded rate somewhere in the 
%                   function code). This is to emulate the observation of 
%                   shrinking colonies in the exp. data

%filter_small -     a boolean value indicating whether objects smaller than
%                   .in_size_min should be removed from the data (emulates 
%                   the filtering step happening on the exp. data). 
%                   Only makes sense if .shrink is true, because otherwise 
%                   all objects stay >= .in_size_min anyway

%alpha_h, mu -      not fully developed: a hazard for single cells to cease
%                   growth. apha_h is the initial hazard for a single cell 
%                   to cease growth, mu is the rate with which this hazard 
%                   changes. This represents one way to implement 
%                   time-dependent changes in effective growth rates. 
%                   Only has effects if the sim_type is "gradual"


%technical/output parameters
%im_dims -          sizes in [x,y]-dimensions of the original images with 
%                   which we want to compare the simulation results. 
%                   Given in pixels
%well_radius -      radius of the circular well in pixels. Should be 
%                   equivalent to the radius of wells in the original data 
%                   and less than half of min(im_dims)
%to_image -         boolean flag indicating whether images should be 
%                   created from the simulation

%to_file -          boolean flag indicating whether the result(s) should be 
%                   stored as a csv-file

%sim_type -         the type of simulation. At the moment only "initial" is 
%                   valid. In this type, growth rates and death state are 
%                   initially set and remain unchanged throughout the 
%                   simulation. Other possible types are "nonspatial", 
%                   where space is not modeled and hence no interaction 
%                   (fusion) between colonies happen, or "gradual", where 
%                   individual cells have a gradually changing probability 
%                   to stop growing (based on .alpha_h and .mu)

%calculate_fused -  boolean flag indicating whether an additional output 
%                   containing fused colonies should be calculated

%calculate_tracked -boolean flag indicating whether an additional output 
%                   containing tracked colonies should be calculated. This 
%                   emulates the tracking happening on experimental data 
%                   to compare the influence of tracking on the resulting
%                   size distributions

%track_max_dist -   the maximum distance between two objects to allow 
%                   mapping during tracking. Only valid if 
%                   .calculate_tracked is true;



%output_file_path - an optional path where the resulting datatable(s) 
%                   should be stored.


%well_idx -         optional: a vector of indices indication which "pixels" 
%                   in a matrix are wells. This is to avoid multiple 
%                   recalcuations of the same values in cases of repetitive
%                   simulation runs (e.g. for parameter fitting)
%seed_idx -         optional: the same concept for indices at which 
%                   colonies are initially seeded. 
%dist_mat -         If seed indices are given, the distance matrix between 
%                   those seeds can be precalculated as well. All three 
%                   options are used in the context of repetitive calls to 
%                   the simulation function e.g. during parameter fitting




%% default values
default.N = 100;
default.PE = 1;
default.N_w = 10;
default.times = 0:3:219;
default.Ds = 0:7;
default.in_size_type = "norm";
default.in_size_range = [90,20];
default.in_size_distr = makedist("Normal",...
    "mu", default.in_size_range(1),...
    "sigma", default.in_size_range(2));
default.in_size_min = 55;
default.doub_time = 24;
default.eta0 = nan;
default.gamma_eta = 0;
default.std_type = "rel";
default.gr_sd = 0.25;
default.gr_sd_f = 0;
default.mit_del = 1;
default.alpha = 0;
default.beta = 0;
default.gamma_dt = 1;
default.alpha_h = 0.01;
default.mu = -log(1/2)/default.doub_time;
default.im_dims = [2240,2240];
default.well_radius = 1000;
default.to_image = false;
default.to_file = false;
default.sim_type = "initial";
default.calculate_fused = false;
default.calculate_tracked = false;
default.track_max_dist = 10;
default.well_idx = nan;%10000:100000;
default.seed_idx = nan;
default.dist_mat = nan;
default.output_file_path = '/home/robinkoch/Documents/SAIVCA/out/simulation/sim_output.csv';
default.shrink = false;
default.filter_small = true;

p = inputParser;
addParameter(p, 'N', default.N);
addParameter(p, 'PE', default.PE);
addParameter(p, 'N_w', default.N_w);
addParameter(p, 'times', default.times);
addParameter(p, 'Ds', default.Ds);
addParameter(p, 'in_size_type', default.in_size_type);
addParameter(p, 'in_size_distr', default.in_size_distr);
addParameter(p, 'alpha', default.alpha);
addParameter(p, 'beta', default.beta);
addParameter(p, 'doub_time', default.doub_time);
addParameter(p, 'eta0', default.eta0);
addParameter(p, 'gamma_eta', default.gamma_eta);
addParameter(p, 'mu', default.mu);
addParameter(p, 'alpha_h', default.alpha_h);
addParameter(p, 'std_type', default.std_type);
addParameter(p, 'gr_sd', default.gr_sd);
addParameter(p, 'gr_sd_f', default.gr_sd_f);
addParameter(p, 'mit_del', default.mit_del);
addParameter(p, 'gamma_dt', default.gamma_dt);
addParameter(p, 'well_radius', default.well_radius);
addParameter(p, 'im_dims', default.im_dims);
addParameter(p, 'in_size_range', default.in_size_range);
addParameter(p, 'in_size_min', default.in_size_min);
addParameter(p, 'to_image', default.to_image);
addParameter(p, 'to_file', default.to_file);
addParameter(p, 'output_file_path', default.output_file_path);
addParameter(p, 'sim_type', default.sim_type);
addParameter(p, 'calculate_fused', default.calculate_fused);
addParameter(p, 'calculate_tracked', default.calculate_tracked);
addParameter(p, 'track_max_dist', default.track_max_dist);
addParameter(p, 'well_idx', default.well_idx);
addParameter(p, 'seed_idx', default.seed_idx);
addParameter(p, 'dist_mat', default.dist_mat);
addParameter(p, 'shrink', default.shrink);
addParameter(p, 'filter_small', default.filter_small);
parse(p, ips);
ips = p.Results;

%% here we go
if length(ips.Ds) ~= length(ips.N_w) %if  there is not a number of wells given for every dose
    ips.N_w = ips.N_w*ones(size(ips.Ds)); %make the number of wells per dose equal to the scalar given in ips
end


%create "image" and "well" for all simulations because they don't change
if ~iscell(ips.seed_idx) && isnan(ips.seed_idx)
%if ~(ips.sim_type=="nonspatial")
    center = floor(ips.im_dims/2);                                                                %detect center of image
    y = 1:ips.im_dims(1);
    x = 1:ips.im_dims(2);
    [xx,yy] = meshgrid(x-center(2), y-center(1));
    R = sqrt(xx.^2+yy.^2);                                                                      %a distance-to-center matrix
    well_idx = find(R<ips.well_radius);                                                             %define parts of the grid (=well) where cells can sit
    ips.well_idx = well_idx(1:10:end);
end



%simulate for each dose
first = true;
for dd = 1:length(ips.Ds)
    disp("Dose:"+ips.Ds(dd));
    [raw_dose_sim, fused_dose_sim] = colony_size_sim_per_dose(ips,dd);
    if first
        first = false;
        raw_sim_output = raw_dose_sim;
        raw_sim_output2 = raw_sim_output;
        fused_sim_output = fused_dose_sim;
        fused_sim_output2 = fused_sim_output;
    else
        raw_sim_output2 = cat(1,raw_sim_output2, raw_dose_sim);
        for f = fields(raw_dose_sim)'
            fn = f{:};
            raw_sim_output.(fn) = cat(1,raw_sim_output.(fn), raw_dose_sim.(fn));
        end
        fused_sim_output2 = cat(1,fused_sim_output2, fused_dose_sim);
        for f = fields(fused_dose_sim)'
            fn = f{:};
            fused_sim_output.(fn) = cat(1,fused_sim_output.(fn), fused_dose_sim.(fn));
        end
    end
end






raw_sim_output.n_cells = colony_count_from_area(raw_sim_output.Area);
raw_sim_output.trackID = raw_sim_output.objID;
raw_sim_output = struct2table(raw_sim_output);



if ips.calculate_fused
    fused_sim_output.n_cells = colony_count_from_area(fused_sim_output.Area);
    fused_sim_output.trackID = fused_sim_output.objID;
    fused_sim_output = struct2table(fused_sim_output);
end

if ips.filter_small
    ix  = raw_sim_output.Area > ips.in_size_min;
    raw_sim_output = raw_sim_output(ix,:);
    if ips.calculate_fused
        ix  = fused_sim_output.Area > ips.in_size_min;
        fused_sim_output = fused_sim_output(ix,:);
    end
end

if ips.calculate_fused && ips.calculate_tracked && ips.sim_type ~= "nonspatial"
    tracked_sim_output = track_colonies(fused_sim_output, ips.track_max_dist, false);
    tracked_sim_output.N_cells = colony_count_from_area(tracked_sim_output.Area);
else
    tracked_sim_output =nan;
end

if ips.to_file
    [fol, fil, ex] = fileparts(ips.output_file_path);
    fil_raw = ['raw_',fil];
    fil_fused = ['fused_',fil];
    fil_tracked = ['tracked_',fil];
    writetable(raw_sim_output, fullfile(fol,[fil_raw, ex]));
    writetable(fused_sim_output, fullfile(fol,[fil_fused, ex]));
    writetable(tracked_sim_output, fullfile(fol,[fil_tracked, ex]));
end
end