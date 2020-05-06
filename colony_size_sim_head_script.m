%% This is the head script from which we start a simulation. Initially, input parameters are defined, then a simulation is performed



%create input structure for model function
inputs = [];
%% experimental setup parameters
%number of objects per well
inputs.N = 385;
%number of wells per dose (can be a scalar or a vector of scalars, one for
%each dose)
inputs.N_w = 12;
%vector of time points (arbitrary units, but need to be coherent with other
%values, usually hours)
inputs.times = 0:3:219;
%vector of doses
inputs.Ds =0:7;

%% "biological" parameters
%plating efficiency: which fraction of colonies stays alive and growing
%without irradiation (growth rate != 0)
inputs.PE = 0.8;
%an optional parameter inspired by the LQM: linear exponential 
%dose-dependent decrease in survival probability 
inputs.alpha = 0;
%an optional parameter inspired by the LQM: quadratic exponential 
%dose-dependent decrease in survival probability 
inputs.beta = 0;
%one way to set growth rates in absence of irradiation: using doubling time 
%(inverse of growth rate) given in hours [h]
inputs.doub_time = 35;
%how does doubling time change with dose? [hGy⁻¹]
inputs.gamma_dt = 1;
%another way: using the growth rate. Both .doub_time and .eta0 indicate the
%mean value in absence of irradiation - given in inverse hours [h⁻¹]
inputs.eta0 = 0.02;
%how does the mean growth rate change per Gy dose? [h⁻¹Gy⁻¹]
inputs.gamma_eta = -0.001;

%an option indicating whether the growth rate standard deviation .gr_sd 
%is to be interpreted as a relative value ("rel") or as an absolute value 
%("abs")
inputs.std_type = "rel";
% dose-independent standard deviation of growth rate. Dimensionless in 
%[0,1] if .std_type is "rel", given in [h⁻¹] if .std_type is "abs"
inputs.gr_sd = 0.32;
%how does the growth rate standard deviation change with dose [Gy⁻¹] or
%[h⁻¹Gy⁻¹] depending on .std_type
inputs.gr_sd_f = 0.015;

%mitotic delay: how much is onset of growth delayed per dose [h⁻¹Gy⁻¹]
inputs.mit_del = 1;

%an option to change the source of initial sizes. If "norm" is given, the
%values in .in_size_range are inpterpreted as mean and std of a normal
%distribution from which to draw sizes. If "distr" is given, sizes are
%drawn from the .in_size_distr, which can eeither be a distribution object
%or a vector of sizes (e.g. from experimental data) from which to draw.
inputs.in_size_type = "norm";
%mean and std of a potential distribution from which to draw initial sizes.
%Only has an effect if .in_size_type is "norm"
inputs.in_size_range = [90,20];
%a distribution or a vector of size values from which to draw. Only has an
%effect if .in_size_type is "distr"
inputs.in_size_distr = makedist('Normal', 'mu', 90,'sigma',20);
%a lower threshold on initial sizes. If sizes below this threshold are 
%drawn as initial sizes, they are redrawn until they are larger than 
%.in_size_min
inputs.in_size_min = 55;
%a boolean flag which indicates whether dead colonies should shrink (with 
%a hardcoded rate somewhere in the function code). This is to emulate the
%observation of shrinking colonies in the exp. data
inputs.shrink = false;
%a boolean value indicating whether objects smaller than .in_size_min
%should be removed from the data (emulates the filtering step happening on
%the exp. data). Only makes sense if .shrink is true, because otherwise all
%objects stay >= .in_size_min anyway
inputs.filter_small = true;



%not fully developed: a hazard for single cells to cease growth. apha_h is 
%the initial hazard for a single cell to cease growth, mu is therate with 
%which this hazard changes. This represents one way to implement 
%time-dependent changes in effective growth rates. Only has effects if the
%sim_type is "gradual"
inputs.alpha_h = 0.01;
inputs.mu = -log(1/2)/inputs.doub_time; %<-log(1/2)/inputs.doub_time> means the hazard halves for every doubling time


%% technical/output parameters
% sizes in [x,y]-dimensions of the original images with which we want to
% compare the simulation results. Given in pixels
inputs.im_dims = [2240,2240];
%radius of the circular well in pixels. Should be equivalent to the radius
%of wells in the original data and less than half of min(im_dims)
inputs.well_radius = 1000;
%boolean flag indicating whether images should be created from the
%simulation
inputs.to_image = false;
%boolean flag indicating whether the result should be stored as a csv-file
inputs.to_file = false;

%the type of simulation. At the moment only "initial" is valid. In this
%type, growth rates and death state are initially set and remain unchanged 
%throughout the simulation. Other possible types are "nonspatial", where 
%space is not modeled and hence no interaction (fusion) between colonies 
%happen, or "gradual", where individual cells have a gradually changing 
%probability to stop growing (based on .alpha_h and .mu)
inputs.sim_type = "initial";
%boolean flag indicating whether an additional output containing fused
%colonies should be calculated
inputs.calculate_fused = true;
%boolean flag indicating whether an additional output containing tracked
%colonies should be calculated. This emulates the tracking happening on
%experimental data to compare the influence of tracking on the resulting
%size distributions
inputs.calculate_tracked = false;
% the maximum distance between two objects to allow mapping during
%tracking. Only valid if .calculate_tracked is true;
inputs.track_max_dist = 10;


%an optional path where the resulting datatable(s) should be stored.
inputs.output_file_path = '/path/to/output/folder/and/file.csv';

%optional: a vector of indices indication which "pixels" in a matrix are
%wells. This is to avoid multiple recalcuations of the same values in cases
%of high-throughput simulation runs (e.g. for parameter fitting)
inputs.well_idx = nan;
%optional: the same concept for indices at which colonies are initially
%seeded. 
inputs.seed_idx = nan;
%If seed indices are given, the distance matrix between those seeds can be
%precalculated as well. All three options are used in the context of
%repetitive calls to the simulation function e.g. during parameter fitting
inputs.dist_mat = nan;



%% Here the magic happens: 
%depending on above parameters, this main function
%runs a simulation and outputs up to three data tables. If
%inputs.calculate_fused is false, the second output is NaN, if 
%inputs.calculate_tracked is false, the third output is NaN

[sim_output_raw, sim_output_fused, sim_output_tracked] = colony_size_sim_main(inputs);


%Further functions post-process the resulting data table(s) to extract
%dose-, well- and time-dependent statistics as well as to visualize the
%distributions of sizes resulting from the simulation. We will deal with
%this later

