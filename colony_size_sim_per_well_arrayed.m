function [raw_out, fused_out] = colony_size_sim_per_well_arrayed(ips,w,dd)
%This is the well-level simulation function, with biological parameters
%depending on dose ips.D(dd)


%check whether a fused version of simulation shall be simulated
calc_fused = ips.calculate_fused;


D = ips.Ds(dd);

%initialize outputs
%rng(1);


%Poisson distributed number of cells per well can be well-approximated by a
%normal distribution if N > 10
if ips.N > 10
    N = round(normrnd(ips.N,sqrt(ips.N)));
else
    N = poissrnd(round(ips.N));
end
fused = zeros(N,1, 'logical');


%put the cells in there
if iscell(ips.seed_idx)     %check whether default seeding positions exist
    seed_idx = ips.seed_idx{dd,w};
    %always choose the first N positions
    %positions = seed_idx(1:N,:);
    %or choose random positions from the predefined ones
    ix = randsample(1:length(seed_idx),N);
    positions = seed_idx(ix,:);
    seed_ys = positions(:,1);
    seed_xs = positions(:,2);
    %if seeding position exists, distance matrix should exist as well
else
    %rng(1);
    seed_idx = randsample(ips.well_idx, N);      %randomly choose a point of seeding for every cell
    [seed_ys, seed_xs] = ind2sub(ips.im_dims, seed_idx);
    positions = [seed_ys,seed_xs];
end
ObjIDs = (1:N)';
MergeIDs = ObjIDs;


if ips.to_image
    x = 1:ips.im_dims(1);
    y = 1:ips.im_dims(2);
    Rs = cell(length(positions),1);
    for i= 1:length(positions)
        pos = positions(i,:);
        [xx,yy] = meshgrid(x-pos(2), y-pos(1));
        Rs{i} = sqrt(xx.^2+yy.^2);
    end
end



%give them a size

%one option: draw from normal distribution
if ips.in_size_type == "norm"
    %rng(1);
    %A_seeds = ips.in_size_range(1) + (ips.in_size_range(2)-ips.in_size_range(1)).*rand(N,1);    %randomly choose a size for every cell (area)
    %A_seeds = normrnd(90,19.5,[N,1]); %randomly choose a size for every cell (area)
    A_seeds = normrnd(ips.in_size_range(1),ips.in_size_range(2),[N,1]); %randomly choose a size for every cell (area)
    smol = A_seeds<ips.in_size_min;
    while any(smol)
        n = sum(smol);
        A_seeds(smol)=normrnd(ips.in_size_range(1),ips.in_size_range(2),[n,1]);
        smol = A_seeds<ips.in_size_min;
    end
%other option: draw from a distribution we deliver with as input option
elseif ips.in_size_type == "distr"
    if class(ips.in_size_distr) == "prob.NormalDistribution"
        A_seeds = random(ips.in_size_distr, [N,1]);
        smol = A_seeds<ips.in_size_min;
        while any(smol)
            n = sum(smol);
            A_seeds(smol)=random(ips.in_size_distr,[n,1]);
            smol = A_seeds<ips.in_size_min;
        end
    else %ips.in_size_distr is already a set of values
        A_seeds = randsample(ips.in_size_distr, N, true);
    end
end




%give them growth rates, doubling times and mitotic delay
if isnan(ips.eta0) %the dose-dependent double-time way
    doub_time = ips.doub_time+ips.gamma_dt*D;
    eta = log(2)/doub_time;
else % the dose-dependent growth-rate way
    eta = ips.eta0 + ips.gamma_eta*D;
end

gr_sd = ips.gr_sd+ips.gr_sd_f*D; %relative width of the growth rate distribution

if ips.std_type == "rel"
    %gr_sd represents the relative standard deviation of growth rates
    growth_rates = normrnd(eta, eta*gr_sd,[N,1]);
else
    %gr_sd represents the absolute standard deviation of growth rates
    growth_rates = normrnd(eta, gr_sd,[N,1]);
end

%first option: every cell with growth rate < 0 is considered dead
quasi_dead = growth_rates<0;


%second option: cells with growth rate < 0 acquire a new growth rate until
%growth rate > 0
%smol = growth_rates<0;
%while any(smol)
%    n = sum(smol);
%    growth_rates(smol) = normrnd(eta, eta*gr_sd, [n,1]);
%    smol = growth_rates<0;
%end

%mitotic delay
md = D*ips.mit_del;

%apply plating efficiency
%rng(1);
deaths = rand(N,1) > ips.PE | quasi_dead;

%irradiate them
if D ~= 0
    %decide on whether a cell is rendered dead after a certain dose, using
    %the LQM
    deaths = deaths | (rand(N,1) > exp(-ips.alpha*D-ips.beta*D^2)); 
end

%try: make dead cells smaller
if ips.shrink
    %A_seeds(deaths) = A_seeds(deaths)./2;
    growth_rates(deaths) = -0.005;
else
     growth_rates(deaths) = 0;
end

r_seeds = sqrt(A_seeds/pi);




%names for output table columns
%varnames_r = {'Area','Centroid_1','Centroid_2','ObjID','Frame','Age','Fused', 'Dead', 'Growth_Rate'};
varnames_r = {'Area','Centroid_1','Centroid_2','objID','frame','age','fused', 'dead', 'growth_rate','wellID','dose'};

if calc_fused
    varnames_f = [varnames_r, 'n_seeds'];
    if ~iscell(ips.dist_mat)
        if isnan(ips.dist_mat)
            dists = pdist2(positions,positions,'euclidean');%a distance matrix between all seeds %maybe this can made faster if distance matrix is calculated beforehand (like positions)
        end
    else
        d = ips.dist_mat{dd,w};
        %as above: either choose first N 
        %dists = d(1:N,1:N);
        %or choose random N
        dists = d(ix,ix);
    end
    glob_fM = zeros(size(dists));
end


%let them grow
%the ugly splitting of all variables into their own vectors happens for
%performance issues. If for every timepoint a table is created and
%concatenated with the rest of the results, the time-consuming table
%concatenation function is called too often. Hence the ugliness.

first = true;
%iterate over time points
for j = 1:length(ips.times)
    t = ips.times(j);
    %disp(t);
    %calculate colony radii
    rs = r_seeds.*sqrt(exp(growth_rates*max(0,t-md)));
    %calculate colony areas
    areas = pi*rs.^2;
    %test whether sum(radius1,radius2) > dists(1,2) -> overlap == fused
    if calc_fused
        
        rr = rs+rs';%pairwise sum of radii
        fM = rr>dists -glob_fM;%matrix containing ones where two colonies are newly fused
        fM = triu(fM,2);%only upper triangle necessary
        glob_fM = glob_fM | fM;
        [xind,yind] = find(fM); %pairs of fused colonies
        for f = 1:length(xind)
            xi = xind(f); %the index of the first colony
            fused(xi) = 1;
            yi = yind(f); %the index of the second colony
            if fused(yi) %if the second colony is already in a fusion state
                ix =  MergeIDs(xi); %we extract the merged ID of the first colony
                MergeIDs(MergeIDs==ix) = MergeIDs(yi); %and set all merge IDs of all potential colonies which are fused with the first one to the merge ID of the second one
            else %if the second is not yet fused
                %ix = MergeIDs(yi); %we extract the merge ID of the second
                %MergeIDs(MergeIDs==ix) = MergeIDs(xi); %and set all merge IDs to the first colonies merge ID (should only be the second colony who has this merge ID since it's not merged yet
                MergeIDs(yi) = MergeIDs(xi); %we add the second colony to the group of colonies with mergeID identical to first colony
                fused(yi) = 1;
            end
        end
    else
        fused = nan*areas;
    end
    frames = j*ones(N,1);
    ages = t*ones(N,1);
    wellIDs =w*ones(N,1);
    doses = D*ones(N,1);
    if calc_fused
        G = findgroups(MergeIDs);
        IDs = unique(G);
        N_seeds = accumarray(G,ones(size(G)));
        Areas = accumarray(G,areas);
        Growth_Rates = accumarray(G,growth_rates)./N_seeds;
        Centroid_1 = accumarray(G,seed_ys)./N_seeds;
        Centroid_2 = accumarray(G,seed_xs)./N_seeds;
        Dead = accumarray(G,deaths)./N_seeds;
        Fused = N_seeds > 1;
        Frames = j*ones(length(IDs),1);
        Ages = t*ones(length(IDs),1);
        WellIDs =w*ones(length(IDs),1);
        Doses = D*ones(length(IDs),1);
    end
    if first
        raw_areas = areas;
        raw_seed_ys = seed_ys;
        raw_seed_xs = seed_xs;
        raw_obj_ids = ObjIDs;
        raw_frames = frames;
        raw_ages = ages;
        raw_wellIDs = wellIDs;
        raw_doses = doses;
        raw_fused = fused;
        raw_deaths = deaths;
        raw_growth_rates = growth_rates;
        if calc_fused
            fused_areas = Areas;
            fused_seed_ys = Centroid_1;
            fused_seed_xs = Centroid_2;
            fused_IDs = IDs;
            fused_frames = Frames;
            fused_ages = Ages;
            fused_fused = Fused;
            fused_dead = Dead;
            fused_growth_rates = Growth_Rates;
            fused_wellIDs = WellIDs;
            fused_doses = Doses;
            fused_nseeds = N_seeds;
        end
        first = false;
    else
        raw_areas = cat(1,raw_areas,areas);
        raw_seed_ys = cat(1,raw_seed_ys,seed_ys);
        raw_seed_xs = cat(1,raw_seed_xs,seed_xs);
        raw_obj_ids = cat(1,raw_obj_ids,ObjIDs);
        raw_frames = cat(1,raw_frames,frames);
        raw_ages = cat(1,raw_ages,ages);
        raw_wellIDs = cat(1,raw_wellIDs,wellIDs);
        raw_doses = cat(1,raw_doses,doses);
        raw_fused = cat(1,raw_fused,fused);
        raw_deaths = cat(1,raw_deaths,deaths);
        raw_growth_rates = cat(1,raw_growth_rates,growth_rates);
        if calc_fused
            fused_areas = cat(1,fused_areas,Areas);
            fused_seed_ys = cat(1,fused_seed_ys,Centroid_1);
            fused_seed_xs = cat(1,fused_seed_xs,Centroid_2);
            fused_IDs = cat(1,fused_IDs,IDs);
            fused_frames = cat(1,fused_frames,Frames);
            fused_ages = cat(1,fused_ages,Ages);
            fused_fused = cat(1,fused_fused,Fused);
            fused_dead = cat(1,fused_dead,Dead);
            fused_growth_rates = cat(1,fused_growth_rates,Growth_Rates);
            fused_wellIDs = cat(1,fused_wellIDs,WellIDs);
            fused_doses = cat(1,fused_doses,Doses);
            fused_nseeds = cat(1,fused_nseeds,N_seeds);
        end
    end
    if ips.to_image
        mask = zeros(ips.im_dims);
        for i= 1:length(positions)
            idx = Rs{i}<rs(i);
            mask(idx) = 1;
        end
    end
end


raw_out = struct();
c = {raw_areas,raw_seed_ys,raw_seed_xs,raw_obj_ids,raw_frames,raw_ages,raw_fused,raw_deaths,raw_growth_rates,raw_wellIDs,raw_doses};
for i = 1:length(c)
    raw_out(:).(varnames_r{i})=c{i};
end



if calc_fused
    c = {fused_areas,fused_seed_ys,fused_seed_xs,fused_IDs,fused_frames,fused_ages,fused_fused,fused_dead,fused_growth_rates,fused_wellIDs,fused_doses,fused_nseeds};
    for i = 1:length(c)
        fused_out.(varnames_f{i})=c{i};
    end
else
    fused_out = struct();
end
