function [datasim] = sim_single_drug_effects(gen_DM_param)
%SIM_DRUG_EFFECTS simulates the growth and cell cycle distribution of a
%cell population under assumptions of dose-dependent effects of a drug on
%cell cycle and apoptotic rates

%% Load the Drug Model (DM) paramters
[kG1S kSG2 kG2M kMG1 papG1 papS papG2 papM Vcc Kcc ncc Vap Kap nap] = gen_DM_param();

%% Define Cell cycle and Apoptosis Model 
cc_names={'G1' 'S' 'G2' 'M'};
%define the untreated transition matrix 
[TM, CTM, AT] = create_TM(kG1S,kSG2,kG2M,kMG1,papG1,papS,papG2,papM);

%% Define the dose-dependent effects of a drug on cell cycle and apoptosis 
%define the logspaced drug dilution used (+ DMSO control)
minDrug=-3;
maxDrug=1;
ndoses=9;
drug=[0 logspace(minDrug,maxDrug,ndoses)];

%calculate dose-dependent drug effect on each cell cycle rate
npoints=length(drug); nparam=length(Vcc);
drug_mod_cc=arrayfun(@rev_hill_eq,repmat(drug,nparam,1),repmat(Vcc,1,npoints),repmat(Kcc,1,npoints),repmat(ncc,1,npoints));

%calculate dose-dependent drug effect on each cell cycle rate
npoints=length(drug); nparam=length(Vap);
drug_mod_ap=arrayfun(@hill_eq,repmat(drug,nparam,1),repmat(Vap,1,npoints),repmat(Kap,1,npoints),repmat(nap,1,npoints));

%% generate the data structure with information about the  model (updated in the relavant sections during model simulation)
datasim.param.cc.k=[kG1S kSG2 kG2M kMG1];
datasim.param.cc.Vcc=Vcc;
datasim.param.cc.Kcc=Kcc;
datasim.param.cc.ncc=ncc;
datasim.param.ap.k=[papG1 papS papG2 papM];
datasim.param.ap.Vap=Vap;
datasim.param.ap.Kap=Kap;
datasim.param.ap.nap=nap;
datasim.param.drug=drug;
datasim.param.phase_names=cc_names;

%% plot dose-dependent drug effects on cell cycle and apoptosis rates as defined
figure('Name','Steady state behaviour');
set(gcf,'Color','White');
%set parameters for plotting 
lw=2; %line width
colors=[0 0 1; 0.8 0.8 0; 1 0 0; 0 0.8 0];%colors
%plot dose-response curves for cell cycle phases on wide drug range
subplot(2,2,1);
n_drug_sampling=100;
drug_sampling=[0 logspace(minDrug-1,maxDrug+1,n_drug_sampling)];
drug_sampling=union(drug_sampling,drug); %add all the measured drugs to the sampled drugs (if missed out by the step size chosen)
%calculate response-curves with higher density points for plotting
npoints=length(drug_sampling); nparam=length(Vcc);
drug_sampling_mod_cc=arrayfun(@rev_hill_eq,repmat(drug_sampling,nparam,1),repmat(Vcc,1,npoints),repmat(Kcc,1,npoints),repmat(ncc,1,npoints));
%calculate response-curves with higher density points for plotting
npoints=length(drug_sampling); nparam=length(Vap);
drug_sampling_mod_ap=arrayfun(@hill_eq,repmat(drug_sampling,nparam,1),repmat(Vap,1,npoints),repmat(Kap,1,npoints),repmat(nap,1,npoints));
%for each cell cycle phase plot the response
for i=1:size(drug_sampling_mod_cc,1)
    plot(drug_sampling,drug_sampling_mod_cc(i,:),'LineWidth',lw,'Color',colors(i,:));
    hold on;
end
%plot dose-response curves for cell cycle phases on sampled drug conc.
for i=1:size(drug_mod_cc,1)
    plot(drug,drug_mod_cc(i,:),'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
    hold on;
end
%plot legend and axis
legend(cc_names);
set(gca,'XScale','log');
ylim([0 max(drug_sampling_mod_cc(:))+0.1]);
xlabel('Drug (\muM)');
ylabel('Rate (h^{-1})');
title('Drug effect on cell cycle progression');

%plot dose-response curves for apoptosis wide drug range
subplot(2,2,2);
%for each cell cycle phase plot the response
for i=1:size(drug_sampling_mod_ap,1)
    plot(drug_sampling,drug_sampling_mod_ap(i,:),'LineWidth',lw,'Color',colors(i,:));
    hold on;
end
%plot dose-response curves for cell cycle phases on sampled drug conc.
for i=1:size(drug_mod_ap,1)
    plot(drug,drug_mod_ap(i,:),'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
    hold on;
end
%plot legend and axis
legend(cc_names);
set(gca,'XScale','log');
ylim([0 max(drug_sampling_mod_ap(:))+0.1]);
xlabel('Drug (\muM)');
ylabel('Rate (h^{-1})');
title('Drug effect on cell-cycle dependent apoptosis');

%% calculate and plot the steady state dose-dependent growth rates and cell cycle distributions
%calculate steady state growth rates and distributions for all sampled drug concentrations (apply @create_TM and @calculate_steady_state)
[TM_drug_sampling] = arrayfun(@create_TM,drug_sampling_mod_cc(1,:),drug_sampling_mod_cc(2,:),drug_sampling_mod_cc(3,:),drug_sampling_mod_cc(4,:),drug_sampling_mod_ap(1,:),drug_sampling_mod_ap(2,:),drug_sampling_mod_ap(3,:),drug_sampling_mod_ap(4,:),'UniformOutput', false);
[ss_growth_rate_drug_sampling, ss_dist_cell_cycle_drug_sampling]=cellfun(@calculate_steady_state,TM_drug_sampling,'UniformOutput', false);
%conver growth_rate from cell array of scalars to a vector and dist_cell_cycle from a cell array of unidemnsional vectors to a matrix
ss_growth_rate_drug_sampling=cell2mat(ss_growth_rate_drug_sampling);
ss_dist_cell_cycle_drug_sampling=cell2mat(ss_dist_cell_cycle_drug_sampling);

%calculate steady state growth rates and distributions for measured drug concentrations (apply @create_TM and @calculate_steady_state)
[TM_drug] = arrayfun(@create_TM,drug_mod_cc(1,:),drug_mod_cc(2,:),drug_mod_cc(3,:),drug_mod_cc(4,:),drug_mod_ap(1,:),drug_mod_ap(2,:),drug_mod_ap(3,:),drug_mod_ap(4,:),'UniformOutput', false);
[ss_growth_rate_drug, ss_dist_cell_cycle_drug]=cellfun(@calculate_steady_state,TM_drug,'UniformOutput', false);
%conver growth_rate from cell array of scalars to a vector and dist_cell_cycle from a cell array of unidemnsional vectors to a matrix
ss_growth_rate_drug=cell2mat(ss_growth_rate_drug);
ss_dist_cell_cycle_drug=cell2mat(ss_dist_cell_cycle_drug);

%plot steady state growth rates
subplot(2,2,3);
plot(drug,ss_growth_rate_drug./ss_growth_rate_drug(1),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on;
for i=1:size(ss_dist_cell_cycle_drug_sampling,1)
    plot(drug,ss_dist_cell_cycle_drug(i,:) .* ss_growth_rate_drug./ss_growth_rate_drug(1),'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
    hold on;
end

%plot the steady state growth rates and scaled fractions for the different drug points
plot(drug_sampling,ss_growth_rate_drug_sampling./ss_growth_rate_drug_sampling(1),'-','LineWidth',lw,'Color','k');
hold on;
for i=1:size(ss_dist_cell_cycle_drug_sampling,1)
    plot(drug_sampling,ss_dist_cell_cycle_drug_sampling(i,:).*ss_growth_rate_drug_sampling./ss_growth_rate_drug_sampling(1),'-','LineWidth',lw,'Color',colors(i,:));
    hold on;
end
set(gca,'XScale','log');
xlabel('Drug (\muM)');
ylabel('\mu_{DMSO}/\mu * fraction population');
title('Drug effect on steady-state growth rate and cell cycle fractions');
legend({'Total ' cc_names{:}}); 
%plot 0 line

%plot steady state cell cycle distributions
subplot(2,2,4);
%plot steady state distribution for drug sampling
for i=1:size(ss_dist_cell_cycle_drug_sampling,1)
    plot(drug_sampling,ss_dist_cell_cycle_drug_sampling(i,:),'-','LineWidth',lw,'Color',colors(i,:));
    hold on;
end
%plot steady state distribution for drug concentrations measured
for i=1:size(ss_dist_cell_cycle_drug,1)
    plot(drug,ss_dist_cell_cycle_drug(i,:),'o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:));
    hold on;
end
set(gca,'XScale','log');
xlabel('Drug (\muM)');
ylabel('Fraction of total population');
title('Drug effect on steady-state cell cycle fractions (normalized)');
legend(cc_names);
ylim([0 1.1]);

%% perform ODE simulation to assess growth dynamics and cell cycle distribution starting from untreated steady state for few observation time-points (12h, etc..)
ncells=2000; %size of untreated population at start 
time_mes=[12 24 48 72 96]; %define the observation times
extratime=4; %defines the extra time to prolong the simulation 
tspan = 0:0.1:(max(time_mes)+extratime); %define the ODE time points results 
tspan=union(tspan, time_mes); %make sure time_mes is in tspan

%calculate the untreated steady state distribution 
[~, ss_dist_cell_cycle_untreated]=calculate_steady_state(TM);
ncell_ss_dist_untreated=ss_dist_cell_cycle_untreated .* ncells;

%ricalculate rates according to the predefiend dose-dependent drug effects (might not be needed but redone to ensure it works even if code moves around)
%calculate dose-dependent drug effect on transitions in each cell cycle rate
npoints=length(drug); nparam=length(Vcc);
drug_mod_cc=arrayfun(@rev_hill_eq,repmat(drug,nparam,1),repmat(Vcc,1,npoints),repmat(Kcc,1,npoints),repmat(ncc,1,npoints));
%calculate dose-dependent drug effect on apoptosis rates for each cell cycle phase
npoints=length(drug); nparam=length(Vap);
drug_mod_ap=arrayfun(@hill_eq,repmat(drug,nparam,1),repmat(Vap,1,npoints),repmat(Kap,1,npoints),repmat(nap,1,npoints));

%for each drug concentration (including d=0, the DMSO control) simulate growth dynamics
for i=1:length(drug_sampling)
    %create the TM matrix for that particular drug dose
    %TM_drug_dose_tmp=create_TM(drug_mod_cc(1,i),drug_mod_cc(2,i),drug_mod_cc(3,i),drug_mod_cc(4,i),drug_mod_ap(1,i),drug_mod_ap(2,i),drug_mod_ap(3,i),drug_mod_ap(4,i));
    TM_drug_dose_tmp=create_TM(drug_sampling_mod_cc(1,i),drug_sampling_mod_cc(2,i),drug_sampling_mod_cc(3,i),drug_sampling_mod_cc(4,i),drug_sampling_mod_ap(1,i),drug_sampling_mod_ap(2,i),drug_sampling_mod_ap(3,i),drug_sampling_mod_ap(4,i));
    %solve the ODE system 
    %[t,y] = ode45(@(t,y)cell_cycle_ODE(t,y,TM_drug_dose_tmp), tspan, ncell_ss_dist_untreated');
    [t,y] = ode45(@(t,y)cell_cycle_ODE(t,y,TM_drug_dose_tmp), tspan, ncell_ss_dist_untreated');
    time_ODE(i,:)=t;
    cellcount_ODE(i,:)=sum(y,2)';
    celldist_ODE(i,:,1:4)=y;   
end    

%update the datasim structure
datasim.ODE.time_mes=time_mes;
datasim.ODE.ncells=ncells;
datasim.ODE.time_ODE=time_ODE;
datasim.ODE.cellcount_ODE=cellcount_ODE;
datasim.ODE.celldist_ODE=celldist_ODE;

%for each drug concentration, plot the corresponding cell growth from the ODE results
figure('Name','ODE dynamics');
set(gcf,'Color','White');
subplot(2,6,[1 7]);
for i=1:length(drug)
    %find the drug in the drug_sampling 
    idrug=find(drug_sampling==drug(i));
    %calculate the color for that drug conc
    col_drug_temp=[(1./length(drug)).*i 0 0];
    %plot the cell count dynamics
    plot(time_ODE(idrug,:),cellcount_ODE(idrug,:),'-','LineWidth',lw,'Color',col_drug_temp); 
    hold on;
end
legend(strread(num2str(drug),'%s'));
for i=1:length(drug)
    %find the drug in the drug_sampling 
    idrug=find(drug_sampling==drug(i));
    %calculate the color for that drug conc
    col_drug_temp=[(1./length(drug)).*i 0 0];
    %plot dots that mark the measurement time points
    [~,itime_mes,~]=intersect(time_ODE(idrug,:),time_mes);
    plot(time_ODE(idrug,itime_mes),cellcount_ODE(idrug,itime_mes),'o','MarkerFaceColor',col_drug_temp,'MarkerEdgeColor',col_drug_temp); 
    hold on;
end
set(gca,'YScale','log');
xlabel('Time (h)');
ylabel('Cell count');
title('Time-course cell growth');

%for each time-point, plot the dose-dependent cell viability (1st subplot) and the cell cycle relative distribution (2nd plot)
for i=1:length(time_mes)
    %find corresponding time vector
    itime=find(time_ODE(1,:)==time_mes(i));
    %1st plot: relative cell viability
    subplot(2,6,i+1);
    rel_cell_viab=cellcount_ODE(:,itime)./cellcount_ODE(1,itime);
    plot(drug_sampling,rel_cell_viab,'-','LineWidth',lw,'Color','k'); 
    hold on;
    xlim([min(drug_sampling) max(drug_sampling)]);
    ylim([0 1.1]);
    set(gca,'XScale','log');
    xlabel('Drug (\muM)');
    ylabel('Population DMSO\\drug');
    title(sprintf('%dh',time_mes(i)));
    legend('Total'); 
    %plot dots for measurement 
    [~,idrug,~]=intersect(drug_sampling,drug);
    plot(drug_sampling(idrug),rel_cell_viab(idrug),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
    %2nd plot: cell cycle distribution 
    subplot(2,6,6+i+1);
    for j=1:size(celldist_ODE,3)
        plot(drug_sampling,celldist_ODE(:,itime,j)./cellcount_ODE(:,itime),'-','LineWidth',lw,'Color',colors(j,:));
        hold on;
    end
    xlim([min(drug_sampling) max(drug_sampling)]);
    ylim([0 1.1]);
    set(gca,'XScale','log');
    xlabel('Drug (\muM)');
    ylabel('Fraction of population');
    legend(cc_names,'Location','northwest');
    %plot dots for measurement
    [~,idrug,~]=intersect(drug_sampling,drug);
    for j=1:size(celldist_ODE,3)
        plot(drug,celldist_ODE(idrug,itime,j)./cellcount_ODE(idrug,itime),'o','MarkerFaceColor',colors(j,:),'MarkerEdgeColor',colors(j,:));
        hold on;
    end
end    

%% simulate each cell dynamics individually and obtain DNA content ssimulation for G1, S, G2, M and apoptotic fraction
%define the limits of the cell cycle phases:
%            G1   S   G2   M
phase_lim=[0   1    2    3   4];
%generate vector that counts apoptotic events
apo_rate=zeros(1,4); %G1 S G2 M
apo_done=0; %counts cells that died
%define DNA content based on cell cycle progression G1=2N, S=from 2 to 4 linear, G2/M=4N
phase_value_DNA=[0 1 2 3 4];
phase_DNA=[2 2 4 4 4];
%define amount of noise to add to DNA measurement
DNA_sigma=normrnd(0.2,0.05); %for normrnd(mu,sigma) or lognrnd(mu,sigma)
%define the pRB content based on cell cycle progression (late G1 to M is pRB)
phase_value_RB=[0 0.7 1 2 3 3.7 4];
phase_RB=[1 1 3 3 3 3 1]; 
%define amount of noise to add to DNA measurement
RB_sigma=normrnd(0.2,0.05); %f
%initiate single cell trackers
single_cells_init=[];
%error model 
%dist_fun=@normrnd;  %normal
dist_fun=@lognrnd; %log normal
%set the bins for DNA 
bins_DNA=0:0.05:6;

%update informations on data sim
datasim.singlecells.time_mes=time_mes;
datasim.singlecells.drug=drug;
datasim.singlecells.ncells=ncells;
datasim.singlecells.DNA_sigma=DNA_sigma;
datasim.singlecells.RB_sigma=RB_sigma;

%populate the single_cells according to the steady state distribution (uniformely spaced in the cell cycle phases)
[~, ss_dist_cell_cycle_DMSP]=calculate_steady_state(TM);
for i=1:length(ss_dist_cell_cycle_DMSP)
    %calculate number of cells in that phase
    ncells_phase=round(ss_dist_cell_cycle_DMSP(i) .* ncells);
    %calculate the uniform spacing of all cells in the phase
    size_phase=(phase_lim(i+1)-phase_lim(i))./ncells_phase;
    cell_phase_pos=phase_lim(i):size_phase:phase_lim(i+1);
    %add the single cells to the vector of single cells
    single_cells_init=[single_cells_init cell_phase_pos];
end

%for each drug concentration simulate single cell progression 
figure('Name','Single cell dynamics');
set(gcf,'Color','White');
%decide if 1 dimensional (DNA) or 2 (DNA-pRB) plot or 3 (DNA-pRB-cell cycle color) plot
dim_plot=1;
for k=1:length(drug)
    %reorganize the rate for the single cell simulation
    kT=[drug_mod_cc(1,k) drug_mod_cc(2,k) drug_mod_cc(3,k) drug_mod_cc(4,k)];
    kA=[drug_mod_ap(1,k) drug_mod_ap(2,k) drug_mod_ap(3,k) drug_mod_ap(4,k)];
    %define and update timesteps depending on rates if necessary
    size_step=0.05;
    %size_step=0.01./min([size_step min(min(abs(kT(kT~=0)))) min(min(abs(kA(kA~=0))))]);
    %size_step=min([size_step 1]);
    %define the time vector to simulate stepwise
    time_steps=0:size_step:(max(time_mes)+extratime);
    %add the timepoint to plot 
    time_steps=union(time_steps, time_mes);
    %simulate single cell progression (perform division and cell death)
    single_cells=single_cells_init;
    cell_pop=[];
    for i=1:length(time_steps)-1
        step=time_steps(i+1)-time_steps(i);
        %use CTM (only cell cycle rates) to calculate actual steps in the phase
        %and AT (only apop rate) to calculate the number of cells that need to die
        for j=1:length(phase_lim)-1
            %cell cycle update
            iphase=find((phase_lim(j)<=single_cells)  & (single_cells <= phase_lim(j+1)));
            single_cells(iphase)=single_cells(iphase) + kT(j) .* step;
            %apoptosis update 
            apo_rate(j)= apo_rate(j) + length(iphase) .* kA(j) .* step; 
            nkill=floor(apo_rate(j));
            %pick nkill random cells in the single cell vector (same phase)
            %ikill = randperm(length(iphase));
            %ikill= ikill(1:min(nkill,length(iphase)));
            ikill =rand_unique(nkill,1,length(iphase));
            %delete the selected cells
            single_cells(iphase(ikill))=[];
            %update kill counter and the rate at which killing should happen
            apo_done=apo_done+nkill;
            apo_rate(j)=apo_rate(j)-nkill;
        end    
        %perform cell divisions for cells that passed M phase
        i_division=single_cells>phase_lim(end);
        %reset them to G1
        single_cells(i_division)=single_cells(i_division)-phase_lim(end);
        %create and identical new copy for each cell
        single_cells=[single_cells single_cells(i_division)];
        cell_pop(end+1)=length(single_cells);
        
        %print progression in simulations
        if (mod(i,1000)==0)
           fprintf('Drug %d(of %d): step %d(of %d), step size %f\n',k,length(drug),i,length(time_steps),size_step);
        end
        %plot cell distributions at measurement times 
        i_time_mes=find(time_mes==time_steps(i));
        if (~isempty(i_time_mes))
            subplot(length(time_mes)+1,length(drug),sub2ind([length(drug) length(time_mes)+1],k,i_time_mes));
            %generate DNA overall distribution
            DNA_dist=dist_fun(interp1(phase_value_DNA,phase_DNA,single_cells),DNA_sigma);
            DNA_dist=log(DNA_dist);
            %generate pRB overall distribution
            RB_dist=dist_fun(interp1(phase_value_RB,phase_RB,single_cells),RB_sigma); 
            RB_dist=log(RB_dist);
            %save in dataset
            datasim.singlecells.DNA_dist{k,i_time_mes}=DNA_dist;
            datasim.singlecells.RB_dist{k,i_time_mes}=RB_dist;
            %plot one-dimensional DNA  or bi-dimensional DNA-pRb distributions
            if (dim_plot==1)
                %plot one dimensional DNA-pRB   
                histogram(DNA_dist,bins_DNA,'FaceColor','k','EdgeColor','none');
                hold on;
            elseif (dim_plot==2)
                %plot two dimensional DNA-pRB   
                dscatter(DNA_dist',RB_dist');
                hold on;
                %plot time on the y-axis if is the first column
                if (k==1)
                    ylabel(sprintf('pRB(%dh)(log)',time_mes(i_time_mes)));
                end 
            elseif (dim_plot==3)
                %do nothing, colors are generated for cell cycle specific phases
            end    
            %plot the different cell phase distributions 
            for j=1:length(phase_lim)-1
                %identify cells in different cell cycle phases
                iphase=find((phase_lim(j)<=single_cells)  & (single_cells <= phase_lim(j+1)));
                DNA_dist_phase=DNA_dist(iphase);
                RB_dist_phase=RB_dist(iphase);
                %save the data
                datasim.singlecells.phase_index{k,i_time_mes,j}=iphase;
                datasim.singlecells.DNA_dist_phase{k,i_time_mes,j}=DNA_dist_phase;
                %DNA_dist=dist_fun(interp1(phase_value_DNA,phase_DNA,single_cells(iphase)),DNA_sigma);
                %DNA_dist=log(DNA_dist);
                %RB_dist=dist_fun(interp1(phase_value_RB,phase_RB,single_cells(iphase)),RB_sigma);
                %RB_dist=log(RB_dist);                
                %plot one-dimensional DNA  or bi-dimensional DNA-pRb distributions
                if (dim_plot==1)
                   %plot one dimensional DNA
                   histogram(DNA_dist_phase,bins_DNA,'FaceAlpha',0.5,'EdgeAlpha',0.5,'FaceColor',colors(j,1:3),'EdgeColor','none');
                   %plot time on the y-axis if is the first column
                   if (k==1)
                       ylabel(sprintf('%dh',time_mes(i_time_mes)));
                   end   
                elseif (dim_plot==2)
                   %do nothing, overall plots is generated from whole cell data above
                elseif (dim_plot==3)
                   %plot two-dimensional with colors depending on cell cycle phase 
                   scatter(DNA_dist_phase',RB_dist_phase',2,colors(j,1:3));
                   hold on;
                   if (k==1)
                    ylabel(sprintf('pRB(%dh)(log)',time_mes(i_time_mes)));
                   end 
                end    
            end
            %set DNA axis limit
            xlim([1 5]);
            %plot DNA content xlabel only if last raw
            if (i_time_mes==length(time_mes))
                xlabel('DNA(log)');
            end
            %plot drug dose only on first row
            if (i_time_mes==1)
                title(sprintf('%suM',num2str(drug(k))));
            end    
        end    
    end
    %plot the cell count behaviour at the end of the plot
    subplot(length(time_mes)+1,length(drug),sub2ind([length(drug) length(time_mes)+1],k,length(time_mes)+1));
    plot(time_steps(1:end-1),cell_pop,'k-','LineWidth',1.5);
    hold on;
    %plot dots at the measurement points 
    [~,iinters,~]=intersect(time_steps,time_mes);
    plot(time_steps(iinters),cell_pop(iinters),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',4);
    set(gca,'YScale','log');
    xlabel('Time(h)'); 
    %plot label only external plots
    if (k==1) 
        ylabel('Cells');
    end
end


end

%function that defines the ODE system for the celly cycle accepts the transition matrix (TM)
function [dydt]=cell_cycle_ODE(t,y,TM)
   %calculate rate of change
   dydt=TM * y;
end

%function that returns n non repeating random number in a range 
function [r]=rand_unique(n,min,max)
  p = min:max;
  m = length(p);
  r = zeros(1,n);
  for k = m:-1:m-n+1
   q = ceil(k*rand);
   r((m+1)-k) = p(q); % <-- n random integers go to 'r'
   p(q) = p(k);
  end
end 
