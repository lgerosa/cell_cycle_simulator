function [] = genSynData()
%GENSYNDATA generates synthetic data for a variety of Drug Models


%% define the type of drug modes that want to simulate
cc_names={'G1' 'S' 'G2' 'M'};
type_DM={};
for i=1:4
    %first define inhibition of only one phase
    type_DM_tmp=[0 0 0 0];
    type_DM_tmp(i)=1;
    %create a single phase block
    type_DM{end+1}=[type_DM_tmp; 0 0 0 0];    
    %create a single phase block with apoptosis in all phases
    type_DM{end+1}=[type_DM_tmp; 1  1  1  1];
    %create a single phase block with apoptosis in one phase
    for j=1:4
        type_DM_tmp_2=[0  0   0   0];
        type_DM_tmp_2(i)=1;
        type_DM{end+1}=[type_DM_tmp; type_DM_tmp_2];    
    end       
end     

%% generate the parameter sets for each case, run simulation and save results
for i=1:length(type_DM)
    %load the type of DM model
    type_DM_tmp=type_DM{i};
    %create a new directory
    dir_name=[cc_names{type_DM_tmp(1,:)==1} '_block_' cc_names{type_DM_tmp(2,:)==1} '_apoptosis'];
    mkdir(['./SynData/' dir_name]);
    %run the same DM multiple time to generate different dynamics
    nRun=3;
    for j=1:nRun
        %create a new subdirectory
        sub_dir_name=['./SynData/' dir_name '/Data' num2str(j) ];
        mkdir(sub_dir_name);
        %run the simulation
        datasim=sim_drug_effects(@()gen_param_set(type_DM_tmp));       
        %save the data and the figures and closes the figures
        save([sub_dir_name '/datasim.mat'], 'datasim');
        f=findobj('type','figure','name','Steady state behaviour');
        figure(f);
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf,'PaperPositionMode','auto')
        saveas(gcf,[sub_dir_name '/Steady state behaviour'],'png');
        close(f);
        f=findobj('type','figure','name','ODE dynamics');
        figure(f);
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf,'PaperPositionMode','auto')
        saveas(gcf,[sub_dir_name '/ODE dynamics'],'png');
        close(f);
        f=findobj('type','figure','name','Single cell dynamics');
        figure(f);
        set(gcf, 'Position', get(0, 'Screensize'));
        set(gcf,'PaperPositionMode','auto')
        saveas(gcf,[sub_dir_name '/Single cell dynamics'],'png');
        close(f);        
    end
    
    
end    




end

%function that generates 
function [kG1S kSG2 kG2M kMG1 papG1 papS papG2 papM Vcc Kcc ncc Vap Kap nap]=gen_param_set(type_DM)

%define the parameters that define the probability from which the kinetic parameters are sampled
muk=[15 15 4 2]; sigmak=[5 5 1 0.5];

%define the untreated baseline
kG1S=1./normrnd(muk(1),sigmak(1));
kSG2=1./normrnd(muk(2),sigmak(2));
kG2M=1./normrnd(muk(3),sigmak(3));
kMG1=1./normrnd(muk(4),sigmak(4));
papG1=0;
papS=0;
papG2=0;
papM=0;

%maximal activity of Vpp is equal to the untreated baseline
Vcc(1,1)=kG1S; Vcc(2,1)=kSG2; Vcc(3,1)=kG2M; Vcc(4,1)=kMG1;
%define the paramters that control the dose-dependen drug effects
KccS=0.00001.*10000 * rand(1,1); nccS=0.1+1.5*rand(1,1);
VapS=0.01+0.3*rand(1,1); KapS=0.00001.*10000 * rand(1,1); napS=0.1+1.5*rand(1,1); 
%sample a single rate for CC and apoptosis and use it for all the defined phases
for i=1:4
    Kcc(i,1)=type_DM(1,i) .* KccS;
    ncc(i,1)=type_DM(1,i) .* nccS;
    Vap(i,1)=(type_DM(2,i) .* VapS);
    if ((sum(type_DM(2,i)))~=0)
         Vap(i,1)= Vap(i,1)./(sum(type_DM(2,i)));
    end     
    Kap(i,1)=type_DM(2,i) .* KapS;
    nap(i,1)=type_DM(2,i) .* napS;
end
%k parameters set to infinite or hill equation to 1 if not necessary
Kcc(Kcc==0)=Inf;
Kap(Kap==0)=Inf;
ncc(ncc==0)=1;
nap(nap==0)=1;

end
