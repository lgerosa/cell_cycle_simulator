function [kG1S kSG2 kG2M kMG1 papG1 papS papG2 papM Vcc Kcc ncc Vap Kap nap] = DM_param_example1_G2_block_S_apoptosis()
%CC_PARAM_EXAMPLE1_G2_BLOCK_S_APOPTOSIS defines a Drug Model (DM parameter set)

%% define uperturbed average time in cell cycle phases and apoptotic rates (in hours)
%defien cell cycle rates
tG1=10; tS=11; tG2=2; tM=2;
%calculate transition rates for cell cycle transition matrix (CTM)
%G1->S   S->G2   G2->M    M->G1
kG1S=1/tG1; kSG2=1/tS; kG2M=1/tG2; kMG1=1/tM;
%define uperturbed average time for apoptosis
taG1=Inf; taS=Inf; taG2=Inf; taM=Inf;
%calculate transition rates for the apoptosis transition matrix (ATM)
%G1->   S->   G2->    M->
papG1=1/taG1; papS=1/taS; papG2=1/taG2; papM=1/taM;

%define the dose-dependent drug effects on Cell Cycle
     %G1   S    G2   M
Vcc=[kG1S; kSG2; kG2M; kMG1]; %maximal rate on rev hill equation 
Kcc=[Inf;  Inf;  0.01;   Inf]; %drug concentration for half maximal effect 
ncc=[1;    1;    1;    1;  ]; %steepness of drug effect 

%define the dose-dependent drug effects on apoptosis
     %G1    S     G2     M
Vap=[papG1; 0.1; papG2;  papM]; %maximal rate on hill equation 
Kap=[Inf;   10;  Inf;     Inf]; %drug concentration for half maximal effect 
nap=[1;     1;    1;     1;  ]; %steepness of drug effect 



end

