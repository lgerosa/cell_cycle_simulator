function [] = sim_two_drugs_effects()
%SIM_TWO_DRUGS_EFFECTS simulates the co-dosing of two drugs that act
%independently on the cell cycle and apopotosis and calculates the level of
%synergy 

%% DEFINE UNTREATED CELLS
%defien cell cycle rates
tG1=13; tS=20; tG2=5; tM=4;
%calculate transition rates for cell cycle transition matrix (CTM)
%G1->S   S->G2   G2->M    M->G1
kG1S=1/tG1; kSG2=1/tS; kG2M=1/tG2; kMG1=1/tM;
kTUnt=[kG1S; kSG2; kG2M; kMG1];
%define uperturbed average time for apoptosis
taG1=Inf; taS=Inf; taG2=Inf; taM=Inf;
%calculate transition rates for the apoptosis transition matrix (ATM)
%G1->   S->   G2->    M->
papG1=1/taG1; papS=1/taS; papG2=1/taG2; papM=1/taM;
pTUnt=[papG1; papS; papG2; papM];

%% DEFINE TREATMENTS
%define drug A
%define the dose-dependent drug effects on Cell Cycle
     %G1   S    G2   M
VccA=[kG1S; kSG2; kG2M; kMG1]; %maximal rate on rev hill equation 
KccA=[0.1; Inf;  0.01;  Inf]; %drug concentration for half maximal effect 
nccA=[1;    1;    0.5;    1;  ]; %steepness of drug effect 
%define the dose-dependent drug effects on apoptosis
     %G1    S     G2     M
VapA=[papG1; papS; papG2;  papM]; %maximal rate on hill equation 
KapA=[Inf;  Inf;  Inf;     Inf]; %drug concentration for half maximal effect 
napA=[1;     1;    1;     1;  ]; %steepness of drug effect 

%define drug B
%define the dose-dependent drug effects on Cell Cycle
     %G1   S    G2   M
VccB=[kG1S; kSG2; kG2M; kMG1]; %maximal rate on rev hill equation 
KccB=[0.0001;  Inf;  Inf;    Inf]; %drug concentration for half maximal effect 
nccB=[1;    1;    1;    1;  ]; %steepness of drug effect 
%define the dose-dependent drug effects on apoptosis

     %G1    S     G2     M
VapB=[1;    1;      1;    1]; %maximal rate on hill equation 
KapB=[0.1;   0.1;  0.1;    0.1; ]; %drug concentration for half maximal effect 
napB=[1;     1;    1;      1;  ]; %steepness of drug effect 

%define the sinergy table
%define the synergy table for cell cycle progression and apoptosis between
%the two drugs

%Synergy in cell cycle control between drug A and B
%       G1      S       G2        M          
STcc = [ 1       1       1        1];
%Synergy in apoptosis between drug A and B
%       G1      S       G2        M   
STap = [ 1       1       1        1];


%% PLOT CELL AND TREATMENT DETAILS

%dilution DrugA
minDrugA=-5;
maxDrugA=1;
ndosesA=100;
ConcA=[0 logspace(minDrugA,maxDrugA,ndosesA)];
%ConcA=[0 linspace(10.^minDrugA,10.^maxDrugA,ndosesA)]; 

%dilution DrugB
minDrugB=-5;
maxDrugB=1;
ndosesB=100;
ConcB=[0 logspace(minDrugB,maxDrugB,ndosesB)];
%ConcB=[0 linspace(10.^minDrugB,10.^maxDrugB,ndosesB)];

%calculate dose-dependent drug effect on each cell cycle rate
nConcA=length(ConcA); nparam=length(VccA);
drugA_mod_cc=arrayfun(@rev_hill_eq,repmat(ConcA,nparam,1),repmat(VccA,1,nConcA),repmat(KccA,1,nConcA),repmat(nccA,1,nConcA));
drugA_mod_ap=arrayfun(@hill_eq,repmat(ConcA,nparam,1),repmat(VapA,1,nConcA),repmat(KapA,1,nConcA),repmat(napA,1,nConcA));
%calculate dose-dependet drug effect on apoptosis
nConcB=length(ConcB); nparam=length(VccB);
drugB_mod_cc=arrayfun(@rev_hill_eq,repmat(ConcB,nparam,1),repmat(VccB,1,nConcB),repmat(KccB,1,nConcB),repmat(nccB,1,nConcB));
drugB_mod_ap=arrayfun(@hill_eq,repmat(ConcB,nparam,1),repmat(VapB,1,nConcB),repmat(KapB,1,nConcB),repmat(napB,1,nConcB));

figure('Name','Steady state behaviour');
set(gcf,'Color','White');
colors=[0 0 1; 0.8 0.8 0; 1 0 0; 0 0.8 0];%colors
cc_names={'G1' 'S' 'G2' 'M'};
lw=2; %line width

%drugA cell cycle
ncol=6;
nrow=4;
subplot(nrow,ncol,1);
for i=1:size(drugA_mod_cc,1)
    plot(ConcA,drugA_mod_cc(i,:),'LineWidth',lw,'Color',colors(i,:));
    hold on;
end
%plot legend and axis
legend(cc_names);
set(gca,'XScale','log');
ylim([0 max(drugA_mod_cc(:))+0.1]);
xlim([min(ConcA) max(ConcA)]);
xlabel('Drug (\muM)');
ylabel('Rate (h^{-1})');
title('Drug A effect on cell cycle');

%drugA apoptosis
subplot(nrow,ncol,7);
for i=1:size(drugA_mod_ap,1)
    plot(ConcA,drugA_mod_ap(i,:),'LineWidth',lw,'Color',colors(i,:));
    hold on;
end
%plot legend and axis
legend(cc_names);
set(gca,'XScale','log');
ylim([0 max(drugA_mod_ap(:))+0.1]);
xlim([min(ConcA) max(ConcA)]);
xlabel('Drug (\muM)');
ylabel('Rate (h^{-1})');
title('Drug A effect on apoptosis');

%drugB cell cycle
subplot(nrow,ncol,2);
for i=1:size(drugB_mod_cc,1)
    plot(ConcB,drugB_mod_cc(i,:),'LineWidth',lw,'Color',colors(i,:));
    hold on;
end
%plot legend and axis
legend(cc_names);
set(gca,'XScale','log');
ylim([0 max(drugA_mod_cc(:))+0.1]);
xlim([min(ConcB) max(ConcB)]);
xlabel('Drug (\muM)');
ylabel('Rate (h^{-1})');
title('Drug B effect on cell cycle');

%drugB apoptosis
subplot(nrow,ncol,8);
for i=1:size(drugB_mod_ap,1)
    plot(ConcB,drugB_mod_ap(i,:),'LineWidth',lw,'Color',colors(i,:));
    hold on;
end
%plot legend and axis
legend(cc_names);
set(gca,'XScale','log');
ylim([0 max(drugB_mod_ap(:))+0.1]);
xlim([min(ConcB) max(ConcB)]);
xlabel('Drug (\muM)');
ylabel('Rate (h^{-1})');
title('Drug B effect on apoptosis');

%% CALCULATE EFFECTS FOR EACH SINGLE DRUG ALSO ONLY CYTOSTATIC AND ONLY CYTOTOXIC AT STEADY STATE

%define the metrix to be used to plot growth inhibition 
f=@GR_Metric;
label_gr='GR value';
%f=@Identity_Metric;
%label_gr='\mu_{DMSO}/\mu_{drug(s)}';

%for drugA
for i=1:length(ConcA)
    %calculate rates modification by each drug
    kA=drugA_mod_cc(:,i);
    pA=drugA_mod_ap(:,i);
    %calculate overall effect
    [TM, CTM, AT] = create_TM(kA(1),kA(2),kA(3),kA(4),pA(1),pA(2),pA(3),pA(4));
    %[AT, ~, ~] = create_TM(kTUnt(1),kTUnt(2),kTUnt(3),kTUnt(4),pA(1),pA(2),pA(3),pA(4));
    ss_gr_A(i)= calculate_steady_state(TM);
    ss_gr_A_static(i)= calculate_steady_state(CTM);
    ss_gr_A_toxic(i)= calculate_steady_state(AT);
    
end    
ss_gr_A_norm=ss_gr_A./ss_gr_A(1);
ss_gr_A_static_norm=ss_gr_A_static./ss_gr_A(1);
ss_gr_A_toxic_norm=ss_gr_A_toxic./ss_gr_A(1);


%for drugB
for i=1:length(ConcB)
    %calculate rates modification by each drug
    kB=drugB_mod_cc(:,i);
    pB=drugB_mod_ap(:,i);
    %calculate overall effect
    [TM, CTM, AT] = create_TM(kB(1),kB(2),kB(3),kB(4),pB(1),pB(2),pB(3),pB(4));
    %[AT, ~, ~] = create_TM(kTUnt(1),kTUnt(2),kTUnt(3),kTUnt(4),pB(1),pB(2),pB(3),pB(4));
    ss_gr_B(i)= calculate_steady_state(TM);
    ss_gr_B_static(i)= calculate_steady_state(CTM);
    ss_gr_B_toxic(i)= calculate_steady_state(AT);
end    
ss_gr_B_norm=ss_gr_B./ss_gr_B(1);
ss_gr_B_static_norm=ss_gr_B_static./ss_gr_B(1);
ss_gr_B_toxic_norm=ss_gr_B_toxic./ss_gr_B(1);

%analyze doses of equal potency (theoretical isobols)
DrugEffectRange=[max([min(ss_gr_A_norm) min(ss_gr_B_norm)]) min([max(ss_gr_A_norm) max(ss_gr_B_norm)])];
nEffectIso=10;
EffectIso=linspace(DrugEffectRange(1), DrugEffectRange(2),nEffectIso);
DrugAEffectIso=interp1(ss_gr_A_norm,ConcA,EffectIso);
DrugBEffectIso=interp1(ss_gr_B_norm,ConcB,EffectIso);

%plot DrugA cytostatic and cytotoxic component
subplot(nrow,ncol,13);
plot(ConcA,f(ss_gr_A_static_norm),'-c','LineWidth',2); 
hold on;
plot(ConcA,f(ss_gr_A_toxic_norm),'-m','LineWidth',2);
legend({'cytostatic' 'cytotoxic'});
xlim([min(ConcA) max(ConcA)]);
set(gca,'XScale','log');
xlabel('Drug(\muM)');
ylabel(label_gr);
title('Drug A static and toxic effects');

%plot DrugA actual and independent effects
subplot(nrow,ncol,19);
plot(ConcA,f(ss_gr_A_norm),'-k','LineWidth',2);
hold on;
plot(ConcA,f(ss_gr_A_static_norm+ss_gr_A_toxic_norm),'--','LineWidth',2,'Color',[0.5 0.5 0.5]);
legend({'actual' 'independent'});
xlim([min(ConcA) max(ConcA)]);
set(gca,'XScale','log');
xlabel('Drug(\muM)');
ylabel(label_gr);
title('Drug A mes and ind effects');
%plot(DrugAEffectIso,f(EffectIso),'o','MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',5);

%plot DrugB cytostatic and cytotoxic component
subplot(nrow,ncol,14);
plot(ConcB,f(ss_gr_B_static_norm),'-c','LineWidth',2); 
hold on;
plot(ConcB,f(ss_gr_B_toxic_norm),'-m','LineWidth',2);
legend({'cytostatic' 'cytotoxic'});
xlim([min(ConcB) max(ConcB)]);
set(gca,'XScale','log');
xlabel('Drug(\muM)');
ylabel(label_gr);
title('Drug A measured and independent');

%plot DrugB actual and independent effects
subplot(nrow,ncol,20);
plot(ConcB,f(ss_gr_B_norm),'-k','LineWidth',2);
hold on;
plot(ConcB,f(ss_gr_B_static_norm+ss_gr_B_toxic_norm),'--','LineWidth',2,'Color',[0.5 0.5 0.5]);
legend({'actual' 'independent'});
xlim([min(ConcB) max(ConcB)]);
set(gca,'XScale','log');
xlabel('Drug(\muM)');
ylabel(label_gr);
title('Drug B mesured and independent');
%plot(DrugBEffectIso,f(EffectIso),'o','MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',5);


%% CALCULATE STEADY-STATE EFFECT ON GROWTH RATE OF THE TWO DRUGS AT DIFFERENT DILUTIONS  

%calculate steady state effects for drug combinations
ss_gr=[];
for i=1:length(ConcA)
    for j=1:length(ConcB)
        %calculate rates modification by each drug
        kA=drugA_mod_cc(:,i);
        pA=drugA_mod_ap(:,i);
        kB=drugB_mod_cc(:,j);
        pB=drugB_mod_ap(:,j);
        %calculate indepenzdence effect
        %kT=(kA+kB)./2;
        %pT=pA+pB;
        %% HOW TO CALCULATE THE INDEPENDENCE AND HOW TO INTRODUCE A PARAMETER TO MODULATE COOPERATIVITY?? THE BIG QUESTION
        kT=1./((1./kTUnt)+((1./kA)-(1./kTUnt))+((1./kB)-(1./kTUnt)));    
        pT=pTUnt+((pA-pTUnt)+(pB-pTUnt));
        %kT=1./((1./kA)+(1./kB));
        %pT=pA+pB;
        %kT=min([kA kB],[],2);
        %pT=max([pA pB],[],2);
        %calculate steady state by solving eingenvector
        [TM, CTM, AT] = create_TM(kT(1),kT(2),kT(3),kT(4),pT(1),pT(2),pT(3),pT(4));
        [ss_growth_rate, ss_dist_cell_cycle] = calculate_steady_state(TM);
        ss_gr(j,i)=ss_growth_rate;
        %calculate steady state growth rate by finding the roots of the 4 degree polynomial        
        %a=; b=; c=; d=; e=;
        %gr=roots([a b c d e]);
        %sol_roots=roots([a b c d e]);
        %gr=sol_roots(imag(sol_roots)==0);       
        %ss_gr(i,j)=gr;
    end
end    

%convert to the right metric
ss_gr=f(ss_gr./ss_gr(1,1));


%% CALCULATE SYNERGY SCORE (ISOBOLGRAM PLOTS)



%IsobolA=((max(ss_gr(1,:))-min(ss_gr(1,:)))./2)+min(ss_gr(1,:));
%IsobolB=((max(ss_gr(:,1))-min(ss_gr(:,1)))./2)+min(ss_gr(:,1));
%DrugAIso=interp1(ss_gr(1,:),ConcA,IsobolA);
%DrugBIso=interp1(ss_gr(:,1),ConcB,IsobolB);

%{

%plot synergy table for CC
subplot(nrow,ncol,5);
heatmap(STcc,cc_names,cc_names,'%2.1f','MinColorValue',-5,'MaxColorValue',5,'ShowAllTicks',1);
colorbar;
title('Drug cooperativity in cell cycle (Not Used)');
xlabel('Drug A');
ylabel('Drug B');

%plot synergy table for AP
subplot(nrow,ncol,9);
heatmap(STcc,cc_names,cc_names,'%2.1f','MinColorValue',-5,'MaxColorValue',5,'ShowAllTicks',1);
colorbar;
title('Drug cooperativity in apoptosis (Not Used)');
xlabel('Drug A');
ylabel('Drug B');
%}

%plot ratio of drug concentrations for same efficiency
%subplot(nrow,ncol,14);
%{
subplot(nrow,ncol,13);
plot(EffectIso,DrugAEffectIso./DrugBEffectIso,'-o','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',5);
xlabel(label_gr);
ylabel('DrugA (\muM) /DrugB (\muM)');
title('Drug Ratio for same effect');
%}

%subplot(nrow,ncol,[ 7 8 11 12 15 16]);
subplot(nrow,ncol,[3 4 5 6 9 10 11 12 15 16 17 18 21 22 23 24]);
%imagesc(ss_gr(end:-1:1,:));
surf(ConcA,ConcB,ss_gr,'EdgeColor','None');
view(2); 
xlabel('Drug A');
ylabel('Drug B');
zlabel('Growth');
set(gca,'XScale','log');
set(gca,'YScale','log');
c=colorbar;
ylabel(c,label_gr);
%get DrugA and DrugB coordinates
%theoIsoDrugA=interp1(ConcA,1:length(ConcA),DrugAIso);
%theoIsoDrugB=interp1(ConcB,length(ConcB):-1:1,DrugBIso);  
%
%plot theoretical linear isobols
ndil=10000;
for i=1:length(DrugAEffectIso)
    A_dil=linspace(ConcA(2),DrugAEffectIso(i),ndil);
    B_dil=DrugBEffectIso(i)-(A_dil.*(DrugBEffectIso(i)./DrugAEffectIso(i)));
    A_dil(B_dil<ConcB(2))=[];
    B_dil(B_dil<ConcB(2))=[];
    hold on;
    plot3(A_dil,B_dil,repmat(EffectIso(i),size(A_dil)),'r--','LineWidth',2);
    
    %line([theoIsoDrugA(i) 0 ],[length(ConcB)+1 theoIsoDrugB(i)],'Color','r','LineWidth',0.5,'LineStyle','--');
end    

%ylabel(c,'GR_{DMSO}/GR_{drugA,drugB}');

%% CALCULATE THE SAME WITH TIME-DEPENDENT MEASUREMENTS (with GR correction and without)

end

%convert to GR value
function [GRvalue] = GR_Metric(mu)
  GRvalue= 2.^mu -1;
end

%identity function
function [y] = Identity_Metric(x)
  y= x;
end