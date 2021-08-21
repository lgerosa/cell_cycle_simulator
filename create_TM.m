function [TM, CTM, AT] = create_TM(kG1S,kSG2,kG2M,kMG1,papG1,papS,papG2,papM)
%CREATE_TM creates a transition matrix (TM) that combines rates of
%transiotion in G1,S,G2,M (kG1,etc..) and rates of apoptosis in
%G1,S,G2,M (taG1, etc..). It returns also the two intermediates matrices
%CTM (only cell cycle transition) and AT (only apoptotic transitions)

%define matrix of state transitions between cell cycle phase
    %G1             S               G2              M
CTM=[-kG1S          0               0               2*kMG1         ;     %G1
     kG1S           -kSG2           0               0              ;     %S
     0              kSG2            -kG2M           0              ;     %G2
     0              0               kG2M            -kMG1         ];     %M  

%Apoptotic matrix: defines the probability of apoptosis in a certain phase  
    %G1             S               G2              M
AT=[ -papG1         0               0               0              ;     %G1
     0              -papS           0               0              ;     %S
     0              0               -papG2          0              ;     %G2
     0              0               0               -papM         ];     %M

%calculate the overall TM 
TM=CTM+AT;

end

