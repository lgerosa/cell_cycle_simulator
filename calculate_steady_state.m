function [ss_growth_rate, ss_dist_cell_cycle] = calculate_steady_state(TM)
%CALCULATE_STEADY_STATE calculates the steady state distribution of a
%number of cells (ncells cycling according to rates in STM (state transition matrix) and
%AM (apoptotic matrix)

%add a line to have conservation of total fractions
TM_const=TM;
TM_const(end+1,:)=[1 1 1 1];
ss_cond=[0;0;0;0;1];
%calculate the steady state growth rate and cell fraction distribution
%using eigenvalues (TM * x = mu * x). The positive eigenvalue is the growth
%rate, the corresponding eigenvector is the fractional distribution in the
%cell cycle
[V,D] = eig(TM); %V: eigenvectors (columns0, D: eigenvalues (diagonal)
E=diag(D);
%find real eigenvalues
E_real=(imag(E)==0); 
%find positive eigenvalues
E_pos=(real(E)>=0); 
%find real eigenvectors
V_real=all(imag(V)==0)'; 
%find eigenvectors with coherent signs
V_sign=sign(real(V));
V_sign_pos=V_sign; V_sign_neg=V_sign;
V_sign_pos(V_sign_pos==0)=1;  V_sign_pos(V_sign_pos==-1)=0; 
V_sign_neg(V_sign_neg==0)=-1; V_sign_neg(V_sign_neg==1)=0; 
V_coherent=(all(V_sign_pos)|all(V_sign_neg))';
%if there is a unique real eigenvalue with a coherent eigenvector (or all are the same value), then that is the stable growth rate
i_sol=intersect(find(E_real==1),find(V_coherent));
if ((~isempty(i_sol)&(length(i_sol)==1))|(sum(E==E(1))==length(E)))
    ss_growth_rate=E(i_sol(1));
    ss_dist_cell_cycle=abs(V(:,i_sol(1)))./sum(abs(V(:,i_sol(1)))); 
else
    %there isn't a single real eigenvalue with coherent eigenvector (either
    %none or multiple with different values), no stable growth rate
    ss_growth_rate=NaN;
    ss_dist_cell_cycle=repmat(NaN,size(TM,1),1);
end   


%{
isol_real=find(imag(E)==0); 
%if there are real solution, then there is a stable growth rate
if (~isempty(isol_real))
    %check if there is a strictly positive growth rate, return it
    isol_pos=find((E>0)& (imag(E)==0)); 
    if (~isempty(isol_pos))
        ss_growth_rate=E(isol_pos);
        ss_dist_cell_cycle=abs(V(:,isol_pos))./sum(abs(V(:,isol_pos)));
    else 
       %otherwise check if there is a negative real eigenvalue with non imaginary vectors all with the same sign
       isol_neg=find((E<0)& (imag(E)==0)); 
       %check eigenvector have all same sign or zero
       V_pos=sign(real(V));       V_neg=sign(real(V));
       V_pos(V_pos==0)=1;         V_neg(V_neg==0)=-1;
       V_pos(V_pos==-1)=0;         V_neg(V_neg==1)=0;
       isol_eig_same_sign=find(all(V_neg,1)|all(V_pos,1));
       isol_neg=intersect(isol_neg,isol_eig_same_sign);
       %check eigenvectors are real
       isol_eig_real=find(sum(abs(imag(V))==0,1)==size(V,1));
       isol_neg=intersect(isol_neg,isol_eig_real);
       if (~isempty(isol_neg))
          %check how many eigenvalue found   
          if(length(isol_neg)==1)
              ss_growth_rate=E(isol_neg);
              ss_dist_cell_cycle=abs(V(:,isol_neg))./sum(abs(V(:,isol_neg)));           
          else
              %check if they are all the same eigenvalue %% not sure is
              %correct!! perhaps should be NaN also here
              E_test=E(isol_neg);
              if (sum(E_test==E_test(1))==length(E_test))
                 ss_growth_rate=E(isol_neg(1));
                 ss_dist_cell_cycle=abs(V(:,isol_neg(1)))./sum(abs(V(:,isol_neg(1))));
              else
                  ss_growth_rate=NaN;
                  ss_dist_cell_cycle=repmat(NaN,size(TM,1),1);
              end 
          end    
       else
           %otherwise all is zero, pick any solution
           ss_growth_rate=E(isol_real(1));
           ss_dist_cell_cycle=abs(V(:,isol_real(1)))./sum(abs(V(:,isol_real(1))));
       end  
    end   
else
    %otherwise set growth rate to NaN because the system does not have a stable gowth rate
    ss_growth_rate=NaN;
    ss_dist_cell_cycle=repmat(NaN,size(TM,1),1);
end  
%}


%% OLD CODE FOR BACKUP
%{
%add a line to have conservation of total fractions
TM_const=TM;
TM_const(end+1,:)=[1 1 1 1];
ss_cond=[0;0;0;0;1];
%calculate the steady state growth rate and cell fraction distribution
%using eigenvalues (TM * x = mu * x). The positive eigenvalue is the growth
%rate, the corresponding eigenvector is the fractional distribution in the
%cell cycle
[V,D] = eig(TM); %V: eigenvectors (columns0, D: eigenvalues (diagonal)
%find, if existing, positive eigenvalues
E=diag(D);
isol=find((E>=0)& (imag(E)==0)); 
%check that isn't the trivial solution of all zeros
if (sum(E==0)==length(E))
   isol=isol(1);
end    
%if there is a positive rational eigenvalue, is the growth rate
if ~isempty(isol)
    ss_growth_rate=E(isol);
    ss_dist_cell_cycle=abs(V(:,isol))./sum(abs(V(:,isol)));
else
    %if there is not a positive rational eigenvalue, search for a negative
    %rational with a coherent eigenvector (all elements same sign), which
    %represents the negative growth rate due to cytotoxycity
    isolE=find((E<=0)& (imag(E)==0)); 
    isolV=find(~std(sign(V)));
    isol=intersect(isolE,isolV);
    if ~isempty(isol)
        ss_growth_rate=E(isol);
         ss_dist_cell_cycle=abs(V(:,isol))./sum(abs(V(:,isol)));
    else
        ss_growth_rate=NaN;
        ss_dist_cell_cycle=repmat(NaN,size(TM,1),1);
    end    
end    
%calculate steady state distribution in cell cycle phases (old way to
%calculate distribution, somehow not exactly similar to the eigenvalue
%method, possibly for numerical reasons)
%ss_dist_cell_cycle_2= linsolve(TM_const,ss_cond);

%}

end

