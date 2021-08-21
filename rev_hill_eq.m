function [m] = rev_hill_eq( x, V, k, n )
%DRUG2EFFECT transforms a drug concentration into a rate modifier using a
%monotonically decreasing hill equation with a maximal value (Vmax), half
%effect concentration (k) and steepness (n)

m=V.*(1-hill_eq(x,1,k,n));


end

