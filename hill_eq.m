function [ y ] = hill_eq( x, V, k, n )
%HILL_EQ is a generic hill equation function

y= V .* ((x.^n)./((x.^n)+(k.^n)));

end

