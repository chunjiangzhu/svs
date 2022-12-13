function [SC] = Schur_Comp(W, index)

degs_site = sum(W, 2);
D = sparse(1:size(W, 1), 1:size(W, 2), degs_site);
L_site = D - W;

n = size(W,1);
temp = 1:1:n;
temp(index) = [];
non_index = temp;

% get SC
L_site_vv = L_site(index, index);
L_site_vvp = L_site(index, non_index);
L_site_vpvp = L_site(non_index, non_index);
L_site_vpv = L_site(non_index, index);
% SC_site = L_site_vv - L_site_vvp * lsqminnorm(L_site_vpvp, L_site_vpv); 
SC_site = L_site_vv - L_site_vvp * inv(L_site_vpvp) * L_site_vpv; 

SC = abs(triu(SC_site,1));
% SC = Graph_site + transpose(Graph_site);