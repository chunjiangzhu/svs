function [result] = check_boundry(W_site_set, s, boundry_index, n)
% Check if non-boundry node only have edges in single site.
% Check if boundry node have edges in different sites

temp = 1:n;
temp(boundry_index) = [];
non_boundry_index = temp;

for i = non_boundry_index
    non_zero_site = 0;
    for j=1:s
        a = full(sum(W_site_set{j}(i,:)));
        if a~=0
            non_zero_site = non_zero_site + 1;
        end
    end
    if non_zero_site > 1 
        fprintf("node %d is not boundy, but have edges in multiple sites\n", i);
    end
end


for i = boundry_index
    non_zero_site = 0;
    for j=1:s
        a = full(sum(W_site_set{j}(i,:)));
        if a~=0
            non_zero_site = non_zero_site + 1;
        end
    end
    if non_zero_site <= 1 
        fprintf("node %d is boundy, but have edges in less than 2 sites\n", i);
    end
end