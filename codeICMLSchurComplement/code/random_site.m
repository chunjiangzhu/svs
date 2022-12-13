function [W_site_set, boundry_index] = random_site(W, s)
% -------------------------------------------------------
% seperate edges into S sites.
% This is randomly seperate version, the boundry size may be huge.
% We provide another version based on coordinates.
% -------------------------------------------------------

W_site_set = cell(s,1);
num_edges = nnz(W);
n = size(W,1);
rand_index = randi([1 s],1,num_edges);
[row,col,v] = find(W);    % Find nonzero index and value for W.

site_dic = containers.Map; % key is row_col index, value is site index.


node_site_dic = containers.Map; % If one edge exists in site s, then 2 end nodes belong to this site also.
                                % If one node belongs to multiple sites, it
                                % is boundry node.
boundry_index = [];                                
for i=1:s
    index_site = find(rand_index == i);
    col_site = col(index_site);
    row_site = row(index_site);
    % Consider that W_site should be symmetric, we only consider upper
    % triangle part. And then make it symmetric.
    W_site = triu(sparse(row_site, col_site, v(index_site), n, n));
    W_site = transpose(W_site) + W_site;
    W_site_set{i} = W_site;
    % get the node site and determine if it is boundry. 
    for j=1:size(row_site,1)
        if row_site(j) < col_site(j)
            if isKey(node_site_dic, num2str(row_site(j)))
                if node_site_dic(num2str(row_site(j))) ~= i
                    boundry_index = [boundry_index, row_site(j)];
                end
            else
                node_site_dic(num2str(row_site(j))) = i;
                
            end
            if isKey(node_site_dic, num2str(col_site(j)))
                if node_site_dic(num2str(col_site(j))) ~= i
                    boundry_index  = [boundry_index, col_site(j)];
                end
            else
                node_site_dic(num2str(col_site(j))) = i;
            end
        end
    end
end