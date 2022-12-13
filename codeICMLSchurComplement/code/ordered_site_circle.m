function [W_site_set, boundry_index] = ordered_site_circle(W, T, s)
% seperate sites based on location (X coordinate).
W_site_set = cell(s,1);
W = triu(W); % only consider the upper half of W.
n = size(W,2);

angles = atan2d(T(2,:), T(1,:));
p = [-180:360/s:180];


% sortT = sort(T(1,:));     % sort according to X coordinates
% p = [-1000000 sortT(1,ceil(n/s):ceil(n/s):n-1) 1000000];

node_site_dic = containers.Map;
% If one edge exists in site s, then 2 end nodes belong to this site also.
% If one node belongs to multiple sites, it is boundry node.

boundry_index = [];
% W_check = zeros(size(W));
W = triu(W);
for i=1:s
    % seperate points into s parts according to angle.
    index_in_range = (find((p(1,i+1)>angles) & (angles>=p(1,i))));
    
    [row_slice,col_slice, v_slice] = find(W(index_in_range,:));   
    row_site = index_in_range(row_slice)';
    col_site = col_slice;
    
    % Consider that W_site should be symmetric, we only consider upper
    % triangle part. And then make it symmetric.
    W_site = triu(sparse(row_site, col_site, v_slice, n, n));
    W_site = transpose(W_site) + W_site;
    W_site_set{i} = W_site;
    
    % get the node site and determine if it is boundry. 
    % If edge belongs to site i, then its 2 end nodes also belong to site i;
    % If node belongs to multiple sites, it is boundry.
    for j=1:size(row_site,1)
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
%     W_check = W_check + W_site;
end
% Check if the union of W_site would equal to W.

% W = W + transpose(W);
% check = isequal(W, W_check);
% [row,col,v] = find(W);
% [row1,col1,v1] = find(W_check);
% row_check = isequal(row, row1);
% col_check = isequal(col, col1);
% fprintf("test W and W_site sum", test)
