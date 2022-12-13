function [ B, IND ] = calIncidence( M )
%calIncidence Calculates edge-vertex incidence matrix from its adjacency matrix.
%
%   Author: 
%   Year  : 


    % number of nodes
    n = size(M, 2);
    m = 0;
    for i=1:n
        for j=i:n
            if M(i,j)~=0
                m = m+1;
            end
        end
    end
    
    indi = [];
    indj = [];
    inds = [];
    IND = zeros(m,2);
    m = 0;
    
    for i=1:n
        for j=i:n
            if M(i,j)~=0
                m = m+1;
                IND(m,1) = i;
                IND(m,2) = j;
                v = sqrt(M(i,j));
                
                indi = [indi,m,m];
                indj = [indj,i,j];
                inds = [inds,v,-v];
            end
        end
    end
    
    B = sparse(indi, indj, inds, m, n);

end