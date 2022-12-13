% Save all H2 and use Julia to sparsify.
% Also save site graphs and use Julia to sparsity (for baseline.)

% for one graph, randomly sample the terminal 5 times for each setting.
clear;

%setting
test_mode = 3
Type = 'MP' % or BB
d_epsilon = 0.3
d_sample_rate = 0.05
d_s = 5

sample_rate_range = [0.05];
s_range = [5];
epsilon_range = [0.3]%0.2:0.1:1.2;

cRange = [5,10,15,20];  % 5-20;      % c control sample rate.
lambdaRange = [0.1, 1, 5, 10, 15]; %[0.05, 0.1, 0.2, 0.5, 0.6, 1, 2,3,4, 5,6,7,8,9 10,11,12,13,14,15, 20]; %[0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30,50,80,100,150,200];  %20:2:64;

if test_mode == 3
    T = dlmread('T_sculpture.txt');
    load('W_sculpture.mat');
elseif test_mode == 1
    T = dlmread('T_gaussian.txt');
%     load('W_gaussian.mat', 'W');
    load('W_gaussian.mat');
    num_clust = 4;
    x_gap = 0;
elseif test_mode == 4
    T = dlmread('T_circle.txt');
    load('W_circle.mat');
    num_clust = 2;
    x_gap = 0;
elseif test_mode == 5
    T = dlmread('T_camera.txt');
    load('W_camera.mat');
    num_clust = 3;
end
m = nnz(W)/2
n = size(W,1)
d = ceil(log(n)/log(2)); 

all_settings = [];
fix_epsilon_settings = [];
for sample_rate = sample_rate_range
    for s = s_range
        for epsilon = epsilon_range
            setting = [sample_rate s epsilon];
            if sample_rate~= d_sample_rate & s ~= d_s & epsilon ~= d_epsilon
                continue
            else
               all_settings = [all_settings; setting];
            end
        end
    end
end

succ_record = [];

data_folder = './experiment_results/';
root_folder = [data_folder 'exp_mode' num2str(test_mode) '/type' num2str(SimGraphType) '_para' num2str(Param) '_sigma' num2str(Sigma) '_seed' num2str(seed) '/'];
mkdir(root_folder);

for i = 1:size(all_settings, 1)
    setting = all_settings(i,:);
    sample_rate = setting(1);
    s = setting(2);
    epsilon = setting(3);
    
    save_filename = [root_folder 'exp' 'rate_' num2str(sample_rate) '_s' num2str(s) '_eps' num2str(epsilon) Type '.csv']
    
    if test_mode == 4
        [W_site_set, boundry_index] = ordered_site_circle(W, T, s);
    else
        [W_site_set, boundry_index] = ordered_site_chunk(W, T, s);
    end
    boundry_index = unique(boundry_index);
    boundry_size = size(boundry_index,2);
    boundry_ratio = 100*boundry_size/n
    fprintf('Boundry node size is %d. %.2f %% of nodes are boundry.\n', boundry_size, boundry_ratio);
    
    
    
    for seed = 1:10
        stream = RandStream('mt19937ar', 'Seed', seed);
        terminal_index = datasample(stream, 1:n,floor(sample_rate*n),'Replace',false);

        % Get the index for TUB and V-T-B;
        t_u_b = unique([boundry_index terminal_index]);
        temp = 1:1:n;
        temp(t_u_b) = [];
        non_tub = temp;

        % Get the Terminal index in t_u_b.
        [tf,idx] = ismember(terminal_index, t_u_b);
        terminal_in_tub = idx;

        fprintf('tub size is %d.\n', size(t_u_b,2));

        % ----------------------------------
        % Schur Complement(SC) 
        % ----------------------------------
        SC = Schur_Comp(W, terminal_index);
        SC_filename = [root_folder 'exp' 'rate_' num2str(sample_rate) '_s' num2str(s) '_eps' num2str(epsilon) Type num2str(seed) '_SC.csv'];
        dlmwrite(SC_filename,full(SC));

        

        terminal_tub_filename = [root_folder 'exp' 'rate_' num2str(sample_rate) '_s' num2str(s) '_eps' num2str(epsilon) Type num2str(seed) '_terminal_in_tub.csv'];
        dlmwrite(terminal_tub_filename,terminal_in_tub);

        terminal_filename = [root_folder 'exp' 'rate_' num2str(sample_rate) '_s' num2str(s) '_eps' num2str(epsilon) Type num2str(seed) '_terminal_index.csv'];
        dlmwrite(terminal_filename,terminal_index);

        cost = 0;
        % H2_set and H3_set are both upper triangular.
        H2_set = cell(s,1);
        H2_union = spalloc(size(t_u_b, 2), size(t_u_b, 2), ceil(m/100));

        for i = 1:s
            W_site = W_site_set{i};
            %     nnz(W_site);
            % compute unnormalized Laplacian
            degs_site = sum(W_site, 2);
            degs_site(degs_site == 0) = 1e-9; %1e-9;
            D_site = sparse(1:size(W_site, 1), 1:size(W_site, 2), degs_site);
            L_site = D_site - W_site;

            % get SC
            L_site_vv = L_site(t_u_b, t_u_b);
            L_site_vvp = L_site(t_u_b, non_tub);
            L_site_vpvp = L_site(non_tub, non_tub);
            L_site_vpv = L_site(non_tub, t_u_b);
            SC_site = L_site_vv - L_site_vvp * inv(L_site_vpvp) * L_site_vpv; 


            % From laplacian to get the graph H2_i.
            SC_site = abs(triu(SC_site,1));
            H2_set{i} = SC_site;
    %         H2_union = H2_union + Graph_site + transpose(Graph_site);

            SC_H2_filename = [root_folder 'exp' 'rate_' num2str(sample_rate) '_s' num2str(s) '_eps' num2str(epsilon) Type num2str(seed) '_SCH2_' num2str(i) '.csv'];
            dlmwrite(SC_H2_filename,full(SC_site));     

            % Save graph site to Julia.
            W_site_filename = [root_folder 'exp' 'rate_' num2str(sample_rate) '_s' num2str(s) '_eps' num2str(epsilon) Type num2str(seed) '_Wsite_' num2str(i) '.csv'];
            W_site = triu(W_site,1);
            [x,y,v] = find(W_site);
            sp_W = [x,y,v];
            dlmwrite(W_site_filename, sp_W);
    %         dlmwrite(W_site_filename,full(W_site));     
        end
        fprintf('Finished build H2 set.\n');
        % clear W_site_set;

        T_tub = T(:,t_u_b);
        W_tub = W(t_u_b,t_u_b);
    end
end
succ_record
    
    