# May 31 version
# randomly sample terminal node 5 times for each settings.
# Get the mean and std.

# using Pkg
# Pkg.add("Laplacians")
# Pkg.add("CSV")

using Laplacians
using CSV
using SparseArrays
using Printf
using LinearAlgebra
using Statistics

type = "MP"
test_mode = 3
SimGraphType = 4
Param = 37
Sigma = 5
seed = 15

adjust_eps = 0


if test_mode==4
    n = 500
elseif test_mode == 3
    n = 6600
elseif test_mode == 1
    n = 2000
end


root_folder = "./experiment_results/"
data_folder = string(root_folder, "exp_mode", test_mode, "/type", SimGraphType, "_para", Param, "_sigma", Sigma, "_seed", seed)
cd(data_folder)

# Find out all settings need to run.
d_epsilon = 0.3
d_sample_rate = 0.05
d_s = 5
sample_rate_range = [0.05]#[0.01, 0.05, 0.10, 0.15]
s_range = [5]#[3, 5, 10, 15]
epsilon_range = 0.2:0.1:1.2#[0.2, 0.3, 0.4, 0.5]

all_settings = [0,0,0]
for sample_rate in sample_rate_range
    for s in s_range
        for epsilon in epsilon_range
            setting = [Float32(sample_rate), Float32(s), Float32(epsilon)]
            if (sample_rate != d_sample_rate) & (s != d_s) & (epsilon != d_epsilon)
                continue
            end
            global all_settings = hcat(all_settings, setting)
        end
    end
end

function schur_comp(A, idx, add_diag = 0)
    non_idx = [i for i in 1:size(A,1) if !(i in idx)]

    degree = sum(A, dims = 2) .+ add_diag
    degree = reshape(degree, size(degree,1))
    D = Diagonal(degree)
    #@show size(D)
    #@show size(A)
    L = D - A

    Lvv = L[idx,idx]
    Lvvp = L[idx,non_idx]
    Lvpvp = L[non_idx,non_idx]
    Lvpv = L[non_idx, idx]

    SC = Lvv - Lvvp * inv(Lvpvp) * Lvpv;
    SC = UpperTriangular(SC)
    SC = broadcast(abs, SC)
    SC[diagind(SC)] .= 0
    SC = sparse(SC)
    SC = SC + SC'
    return SC
end

function baseline(SC, file_prefix, s, epsilon, terminal_idx, graph_content,n, adjust_eps = 0)
    cost = 0
    for j in collect(1:s)
        Wsite_file = string(file_prefix, graph_content, Int(j), ".csv")
        W_site = CSV.read(Wsite_file, header=false)
        W_site = Matrix(W_site)
        if graph_content == "_SCH2_"
            W_site = W_site + transpose(W_site)
        else
            x = W_site[:,1]
            y = W_site[:,2]
            v = W_site[:,3]
            W_site = sparse(x,y,v,n,n)
            W_site = W_site + W_site'
            W_site = Matrix(W_site)
        end

        # eliminate the disjoint nodes first.
        row_sum = sum(W_site, dims=2)
        zero_row_idx = [i for i=1:size(W_site,1) if row_sum[i]<1e-4]
        non_zero_idx = [i for i=1:size(W_site,1) if row_sum[i]>=1e-4]
        W_site_select = W_site[non_zero_idx, non_zero_idx]

        # sparsify
        W_site_select = sparse(W_site_select)
        #@show nnz(W_site_select)/2
        W_site_sp = sparsify(W_site_select, ep=epsilon+adjust_eps)

        # map the sparsified graph to original size.
        W_site_back = zeros(size(W_site))
        W_site_back[non_zero_idx,non_zero_idx] = Matrix(W_site_sp)
        # SC_H2_back = sparse(SC_H2_back)

        cost = cost + nnz(W_site_sp)/2
        if j == 1
            global W = W_site_back
        else
            global W = W + W_site_back
        end
    end

    # @show SC_H2_union
    if graph_content == "_SCH2_"
        add_diag = 1e-6
    else
        add_diag = 1e-6
    end
    # edges = nnz(sparse(W))/2
    SC_H3 = schur_comp(W, terminal_idx, add_diag)
    edges = nnz(SC_H3)/2

    eps = 100
    if !isConnected(SC_H3)
        @printf(graph_content, "SC_H3 not connected.\n")

    else
        for i in [1,2,3]
            try_eps = approxQual(SC, SC_H3)
            if try_eps < eps
                eps = try_eps
            end
        end

    end
    return cost, eps, edges
end



@show all_settings
for i in collect(2:size(all_settings,2))
    sample_rate = all_settings[1, i]
    s = Int(all_settings[2, i])
    epsilon = all_settings[3, i]

    # if sample_rate == 0
    #     continue
    # end

    setting_edge_num_1 = []
    setting_edge_num_2 = []
    setting_edge_num_3 = []

    setting_cost = []
    setting_cost_base = []

    setting_eps = []
    setting_eps_base = []

    for j in 1:10
        file_prefix = string("exprate_", sample_rate, "_s", s, "_eps", 0.3, type, j)
        # @show file_prefix
        SC_file = string(file_prefix, "_SC.csv")
        SC = CSV.read(SC_file, header=false)
        SC = sparse(Matrix(SC))
        SC = SC + SC'
        if !isConnected(SC)
            @printf("For setting %f %d %f, SC not connected.\n", sample_rate, s, epsilon)
            @printf("Skip %d.\n", j)
            continue
        end

        edge_num_1 = nnz(SC)/2

        # get the index of terminal, need use them to calculate SC.
        terminal_in_tub_file = string(file_prefix, "_terminal_in_tub.csv")
        terminal_in_tub = Array(CSV.read(terminal_in_tub_file, header=false))
        terminal_in_tub = reshape(terminal_in_tub, size(terminal_in_tub,2))

        terminal_index_file = string(file_prefix, "_terminal_index.csv")
        terminal_index = Array(CSV.read(terminal_index_file, header=false))
        terminal_index = reshape(terminal_index, size(terminal_index,2))

        # our method and baseline method
        #cost, eps = g(SC, file_prefix,s, epsilon, terminal_in_tub)
        cost, eps, edge_num_2 = baseline(SC, file_prefix,s, epsilon, terminal_in_tub, "_SCH2_",n, adjust_eps)
        cost_base, eps_base, edge_num_3 = baseline(SC, file_prefix,s, epsilon, terminal_index, "_Wsite_",n, adjust_eps)

        # print out results
        # @printf("In setting sample_rate:%f s:%d epsilon:%f, has SC edge_num: %d\n", sample_rate, s, epsilon, edge_num_1)
        # @printf("SC method get cost %d, eps %f, edge_num %d\n", cost, eps, edge_num_2)
        # @printf("Baseline get cost %d, eps %f, edge_num %d\n", cost_base, eps_base, edge_num_3)

        setting_edge_num_1 = append!(setting_edge_num_1, edge_num_1)
        setting_edge_num_2 = append!(setting_edge_num_2, edge_num_2)
        setting_edge_num_3 = append!(setting_edge_num_3, edge_num_3)

        setting_cost = append!(setting_cost, cost)
        setting_cost_base = append!(setting_cost_base, cost_base)

        setting_eps = append!(setting_eps, eps)
        setting_eps_base = append!(setting_eps_base, eps_base)

        # @show size(setting_edge_num_1)[1]
        if size(setting_edge_num_1)[1] == 5
            break
        end
    end
    @printf("-----------------------------------------------------\n")
    @printf("In setting sample_rate:%f s:%d epsilon:%f\n", sample_rate, s, epsilon)


    @printf("The Cost for LocalSS are: ")
    for j = 1:5
        @printf("%d, ", setting_cost_base[j])
    end
    @printf(". The mean: %f, std: %f\n", mean(setting_cost_base), std(setting_cost_base))

    @printf("The Cost for LocalSC are: ")
    for j = 1:5
        @printf("%d, ", setting_cost[j])
    end
    @printf(". The mean: %f, std: %f\n", mean(setting_cost), std(setting_cost))


    @printf("The Quality for LocalSS are: ")
    for j = 1:5
        @printf("%f, ", setting_eps_base[j])
    end
    @printf(". The mean: %f, std: %f\n", mean(setting_eps_base), std(setting_eps_base))

    @printf("The Quality for LocalSC are: ")
    for j = 1:5
        @printf("%f, ", setting_eps[j])
    end
    @printf(". The mean: %f, std: %f\n", mean(setting_eps), std(setting_eps))




    @printf("The edges_num for original SC are: ")
    for j = 1:5
        @printf("%d, ", setting_edge_num_1[j])
    end
    @printf(". The mean: %f, std: %f\n", mean(setting_edge_num_1), std(setting_edge_num_1))

    @printf("The edges_num for LocalSS are: ")
    for j = 1:5
        @printf("%d, ", setting_edge_num_3[j])
    end
    @printf(". The mean: %f, std: %f\n", mean(setting_edge_num_3), std(setting_edge_num_3))

    @printf("The edges_num for LocalSC are: ")
    for j = 1:5
        @printf("%d, ", setting_edge_num_2[j])
    end
    @printf(". The mean: %f, std: %f\n", mean(setting_edge_num_2), std(setting_edge_num_2))

    # @printf("-----------------------------------------------------\n")
end
