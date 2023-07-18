
using GraphRecipes, Plots
using Test
using CairoMakie

@testset "2 modules" begin
    A1 = [0 1 0 1; 0 0 0 1; 0 1 0 0; 0 0 1 0]
    A2 = [0 0 1 1 0; 1 0 0 1 1; 1 1 0 0 0; 0 1 0 0 1; 1 0 1 1 0]
    A1 = convert_int_to_symbol(A1)
    A2 = convert_int_to_symbol(A2)
    
    As = [A1, A2]
    num_connections = [0 4; 4 0]
    inter_modular_specifier = Dict("num_connections"=>num_connections)
    A3 = connect_modules(As, inter_modular_specifier)
    @test A3[1:4, 1:4] == A1 
    @test A3[5:end, 5:end] == A2

    l1=count_connections(A3[1:4, 5:end])
    l2=count_connections(A3[5:end, 1:4])
    @test l1+l2 == sum(num_connections)
    fig, axs = plot_adjacency_matrix(A3, As); fig
end


@testset "3 modules" begin
    A1 = [0 1 0 1; 0 0 0 1; 0 1 0 0; 0 0 1 0]
    A2 = [0 0 1 1 0; 1 0 0 1 1; 1 1 0 0 0; 0 1 0 0 1; 1 0 1 1 0]
    A3 = [0 1 0 0 0; 0 0 1 0 1; 0 1 0 0 1; 1 1 0 0 1; 1 1 1 1 0]
    A1 = convert_int_to_symbol(A1)
    A2 = convert_int_to_symbol(A2)
    A3 = convert_int_to_symbol(A3)

    As = [A1, A2, A3]
    num_connections = [0 4 2; 4 0 1; 0 0 0]
    
    inter_modular_specifier = Dict("num_connections"=>num_connections)
    Afull = connect_modules(As, inter_modular_specifier)
    @test  Afull[1:4, 1:4] == A1 
    @test  Afull[5:9, 5:9] == A2
    @test  Afull[10:14, 10:14] == A3
    
    idxs_modules = [1:4, 5:9, 10:14]
    ls = zeros(Int64, (3,3))
    for i=1:3 
        for j=1:3
            if i == j continue  end
            ls[i,j]=count_connections(Afull[idxs_modules[i], idxs_modules[j]])
        end
    end
    @test ls == num_connections
    ls 
    num_connections
    

    fig, axs = plot_adjacency_matrix(Afull, As); fig
end

@testset "2 modules I and E" begin
    A1 = [0 :E 0 :E 0; 0 0 0 :I 0; 0 :I 0 0 0; 0 0 :E 0 0; :I :I :I 0 0]
    A2 = [0 :E 0 :E 0; 0 0 0 :I 0; 0 :I 0 0 0; 0 0 :E 0 0; :E :E :E 0 0]
    A1 = convert_int_to_symbol(A1)
    A2 = convert_int_to_symbol(A2)
    
    As = [A1, A2]
    num_connections = [0 2; 2 0]
    inter_modular_specifier = Dict("num_connections"=>num_connections)
    A3 = connect_modules(As, inter_modular_specifier, _decide_source_and_targets_only_E_to_I)
    @test A3[1:5, 1:5] == A1 
    @test A3[6:end, 6:end] == A2
    @test unique(A3[1:5, 6:end]) == [:O, :E]
    @test unique(A3[6:end, 1:5]) == [:O, :E]

    l1=count_connections(A3[1:5, 6:end])
    l2=count_connections(A3[6:end, 1:5])
    @test l1+l2 == sum(num_connections)
    fig, axs = plot_adjacency_matrix(A3, As); #white is exc, black is inh
    fig

    (N,N) = size(A3)
    node_names = 1:N
    node_types = identify_node_types(A3)
    node_color = replace(node_types, :E=>"red", :I=>"blue") 
    edgelabels = replace(A3, :O=>nothing, :E=>"E", :I=>"I")
    _A3 = convert_symbol_to_int(A3)
    graphplot(_A3, names=node_names, edgelabel=edgelabels, nodeshape=:rect, nodecolor=node_color)
end

function plot_graph(A)
    (N,N) = size(A)
    node_names = 1:N
    node_types = identify_node_types(A)
    node_color = replace(node_types, :E=>"red", :I=>"blue") 
    edgelabels = replace(A, :O=>nothing, :E=>"E", :I=>"I")
    _A = convert_symbol_to_int(A)
    graphplot(_A, names=node_names, edgelabel=edgelabels, nodeshape=:rect, nodecolor=node_color)
end

@testset "conversions" begin 

    A3 = [:O :E :O :E :O :O :O :E :O :O; :O :O :O :I :O :O :O :O :O :O; :O :I :O :O :O :O :O :O :O :O; :O :O :E :O :O :O :O :E :O :O; :I :I :I :O :O :O :O :O :O :O; :O :E :O :O :O :O :E :O :E :O; :O :O :O :O :O :O :O :O :I :O; :O :O :O :O :O :O :I :O :O :O; :O :E :O :O :O :O :O :E :O :O; :O :O :O :O :O :E :E :E :O :O]
    connsl = adjmat_to_connsl(A3)
end


@testset "E-E && interface" begin 
    dfull =  Dict( #Basic test of the effect of a single neuron sending E and I to other neurons
    #sim details
    "ttrans" => (0, "SP4"),
    "tend" => (600, "SP4"),
    "Δt" => (0.1, "SP2"),
    # "solver" => ("autotsi5trbdf2", "SP2"),
    "solver" => ("tsit5", "SP2"),
    #unit dynamicssol[1,
    "unitm" => ("FHN", "SP1"),
    "I" => ([[0.0]], "SP2"),
    "a" => ([[-0.5]],"SP2"),
    "b" => ([[0.1]],"SP2"),
    "c" => ([[0.2]],"SP2"),
    "d" => ([[0.1]],"SP2"),
    #coupling
    "ϵ" => ([0.0, 5.0, 7.5, 10.0], "CP"),
    # "gmax_exc" => (+0.01, "SP2"), #DEFAULT
    "gmax_exc" => (+0.05, "SP2"),
    "gmax_inh" => (0.05, "SP2"),
    "E_inh" =>  (-1.5, "SP2"),
    "E_exc" =>  (1.9, "SP2"),
    "τs_exc" => (2.0, "SP2"),  #DEFAULT
    "τs_inh" => (1.0, "SP2"), #DEFAULT
    "numstoredspiketimes" => (5, "SP2"),
    "Vth" => (0.5, "SP2"),
    #topology
    "N" => (10, "SP1"),
    "Nmod" => (5, "SP1"),
    "k" => (2, "SP1"),
    "topseed" => (1, "SP4"),
    "topm" => ("modular", "SP0"),
    "modules_topologies" => ([["ER", "ER"]], "SP3"),
    "intermodseed" => (1, "SP4"),
    "num_connections" => ([[0 2; 3 0]], "SP4"),
    "connection_decider" =>  (_decide_source_and_targets_uniform, "SP1"),
    "EIseed" =>  (0, "SP1"),
    "types" => ([[:E, :I]], "SP3"),
    "probabilities" => ([[1.0, 0.0]], "SP3"),
    #ics
    "ictype" =>  ("uniform", "SP1"),
    "icmin" =>  (-1.5, "IG"),
    "icmax" =>  (+1.5, "IG"),
    "icseed" => (1, "IG"),
);
pvals_all, ptypes = separatedicts(dfull);
pvals_list = dict_list(pvals_all);
pvals = pvals_list[1]
connsl = get_connsl(pvals)
A = connsl_to_adjmat(connsl)

plot_graph(A)
end