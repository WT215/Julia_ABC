function Data_sim(;Param_Space,BETA_vec)

    num_cells=size(BETA_vec)[1]
    aa=rand.(Beta.(Param_Space[:,1],Param_Space[:,2]),num_cells)
    aa=transpose(reduce(hcat,aa))
    aa2=rand.(Poisson.(aa.*Param_Space[:,3]))
    sim_dat=transpose(rand.(Binomial.(transpose(aa2),BETA_vec)))
    return sim_dat
end

function Data_sim_sum2(;Param_Space,BETA_vec)

    num_cells=size(BETA_vec)[1]
    aa=rand.(Beta.(Param_Space[:,1],Param_Space[:,2]),num_cells)
    aa=transpose(reduce(hcat,aa))
    aa2=rand.(Poisson.(aa.*Param_Space[:,3]))

    aa_2=rand.(Beta.(Param_Space[:,1],Param_Space[:,2]),num_cells)
    aa_2=transpose(reduce(hcat,aa_2))
    aa2_2=rand.(Poisson.(aa_2.*Param_Space[:,3]))
    aa_sum=aa2 .+ aa2_2
    
    sim_dat=transpose(rand.(Binomial.(transpose(aa_sum),BETA_vec)))
    
    return sim_dat
end




function ABC_sim(;num_samples=1)
    # sim_kon=10 .^ rand(Uniform(-2, 1),num_samples)
    # sim_inter=rand(Uniform(0, 1),num_samples)
    # sim_koff=sim_kon ./ sim_inter .- sim_kon
    # mRNA=exp.(rand(Uniform(log(5), log(200)),num_samples)) ./ 0.1
    # sim_ksyn=mRNA .* (sim_kon .+ sim_koff) ./ sim_kon


    #version 2
    # sim_kon=10 .^ rand(Uniform(-2, 1),num_samples)
    # sim_koff=10 .^ rand(Uniform(-2, 1),num_samples)
    # mRNA=10 .^ (rand(Uniform(0, 3),num_samples)) 
    # sim_inter=10 .^ rand(Uniform(0, 1.5),num_samples)
    # sim_ksyn=(sim_inter .- 1) .* (sim_kon .+ sim_koff) .* (sim_kon .+ sim_koff .+ 1) ./ sim_koff

    #version3
    # sim_kon=10 .^ rand(Uniform(-2, 1),num_samples)
    # sim_inter=rand(Uniform(0, 0.5),num_samples)
    # sim_koff=sim_kon ./ sim_inter .- sim_kon
    # sim_inter=10 .^ rand(Uniform(0.0005, 1.5),num_samples)
    # sim_ksyn=(sim_inter .- 1) .* (sim_kon .+ sim_koff) .* (sim_kon .+ sim_koff .+ 1) ./ sim_koff

    #version4
    # sim_kon=10 .^ rand(Uniform(-2, 1),num_samples)
    # qr=sim_kon ./ (0.01 .+sim_kon)
    # sim_inter=rand.(Uniform.(0,qr ),1)
    # sim_inter=reduce(vcat,sim_inter)
    # sim_koff=sim_kon ./ sim_inter .- sim_kon
    # sim_fano=10 .^ rand(Uniform(log10(1.001), log10(30)),num_samples)
    # sim_ksyn=(sim_fano .- 1) .* (sim_kon .+ sim_koff) .* (sim_kon .+ sim_koff .+ 1) ./ sim_koff

    #version5
    # sim_kon=10 .^ rand(Uniform(-2, 2),num_samples)
    # sim_koff=10 .^ rand(Uniform(-2, 2),num_samples)
    # sim_fano=10 .^ rand(Uniform(log10(1.001), log10(30)),num_samples)
    # sim_ksyn=(sim_fano .- 1) .* (sim_kon .+ sim_koff) .* (sim_kon .+ sim_koff .+ 1) ./ sim_koff

    #version6
    sim_kon=10 .^ rand(Uniform(log10(0.01), log10(100)),num_samples)
    #sim_inter=10 .^ rand(Uniform(log10(0.01), log10(100)),num_samples)
    sim_inter=abs.(rand(Normal(0.05, 0.5),num_samples))
    sim_koff=sim_kon ./ sim_inter
    sim_fano=10 .^ rand(Uniform(log10(1.001), log10(30)),num_samples)
    sim_ksyn=(sim_fano .- 1) .* (sim_kon .+ sim_koff) .* (sim_kon .+ sim_koff .+ 1) ./ sim_koff

    new_param_space=[sim_kon sim_koff sim_ksyn]


    return new_param_space
end


function ABC_sim_KP(;population_dat,num_samples=1)

    leftt=minimum(population_dat[:,1],dims=1)
    rightt=maximum(population_dat[:,1],dims=1)

    leftt_inter=minimum(population_dat[:,1] ./population_dat[:,2],dims=1)
    rightt_inter=maximum(population_dat[:,1]./population_dat[:,2],dims=1)
    fanotemp=population_dat[:,3] .* population_dat[:,2] ./ (population_dat[:,1].+population_dat[:,2])./ (population_dat[:,1].+population_dat[:,2] .+1) .+1 
    leftt_fanotemp=minimum(fanotemp)
    rightt_fanotemp=maximum(fanotemp)

    sim_kon=10 .^ rand(Uniform(log10.(leftt)[1], log10.(rightt)[1]),num_samples)
    sim_inter=10 .^ rand(Uniform(log10.(leftt_inter)[1], log10.(rightt_inter)[1]),num_samples)
    sim_koff=sim_kon ./ sim_inter
    sim_fano=10 .^ rand(Uniform(log10(leftt_fanotemp), log10(rightt_fanotemp)),num_samples)
    sim_ksyn=(sim_fano .- 1) .* (sim_kon .+ sim_koff) .* (sim_kon .+ sim_koff .+ 1) ./ sim_koff

    new_param_space=[sim_kon sim_koff sim_ksyn]


    return new_param_space
end


function param_check(;params)
    sim_kon=params[1,1]
    sim_koff=params[1,2]
    sim_ksyn=params[1,3]

    variancee=sim_ksyn .* sim_kon ./ (sim_kon .+ sim_koff) .+ sim_ksyn.^2 .*sim_kon .*sim_koff ./( (sim_kon .+sim_koff) .^2 .* (sim_kon.+sim_koff .+1)) 
    meann=sim_ksyn .* sim_kon ./ (sim_kon .+ sim_koff)
    vmratio=variancee ./ meann
    return(vmratio)
end

















