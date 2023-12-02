function ABC_OneParam_test(sim_vec, true_vec)
    freq_1=freqtable(true_vec)
    freq_2= freqtable(sim_vec)
    lenn=maximum([size(freq_1)[1],size(freq_2)[1]])
    x_1=zeros(lenn)
    x_2=zeros(lenn)
    x_1[1:size(freq_1)[1]]=freq_1
    x_2[1:size(freq_2)[1]]=freq_2
    
    x_1=x_1[:,1] ./ sum(x_1[:,1])
    x_2=x_2[:,1] ./ sum(x_2[:,1])
    out_metric=hellinger2(x_1,x_2)
    return out_metric
end


function ABC_OneParam_test_bins(sim_vec,true_vec)
    tempcom=vcat(true_vec,sim_vec)
    ncells=length(true_vec)


    #basevec=tempcom
    basevec=true_vec

    h=2*iqr(basevec)/(ncells^(1/3))
    extt=extrema(basevec)
    #Freedman–Diaconis rule
    if floor((extt[2]-extt[1])/h)== Inf
        #nbins=Int(size(true_vec)[1]/10)
        nbins=Int(ceil(size(true_vec)[1]^(1/2)))
    else
        nbins=Int(ceil((extt[2]-extt[1])/h))
    end

    #nbins=Int(size(true_vec)[1]/10)
    bins = range(extt[1],stop=extt[2], length=nbins)
    
    h_1 = fit(Histogram, true_vec,bins)
    x_1=h_1.weights
    
    #sim_vec=sim_fun(Param_Space=reshape(TRUE_params[geneind,:],1,3),BETA_vec=BETA_vec)[1,:]
    h_2 = fit(Histogram, sim_vec,bins)
    x_2=h_2.weights
    
    #Hellinger distance
    out_metric=hellinger2(x_1,x_2)

    #Wasserstain distance
    # C = pairwise(SqEuclidean(), x_2', x_1');
    # C=convert(Array{Float64}, C)
    # p2=x_2./sum(x_2)
    # p1=x_1./sum(x_1)
    # out_metric=Float64(1.)
    # try 
    #     out_metric=OptimalTransport.emd2(p2,p1, C,Tulip.Optimizer())
    # catch e
    #     out_metric=Inf
    # end


    return out_metric

end





function ABC_test(sim_dat, true_vec)

    store_vec=zeros(size(sim_dat)[1])
    for i in 1:size(sim_dat)[1]
        #print(i)
        #change to binning
        #temp=ABC_OneParam_test(sim_dat[i,:], true_vec)
        temp=ABC_OneParam_test_bins(sim_dat[i,:], true_vec)
        store_vec[i]=temp
    end
    return store_vec
end

#not use
function ABC_OneGene_test(;Data,BETA_vec,gene_ind=1, loopsnum=2000,threshold=0.2)

    aa=@distributed (vcat) for i in 1:loopsnum
    
        new_param_space=ABC_sim(num_samples=1)
    
        input_obs=Data[gene_ind,:]
        sim_vec=rand.(Binomial.(rbp_self(cell_numbers,new_param_space[1,1],new_param_space[1,2],new_param_space[1,3]),BETA_vec))
        tout=ABC_OneParam_test(sim_vec,input_obs)
        ff=[new_param_space tout]
    
        # if tout<0.2
        #      ff=[new_param_space tout]
        #      #push!(store_dat,ff)
        #      #store_dat=vcat(store_dat,ff)
        # end
        ff
    end
    
    aa2=aa[findall(aa[:,4] .<threshold),:]
    return(aa2)

end



@everywhere function param_to_mean(;x,meanbeta)
    tempp=x[3]*(x[1]/(x[1]+x[2])) *meanbeta
    return (tempp)
end





