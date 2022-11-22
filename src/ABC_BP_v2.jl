
function ABC_initialization_v2(; Data,BETA_vec,num_boot=100,threshold_on_distance=0.95,mode_allele=true,verbose=false)
    scaledData= transpose(transpose(Data) ./ BETA_vec)
    n_genes=size(Data)[1]
    n_cells=size(Data)[2]
    ww=fill(1/size(Data)[2], size(Data)[2])

    if mode_allele==true
        #print("Running ABC based on allele counts.")
        MME_fun=MomentInference_allele
        sim_fun=Data_sim
    else
        #print("Running ABC based on UMI counts.")
        MME_fun=MomentInference_UMI
        sim_fun=Data_sim_sum2
    end

    MME_ori=MME_fun(Data=scaledData,BETA_vec=nothing)

    M1_boot=SharedArray{Float64,2}(n_genes, num_boot)
    M2_boot=SharedArray{Float64,2}(n_genes, num_boot)
    M3_boot=SharedArray{Float64,2}(n_genes, num_boot)
    Dist_boot=SharedArray{Float64,2}(n_genes, num_boot)

    if verbose==true
        progress = Progress(num_boot, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
        channel = RemoteChannel(()->Channel{Bool}(num_boot), 1)

        @sync begin
            # this task prints the progress bar
            @async while take!(channel)
                next!(progress)
            end

            # this task does the computation
            @async begin
                @distributed vcat for boot_ind in range(1, stop = num_boot)
                    colindex=wsample(range(1,stop=size(Data)[2]),ww,size(Data)[2],replace=true)
                    bootdat=Data[:,colindex]
                    scaledbootdat=scaledData[:,colindex]       
                    
                    dist_genes=hellinger2_data(Data, bootdat)
                    MME_out=MME_fun(Data=scaledbootdat,BETA_vec=nothing)
            
                    M1_boot[:,boot_ind]=MME_out[:,1]
                    M2_boot[:,boot_ind]=MME_out[:,2]
                    M3_boot[:,boot_ind]=MME_out[:,3]
                    Dist_boot[:,boot_ind]=dist_genes
                    put!(channel, true)
                end # end of each gene (end of distributed)
            put!(channel, false) # this tells the printing task to finish
            end
        end #end of begin for progress bar
    else
        @distributed vcat for boot_ind in range(1, stop = num_boot)
            colindex=wsample(range(1,stop=size(Data)[2]),ww,size(Data)[2],replace=true)
            bootdat=Data[:,colindex]
            scaledbootdat=scaledData[:,colindex]       
            
            dist_genes=hellinger2_data(Data, bootdat)
            MME_out=MME_fun(Data=scaledbootdat,BETA_vec=nothing)
    
            M1_boot[:,boot_ind]=MME_out[:,1]
            M2_boot[:,boot_ind]=MME_out[:,2]
            M3_boot[:,boot_ind]=MME_out[:,3]
            Dist_boot[:,boot_ind]=dist_genes
            #put!(channel, true)
        end # end of each gene (end of distributed)

    end

    interval_1=[mapslices(X -> quantile_fun(X, 0.025), M1_boot, dims=(2)) mapslices(X -> quantile_fun(X, 0.975), M1_boot, dims=(2))]
    interval_2=[mapslices(X -> quantile_fun(X, 0.025), M2_boot, dims=(2)) mapslices(X -> quantile_fun(X, 0.975), M2_boot, dims=(2))]
    interval_3=[mapslices(X -> quantile_fun(X, 0.025), M3_boot, dims=(2)) mapslices(X -> quantile_fun(X, 0.975), M3_boot, dims=(2))]

    threshold_data=mapslices(X -> quantile(X, threshold_on_distance), Dist_boot, dims=(2))

    return threshold_data, interval_1, interval_2, interval_3, M1_boot, M2_boot, M3_boot
end



function ABC_initialization_v3(; Data,BETA_vec,num_boot=100,threshold_on_distance=0.95,mode_allele=true,verbose=false)
    scaledData= transpose(transpose(Data) ./ BETA_vec)
    n_genes=size(Data)[1]
    n_cells=size(Data)[2]
    ww=fill(1/size(Data)[2], size(Data)[2])

    if mode_allele==true
        #print("Running ABC based on allele counts.")
        MME_fun=MomentInference_allele
        sim_fun=Data_sim
    else
        #print("Running ABC based on UMI counts.")
        MME_fun=MomentInference_UMI
        sim_fun=Data_sim_sum2
    end

    MME_ori=MME_fun(Data=scaledData,BETA_vec=nothing)

    M1_boot=SharedArray{Float64,2}(n_genes, num_boot)
    M2_boot=SharedArray{Float64,2}(n_genes, num_boot)
    M3_boot=SharedArray{Float64,2}(n_genes, num_boot)
    Dist_boot=SharedArray{Float64,2}(n_genes, num_boot)

    if verbose==true
        progress = Progress(num_boot, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
        channel = RemoteChannel(()->Channel{Bool}(num_boot), 1)

        @sync begin
            # this task prints the progress bar
            @async while take!(channel)
                next!(progress)
            end

            # this task does the computation
            @async begin
                @distributed vcat for boot_ind in range(1, stop = num_boot)
                    colindex=wsample(range(1,stop=size(Data)[2]),ww,size(Data)[2],replace=true)
                    bootdat=Data[:,colindex]
                    scaledbootdat=scaledData[:,colindex]       
                    
                    dist_genes=hellinger2_data(Data, bootdat)
                    MME_out=MME_fun(Data=scaledbootdat,BETA_vec=nothing)
            
                    M1_boot[:,boot_ind]=MME_out[:,1]
                    M2_boot[:,boot_ind]=MME_out[:,2]
                    M3_boot[:,boot_ind]=MME_out[:,3]
                    Dist_boot[:,boot_ind]=dist_genes
                    put!(channel, true)
                end # end of each gene (end of distributed)
            put!(channel, false) # this tells the printing task to finish
            end
        end #end of begin for progress bar
    else
        @distributed vcat for boot_ind in range(1, stop = num_boot)
            colindex=wsample(range(1,stop=size(Data)[2]),ww,size(Data)[2],replace=true)
            bootdat=Data[:,colindex]
            scaledbootdat=scaledData[:,colindex]       
            
            dist_genes=hellinger2_data(Data, bootdat)
            MME_out=MME_fun(Data=scaledbootdat,BETA_vec=nothing)
    
            M1_boot[:,boot_ind]=MME_out[:,1]
            M2_boot[:,boot_ind]=MME_out[:,2]
            M3_boot[:,boot_ind]=MME_out[:,3]
            Dist_boot[:,boot_ind]=dist_genes
            #put!(channel, true)
        end # end of each gene (end of distributed)

    end

    # interval_1=[mapslices(X -> quantile_fun(X, 0.05), M1_boot, dims=(2)) mapslices(X -> quantile_fun(X, 0.95), M1_boot, dims=(2))]
    # interval_2=[mapslices(X -> quantile_fun(X, 0.05), M2_boot, dims=(2)) mapslices(X -> quantile_fun(X, 0.95), M2_boot, dims=(2))]
    # interval_3=[mapslices(X -> quantile_fun(X, 0.05), M3_boot, dims=(2)) mapslices(X -> quantile_fun(X, 0.95), M3_boot, dims=(2))]
    # interval_1=[minimum(M1_boot,dims=2) maximum(M1_boot,dims=2)]
    # interval_2=[minimum(M2_boot,dims=2) maximum(M2_boot,dims=2)]
    # interval_3=[minimum(M3_boot,dims=2) maximum(M3_boot,dims=2)]
    
    interval_1=[mapslices(X -> minimum_fun(X), M1_boot, dims=(2)) mapslices(X -> maximum_fun(X), M1_boot, dims=(2))]
    interval_2=[mapslices(X -> minimum_fun(X), M2_boot, dims=(2)) mapslices(X -> maximum_fun(X), M2_boot, dims=(2))]
    interval_3=[mapslices(X -> minimum_fun(X), M3_boot, dims=(2)) mapslices(X -> maximum_fun(X), M3_boot, dims=(2))]


    threshold_data=mapslices(X -> quantile(X, threshold_on_distance), Dist_boot, dims=(2))

    return threshold_data, interval_1, interval_2, interval_3, M1_boot, M2_boot, M3_boot
end


function extract_pointInterval_v3(;MME_ori, interval_1, interval_2, interval_3, M1_boot, M2_boot,M3_boot,gene_ind,particles=1000,mode_allele=true)

    if  ((interval_1[gene_ind,1] == interval_1[gene_ind,2])  ||  (interval_2[gene_ind,1] == interval_2[gene_ind,2])  ||  (interval_3[gene_ind,1] == interval_3[gene_ind,2] ))

        new_param_space=[]
    else


        qq1_interval = (interval_1[gene_ind,:]);
        qq2_interval = (interval_2[gene_ind,:]);
        qq3_interval = (interval_3[gene_ind,:]);

        if  ((qq1_interval[1] >= qq1_interval[2])  || (qq2_interval[1] >= qq2_interval[2])  ||  (qq3_interval[1] >= qq3_interval[2]) || isnan(qq1_interval[1]) || isnan(qq1_interval[2]) || isnan(qq2_interval[1]) || isnan(qq2_interval[2]) || isnan(qq3_interval[1]) || isnan(qq3_interval[2])  || qq1_interval[1]<0 ||   qq2_interval[1]<0  || qq3_interval[1]<0   || qq1_interval[2]<0 ||   qq2_interval[2]<0  || qq3_interval[2]<0)

            new_param_space=[]

        else
        
            qq1=rand.(Uniform.(log10(qq1_interval[1]/10),        log10(qq1_interval[2]*10)),particles)
            qq2=rand.(Uniform.(log10(qq2_interval[1]/10),        log10(qq2_interval[2]*10)),particles)
            qq3=rand.(Uniform.(log10(qq3_interval[1]/10),        log10(qq3_interval[2]*10)),particles)

            new_param_space=10 .^ [qq1 qq2 qq3]
    

            
        end
    end
    return new_param_space
end




function extract_pointInterval_fano(;Data,MME_ori, interval_1, interval_2, interval_3, M1_boot, M2_boot,M3_boot,gene_ind,particles=1000,mode_allele=true)

    if  ((interval_1[gene_ind,1] == interval_1[gene_ind,2])  ||  (interval_2[gene_ind,1] == interval_2[gene_ind,2])  ||  (interval_3[gene_ind,1] == interval_3[gene_ind,2] ))

        new_param_space=[]
    else


        qq1_interval = (interval_1[gene_ind,:]);
        qq2_interval = (interval_2[gene_ind,:]);
        qq3_interval = (interval_3[gene_ind,:]);

        sim_ratio=(MME_ori[:,2] ./MME_ori[:,1])[gene_ind]
        sim_fano=(var(Data,dims=2) ./ mean(Data,dims=2))[gene_ind]

        if  ((qq1_interval[1] >= qq1_interval[2])  || (qq2_interval[1] >= qq2_interval[2])  ||  (qq3_interval[1] >= qq3_interval[2]) || isnan(qq1_interval[1]) || isnan(qq1_interval[2]) || isnan(qq2_interval[1]) || isnan(qq2_interval[2]) || isnan(qq3_interval[1]) || isnan(qq3_interval[2])  || qq1_interval[1]<0 ||   qq2_interval[1]<0  || qq3_interval[1]<0   || qq1_interval[2]<0 ||   qq2_interval[2]<0  || qq3_interval[2]<0) || isnan(sim_ratio)

            new_param_space=[]

        else

            #(sim_fano .- 1) .* (qq1_interval .+ qq2_interval) .* (qq1_interval .+ qq2_interval .+ 1) ./ qq2_interval
        


            qq1=rand.(Uniform.(log10(qq1_interval[1]/10),        log10(qq1_interval[2]*10)),particles)
            #qq2=rand.(Uniform.(log10(qq2_interval[1]/10),        log10(qq2_interval[2]*10)),particles)
            qq2=rand.(Uniform.(log10(qq1_interval[1] .*sim_ratio),        log10(qq1_interval[2] .*sim_ratio)),particles)
            qq3=rand.(Uniform.(log10(qq3_interval[1]/10),        log10(qq3_interval[2]*10)),particles)



   

            # ii=109
            # (TRUE_params[:,3] .* TRUE_params[:,2] ./ (TRUE_params[:,1] .+ TRUE_params[:,2] .+1) ./ (TRUE_params[:,1] .+ TRUE_params[:,2]))[ii]
            # (var(Data,dims=2) ./ mean(Data,dims=2))[ii]
            # (var(scaledData,dims=2) ./ mean(scaledData,dims=2))[ii]

            # ii=2
            # (TRUE_params[:,2] ./TRUE_params[:,1])[ii]
            # (MME_out[:,2] ./MME_out[:,1])[ii]
            # scatter((TRUE_params[:,2] ./TRUE_params[:,1]),(MME_out[:,2] ./MME_out[:,1]))
            # Plots.abline!(1, 0, line=:dash)



            new_param_space=10 .^ [qq1 qq2 qq3]
            #new_param_space[:,2]=new_param_space[:,1] .* sim_ratio

            sim_kon=new_param_space[:,1]
            sim_koff=new_param_space[:,2]

            new_param_space[:,3]=(sim_fano .- 1) .* (sim_kon .+ sim_koff) .* (sim_kon .+ sim_koff .+ 1) ./ sim_koff

            
        end
    end
    return new_param_space
end






function ABC_BP_v2(; Data,BETA_vec, mode_allele=true, particles=10000,threshold_on_distance=0.95,verbose=false)
    n_genes=size(Data)[1]
    n_cells=size(Data)[2]

    #calculate threshold
    #threshold_Data=ABC_threshold(Data=Data, particles=1000,threshold_on_distance=threshold_on_distance,verbose=verbose)
    threshold_Data, interval_1, interval_2, interval_3, M1_boot, M2_boot, M3_boot= ABC_initialization_v2(Data=Data,BETA_vec=BETA_vec,num_boot=1000,threshold_on_distance=threshold_on_distance,mode_allele=mode_allele,verbose=verbose)
    

    if mode_allele==true
        print("Running ABC based on allele counts.")
        #MME_out=MomentInference_data(Data=Data,BETA_vec=BETA_vec)
        MME_out=MomentInference_allele(Data=Data,BETA_vec=BETA_vec)
        sim_fun=Data_sim
    else
        print("Running ABC based on UMI counts.")
        #MME_out=MomentInference_sum2alleles_data(Data=Data,BETA_vec=BETA_vec)
        MME_out=MomentInference_UMI(Data=Data,BETA_vec=BETA_vec)
        sim_fun=Data_sim_sum2
    end
    @everywhere MME_out=$MME_out
    @everywhere particles=$particles

    store_dat_all_median_3 =SharedArray{Float64,2}(n_genes, 6)


    if verbose==true
        progress = Progress(n_genes, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
        channel = RemoteChannel(()->Channel{Bool}(n_genes), 1)


        @sync begin
            # this task prints the progress bar
            @async while take!(channel)
                next!(progress)
            end
        
            # this task does the computation
            @async begin
                @distributed vcat for gene_ind in range(1, stop = n_genes)
                    #@distributed vcat for gene_ind in range(1, stop = size(Data)[1])
                    #for gene_ind in range(1, stop = size(Data)[1])
                        #thresholdd=threshold_Data[gene_ind,1]
                        cellsused=findall(.!ismissing.(Data[gene_ind,:]))
                    
                        #new_param_space=extract_pointInterval_v3(;MME_ori=MME_out,interval_1= interval_1,interval_2= interval_2,interval_3= interval_3,M1_boot=M1_boot, M2_boot=M2_boot,M3_boot=M3_boot,gene_ind=gene_ind,particles=particles,mode_allele=mode_allele)
                        new_param_space=extract_pointInterval_fano(;Data=Data,MME_ori=MME_out,interval_1= interval_1,interval_2= interval_2,interval_3= interval_3,M1_boot=M1_boot, M2_boot=M2_boot,M3_boot=M3_boot,gene_ind=gene_ind,particles=particles,mode_allele=mode_allele)


                        if size(new_param_space)[1]==0
                            preselect=[]
                        else
                            preselect=findall( (new_param_space[:,1] .>0) .& (new_param_space[:,2] .>0) .& (new_param_space[:,3] .>0))
                        end

                        
                        if length(preselect)>1000
                            ww=fill(1/length(preselect), length(preselect))
                            ranind=wsample(range(1,stop=length(preselect)),ww,1000,replace=false)
                            preselect=preselect[ranind]
                        end

                        if length(preselect)==0
                            store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
                        else
                            new_param_space=new_param_space[preselect,:]
                            temp_dat=sim_fun(Param_Space=new_param_space,BETA_vec=BETA_vec[cellsused])
                    
                            ff=ABC_test(temp_dat, Data[gene_ind,:])
                            thresholdd=quantile(ff,0.01)
                            particles_selected=findall( (ff .<= thresholdd) )


                            population_dat_3=new_param_space[particles_selected,:]


                            accepted_particles=size(population_dat_3)[1]
                            if accepted_particles==0
                            final_estimates=[NaN NaN NaN]
                            else
                            final_estimates=extract_ml_particle_wt(transpose(population_dat_3))
                            end
        
                            store_dat_all_median_3[gene_ind,:] = [final_estimates[1] final_estimates[2] final_estimates[3] thresholdd accepted_particles particles] 

                        end
                        put!(channel, true)
                    end # end of each gene (end of distributed)
                put!(channel, false) # this tells the printing task to finish
            end
        end #end of begin for progress bar

    else

        @distributed vcat for gene_ind in range(1, stop = n_genes)
            #@distributed vcat for gene_ind in range(1, stop = size(Data)[1])
            #for gene_ind in range(1, stop = size(Data)[1])
                #thresholdd=threshold_Data[gene_ind,1]
                print(string(gene_ind ))

                #gene_ind =193

                cellsused=findall(.!ismissing.(Data[gene_ind,:]))
            
                #new_param_space=extract_pointInterval_v3(;MME_ori=MME_out,interval_1= interval_1,interval_2= interval_2,interval_3= interval_3,M1_boot=M1_boot, M2_boot=M2_boot,M3_boot=M3_boot,gene_ind=gene_ind,particles=particles,mode_allele=mode_allele)
                new_param_space=extract_pointInterval_fano(;Data=Data,MME_ori=MME_out,interval_1= interval_1,interval_2= interval_2,interval_3= interval_3,M1_boot=M1_boot, M2_boot=M2_boot,M3_boot=M3_boot,gene_ind=gene_ind,particles=particles,mode_allele=mode_allele)


                if size(new_param_space)[1]==0
                    preselect=[]
                else
                    preselect=findall( (new_param_space[:,1] .>0) .& (new_param_space[:,2] .>0) .& (new_param_space[:,3] .>0))
                end

                
                if length(preselect)>1000
                    ww=fill(1/length(preselect), length(preselect))
                    ranind=wsample(range(1,stop=length(preselect)),ww,1000,replace=false)
                    preselect=preselect[ranind]
                end

                if length(preselect)==0
                    store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
                else
                    new_param_space=new_param_space[preselect,:]
                    temp_dat=sim_fun(Param_Space=new_param_space,BETA_vec=BETA_vec[cellsused])
            
                    ff=ABC_test(temp_dat, Data[gene_ind,:])

                    
                   # CSV.write(string("D:/RNAseqProject/Julia_allele/ABC_debug/ff_emd.csv"),DataFrame(reshape(ff,size(ff)[1],1))) 


                    thresholdd=quantile(ff,0.01)
                    particles_selected=findall( (ff .<= thresholdd) )


                    population_dat_3=new_param_space[particles_selected,:]


                    accepted_particles=size(population_dat_3)[1]
                    if accepted_particles==0
                    final_estimates=[NaN NaN NaN]
                    else
                    final_estimates=extract_ml_particle_wt(transpose(population_dat_3))
                    end

                    store_dat_all_median_3[gene_ind,:] = [final_estimates[1] final_estimates[2] final_estimates[3] thresholdd accepted_particles particles] 

                end
                #put!(channel, true)
            end # end of each gene (end of distributed)
    end# end of verbose ifelse


    outputdata=rename!(DataFrame(store_dat_all_median_3),[:kon,:koff,:ksyn,:threshold,:accepted_particles,:particles])

    return (outputdata)
end








function ABC_BP_base_interval(; Data,BETA_vec, mode_allele=true, particles=10000,threshold_on_distance=0.95,verbose=false)
    n_genes=size(Data)[1]
    n_cells=size(Data)[2]

    #calculate threshold
    #threshold_Data=ABC_threshold(Data=Data, particles=1000,threshold_on_distance=threshold_on_distance,verbose=verbose)
    threshold_Data, interval_1, interval_2, interval_3= ABC_initialization(Data=Data,BETA_vec=BETA_vec,num_boot=1000,threshold_on_distance=threshold_on_distance,mode_allele=mode_allele,verbose=verbose)


    
    if mode_allele==true
        print("Running ABC based on allele counts.")
        #MME_out=MomentInference_data(Data=Data,BETA_vec=BETA_vec)
        MME_out=MomentInference_allele(Data=Data,BETA_vec=BETA_vec)
        sim_fun=Data_sim
    else
        print("Running ABC based on UMI counts.")
        #MME_out=MomentInference_sum2alleles_data(Data=Data,BETA_vec=BETA_vec)
        MME_out=MomentInference_UMI(Data=Data,BETA_vec=BETA_vec)
        sim_fun=Data_sim_sum2
    end
    

    #simulate parameter space using the same way as that in simulation of input data
    #Random.seed!(random_seed)
    tempspace=ABC_sim(num_samples=particles)
    @everywhere new_param_space=$tempspace

    Moments_inferred=kinetic_to_moments(kon=new_param_space[:,1],koff=new_param_space[:,2],ksyn=new_param_space[:,3],mode_allele=mode_allele)
    @everywhere Moments_inferred=$Moments_inferred


    store_dat_all_median_3 =SharedArray{Float64,2}(n_genes, 15)
    # int_1 =SharedArray{Float64,2}(n_genes, 2)
    # int_2 =SharedArray{Float64,2}(n_genes, 2)
    # int_3 =SharedArray{Float64,2}(n_genes, 2)

    if verbose==true
        progress = Progress(n_genes, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
        channel = RemoteChannel(()->Channel{Bool}(n_genes), 1)


        @sync begin
            # this task prints the progress bar
            @async while take!(channel)
                next!(progress)
            end
        
            # this task does the computation
            @async begin
                @distributed vcat for gene_ind in range(1, stop = n_genes)
                    #@distributed vcat for gene_ind in range(1, stop = size(Data)[1])

                        #Moments_inferred=kinetic_to_moments(new_param_space[:,1],new_param_space[:,2],new_param_space[:,3])
                        preselect=findall( (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) )
         
                         if length(preselect)>10000
                             ww=fill(1/length(preselect), length(preselect))
                             ranind=wsample(range(1,stop=length(preselect)),ww,10000,replace=false)
                             preselect=preselect[ranind]
                         end
         
         
                         if length(preselect)==0
                            #store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
                             store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]
                         else
                            new_param_space_used=new_param_space[preselect,:]
                            temp_dat=sim_fun(Param_Space=new_param_space_used,BETA_vec=BETA_vec)
                            ff=ABC_test(temp_dat, Data[gene_ind,:])
                            thresholdd=quantile(ff,0.05)
                            population_dat_3=new_param_space_used[findall(ff .<= thresholdd),:]
                            accepted_particles=size(population_dat_3)[1]
                            final_estimates=extract_ml_particle_wt(transpose(population_dat_3))
                            rr=population_dat_3[:,3] ./population_dat_3[:,2]



                            store_dat_all_median_3[gene_ind,:] = [final_estimates[1] final_estimates[2] final_estimates[3] median(rr) thresholdd accepted_particles particles quantile_fun(population_dat_3[:,1], 0.025) quantile_fun(population_dat_3[:,1], 0.975) quantile_fun(population_dat_3[:,2], 0.025) quantile_fun(population_dat_3[:,2], 0.975) quantile_fun(population_dat_3[:,3], 0.025) quantile_fun(population_dat_3[:,3], 0.975) quantile_fun(rr, 0.025) quantile_fun(rr, 0.975)] 
                            # int_1[gene_ind,:]=[quantile_fun(population_dat_3[:,1], 0.025) quantile_fun(population_dat_3[:,1], 0.975)]
                            # int_2[gene_ind,:]=[quantile_fun(population_dat_3[:,2], 0.025) quantile_fun(population_dat_3[:,2], 0.975)]
                            # int_3[gene_ind,:]=[quantile_fun(population_dat_3[:,3], 0.025) quantile_fun(population_dat_3[:,3], 0.975)]
                         end
       
                        put!(channel, true)
                    end # end of each gene (end of distributed)
                put!(channel, false) # this tells the printing task to finish
            end
        end #end of begin for progress bar

    else
        @distributed vcat for gene_ind in range(1, stop = n_genes)
            #@distributed vcat for gene_ind in range(1, stop = size(Data)[1])
                print(string(gene_ind , "\r"))

                preselect=findall( (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) )
 
                 if length(preselect)>10000
                     ww=fill(1/length(preselect), length(preselect))
                     ranind=wsample(range(1,stop=length(preselect)),ww,10000,replace=false)
                     preselect=preselect[ranind]
                 end
 
 
                 if length(preselect)==0
                     #store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
                     store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]
                 else
                    new_param_space_used=new_param_space[preselect,:]
                    temp_dat=sim_fun(Param_Space=new_param_space_used,BETA_vec=BETA_vec)
                    ff=ABC_test(temp_dat, Data[gene_ind,:])
                    thresholdd=quantile(ff,0.05)
                    population_dat_3=new_param_space_used[findall(ff .<= thresholdd),:]
                    accepted_particles=size(population_dat_3)[1]
                    final_estimates=extract_ml_particle_wt(transpose(population_dat_3))
                    rr=population_dat_3[:,3] ./population_dat_3[:,2]



                    store_dat_all_median_3[gene_ind,:] = [final_estimates[1] final_estimates[2] final_estimates[3] median(rr) thresholdd accepted_particles particles quantile_fun(population_dat_3[:,1], 0.025) quantile_fun(population_dat_3[:,1], 0.975) quantile_fun(population_dat_3[:,2], 0.025) quantile_fun(population_dat_3[:,2], 0.975) quantile_fun(population_dat_3[:,3], 0.025) quantile_fun(population_dat_3[:,3], 0.975) quantile_fun(rr, 0.025) quantile_fun(rr, 0.975)] 
                    # int_1[gene_ind,:]=[quantile_fun(population_dat_3[:,1], 0.025) quantile_fun(population_dat_3[:,1], 0.975)]
                    # int_2[gene_ind,:]=[quantile_fun(population_dat_3[:,2], 0.025) quantile_fun(population_dat_3[:,2], 0.975)]
                    # int_3[gene_ind,:]=[quantile_fun(population_dat_3[:,3], 0.025) quantile_fun(population_dat_3[:,3], 0.975)]

                 end
            end # end of each gene (end of distributed)

    end # end of verbose ifelse

    
    #outputdata=rename!(DataFrame(store_dat_all_median_3),[:kon,:koff,:ksyn,:threshold,:accepted_particles,:particles])
    outputdata=rename!(DataFrame(store_dat_all_median_3),[:kon,:koff,:ksyn,:bs,:threshold,:accepted_particles,:particles,:lower_kon,:upper_kon,:lower_koff,:upper_koff,:lower_ksyn,:upper_ksyn ,:lower_bs,:upper_bs])

    #return outputdata,int_1,int_2,int_3
    return outputdata
end