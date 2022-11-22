function ABC_initialization(; Data,BETA_vec,num_boot=100,threshold_on_distance=0.95,mode_allele=true,verbose=false)
    scaledData= transpose(transpose(Data) ./ BETA_vec)
    n_genes=size(Data)[1]
    n_cells=size(Data)[2]
    ww=fill(1/size(Data)[2], size(Data)[2])

    if mode_allele==true
        m1_ori = mean(scaledData,dims=2);
        m2_ori = sum(scaledData .* (scaledData .- 1.),dims=2) ./n_cells;
        m3_ori = sum(scaledData .* (scaledData .- 1.).* (scaledData .- 2.),dims=2) ./n_cells;
    else
        M1 = mean(scaledData,dims=2);
        M2 = sum(scaledData .* (scaledData .- 1.),dims=2) ./n_cells;
        M3 = sum(scaledData .* (scaledData .- 1.) .* (scaledData .- 2.),dims=2) ./n_cells;
        m1_ori = M1 ./2;
        m2_ori = M2 ./2 .- m1_ori .*m1_ori;
        m3_ori = M3 ./2 .- 3 .* m1_ori .* m2_ori;
    
    end

    #num_boot=10
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
                    #if mode_allele==true
                    m1 = mean(scaledbootdat,dims=2);
                    m2 = sum(scaledbootdat .* (scaledbootdat .- 1.),dims=2) ./n_cells;
                    m3 = sum(scaledbootdat .* (scaledbootdat .- 1.).* (scaledbootdat .- 2.),dims=2) ./n_cells;
                    # else
                    #     M1 = mean(scaledbootdat,dims=2);
                    #     M2 = sum(scaledbootdat .* (scaledbootdat .- 1.),dims=2) ./n_cells;
                    #     M3 = sum(scaledbootdat .* (scaledbootdat .- 1.) .* (scaledbootdat .- 2.),dims=2) ./n_cells;
                    #     m1 = M1 ./2;
                    #     m2 = M2 ./2 .- m1 .*m1;
                    #     m3 = M3 ./2 .- 3 .* m1 .* m2;
                    # end
            
                    M1_boot[:,boot_ind]=m1
                    M2_boot[:,boot_ind]=m2
                    M3_boot[:,boot_ind]=m3
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
            #if mode_allele==true
                m1 = mean(scaledbootdat,dims=2);
                m2 = sum(scaledbootdat .* (scaledbootdat .- 1.),dims=2) ./n_cells;
                m3 = sum(scaledbootdat .* (scaledbootdat .- 1.).* (scaledbootdat .- 2.),dims=2) ./n_cells;
            # else
            #     M1 = mean(scaledbootdat,dims=2);
            #     M2 = sum(scaledbootdat .* (scaledbootdat .- 1.),dims=2) ./n_cells;
            #     M3 = sum(scaledbootdat .* (scaledbootdat .- 1.) .* (scaledbootdat .- 2.),dims=2) ./n_cells;
            #     m1 = M1 ./2;
            #     m2 = M2 ./2 .- m1 .*m1;
            #     m3 = M3 ./2 .- 3 .* m1 .* m2;
            # end
    
            M1_boot[:,boot_ind]=m1
            M2_boot[:,boot_ind]=m2
            M3_boot[:,boot_ind]=m3
            Dist_boot[:,boot_ind]=dist_genes
            #put!(channel, true)
        end # end of each gene (end of distributed)

    end

    # interval_1=[(m1_ori .-2 .*std(M1_boot,dims=2)) (m1_ori .+ 2 .*std(M1_boot,dims=2))]
    # interval_2=[(m2_ori .-2 .*std(M2_boot,dims=2)) (m2_ori .+ 2 .*std(M2_boot,dims=2))]
    # interval_3=[(m3_ori .-2 .*std(M3_boot,dims=2)) (m3_ori .+ 2 .*std(M3_boot,dims=2))]

    interval_1=[mapslices(X -> quantile(X, 0.05), M1_boot, dims=(2)) mapslices(X -> quantile(X, 0.95), M1_boot, dims=(2))]
    interval_2=[mapslices(X -> quantile(X, 0.05), M2_boot, dims=(2)) mapslices(X -> quantile(X, 0.95), M2_boot, dims=(2))]
    interval_3=[mapslices(X -> quantile(X, 0.05), M3_boot, dims=(2)) mapslices(X -> quantile(X, 0.95), M3_boot, dims=(2))]
    
    #std(M3_boot,dims=2)[gene_ind]

    threshold_data=mapslices(X -> quantile(X, threshold_on_distance), Dist_boot, dims=(2))

    return threshold_data, interval_1, interval_2, interval_3
end

#threshold_data, interval_1, interval_2, interval_3= ABC_initialization(Data=Data,BETA_vec=BETA_vec,num_boot=100,threshold_on_distance=0.95,mode_allele=true,verbose=false)



function ABC_MME(; Data,BETA_vec,num_boot=100,mode_allele=true,verbose=false)
    scaledData= transpose(transpose(Data) ./ BETA_vec)
    n_genes=size(Data)[1]
    n_cells=size(Data)[2]
    ww=fill(1/size(Data)[2], size(Data)[2])

    if mode_allele==true
        m1_ori = mean(scaledData,dims=2);
        m2_ori = sum(scaledData .* (scaledData .- 1.),dims=2) ./n_cells;
        m3_ori = sum(scaledData .* (scaledData .- 1.).* (scaledData .- 2.),dims=2) ./n_cells;
    else
        M1 = mean(scaledData,dims=2);
        M2 = sum(scaledData .* (scaledData .- 1.),dims=2) ./n_cells;
        M3 = sum(scaledData .* (scaledData .- 1.) .* (scaledData .- 2.),dims=2) ./n_cells;
        m1_ori = M1 ./2;
        m2_ori = M2 ./2 .- m1_ori .*m1_ori;
        m3_ori = M3 ./2 .- 3 .* m1_ori .* m2_ori;
    
    end

    #num_boot=10
    kon_boot=SharedArray{Float64,2}(n_genes, num_boot)
    koff_boot=SharedArray{Float64,2}(n_genes, num_boot)
    ksyn_boot=SharedArray{Float64,2}(n_genes, num_boot)
    #Dist_boot=SharedArray{Float64,2}(n_genes, num_boot)

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
                    
                    #dist_genes=hellinger2_data(Data, bootdat)
                    if mode_allele==true
                        m1 = mean(scaledbootdat,dims=2);
                        m2 = sum(scaledbootdat .* (scaledbootdat .- 1.),dims=2) ./n_cells;
                        m3 = sum(scaledbootdat .* (scaledbootdat .- 1.).* (scaledbootdat .- 2.),dims=2) ./n_cells;
                    else
                        M1 = mean(scaledbootdat,dims=2);
                        M2 = sum(scaledbootdat .* (scaledbootdat .- 1.),dims=2) ./n_cells;
                        M3 = sum(scaledbootdat .* (scaledbootdat .- 1.) .* (scaledbootdat .- 2.),dims=2) ./n_cells;
                        m1 = M1 ./2;
                        m2 = M2 ./2 .- m1 .*m1;
                        m3 = M3 ./2 .- 3 .* m1 .* m2;
                    end

                    r1=m1;
                    r2=m2 ./m1;
                    r3=m3 ./m2;
                    lambda_est = (2 .*r1 .*(r3 .-r2)) ./(r1 .*r2-2 .*r1 .*r3 + r2 .*r3);
                    mu_est = (2 .*(r3 .-r2) .*(r1 .-r3) .*(r2 .-r1)) ./((r1 .*r2  .- 2 .*r1 .*r3 .+ r2 .*r3) .*(r1 .-2 .*r2 .+r3));
                    v_est = (2 .*r1 .*r3 - r1 .*r2 .- r2 .*r3) ./(r1 .- 2 .*r2 .+ r3);


                    dropp=findall((lambda_est .<0) .| (mu_est .<0) .| (v_est .<0))
                    lambda_est[dropp] .=NaN
                    mu_est[dropp] .=NaN
                    v_est[dropp] .=NaN

            
                    kon_boot[:,boot_ind]=lambda_est
                    koff_boot[:,boot_ind]=mu_est
                    ksyn_boot[:,boot_ind]=v_est
                    #Dist_boot[:,boot_ind]=dist_genes
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
            
            #dist_genes=hellinger2_data(Data, bootdat)
            if mode_allele==true
                m1 = mean(scaledbootdat,dims=2);
                m2 = sum(scaledbootdat .* (scaledbootdat .- 1.),dims=2) ./n_cells;
                m3 = sum(scaledbootdat .* (scaledbootdat .- 1.).* (scaledbootdat .- 2.),dims=2) ./n_cells;
            else
                M1 = mean(scaledbootdat,dims=2);
                M2 = sum(scaledbootdat .* (scaledbootdat .- 1.),dims=2) ./n_cells;
                M3 = sum(scaledbootdat .* (scaledbootdat .- 1.) .* (scaledbootdat .- 2.),dims=2) ./n_cells;
                m1 = M1 ./2;
                m2 = M2 ./2 .- m1 .*m1;
                m3 = M3 ./2 .- 3 .* m1 .* m2;
            end

            r1=m1;
            r2=m2 ./m1;
            r3=m3 ./m2;
            lambda_est = (2 .*r1 .*(r3 .-r2)) ./(r1 .*r2-2 .*r1 .*r3 + r2 .*r3);
            mu_est = (2 .*(r3 .-r2) .*(r1 .-r3) .*(r2 .-r1)) ./((r1 .*r2  .- 2 .*r1 .*r3 .+ r2 .*r3) .*(r1 .-2 .*r2 .+r3));
            v_est = (2 .*r1 .*r3 - r1 .*r2 .- r2 .*r3) ./(r1 .- 2 .*r2 .+ r3);


            dropp=findall((lambda_est .<0) .| (mu_est .<0) .| (v_est .<0))
            lambda_est[dropp] .=NaN
            mu_est[dropp] .=NaN
            v_est[dropp] .=NaN

    
            kon_boot[:,boot_ind]=lambda_est
            koff_boot[:,boot_ind]=mu_est
            ksyn_boot[:,boot_ind]=v_est
            #Dist_boot[:,boot_ind]=dist_genes
            #put!(channel, true)
        end # end of each gene (end of distributed)

    end


    # interval_1=[mapslices(X -> quantile(X, 0.05), M1_boot, dims=(2)) mapslices(X -> quantile(X, 0.95), M1_boot, dims=(2))]
    # interval_2=[mapslices(X -> quantile(X, 0.05), M2_boot, dims=(2)) mapslices(X -> quantile(X, 0.95), M2_boot, dims=(2))]
    # interval_3=[mapslices(X -> quantile(X, 0.05), M3_boot, dims=(2)) mapslices(X -> quantile(X, 0.95), M3_boot, dims=(2))]

    point_kon=mapslices(X -> median_fun(X), kon_boot, dims=(2))
    point_koff=mapslices(X -> median_fun(X), koff_boot, dims=(2))
    point_ksyn=mapslices(X -> median_fun(X), ksyn_boot, dims=(2))

    #threshold_data=mapslices(X -> quantile(X, threshold_on_distance), Dist_boot, dims=(2))
    store_dat_all_median_3=[point_kon point_koff point_ksyn]
    outputdata=rename!(DataFrame(store_dat_all_median_3),[:kon,:koff,:ksyn])
    return  outputdata
end





function ABC_threshold(; Data, particles=1000,threshold_on_distance=0.95,verbose=false)

    n_genes=size(Data)[1]
    n_cells=size(Data)[2]

    threshold_data =SharedArray{Float64,2}(n_genes, 1)


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
                    p=Data[gene_ind,:]
                    ww=fill(1/size(p)[1], size(p)[1])
                    aa=Array{Float64,2}(undef,particles,1)
                    for i in 1:size(aa)[1]
                        q=wsample(p,ww,size(p)[1],replace=true)
                        tt=ABC_OneParam_test(p, q)
                        aa[i,1]=tt        
                    end
                    thresholdd=quantile(aa[:,1],threshold_on_distance)


                    threshold_data[gene_ind,1]=thresholdd

                    put!(channel, true)


                    end # end of each gene (end of distributed)
                put!(channel, false) # this tells the printing task to finish
            end
        end #end of begin for progress bar


    else #else not show progress bar
        @distributed vcat for gene_ind in range(1, stop = n_genes)
            p=Data[gene_ind,:]
            ww=fill(1/size(p)[1], size(p)[1])
            aa=Array{Float64,2}(undef,particles,1)
            for i in 1:size(aa)[1]
                q=wsample(p,ww,size(p)[1],replace=true)
                tt=ABC_OneParam_test(p, q)
                aa[i,1]=tt        
            end
            thresholdd=quantile(aa[:,1],threshold_on_distance)
            threshold_data[gene_ind,1]=thresholdd
            end # end of each gene (end of distributed)

    end #end verbose ifelse

    
    return(threshold_data)
end


#threshold_Data=ABC_threshold(Data=Data[1:10,:], particles=10000,verbose=false)


function ABC_BP(; Data,BETA_vec, mode_allele=true, particles=10000,threshold_on_distance=0.95,verbose=false)
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
                    
                        #new_param_space=10 .^ [qq1[gene_ind] qq2[gene_ind] qq3[gene_ind]]
                        new_param_space=extract_pointMME(; MME_out= MME_out,gene_ind=gene_ind,particles=particles)
                        Moments_inferred=kinetic_to_moments(kon=new_param_space[:,1],koff=new_param_space[:,2],ksyn=new_param_space[:,3],mode_allele=mode_allele)

                        preselect=findall( (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) )
                        #preselect=findall( (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) .& (interval_2[gene_ind,1] .<Moments_inferred[:,2].<  interval_2[gene_ind,2]) .& (interval_3[gene_ind,1] .<Moments_inferred[:,3].<  interval_3[gene_ind,2]))
                        
                        if length(preselect)>10000
                            ww=fill(1/length(preselect), length(preselect))
                            ranind=wsample(range(1,stop=length(preselect)),ww,10000,replace=false)
                            preselect=preselect[ranind]
                        end

                        if length(preselect)==0
                            store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
                        else
                            new_param_space=new_param_space[preselect,:]
                            temp_dat=sim_fun(Param_Space=new_param_space,BETA_vec=BETA_vec[cellsused])
                    
                            ff=ABC_test(temp_dat, Data[gene_ind,:])
                            thresholdd=quantile(ff,0.05)
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
                #print(string(gene_ind , "\r"))
                cellsused=findall(.!ismissing.(Data[gene_ind,:]))
            
                #new_param_space=10 .^ [qq1[gene_ind] qq2[gene_ind] qq3[gene_ind]]
                new_param_space=extract_pointMME(; MME_out= MME_out,gene_ind=gene_ind,particles=particles)
                Moments_inferred=kinetic_to_moments(kon=new_param_space[:,1],koff=new_param_space[:,2],ksyn=new_param_space[:,3],mode_allele=mode_allele)


               preselect=findall( (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) )
               #preselect=findall( (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) .& (interval_2[gene_ind,1] .<Moments_inferred[:,2].<  interval_2[gene_ind,2]) .& (interval_3[gene_ind,1] .<Moments_inferred[:,3].<  interval_3[gene_ind,2]))

                if length(preselect)>10000
                    ww=fill(1/length(preselect), length(preselect))
                    ranind=wsample(range(1,stop=length(preselect)),ww,10000,replace=false)
                    preselect=preselect[ranind]
                end


                if length(preselect)==0
                    store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
                else
                    new_param_space=new_param_space[preselect,:]
                    temp_dat=sim_fun(Param_Space=new_param_space,BETA_vec=BETA_vec[cellsused])
            
                    ff=ABC_test(temp_dat, Data[gene_ind,:])
                    thresholdd=quantile(ff,0.05)
		            particles_selected=findall( (ff .<= thresholdd))


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


function ABC_BP_base(; Data,BETA_vec, mode_allele=true, particles=10000,threshold_on_distance=0.95,verbose=false)
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

                        #Moments_inferred=kinetic_to_moments(new_param_space[:,1],new_param_space[:,2],new_param_space[:,3])
                        preselect=findall( (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) )
         
                         if length(preselect)>10000
                             ww=fill(1/length(preselect), length(preselect))
                             ranind=wsample(range(1,stop=length(preselect)),ww,10000,replace=false)
                             preselect=preselect[ranind]
                         end
         
         
                         if length(preselect)==0
                             store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
                         else
                            new_param_space_used=new_param_space[preselect,:]
                            temp_dat=sim_fun(Param_Space=new_param_space_used,BETA_vec=BETA_vec)
                            ff=ABC_test(temp_dat, Data[gene_ind,:])
                            thresholdd=quantile(ff,0.05)
                            population_dat_3=new_param_space_used[findall(ff .<= thresholdd),:]
                            accepted_particles=size(population_dat_3)[1]
                            final_estimates=extract_ml_particle_wt(transpose(population_dat_3))
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
                print(string(gene_ind , "\r"))

                preselect=findall( (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) )
 
                 if length(preselect)>10000
                     ww=fill(1/length(preselect), length(preselect))
                     ranind=wsample(range(1,stop=length(preselect)),ww,10000,replace=false)
                     preselect=preselect[ranind]
                 end
 
 
                 if length(preselect)==0
                     store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
                 else
                    new_param_space_used=new_param_space[preselect,:]
                    temp_dat=sim_fun(Param_Space=new_param_space_used,BETA_vec=BETA_vec)
                    ff=ABC_test(temp_dat, Data[gene_ind,:])
                    thresholdd=quantile(ff,0.05)
                    population_dat_3=new_param_space_used[findall(ff .<= thresholdd),:]
                    accepted_particles=size(population_dat_3)[1]
                    final_estimates=extract_ml_particle_wt(transpose(population_dat_3))
                    store_dat_all_median_3[gene_ind,:] = [final_estimates[1] final_estimates[2] final_estimates[3] thresholdd accepted_particles particles] 
                 end
            end # end of each gene (end of distributed)

    end # end of verbose ifelse

    

    outputdata=rename!(DataFrame(store_dat_all_median_3),[:kon,:koff,:ksyn,:threshold,:accepted_particles,:particles])

    return (outputdata)
end


#ouuput=ABC_BP(Data=raw_dat[1:10,:],BETA_vec=raw_beta, mode_allele=true, particles=10000,threshold_on_distance=0.95,verbose=false)
#ouuput=ABC_BP_base(Data=raw_dat[1:10,:],BETA_vec=raw_beta, mode_allele=true, particles=10000,threshold_on_distance=0.95,verbose=false)

#ouuput=ABC_BP(Data=raw_dat[1:10,:],BETA_vec=raw_beta, mode_allele=false, particles=10000,threshold_on_distance=0.95,verbose=false)
#ouuput=ABC_BP_base(Data=raw_dat[1:10,:],BETA_vec=raw_beta, mode_allele=false, particles=10000,threshold_on_distance=0.95,verbose=false)


# rename!(DataFrame(store_dat_all_median_3),[:kon,:koff,:ksyn,:threshold,:accepted_particles,:particles])
# sdf=["sdf1", "sdf2", "sdf3", "sdf4", "sdf5", "sdf6"]
# names!(transpose(ouuput[1:6,:]), Symbol.(sdf))



function flat_priors(;MME_out,particles=1000)
    
    x1=MME_out[:,1]
    x2=MME_out[:,2]
    x3=MME_out[:,3]


    x1[findall(x1 .<=0)] .=NaN
    x2[findall(x2 .<=0)] .=NaN
    x3[findall(x3 .<=0)] .=NaN

    lower_x1=log10.(x1./10)
    upper_x1=log10.(x1.*10)

    lower_x1[isnan.(lower_x1)].=log10(1e-2) 
    upper_x1[isnan.(upper_x1)].=log10(1e2) 

    lower_x2=log10.(x2./10)
    upper_x2=log10.(x2.*10)

    lower_x2[isnan.(lower_x2)].=log10(1e-3) 
    upper_x2[isnan.(upper_x2)].=log10(1e6) 

    lower_x3=log10.(x3./10)
    upper_x3=log10.(x3.*10)

    lower_x3[isnan.(lower_x3)].=log10(1e-3) 
    upper_x3[isnan.(upper_x3)].=log10(1e6) 


    qq1=rand.(Uniform.(lower_x1, upper_x1),particles)
    qq2=rand.(Uniform.(lower_x2, upper_x2),particles)
    qq3=rand.(Uniform.(lower_x3, upper_x3),particles)
    boundary_1=10 .^[lower_x1 upper_x1]
    boundary_2=10 .^[lower_x2 upper_x2]
    boundary_3=10 .^[lower_x3 upper_x3]

    return qq1, qq2, qq3 , boundary_1, boundary_2, boundary_3

end





# for gene_ind in range(1,stop=7000)
#     print(gene_ind)
#     extract_pointMME(; MME_out,gene_ind,particles=1000)
# end


function extract_pointMME(; MME_out,gene_ind,particles=1000)
    x1=MME_out[gene_ind,1]
    x2=MME_out[gene_ind,2]
    x3=MME_out[gene_ind,3]
    
    if (x1<=0 || isnan(x1) || isinf(x1))
        x1=NaN
        lower_x1=log10(1e-2) 
        upper_x1=log10(1e2) 

    else
        lower_x1=log10.(x1./10)
        upper_x1=log10.(x1.*10)

    end

    if (x2<=0  || isnan(x2) || isinf(x2))
        x2=NaN
        lower_x2=log10(1e-3) 
        upper_x2=log10(1e6) 
    else
        lower_x2=log10.(x2./10)
        upper_x2=log10.(x2.*10)
    end
    if (x3<=0  || isnan(x3) || isinf(x3))
        x3=NaN
        lower_x3=log10(1e-3) 
        upper_x3=log10(1e6) 
    else
        lower_x3=log10.(x3./10)
        upper_x3=log10.(x3.*10)
    end

    qq1=rand.(Uniform.(lower_x1, upper_x1),particles)
    qq2=rand.(Uniform.(lower_x2, upper_x2),particles)
    qq3=rand.(Uniform.(lower_x3, upper_x3),particles)
    new_param_space=10 .^[qq1 qq2 qq3]
    return new_param_space
end


# function extract_pointMME(; MME_out,gene_ind,particles=1000)
#     x1=MME_out[gene_ind,1]
#     x2=MME_out[gene_ind,2]
#     x3=MME_out[gene_ind,3]
    
#     # if (x1<=0 || isnan(x1) || isinf(x1))
#     #     x1=NaN
#         lower_x1=log10(1e-2) 
#         upper_x1=log10(1e2) 

#     # else
#     #     lower_x1=log10.(x1./10)
#     #     upper_x1=log10.(x1.*10)

#     # end

#     # if (x2<=0  || isnan(x2) || isinf(x2))
#     #     x2=NaN
#         lower_x2=log10(1e-3) 
#         upper_x2=log10(1e6) 
#     # else
#     #     lower_x2=log10.(x2./10)
#     #     upper_x2=log10.(x2.*10)
#     # end
#     # if (x3<=0  || isnan(x3) || isinf(x3))
#     #     x3=NaN
#         lower_x3=log10(1e-3) 
#         upper_x3=log10(1e6) 
#     # else
#     #     lower_x3=log10.(x3./10)
#     #     upper_x3=log10.(x3.*10)
#     # end

#     qq1=rand.(Uniform.(lower_x1, upper_x1),particles)
#     qq2=rand.(Uniform.(lower_x2, upper_x2),particles)
#     qq3=rand.(Uniform.(lower_x3, upper_x3),particles)
#     new_param_space=10 .^[qq1 qq2 qq3]
#     return new_param_space
# end





function ABC_BP_moments(; Data,BETA_vec, mode_allele=true, particles=10000,threshold_on_distance=0.95,verbose=false)
    n_genes=size(Data)[1]
    n_cells=size(Data)[2]

    #calculate threshold
    #threshold_Data=ABC_threshold(Data=Data, particles=1000,threshold_on_distance=threshold_on_distance,verbose=verbose)
    threshold_Data, interval_1, interval_2, interval_3= ABC_initialization(Data=Data,BETA_vec=BETA_vec,num_boot=100,threshold_on_distance=threshold_on_distance,mode_allele=mode_allele,verbose=verbose)


    if mode_allele==true
        print("Running ABC based on allele counts.")
        #MME_out=MomentInference_data(Data=Data,BETA_vec=BETA_vec)
        MME_out=MomentInference_allele(Data=Data,BETA_vec=BETA_vec)
    else
        print("Running ABC based on UMI counts.")
        #MME_out=MomentInference_sum2alleles_data(Data=Data,BETA_vec=BETA_vec)
        MME_out=MomentInference_UMI(Data=Data,BETA_vec=BETA_vec)
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
                        thresholdd=threshold_Data[gene_ind,1]
                    
                        #new_param_space=10 .^ [qq1[gene_ind] qq2[gene_ind] qq3[gene_ind]]
                        new_param_space=extract_pointMME(; MME_out= MME_out,gene_ind=gene_ind,particles=particles)
                
                        # if mode_allele==true
                        #     temp_dat=Data_sim(Param_Space=new_param_space,BETA_vec=BETA_vec)
                        # else
                        #     temp_dat=Data_sim_sum2(Param_Space=new_param_space,BETA_vec=BETA_vec)
                        # end
                
                        # ff=ABC_test(temp_dat, Data[gene_ind,:])
                        #thresholdd=quantile(ff,0.05)
                        Moments_inferred=kinetic_to_moments(kon=new_param_space[:,1],koff=new_param_space[:,2],ksyn=new_param_space[:,3],mode_allele=mode_allele)
                        #ff=ABC_test(temp_dat, Data[gene_ind,:])
                        #particles_selected=findall( (ff .<= thresholdd)   .& (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) .& (interval_2[gene_ind,1] .<Moments_inferred[:,2].<  interval_2[gene_ind,2]) .& (interval_3[gene_ind,1] .<Moments_inferred[:,3].<  interval_3[gene_ind,2]))
                        particles_selected=findall( (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) .& (interval_2[gene_ind,1] .<Moments_inferred[:,2].<  interval_2[gene_ind,2]) .& (interval_3[gene_ind,1] .<Moments_inferred[:,3].<  interval_3[gene_ind,2]))
                    
                        #if sum(prop_vec)==0
                        if size( particles_selected)[1]==0
                    
                            store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
                    
                        else
                
                            population_dat_3=new_param_space[particles_selected,:]
                    
                    
                            if size(population_dat_3)[1]==0
                                store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
                    
                            else
                                accepted_particles=size(population_dat_3)[1]
                                extract_anthony=extract_ml_particle_wt(transpose(population_dat_3))
                                store_dat_all_median_3[gene_ind,:] = [extract_anthony[1] extract_anthony[2] extract_anthony[3] thresholdd accepted_particles particles] 
                
                            end
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
                thresholdd=threshold_Data[gene_ind,1]
                print(gene_ind)
            
                #new_param_space=10 .^ [qq1[gene_ind] qq2[gene_ind] qq3[gene_ind]]
                new_param_space=extract_pointMME(; MME_out= MME_out,gene_ind=gene_ind,particles=particles)
        
                # if mode_allele==true
                #     temp_dat=Data_sim(Param_Space=new_param_space,BETA_vec=BETA_vec)
                # else
                #     temp_dat=Data_sim_sum2(Param_Space=new_param_space,BETA_vec=BETA_vec)
                # end
        
                # ff=ABC_test(temp_dat, Data[gene_ind,:])
                #thresholdd=quantile(ff,0.05)
            
                #if sum(prop_vec)==0
                Moments_inferred=kinetic_to_moments(new_param_space[:,1],new_param_space[:,2],new_param_space[:,3])
                #ff=ABC_test(temp_dat, Data[gene_ind,:])
                #particles_selected=findall( (ff .<= thresholdd)   .& (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) .& (interval_2[gene_ind,1] .<Moments_inferred[:,2].<  interval_2[gene_ind,2]) .& (interval_3[gene_ind,1] .<Moments_inferred[:,3].<  interval_3[gene_ind,2]))
                particles_selected=findall( (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) .& (interval_2[gene_ind,1] .<Moments_inferred[:,2].<  interval_2[gene_ind,2]) .& (interval_3[gene_ind,1] .<Moments_inferred[:,3].<  interval_3[gene_ind,2]))


                if size( particles_selected)[1]==0
            
                    store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
            
                else
        
                    population_dat_3=new_param_space[particles_selected,:]
            
            
                    if size(population_dat_3)[1]==0
                        store_dat_all_median_3[gene_ind,:]=[NaN NaN NaN NaN NaN NaN]
            
                    else
                        accepted_particles=size(population_dat_3)[1]
                        extract_anthony=extract_ml_particle_wt(transpose(population_dat_3))
                        store_dat_all_median_3[gene_ind,:] = [extract_anthony[1] extract_anthony[2] extract_anthony[3] thresholdd accepted_particles particles] 
        
                    end
                end
            end # end of each gene (end of distributed)

    end# end of verbose ifelse


    outputdata=rename!(DataFrame(store_dat_all_median_3),[:kon,:koff,:ksyn,:threshold,:accepted_particles,:particles])

    return (outputdata)
end

function preselection(;Moments_inferred,interval_1,interval_2, interval_3,gene_ind)
    preselect=findall( (interval_1[gene_ind,1] .<Moments_inferred[:,1].<  interval_1[gene_ind,2]) .& (interval_2[gene_ind,1] .<Moments_inferred[:,2].<  interval_2[gene_ind,2]) .& (interval_3[gene_ind,1] .<Moments_inferred[:,3].<  interval_3[gene_ind,2]))
    return preselect
end

function while_fun(;Data,BETA_vec,gene_ind,population_data_initial,interval_1,interval_2, interval_3,threshold,particles=1000,while_thres=20, mode_allele=true,verbose=false)
    population_data=population_data_initial

    if size(population_data)[1]<=5
        final_estimates=extract_ml_particle_wt(transpose(population_data))
        population_data=population_data

    else
        mu=mean(population_data,dims=1)
        sd=std(population_data,dims=1)
        counttt=1
        while  sum(sd .< 0.5)!=3
            # if verbose
                 print(counttt)
            # end
            
            mu=mean(population_data,dims=1)
            sd=std(population_data,dims=1)
            minn=minimum(population_data,dims=1)
            maxx=maximum(population_data,dims=1)

            upper_x1=log10(mu[1]+sd[1])
            upper_x2=log10(mu[2]+sd[2])
            upper_x3=log10(mu[3]+sd[3])
            if (mu[1]-sd[1])<=0
                lower_x1=log10(mu[1]/2)
            else
                lower_x1=log10(mu[1]-sd[1])
            end
            if  (mu[2]-sd[2])<=0
                lower_x2=log10(mu[2]/2)
            else
                lower_x2=log10(mu[2]-sd[2])
            end
            if  (mu[3]-sd[3])<=0
                lower_x3=log10(mu[3]/2)
            else
                lower_x3=log10(mu[3]-sd[3])
            end
        
            qq1=rand.(Uniform.(lower_x1, upper_x1),particles)
            qq2=rand.(Uniform.(lower_x2, upper_x2),particles)
            qq3=rand.(Uniform.(lower_x3, upper_x3),particles)
            space2=10 .^ [qq1 qq2 qq3]
            Moments_temp=kinetic_to_moments( space2[:,1], space2[:,2], space2[:,3])
            preselect=preselection(Moments_inferred=Moments_temp,interval_1=interval_1,interval_2=interval_2, interval_3=interval_3,gene_ind=gene_ind)
            space2=space2[preselect,:]
            if size(space2)[1]==0
                final_estimates=[NaN NaN NaN]
                population_data=population_data
                break

            else
                if mode_allele==true
                    abc_dat=Data_sim(Param_Space=space2,BETA_vec=BETA_vec)
                else
                    abc_dat=Data_sim_sum2(Param_Space=space2,BETA_vec=BETA_vec)
                end
                abc_dist=ABC_test(abc_dat, Data[gene_ind,:])
                population_data=space2[findall(abc_dist .<= threshold),:]
 
                counttt=counttt+1
                if counttt==while_thres
                    break
                end 
            end
            sd=std(space2,dims=1)
    
           
        end# end while
        if size(population_data)[1]==0
            final_estimates=[NaN NaN NaN]
        else
            final_estimates=extract_ml_particle_wt(transpose(population_data))
        end
        
    end

    return final_estimates,population_data

   
end





function extract_pointInterval(; interval_1, interval_2, interval_3,gene_ind,particles=1000,mode_allele=true)

    if  ((interval_1[gene_ind,1] == interval_1[gene_ind,2])  ||  (interval_2[gene_ind,1] == interval_2[gene_ind,2])  ||  (interval_3[gene_ind,1] == interval_3[gene_ind,2] ))

        new_param_space=[]
    else
        lower_x1=(interval_1[gene_ind,1]) 
        upper_x1=(interval_1[gene_ind,2]) 
    
        lower_x2=interval_2[gene_ind,1]
        upper_x2=interval_2[gene_ind,2]
    
        lower_x3=interval_3[gene_ind,1] 
        upper_x3=interval_3[gene_ind,2]
    
        m1=rand.(Uniform.(lower_x1, upper_x1),particles)
        m2=rand.(Uniform.(lower_x2, upper_x2),particles)
        m3=rand.(Uniform.(lower_x3, upper_x3),particles)
        if mode_allele==true
            r1=m1;
            r2=m2 ./m1;
            r3=m3 ./m2;
        else
            M1 = m1 ./2;
            M2 = m2 ./2 .- M1 .*M1;
            M3 = m3 ./2 .- 3 .* M1 .* M2;
            r1=M1;
            r2=M2 ./M1;
            r3=M3 ./M2;
        end
        qq1 = (2 .*r1 .*(r3 .-r2)) ./(r1 .*r2-2 .*r1 .*r3 + r2 .*r3);
        qq2 = (2 .*(r3 .-r2) .*(r1 .-r3) .*(r2 .-r1)) ./((r1 .*r2  .- 2 .*r1 .*r3 .+ r2 .*r3) .*(r1 .-2 .*r2 .+r3));
        qq3 = (2 .*r1 .*r3 - r1 .*r2 .- r2 .*r3) ./(r1 .- 2 .*r2 .+ r3);
    
        new_param_space=[qq1 qq2 qq3]
    end

   
    return new_param_space
end
function extract_pointInterval_v2(; interval_1, interval_2, interval_3,gene_ind,particles=1000,mode_allele=true)

    if  ((interval_1[gene_ind,1] == interval_1[gene_ind,2])  ||  (interval_2[gene_ind,1] == interval_2[gene_ind,2])  ||  (interval_3[gene_ind,1] == interval_3[gene_ind,2] ))

        new_param_space=[]
    else


        if mode_allele==true
            r1=interval_1[gene_ind,:];
            r2=interval_2[gene_ind,:] ./interval_1[gene_ind,:];
            r3=interval_3[gene_ind,:] ./interval_2[gene_ind,:];
        else
            M1 = interval_1[gene_ind,:] ./2;
            M2 = interval_2[gene_ind,:] ./2 .- M1 .*M1;
            M3 = interval_3[gene_ind,:] ./2 .- 3 .* M1 .* M2;
            r1=M1;
            r2=M2 ./M1;
            r3=M3 ./M2;
        end
        qq1_interval = (2 .*r1 .*(r3 .-r2)) ./(r1 .*r2-2 .*r1 .*r3 + r2 .*r3);
        qq2_interval = (2 .*(r3 .-r2) .*(r1 .-r3) .*(r2 .-r1)) ./((r1 .*r2  .- 2 .*r1 .*r3 .+ r2 .*r3) .*(r1 .-2 .*r2 .+r3));
        qq3_interval = (2 .*r1 .*r3 - r1 .*r2 .- r2 .*r3) ./(r1 .- 2 .*r2 .+ r3);

        if  ((qq1_interval[1] >= qq1_interval[2])  || (qq2_interval[1] >= qq2_interval[2])  ||  (qq3_interval[1] >= qq3_interval[2]) || isnan(qq1_interval[1]) || isnan(qq1_interval[2]) || isnan(qq2_interval[1]) || isnan(qq2_interval[2]) || isnan(qq3_interval[1]) || isnan(qq3_interval[2])  || qq1_interval[1]<0 ||   qq2_interval[1]<0  || qq3_interval[1]<0   || qq1_interval[2]<0 ||   qq2_interval[2]<0  || qq3_interval[2]<0)

            new_param_space=[]

        else
            qq1=rand.(Uniform.(qq1_interval[1],qq1_interval[2]),particles)
            qq2=rand.(Uniform.(qq2_interval[1],qq2_interval[2]),particles)
            qq3=rand.(Uniform.(qq3_interval[1],qq3_interval[2]),particles)
            new_param_space=[qq1 qq2 qq3]
        end
    

    end

   
    return new_param_space
end



