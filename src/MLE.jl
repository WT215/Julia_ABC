
function Raj_PB(n,ract,rdeact,ron,rdeg=1;nemoaccur=64)

    CC = Nemo.ComplexField(nemoaccur)
    # n=4
    # ract=87.9744
    # rdeact=27.1511
    # ron=803.149
    # ract=500
    # rdeact=100
    # ron=1000
    # rdeg=1
  
      x1=CC(ract/rdeg)
      x2=CC(rdeact/rdeg)
      x3= broadcast(CC,ron ./ rdeg)
      n=broadcast(CC,n)
      #rdeg=CC(rdeg)
      part1=Nemo.gamma(x1+x2) .* x3.^n .* broadcast(Nemo.gamma,x1 .+n) ./ Nemo.gamma(x1) ./broadcast(Nemo.gamma,n.+1) ./broadcast(Nemo.gamma,x1.+x2.+n)
      part2=Nemo.hyp1f1.(x1 .+ n, x1.+x2.+n, -x3)
      outtt=part1 .* part2
  
      return convert.(Float64,real(outtt))
  end
  
  #broadcast(Nemo.gamma,CC(ract/rdeg).+CC(rdeact/rdeg).+broadcast(CC,n))
  
  #gamma(x1.+x2.+n)
  
  #Nemo.hyp1f1(CC(1.0),CC(201.0),CC(-0.01))
  #GSL.hypergeom(1.0,201.0,-0.01)
  
  
  function Raj_PB_raw(n,ract,rdeact,ron,rdeg=1)
      x1=(ract/rdeg)
      x2=(rdeact/rdeg)
      x3= ron ./ rdeg
      part1=gamma.(x1+x2) .* x3.^n .* gamma.(x1 .+n) ./ gamma.(x1) ./gamma.(n.+1) ./ gamma.(x1.+x2.+n)
      part2=hypergeom.(x1 .+ n, x1.+x2.+n, -x3)
      outtt=part1 .* part2
      return outtt
  end
  
  #= function Raj_PB_py(n,ract,rdeact,ron,rdeg=1)
  
  #mu ron
  #lambda x1: ract
  #gam x2 rdeact
    x = sympy.Symbol("x") # number
  
    lambda = sympy.Symbol("lambda")
    gam = sympy.Symbol("gamma")
  
    mu = sympy.Symbol("mu")
  
    s = sympy.Symbol("s") # size
    
    zz=sympy.gamma(lambda+x)/(sympy.gamma(x+1)*sympy.gamma(lambda+gam+x)) * sympy.gamma(lambda+gam)/sympy.gamma(lambda) * (s*mu)^x * sympy.hyper([lambda+x], [lambda+gam+x], -s*mu)
    zz=zz.subs(lambda,ract).subs(mu,ron).subs(gam,rdeact).subs(s,1)
  
  
    
    #convert(Float64,zz.subs(x,4))
  
      return convert(Float64,zz.subs(x,n))
  end =#
  
  
  #Raj_PB_py([1,5,6],4,5,[2,3,5],1)
  
function Raj_PB_py_final(n,ract,rdeact,ron,rdeg=1)

  tempp=pmap((x,y)->Raj_PB_py(x,ract,rdeact,y,1),n,ron)
  #tempp=map((x,y)->Raj_PB_py(x,ract,rdeact,y,1),n,ron)
  
      return tempp
end

# Raj_PB(1:4,500,100,4:7,1)
# Raj_PB_py_final(1:4,500,100,4:7,1)
#pmap((x,y)->Raj_PB_py(x,4,5,y,1),[1,5,6],[2,3,5])
#map((x,y)->Raj_PB_py(x,4,5,y,1),[1,5,6],[2,3,5])

# Raj_PB_py(4,500,100,1000,1)
# Raj_PB(1:4,500,100,4:7,1)
# pmap(x->Raj_PB_py(x,500,100,1000,1): x,1:4)
# pmap((x,y)->Raj_PB_py(x,500,100,y,1),1:4,4:7)


# Raj_PB_py([4,6],500,100,1000 .*[0.5,0.8],1)

#Raj_PB(15,0.1,0.8,0.31)
#Raj_PB([15,12,1],0.1,0.8,[0.5,0.01,0.6])

# vals=[15,12,1]
# betaa=[0.5,0.01,0.6]
# x=[1,0.5,0.01]




function loglikelihood_Raj(x; vals,BETA_vec,nemoaccur=256)
  kon = x[1]
  koff = x[2]
  ksyn = x[3]

  tempval=  Raj_PB_raw(vals,kon,koff,ksyn *BETA_vec)
  #tempval=  Raj_PB(vals,kon,koff,ksyn*BETA_vec,nemoaccur=nemoaccur)

  if sum(isnan.(tempval))==length(tempval)
    outtt=Inf
    return outtt
  else
    if sum(isnan.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isnan.(tempval))]
    end
    if sum(isinf.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isinf.(tempval))]
    end
    tempval.=tempval .+ 1.0e-10
    tempval = tempval[findall(tempval .> 0)]
    

    outtt=-sum(log.(tempval))
    return outtt
  #return convert.(Float64,real(outtt))
  end
end


###analytic form of BP distribution
function MLE_LBFGS(;Data,BETA_vec, mode_allele=true,lxa = [1e-3,1e-3,1e-3],uxa = [100.,1e10,1e10])

  store_dat_all=SharedArray{Float64,2}(size(Data)[1], 6)

  @distributed vcat for Gene_ind in range(1, stop = size(Data)[1])
      print(Gene_ind)

      vals=Data[Gene_ind,:]
      if mode_allele==true
        init_guessa=MomentInference_RCPP(vals./BETA_vec )
      else
        init_guessa=MomentInference_sum2alleles(vals=vals./BETA_vec )
      
      end


      if (init_guessa[1]<0) | (init_guessa[2]<0) | (init_guessa[3]<0)
          init_guessa=[10,10,10]
      end

      if isnan(init_guessa[1]) | isnan(init_guessa[2]) | isnan(init_guessa[3]) | (sum(vals)==0)  
          init_guessa=[10,10,10]
      end

      try 
          function min_function(x)
              #xx=loglikelihood_larssonBETA(x; vals=vals,BETA_vec=BETA_vec)
              xx=loglikelihood_Raj(x; vals=vals,BETA_vec=BETA_vec,nemoaccur=64)
          return xx
          end
          result = optimize(min_function,lxa ,uxa , init_guessa, Fminbox(LBFGS()), Optim.Options(g_tol = 1e-8,iterations = 10))
          push!(result.minimizer,  result.minimizer[3]/result.minimizer[2] , Optim.minimum(result),Optim.converged(result))
          store_dat_all[Gene_ind,:]=result.minimizer
      catch e
          store_dat_all[Gene_ind,:]=[NaN,NaN, NaN,NaN, NaN,NaN]
      end

  end
  return store_dat_all

end




###larsson's Code
function Pf(at,m)
  if maximum(m)<1e6
      X=Poisson.(m)
      return(broadcast(pdf,X,at'))
  else
      X=Normal.(m,sqrt.(m))
      return(broadcast(pdf,X,at'))
  end
end

function dBP(at,alpha,bet,lam)
  x,w = gaussjacobi(50,bet - 1, alpha - 1)
  gs = sum(w .*Pf(at, lam.*(1 .+x) ./2),dims=1)
  prob = 1 ./beta(alpha, bet) .*2^(-alpha-bet +1.0) .*gs
  return prob

end

function dBP_vec(at,alpha,bet,lam)
  x,w = gaussjacobi(50,bet - 1, alpha - 1)
  gs_vec=sum(w' .*  Pf.(at, (lam  .* (1 .+x)' ./2)),dims=2)
  prob_vec = 1 ./beta(alpha, bet) .*2^(-alpha-bet +1.0) .*gs_vec
  return prob_vec

end

# sum(w .*Pf(at, lam.*BETA_vec[1].*(1 .+x) ./2),dims=1)
# sum(w' .*  Pf.(vals, (lam .* BETA_vec  .* (1 .+x)' ./2)),dims=2)[1]

# Pf.(at, (lam  .* (1 .+x)' ./2))
# Pf.(vals, (lam .* BETA_vec  .* (1 .+x)' ./2))[1,:]
  



function loglikelihood_larsson(x ; vals)
  kon = x[1]
  koff = x[2]
  ksyn = x[3]
  tempval=  dBP(vals,kon,koff,ksyn)
  if sum(isnan.(tempval))==length(tempval)
    outtt=Inf
    return outtt
  else
    if sum(isnan.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isnan.(tempval))]
    end
    if sum(isinf.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isinf.(tempval))]
    end
    tempval.=tempval .+ 1.0e-10

    outtt=-sum(log.(tempval))
    return outtt
  end
end
  
function loglikelihood_larssonBETA(x ; vals, BETA_vec)
  kon = x[1]
  koff = x[2]
  ksyn = x[3]
  tempval=repeat([0.],length(vals))
  for i in 1:length(vals)
    tempval[i] =  dBP(vals[i],kon,koff,ksyn*BETA_vec[i])[1]
  end
  if sum(isnan.(tempval))==length(tempval)
    outtt=Inf
    return outtt
  else
    if sum(isnan.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isnan.(tempval))]
    end
    if sum(isinf.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isinf.(tempval))]
    end
    tempval.=tempval .+ 1.0e-10

    outtt=-sum(log.(tempval))
    return outtt
  end
end
  
  function MLE_LBFGS_larsson(;Data,lxa = [1e-3,1e-3,1e-3],uxa = [100.,1e10,1e10])
  
    store_dat_all=SharedArray{Float64,2}(size(Data)[1], 6)
  
    @distributed vcat for Gene_ind in range(1, stop = size(Data)[1])
        print(Gene_ind)
  
        vals=Data[Gene_ind,:]
        init_guessa=MomentInference_RCPP(vals )
  
  
        if (init_guessa[1]<0) | (init_guessa[2]<0) | (init_guessa[3]<0)
            init_guessa=[10,10,10]
        end
  
        if isnan(init_guessa[1]) | isnan(init_guessa[2]) | isnan(init_guessa[3]) | (sum(vals)==0) 
          init_guessa=[10,10,10]
        end
  
        try 
            function min_function(x)
                xx=loglikelihood_larsson(x; vals=vals)
            return xx
            end
            result = optimize(min_function,lxa ,uxa , init_guessa, Fminbox(LBFGS()), Optim.Options(g_tol = 1e-8,iterations = 10))
            #result = optimize(min_function,lxa ,uxa , init_guessa, Fminbox(LBFGS()))
            #store_dat_all[Gene_ind,:]=[result.minimizer Optim.minimum(result) Optim.converged(result)]
            push!(result.minimizer,  result.minimizer[3]/result.minimizer[2] , Optim.minimum(result),Optim.converged(result))
            store_dat_all[Gene_ind,:]=result.minimizer
        catch e
            store_dat_all[Gene_ind,:]=[NaN,NaN, NaN,NaN, NaN,NaN]
        end
  
    end
    return store_dat_all
  
  end
  
  function MLE_LBFGS_larssonBETA(;Data,BETA_vec, mode_allele=true,lxa = [1e-3,1e-3,1e-3],uxa = [100.,1e10,1e10])
  
    store_dat_all=SharedArray{Float64,2}(size(Data)[1], 6)
  
    @distributed vcat for Gene_ind in range(1, stop = size(Data)[1])
        print(Gene_ind)
  
        vals=Data[Gene_ind,:]
        if mode_allele==true
          init_guessa=MomentInference_RCPP(vals./BETA_vec )
        else
          init_guessa=MomentInference_sum2alleles(vals=vals./BETA_vec )
        
        end
  
  
        if (init_guessa[1]<0) | (init_guessa[2]<0) | (init_guessa[3]<0)
            init_guessa=[10,10,10]
        end
  
        if isnan(init_guessa[1]) | isnan(init_guessa[2]) | isnan(init_guessa[3]) | (sum(vals)==0)  
            init_guessa=[10,10,10]
        end
  
        try 
            function min_function(x)
                xx=loglikelihood_larssonBETA(x; vals=vals,BETA_vec=BETA_vec)
            return xx
            end
            result = optimize(min_function,lxa ,uxa , init_guessa, Fminbox(LBFGS()), Optim.Options(g_tol = 1e-8,iterations = 10))
            #result = optimize(min_function,lxa ,uxa , init_guessa, Fminbox(LBFGS()))
            #store_dat_all[Gene_ind,:]=[result.minimizer Optim.minimum(result) Optim.converged(result)]
            push!(result.minimizer,  result.minimizer[3]/result.minimizer[2] , Optim.minimum(result),Optim.converged(result))
            store_dat_all[Gene_ind,:]=result.minimizer
        catch e
            store_dat_all[Gene_ind,:]=[NaN,NaN, NaN,NaN, NaN,NaN]
        end
  
    end
    return store_dat_all
  
  end


#integration

# at=vals[1]
# i=1
# kon=alpha
# koff=bet
# ksyn=lam

# i=3
# aa=dBP_vec(vals,kon,koff,ksyn .*BETA_vec)[i]
# bb=dBP(vals[i],kon,koff,ksyn*BETA_vec[i])[1]
# isapprox( aa, bb)

# dBP_vec(vals,kon,koff,ksyn .*BETA_vec)[i]==dBP(vals[i],kon,koff,ksyn*BETA_vec[i])[1]




# aa=[1,2,3]
# bb=[[1,2,3] [1,2,3] [1,2,3]]
# aa' .* bb



function loglikelihood_all(x ; vals, BETA_vec, density_fun="integration")
  kon  = x[1]
  koff = x[2]
  ksyn = x[3]

  if density_fun=="integration"
    tempval=repeat([0.],length(vals))
    for i in 1:length(vals)
      tempval[i] =  dBP(vals[i],kon,koff,ksyn*BETA_vec[i])[1]
    end
  elseif density_fun=="integration_vec"
    tempval=dBP_vec(vals,kon,koff,ksyn*BETA_vec)
  else
    tempval=Raj_PB_raw(vals,kon,koff,ksyn *BETA_vec)
  end


  if sum(isnan.(tempval))==length(tempval)
    outtt=Inf
    return outtt
  else
    if sum(isnan.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isnan.(tempval))]
    end
    if sum(isinf.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isinf.(tempval))]
    end
    tempval.=tempval .+ 1.0e-10

    outtt=-sum(log.(tempval))
    return outtt
  end
end



function loglikelihood_umi(x ; vals, BETA_vec, density_fun="integration")
  kon  = x[1]
  koff = x[2]
  ksyn = x[3]

  tempval=repeat([0.],length(vals))

  
  for i in 1:length(vals)
    obs=vals[i]
    bett=BETA_vec[i]

    mult_sum=0.
    if obs%2==0
      for j in 0:obs/2
        p1=dBP_vec(j,kon,koff,ksyn*bett)
        p2=dBP_vec(obs-j,kon,koff,ksyn*bett)
        mult=p1*p2
        if j!=(obs-j)
          mult_sum=mult_sum+mult[1,1]*2
        else
          mult_sum=mult_sum+mult[1,1]
        end
      end#end of even number
    else
      for j in 0:floor(obs/2)
        p1=dBP_vec(j,kon,koff,ksyn*bett)
        p2=dBP_vec(obs-j,kon,koff,ksyn*bett)
        mult=p1*p2
        mult_sum=mult_sum+mult[1,1]*2
      end#end of odd number
    end #end of one cell
    tempval[i]=mult_sum
  end


  # mult_sum2=0.
  # for j in 0:obs
  #   p1=dBP_vec(j,kon,koff,ksyn*bett)
  #   p2=dBP_vec(obs-j,kon,koff,ksyn*bett)
  #   mult=p1*p2
  #   mult_sum2=mult_sum2+mult[1,1]
  # end

  if sum(isnan.(tempval))==length(tempval)
    outtt=Inf
    return outtt
  else
    if sum(isnan.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isnan.(tempval))]
    end
    if sum(isinf.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isinf.(tempval))]
    end
    tempval.=tempval .+ 1.0e-10

    outtt=-sum(log.(tempval))
    return outtt
  end
end


function loglikelihood_umivec(x ; vals, BETA_vec, density_fun="integration")
  kon  = x[1]
  koff = x[2]
  ksyn = x[3]

  tempval=repeat([0.],length(vals))

  for i in 1:length(vals)
    obs=vals[i]
  
    if obs%2==0
      valle=range(0,(obs/2-1))
      bett_vec=repeat([ksyn*BETA_vec[i]],length(valle))
      mult1=dBP_vec(valle,kon,koff,bett_vec)
      mult2=dBP_vec((obs .- valle),kon,koff,bett_vec)
      mult_sum=sum(mult1 .* mult2) *2 +dBP_vec(obs/2,kon,koff,ksyn*BETA_vec[i])[1,1]^2

    else
      valle=range(0,floor(obs/2))
      bett_vec=repeat([ksyn*BETA_vec[i]],Int(length(valle)))
      mult1=dBP_vec(valle,kon,koff,bett_vec)
      mult2=dBP_vec((obs .- valle),kon,koff,bett_vec)

      mult_sum=sum(mult1 .* mult2) *2
    end #end of one cell


    tempval[i]=mult_sum
  end

  if sum(isnan.(tempval))==length(tempval)
    outtt=Inf
    return outtt
  else
    if sum(isnan.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isnan.(tempval))]
    end
    if sum(isinf.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isinf.(tempval))]
    end
    tempval.=tempval .+ 1.0e-10

    outtt=-sum(log.(tempval))
    return outtt
  end
end






function MLE_LBFGS_all(;Data,BETA_vec, mode_allele=true, density_fun="integration",lxa = [1e-3,1e-3,1e-3],uxa = [100.,1e10,1e10])

  store_dat_all=SharedArray{Float64,2}(size(Data)[1], 6)

  @distributed vcat for Gene_ind in range(1, stop = size(Data)[1])
      print(Gene_ind)

      vals=Data[Gene_ind,:]
      if mode_allele==true
        init_guessa=MomentInference_RCPP(vals./BETA_vec )
        likelihood_fun=loglikelihood_all
      else
        init_guessa=MomentInference_sum2alleles(vals=vals./BETA_vec )
        likelihood_fun=loglikelihood_umi
      
      end

      if (init_guessa[1]<0) | (init_guessa[2]<0) | (init_guessa[3]<0)
          init_guessa=[10,10,10]
      end

      if isnan(init_guessa[1]) | isnan(init_guessa[2]) | isnan(init_guessa[3]) | (sum(vals)==0)  
          init_guessa=[10,10,10]
      end

      try 
          function min_function(x)
              xx=likelihood_fun(x; vals=vals,BETA_vec=BETA_vec,density_fun=density_fun)
          return xx
          end
          result = optimize(min_function,lxa ,uxa , init_guessa, Fminbox(LBFGS()), Optim.Options(g_tol = 1e-8,iterations = 10))
          push!(result.minimizer,  result.minimizer[3]/result.minimizer[2] , Optim.minimum(result),Optim.converged(result))
          store_dat_all[Gene_ind,:]=result.minimizer
      catch e
          store_dat_all[Gene_ind,:]=[NaN,NaN, NaN,NaN, NaN,NaN]
      end

  end
  return store_dat_all

end



function MLE_LBFGS_bounded(;Data,BETA_vec, mode_allele=true, density_fun="integration")
  threshold_Data, interval_1, interval_2, interval_3, M1_boot, M2_boot, M3_boot= ABC_initialization_v2(Data=Data,BETA_vec=BETA_vec,num_boot=1000,threshold_on_distance=0.95,mode_allele=mode_allele,verbose=true)


  dat_lower=[interval_1[:,1] interval_2[:,1] interval_3[:,1]]./10
  dat_upper=[interval_1[:,2] interval_2[:,2] interval_3[:,2]].*10

  @everywhere dat_lower=$dat_lower
  @everywhere dat_upper=$dat_upper

  store_dat_all=SharedArray{Float64,2}(size(Data)[1], 6)

  @distributed vcat for Gene_ind in range(1, stop = size(Data)[1])
      print(Gene_ind)

      vals=Data[Gene_ind,:]
      if mode_allele==true
        init_guessa=MomentInference_RCPP(vals./BETA_vec )
      else
        init_guessa=MomentInference_sum2alleles(vals=vals./BETA_vec )
      
      end

      if (init_guessa[1]<0) | (init_guessa[2]<0) | (init_guessa[3]<0)
          init_guessa=[10,10,10]
      end

      if isnan(init_guessa[1]) | isnan(init_guessa[2]) | isnan(init_guessa[3]) | (sum(vals)==0)  
          init_guessa=[10,10,10]
      end

      try 
          function min_function(x)
              xx=loglikelihood_all(x; vals=vals,BETA_vec=BETA_vec,density_fun=density_fun)
          return xx
          end
          lxa=dat_lower[Gene_ind,:]
          uxa=dat_upper[Gene_ind,:]

          result = optimize(min_function,lxa ,uxa , init_guessa, Fminbox(LBFGS()), Optim.Options(g_tol = 1e-8,iterations = 10))
          push!(result.minimizer,  result.minimizer[3]/result.minimizer[2] , Optim.minimum(result),Optim.converged(result))
          store_dat_all[Gene_ind,:]=result.minimizer
      catch e
          store_dat_all[Gene_ind,:]=[NaN,NaN, NaN,NaN, NaN,NaN]
      end

  end
  return store_dat_all

end



function loglikelihood_all_fano(x ; vals, BETA_vec, density_fun="integration")
  kon  = x[1]
  koff = x[2]
  ksyn = x[3]
  fano=var(vals)/mean(vals)
  sim_fano=ksyn*koff/(kon+koff)/(kon+koff+1)

  if density_fun=="integration"
    tempval=repeat([0.],length(vals))
    for i in 1:length(vals)
      tempval[i] =  dBP(vals[i],kon,koff,ksyn*BETA_vec[i])[1]
    end
  elseif density_fun=="integration_vec"
    tempval=dBP_vec(vals,kon,koff,ksyn*BETA_vec)
  else
    tempval=Raj_PB_raw(vals,kon,koff,ksyn *BETA_vec)
  end


  if sum(isnan.(tempval))==length(tempval)
    outtt=Inf
    return outtt
  else
    if sum(isnan.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isnan.(tempval))]
    end
    if sum(isinf.(tempval))>0
        tempval=tempval[.!convert(Array{Bool},isinf.(tempval))]
    end
    tempval.=tempval .+ 1.0e-10

    outtt=-sum(log.(tempval))+abs(fano-sim_fano)
    return outtt
  end
end

function MLE_LBFGS_bounded2(;Data,BETA_vec, mode_allele=true, density_fun="integration")
  threshold_Data, interval_1, interval_2, interval_3, M1_boot, M2_boot, M3_boot= ABC_initialization_v2(Data=Data,BETA_vec=BETA_vec,num_boot=1000,threshold_on_distance=0.95,mode_allele=mode_allele,verbose=false)
  
  if mode_allele==true
    MME_out=MomentInference_allele(Data=Data,BETA_vec=BETA_vec)
    sim_fun=Data_sim
    likelihood_fun=loglikelihood_all
  else
    MME_out=MomentInference_UMI(Data=Data,BETA_vec=BETA_vec)
    sim_fun=Data_sim_sum2
    likelihood_fun=loglikelihood_umivec
  end


  store_dat_all=SharedArray{Float64,2}(size(Data)[1], 6)

  @distributed vcat for Gene_ind in range(1, stop = size(Data)[1])
      print(Gene_ind)

      vals=Data[Gene_ind,:]
      cellsused=findall(.!ismissing.(Data[Gene_ind,:]))

      MME_init=MME_out[Gene_ind,:]
            
      new_param_space=extract_pointInterval_v3(;MME_ori=MME_out,interval_1= interval_1,interval_2= interval_2,interval_3= interval_3,M1_boot=M1_boot, M2_boot=M2_boot,M3_boot=M3_boot,gene_ind=Gene_ind,particles=10000,mode_allele=mode_allele)

      if size(new_param_space)[1]==0
          preselect=[]
          init_guessa=[10.,10.,10.]
          lxa=[1e-3,1e-3,1e-3]
          uxa = [100.,1e10,1e10]
      else
          preselect=findall( (new_param_space[:,1] .>0) .& (new_param_space[:,2] .>0) .& (new_param_space[:,3] .>0))
      end

      
      if length(preselect)>1000
          ww=fill(1/length(preselect), length(preselect))
          ranind=wsample(range(1,stop=length(preselect)),ww,1000,replace=false)
          preselect=preselect[ranind]
      end

      if length(preselect)==0
        init_guessa=[10.,10.,10.]
        lxa=[1e-3,1e-3,1e-3]
        uxa = [100.,1e10,1e10]
      else
          new_param_space=new_param_space[preselect,:]
          temp_dat=sim_fun(Param_Space=new_param_space,BETA_vec=BETA_vec[cellsused])
          ff=ABC_test(temp_dat, Data[Gene_ind,:])
          thresholdd=quantile(ff,0.01)
          particles_selected=findall( (ff .<= thresholdd) )
          population_dat_3=new_param_space[particles_selected,:]
          lxa=minimum(population_dat_3,dims=1)[1,:]
          uxa=maximum(population_dat_3,dims=1)[1,:]
          init_guessa=extract_ml_particle_wt(transpose(population_dat_3))
      end

      try 
          function min_function(x)
              xx=likelihood_fun(x; vals=vals,BETA_vec=BETA_vec,density_fun=density_fun)
              #xx=loglikelihood_all_fano(x; vals=vals,BETA_vec=BETA_vec,density_fun=density_fun)
          return xx
          end


          result = optimize(min_function,lxa ,uxa , init_guessa, Fminbox(LBFGS()), Optim.Options(g_tol = 1e-8,iterations = 10))
          #aa=[result.minimizer[:,1]  result.minimizer[3]/result.minimizer[2]  Optim.minimum(result) Optim.converged(result)]
          aa=result.minimizer[:,1]
          push!(aa, result.minimizer[3]/result.minimizer[2] , Optim.minimum(result), Optim.converged(result))

          store_dat_all[Gene_ind,:]=aa
      catch e
          store_dat_all[Gene_ind,:]=[NaN,NaN, NaN,NaN, NaN,NaN]
      end

  end
  return store_dat_all

end






function MLE_LBFGS_bounded3(;Data,BETA_vec, mode_allele=true, density_fun="integration")
  threshold_Data, interval_1, interval_2, interval_3, M1_boot, M2_boot, M3_boot= ABC_initialization_v2(Data=Data,BETA_vec=BETA_vec,num_boot=1000,threshold_on_distance=0.95,mode_allele=mode_allele,verbose=false)
  
  
  if mode_allele==true
    MME_out=MomentInference_allele(Data=Data,BETA_vec=BETA_vec)
    sim_fun=Data_sim
  else
      MME_out=MomentInference_UMI(Data=Data,BETA_vec=BETA_vec)
      sim_fun=Data_sim_sum2
  end


  store_dat_all=SharedArray{Float64,2}(size(Data)[1], 6)

  # interval_1=[mapslices(X -> minimum_fun(X), M1_boot, dims=(2)) mapslices(X -> maximum_fun(X), M1_boot, dims=(2))]
  # interval_2=[mapslices(X -> minimum_fun(X), M2_boot, dims=(2)) mapslices(X -> maximum_fun(X), M2_boot, dims=(2))]
  # interval_3=[mapslices(X -> minimum_fun(X), M3_boot, dims=(2)) mapslices(X -> maximum_fun(X), M3_boot, dims=(2))]

  
  interval_1=[mapslices(X -> quantile_fun(X, 0.05), M1_boot, dims=(2)) mapslices(X -> quantile_fun(X, 0.95), M1_boot, dims=(2))]
  interval_2=[mapslices(X -> quantile_fun(X, 0.05), M2_boot, dims=(2)) mapslices(X -> quantile_fun(X, 0.95), M2_boot, dims=(2))]
  interval_3=[mapslices(X -> quantile_fun(X, 0.05), M3_boot, dims=(2)) mapslices(X -> quantile_fun(X, 0.95), M3_boot, dims=(2))]

  @distributed vcat for Gene_ind in range(1, stop = size(Data)[1])
      print(Gene_ind)

      vals=Data[Gene_ind,:]
      cellsused=findall(.!ismissing.(Data[Gene_ind,:]))
      temppp=[interval_1[Gene_ind,:] interval_2[Gene_ind,:] interval_3[Gene_ind,:]]
      lxa=temppp[1,:] 
      uxa=temppp[2,:] 
      MME_init=MME_out[Gene_ind,:]
      
      

      # if (sum(isnan.(MME_init))==0) & (sum(MME_init .>0)==3) & (sum(isnan.(lxa))==0) & (sum(isnan.(uxa))==0)
      if (sum((MME_init - lxa) .>0)==3) & (sum((uxa- MME_init) .>0)==3)
        init_guessa=MME_init
        lxa=temppp[1,:] ./10
        uxa=temppp[2,:] .*10

      else
        new_param_space=extract_pointInterval_v3(;MME_ori=MME_out,interval_1= interval_1,interval_2= interval_2,interval_3= interval_3,M1_boot=M1_boot, M2_boot=M2_boot,M3_boot=M3_boot,gene_ind=Gene_ind,particles=100000,mode_allele=mode_allele)

        if size(new_param_space)[1]==0
            preselect=[]
            init_guessa=[10.,10.,10.]
            lxa=[1e-3,1e-3,1e-3]
            uxa = [100.,1e10,1e10]
        else
            preselect=findall( (new_param_space[:,1] .>0) .& (new_param_space[:,2] .>0) .& (new_param_space[:,3] .>0))
        end
  

        if length(preselect)>10000
            ww=fill(1/length(preselect), length(preselect))
            ranind=wsample(range(1,stop=length(preselect)),ww,10000,replace=false)
            preselect=preselect[ranind]
        end

        if length(preselect)==0
          init_guessa=[10.,10.,10.]
          lxa=[1e-3,1e-3,1e-3]
          uxa = [100.,1e10,1e10]
        else
          new_param_space=new_param_space[preselect,:]
          temp_dat=sim_fun(Param_Space=new_param_space,BETA_vec=BETA_vec[cellsused])
          ff=ABC_test(temp_dat, Data[Gene_ind,:])
          thresholdd=quantile(ff,0.01)
          particles_selected=findall( (ff .<= thresholdd) )
          population_dat_3=new_param_space[particles_selected,:]
          lxa=minimum(population_dat_3,dims=1)[1,:]
          uxa=maximum(population_dat_3,dims=1)[1,:]
          init_guessa=extract_ml_particle_wt(transpose(population_dat_3))
        end


      end # end of setting initial guess

      try 
          function min_function(x)
              xx=loglikelihood_all(x; vals=vals,BETA_vec=BETA_vec,density_fun=density_fun)
              #xx=loglikelihood_all_fano(x; vals=vals,BETA_vec=BETA_vec,density_fun=density_fun)
          return xx
          end
          result = optimize(min_function,lxa ,uxa , init_guessa, Fminbox(LBFGS()), Optim.Options(g_tol = 1e-8,iterations = 10))
          #aa=[result.minimizer[:,1]  result.minimizer[3]/result.minimizer[2]  Optim.minimum(result) Optim.converged(result)]
          aa=result.minimizer[:,1]
          push!(aa, result.minimizer[3]/result.minimizer[2] , Optim.minimum(result), Optim.converged(result))

          store_dat_all[Gene_ind,:]=aa
      catch e
          store_dat_all[Gene_ind,:]=[NaN,NaN, NaN,NaN, NaN,NaN]
      end

  end
  return store_dat_all

end