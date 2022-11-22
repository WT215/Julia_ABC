
#= Pkg.add("RDatasets")
Pkg.add("JuliaDB")
Pkg.add("CodecBzip2")
Pkg.add("SpecialFunctions")
Pkg.add("Distributed")
Pkg.add("ProgressMeter")
Pkg.add("SharedArrays")
Pkg.add("Random")
Pkg.add("Distributions")
Pkg.add("Optim")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("FreqTables")
Pkg.add("Plots")
Pkg.add("SparseArrays")
Pkg.add("MatrixMarket")
Pkg.add("StatsBase")
Pkg.add("GSL")
Pkg.add("FastGaussQuadrature")
Pkg.add("Nemo")
Pkg.add("LinearAlgebra")
Pkg.add("CovarianceEstimation")
Pkg.add("Distributions")
Pkg.add("LsqFit")
Pkg.add("DelimitedFiles")
Pkg.add("Sundials")
Pkg.add("Dierckx")
Pkg.add("DataStructures")
Pkg.add("JLD2")
Pkg.add("DiffEqBase")
Pkg.add("KernelDensity")
Pkg.add("FileIO")
Pkg.add("JSON")
Pkg.add("Distances")
Pkg.add("CovarianceEstimation")
Pkg.add("OptimalTransport")
Pkg.add("Tulip")
 =#




 
#using BlackBoxOptim
#using StatsKit
#using PyCall
#using Conda
#using BlackBoxOptim: num_func_evals
#using Cubature
#import Nemo
#PyCall.pyimport_conda("sympy", "sympy")
#sympy = PyCall.pyimport("sympy")

#include("D:/RNAseqProject/Julia_allele/ABC_wt/ABC_metrics.jl")
#include("D:/RNAseqProject/Julia_allele/ABC_wt/ABC_test.jl")
#include("D:/RNAseqProject/Julia_allele/ABC_wt/ABC_sim.jl")
#include("D:/RNAseqProject/Julia_allele/ABC_wt/ABC_BP.jl")
#include("D:/RNAseqProject/Julia_allele/ABC_wt/ABC_BP_v2.jl")
#include("D:/RNAseqProject/Julia_allele/ABC_wt/MLE.jl")

# include("/rdsgpfs/general/user/wt215/home/wt_ABC/ABC_wt/ABC_metrics.jl")
# include("/rdsgpfs/general/user/wt215/home/wt_ABC/ABC_wt/ABC_test.jl")
# include("/rdsgpfs/general/user/wt215/home/wt_ABC/ABC_wt/ABC_sim.jl")
# include("/rdsgpfs/general/user/wt215/home/wt_ABC/ABC_wt/ABC_BP.jl")
# include("/rdsgpfs/general/user/wt215/home/wt_ABC/ABC_wt/ABC_BP_v2.jl")
# include("/rdsgpfs/general/user/wt215/home/wt_ABC/ABC_wt/MLE.jl")

# hellinger2([1.,2.,3.],[3.,2.,3.])



function minimum_fun(X)
  et=findall(.!isnan.(X))
  if length(et)>0
      return minimum(X[findall(.!isnan.(X))])
  else
      return NaN
  end
end

function maximum_fun(X)
  et=findall(.!isnan.(X))
  if length(et)>0
      return maximum(X[findall(.!isnan.(X))])
  else
      return NaN
  end
end


function median_fun(X)
  et=findall(.!isnan.(X))
  if length(et)>0
      return median(X[findall(.!isnan.(X))])
  else
      return NaN
  end
end

function quantile_fun(X,quantile_num)
  et=findall(.!isnan.(X))
  if length(et)>0
      return quantile(X[findall(.!isnan.(X))],quantile_num)
  else
      return NaN
  end
end

function std_fun(X)
  et=findall(.!isnan.(X))
  if length(et)>0
      return std(X[findall(.!isnan.(X))])
  else
      return NaN
  end
end



function rbp_self(sizee,alpha, beta, c)
    simout=rand.(Poisson.(rand(Beta(alpha,beta),sizee).*c))
    return simout
end



function DownSampling(;Data,BETA_vec)
  Counts_downsampling=Array{Int64,2}(undef, size(Data))

  for i in 1:size(Data)[1]
      for j in 1:size(Data)[2]
          d=Binomial(Int(Data[i,j]),BETA_vec[j])
          Counts_downsampling[i,j]= rand(d, 1)[1]
      end
  end

  return(Counts_downsampling)
end



function SimData(;Param_Space,BETA_vec)
  gene_num=size(Param_Space)[1]
  cell_numbers=length(BETA_vec)

  Dat_real=zeros(gene_num,cell_numbers)
  for ii=1:gene_num
      tempvec=rbp_self(cell_numbers,LAR_true_new[ii,1],LAR_true_new[ii,2],LAR_true_new[ii,3])
      Dat_real[ii,:]=tempvec
  end
  Dat_down=DownSampling(Data=Dat_real,BETA_vec=BETA_vec)

return(Dat_down)
end



function empirBETA(;Data,meanBeta=0.06)
  BETA_vec = sum(Data,dims=1)./mean(sum(Data,dims=1)).*meanBeta
  BETA_vec[BETA_vec.<=0].=minimum(BETA_vec[BETA_vec.>0])
  BETA_vec[BETA_vec.>=1].=maximum(BETA_vec[BETA_vec.<1])
  BETA_vec=BETA_vec[1,:]
  return BETA_vec
end




function cpp_gmRNA_switch(n, r_act, r_deact,r_on, r_degr) 

  if(!isa(n, Int)) 
    return (0);
  end

  res=Array{Int,1}(undef,n);
  t0 = 0.
  tmax = 20. / r_degr 
  x0=[1.,0.,0.]
  x=Array{Float64,1}(undef,3);

  for i = 1:n 
    tx = t0;
    x = x0;
    #step 1
    r_act1 = r_act * x[1];
    r_act2 = r_deact * x[2];
    r_act3 = r_on * x[2];
    r_act4 = r_degr * x[3];
    r_actx = r_act1 + r_act2 + r_act3 + r_act4;

    # step 2
    dd=Exponential(r_actx)
    #tau_vec = rexp(1, r_actx);
    tau_vec=rand(dd,1)[1];
    tau = tau_vec[1];
    tau_stern = minimum([tau, tmax - tx]);
    tx = tx+ tau_stern;
    dd_u=Uniform(0,1)
    while tx < tmax
      # step 3
      u_vec = rand(dd_u,1);
      u = u_vec[1];
      if u <= r_act1/r_actx
        x[1]=x[1]-1
        x[2]=x[2]+1
        break;

      elseif (u <= (r_act1+r_act2)/r_actx)
        x[1]=x[1]+1;
        x[2]=x[2]-1;
        break;
      elseif (u <= (r_act1+r_act2+r_act3)/r_actx)
        x[3]=x[3]+1;
        break;
      else
        x[3]=x[3]-1;
        break;
      end
      # step 5
      r_act1 = r_act * x[1];
      r_act2 = r_deact * x[2];
      r_act3 = r_on * x[2];
      r_act4 = r_degr * x[3];
      r_actx = r_act1 + r_act2 + r_act3 + r_act4;

      # step 6
      dd3=Exponential(r_actx)
      tau_vec = rand(dd3,1)[1];
      tau = tau_vec[1];
      tau_stern = minimum([tau, tmax - tx]);
      tx = tx+ tau_stern;
    end # end while
      res[i] = x[3];
    end#end for

  return res;
end


function kinetic_to_moments(;kon,koff,ksyn,mode_allele=true)
  if mode_allele==true
      r1=kon .* ksyn ./ (kon  .+koff)
      r2=(kon .*ksyn .+ksyn) ./(kon .+koff .+1)
      r3=(kon.*ksyn .+ 2 .*ksyn) ./(kon .+koff .+2)
      M1=r1
      M2=r1 .*r2
      M3=r1 .*r2 .*r3
  else
      r1=kon .* ksyn ./ (kon  .+koff)
      r2=(kon .*ksyn .+ksyn) ./(kon .+koff .+1)
      r3=(kon.*ksyn .+ 2 .*ksyn) ./(kon .+koff .+2)

      # m1=r1
      # m2=r1 .*r2
      # m3=r1 .*r2 .*r3
      # M1=2 .* m1
      # M2=(m2 .+ m1 .* m1 ) .*2
      # M3= (m3 .+ 3 .* m1 .*m2) .*2

      M1=r1
      M2=r1 .*r2
      M3=r1 .*r2 .*r3
  end

  return [M1 M2 M3]
  
end


function MomentInference_data(;Data,BETA_vec)
  Data_used=transpose(transpose(Data) ./ BETA_vec)
  aa=@distributed (hcat) for i in 1:size(Data_used)[1]
  kk=MomentInference_RCPP(Data_used[i,findall(.!ismissing.(Data_used[i,:]))])

  kk

  end
return transpose(aa)

end


function MomentInference_allele(;Data,BETA_vec=nothing)

  if BETA_vec!=nothing
    Data=transpose(transpose(Data) ./ BETA_vec)
  end

  lenn=size(Data)[2]
  n_genes=size(Data)[1]
  xx=Array{Float64}(undef, 3) 

  m1 = mean(Data,dims=2);
  m2 = sum(Data .* (Data .- 1.),dims=2) ./lenn;
  m3 =   sum(Data .* (Data .- 1.).* (Data .- 2.),dims=2) ./lenn;

  MME_out=Array{Float64,2}(undef,n_genes, 3) 



  conten=findall(sum(Data,dims=2) .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[conten, :] .=fill(NaN, (lx,3))
  end

  conten=findall(m1 .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end

  conten=findall(m2 .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end


  r1=m1;
  r2=m2 ./m1;
  r3=m3 ./m2;


  conten=findall( (r1 .*r2 .-2 .*r1 .*r3 .+ r2 .*r3)  .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end


  conten=findall( ((r1.*r2 .- 2 .*r1 .*r3 .+ r2 .*r3) .*(r1 .-2 .*r2 .+r3))  .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end

  conten=findall( (r1 .- 2 .*r2 .+ r3)  .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end

  lambda_est = (2 .*r1 .*(r3 .-r2)) ./(r1 .*r2-2 .*r1 .*r3 + r2 .*r3);
  mu_est = (2 .*(r3 .-r2) .*(r1 .-r3) .*(r2 .-r1)) ./((r1 .*r2  .- 2 .*r1 .*r3 .+ r2 .*r3) .*(r1 .-2 .*r2 .+r3));
  v_est = (2 .*r1 .*r3 - r1 .*r2 .- r2 .*r3) ./(r1 .- 2 .*r2 .+ r3);
  MME_out[:,1]=lambda_est;
  MME_out[:,2]=mu_est;
  MME_out[:,3]=v_est;


  conten=findall( (MME_out[:,1] .<0) .| (MME_out[:,2] .<0) .| (MME_out[:,3] .<0))
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end


  return MME_out

end



function MomentInference_UMI(; Data,BETA_vec=nothing)

  if BETA_vec!=nothing
      Data=transpose(transpose(Data) ./ BETA_vec)
  end
  lenn=size(Data)[2]
  n_genes=size(Data)[1]


  M1 = mean(Data,dims=2);
  M2 = sum(Data .* (Data .- 1.),dims=2) ./lenn;
  M3 = sum(Data .* (Data .- 1.) .* (Data .- 2.),dims=2) ./lenn;
  m1 = M1 ./2;
  m2 = M2 ./2 .- m1 .*m1;
  m3 = M3 ./2 .- 3 .* m1 .* m2;
  MME_out=Array{Float64,2}(undef,n_genes, 3) 


  conten=findall(sum(Data,dims=2) .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end
  conten=findall(M1 .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end

  conten=findall(M2 .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end


  r1=m1;
  r2=m2 ./m1;
  r3=m3 ./m2;

  conten=findall( (r1 .*r2 .-2 .*r1 .*r3 .+ r2 .*r3)  .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end


  conten=findall( ((r1.*r2 .- 2 .*r1 .*r3 .+ r2 .*r3) .*(r1 .-2 .*r2 .+r3))  .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end

  conten=findall( (r1 .- 2 .*r2 .+ r3)  .==0)
  conten=getindex.(conten,1)
  lx=length(conten)
  if  lx!=0
    MME_out[ conten,:] .=fill(NaN, (lx,3))
  end


  lambda_est = (2 .*r1 .*(r3 .-r2)) ./(r1 .*r2-2 .*r1 .*r3 + r2 .*r3);
  mu_est = (2 .*(r3 .-r2) .*(r1 .-r3) .*(r2 .-r1)) ./((r1 .*r2  .- 2 .*r1 .*r3 .+ r2 .*r3) .*(r1 .-2 .*r2 .+r3));
  v_est = (2 .*r1 .*r3 - r1 .*r2 .- r2 .*r3) ./(r1 .- 2 .*r2 .+ r3);
  MME_out[:,1]=lambda_est;
  MME_out[:,2]=mu_est;
  MME_out[:,3]=v_est;

  return  MME_out

end



function MomentInference_RCPP(vals)
    lenn=length(vals)
    xx=Array{Float64}(undef, 3) 

    m1 = mean(vals);
    m2 = sum(vals .* (vals .- 1.))/lenn;
    m3 = sum(vals .* (vals .- 1.) .* (vals .- 2.))/lenn;

    if sum(vals) == 0
        return [NaN,NaN, NaN];
    end
    
    if m1 == 0
        return [NaN,NaN, NaN];
    end
    if m2 == 0
        return [NaN,NaN, NaN];
    end


    r1=m1;
    r2=m2/m1;
    r3=m3/m2;
    if (r1*r2-2*r1*r3 + r2*r3) == 0
        return [NaN,NaN, NaN];
    end

    if ((r1*r2 - 2*r1*r3 + r2*r3)*(r1-2*r2+r3)) == 0
    return [NaN,NaN, NaN];
    end


    if (r1 - 2*r2 + r3) == 0
        return [NaN,NaN, NaN];
    end

    lambda_est = (2*r1*(r3-r2))/(r1*r2-2*r1*r3 + r2*r3);
    mu_est = (2*(r3-r2)*(r1-r3)*(r2-r1))/((r1*r2 - 2*r1*r3 + r2*r3)*(r1-2*r2+r3));
    v_est = (2*r1*r3 - r1*r2 - r2*r3)/(r1 - 2*r2 + r3);
    xx[1]=lambda_est;
    xx[2]=mu_est;
    xx[3]=v_est;

    if (xx[1]<0) | (xx[2]<0) | (xx[3]<0)
      return [NaN,NaN, NaN];
    else
      return xx;
    end

end

#MomentInference_RCPP(vals)



#MomentInference_sum2alleles(vals=vals)
function MomentInference_sum2alleles(; vals)
  lenn=size(vals)[1]
  xx=Array{Float64}(undef, 3) 

  M1 = mean(vals);
  M2 = sum(vals .* (vals .- 1)) ./lenn;
  M3 = sum(vals .* (vals .- 1) .* (vals .- 2)) ./lenn;
  m1 = M1 ./2;
  m2 = M2 ./2 .- m1 .*m1;
  m3 = M3 ./2 .- 3 .* m1 .* m2;


  if sum(vals) == 0
    xx[1]=NaN;
    xx[2]=NaN;
    xx[3]=NaN;
    return xx;
  end
  if m1 == 0
    xx[1]=NaN;
    xx[2]=NaN;
    xx[3]=NaN;
    return xx;
  end
  if m2 == 0
    xx[1]=NaN;
    xx[2]=NaN;
    xx[3]=NaN;
    return xx;
  end
  r1=m1;
  r2=m2/m1;
  r3=m3/m2;

  if (r1*r2-2*r1*r3 + r2*r3) == 0
    xx[1]=NaN;
    xx[2]=NaN;
    xx[3]=NaN;
    return xx;
  end
  if ((r1*r2 - 2*r1*r3 + r2*r3)*(r1-2*r2+r3)) == 0
    xx[1]=NaN;
    xx[2]=NaN;
    xx[3]=NaN;
    return xx;
  end
  if (r1 - 2*r2 + r3) == 0
    xx[1]=NaN;
    xx[2]=NaN;
    xx[3]=NaN;
    return xx;
  end


  lambda_est = (2*r1*(r3-r2))/(r1*r2-2*r1*r3 + r2*r3);
  mu_est = (2*(r3-r2)*(r1-r3)*(r2-r1))/((r1*r2 - 2*r1*r3 + r2*r3)*(r1-2*r2+r3));
  v_est = (2*r1*r3 - r1*r2 - r2*r3)/(r1 - 2*r2 + r3);

  xx[1]=lambda_est;
  xx[2]=mu_est;
  xx[3]=v_est;

  return xx


end

function MomentInference_sum2alleles_data(;Data,BETA_vec)
  Data_used=transpose(transpose(Data) ./ BETA_vec)
  aa=@distributed (hcat) for i in 1:size(Data_used)[1]
  kk=MomentInference_sum2alleles(vals=Data_used[i,findall(.!ismissing.(Data_used[i,:]))])

  kk

  end
return transpose(aa)

end

