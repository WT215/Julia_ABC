function hellinger2(p, q)
    return  sqrt(sum((sqrt.(p)-sqrt.(q)).^2)) / sqrt(2.)
end

function hellinger2_data(p, q)
  nrows=size(p)[1]
  Dist_mat=Array{Float64,2}(undef,nrows,1)
  for i in range(1,stop=nrows)
    Dist_mat[i,1]= ABC_OneParam_test(p[i,:],q[i,:])
  end
  #Dist_mat=sqrt.(sum((sqrt.(p)-sqrt.(q)) .^2,dims=2)) ./ sqrt(2.)
  return  Dist_mat
end



function mode(dens::UnivariateKDE)
    ind=findmax(dens.density)[2]
    dens.x[ind]
end

function mode_univar(samp)
  md=Array{Float64}(undef,size(samp)[1])
  for i in 1:size(samp)[1]
    dens=kde(samp[i,:])
    md[i]=mode(dens)
  end
  return(md)
end


# function extract_ml_particle_wt(samp)
#     temp=samp
#     md=mode_univar(temp)
#     sepmat=broadcast(-,temp,md)
#     for i in 1:size(sepmat)[1]
#       sepmat[i,:]=sepmat[i,:]./std(sepmat[i,:])
#     end
#     dists=sqrt.(sum(sepmat.^2,dims=1))
#     ind=findmin(dists)
#     #return((temp[:,ind[2][2]],ind[1]))
#     return temp[:,ind[2][2]]
# end

function extract_ml_particle_wt(samp)
  temp=samp
  aa=median(temp,dims=2)
  return aa
end


