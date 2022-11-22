# Julia_ABC

Source all `.jl` files and `using` all packages listed in `Julia_ABC.jl`, then run `ABC_BP_v2`.

## Quick start
```julia
#A toy scRNA-seq data
Daa=[1:10 1:10 3:12 4:13] 
Bea=[0.1,0.5,0.3,0.4]
outt=ABC_BP_v2(Data=Daa,BETA_vec=Bea, mode_allele=true, particles=1000,threshold_on_distance=0.95,verbose=false)
    
```
