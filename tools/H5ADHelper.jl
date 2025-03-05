# A small module for Write results back to h5ad

module H5Adh

using HDF5

export insertddh

"Replace data if the key exists"
function writeobs(fid,k,x)
    h = "obs/"*k
    haskey(fid, h) && delete_object(fid,h) 
    fid[h] = x
    orders = attrs(fid["obs"])["column-order"]
    if k ∉ orders
        attrs(fid["obs"])["column-order"] = [orders;k]
    end
end

function writeobsm(fid,k,x)
    h = "obsm/"*k
    arr = permutedims(cat(x...,dims=3),(2,1,3))
    haskey(fid, h) && delete_object(fid,h) 
    fid[h] = arr
end

# Insert vertex-level ddHodge results to .obs of anndata
function insertddh(adfile,ddh;keys=[:u,:div,:rot,:vgrass],prefix="ddh_")
    h5open(adfile,"r+") do fd
        for k in keys
            writeobs(fd,"prefix$(k)", getfield(ddh,k))
        end
    end
end
end # module

