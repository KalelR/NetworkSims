function construct_z(x, y, zfunc)
    z = zeros(length(x), length(y))
    # z = zeros(length(x), length(y)+1)
    for idxrow in 1:size(z, 1)
        for idxcol in 1:size(z, 2)
            z[:, idxcol] = zfunc(x, y)
        end 
    end 
    # z[:, end] = z[:, end-1] .+ (z[:, end-1] .- z[:, end-2])
    return z 
end

function Makie.surface(x, y, f::Function; kwargs...)
    zs = construct_z(x, y, f)
    surface(x, y, zs; kwargs...)
end

function Makie.surface!(ax, x, y, f::Function; kwargs...)
    zs = construct_z(x, y, f)
    surface!(ax, x, y, zs; kwargs...)
end