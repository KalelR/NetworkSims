function diffusive_coupling_all!(u, t, xs, ys, coup_params)
    @inbounds begin
        @unpack ϵ, adjl, Icoupx, Icoupy = coup_params

        Icoupx = get_tmp(Icoupx, u)
        Icoupy = get_tmp(Icoupy, u)

        fill!(Icoupx, 0.0); fill!(Icoupy, 0.0); 
        for i in eachindex(xs)
            if isempty(adjl[i])  continue end
            for j ∈ adjl[i]
                Icoupx[i] += xs[j] - xs[i]
                Icoupy[i] += ys[j] - ys[i]
            end
        end
        @. Icoupx *= ϵ
        @. Icoupy *= ϵ
    end
end