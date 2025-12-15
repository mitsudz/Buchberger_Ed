" Computes a (not necessarily reduced) Grobner basis of I = <g | g in `G`>."
function buchberger(G::Vector{T})::Vector{T} where {T <: AbstractPolynomial}
	G = deepcopy(G)
    while true 
        n = length(G)
        
        # Apply Buchberger criterion to extend G from G'
        for i in 1:n-1
            for j in 2:n
                S = spoly(G[i], G[j])
                r = polyrem(S, G)
                if r != 0
                    # Slight modification to original algorithm - we divide by the updated G in the same iteration
                    push!(G, r)
                end
            end
        end
        
        # Criterion satisfied if no extensions to G
        if length(G) == n
            break
        end
    end

	return G
end

#=
We can optimise this initially by including conditions to ignore S-polynomials
from the end of chapter 2, and we can also auto-reduce the intermediate basis
as we go to avoid too many extraneous computations.
=#

