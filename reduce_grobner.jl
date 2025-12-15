"""
Reduce a Grobner basis G to the unique minimal grobner basis for a given
ideal <G> for a particular monomial ordering.
"""
function reduce_gb(G::Vector{T})::Vector{T} where {T <: AbstractPolynomial}
	G = filter(p -> !iszero(p), G) # Remove 0 polynomials

	# Make all monic
	map!( p -> p / leading_coefficient(p), G, G)
	length(G) == 1 && return G # Short circuit edge case

	# Repeatedly reduce G until it is stable
	while true
		changed = false

        # gi -> gi ÷ (G / gi)
        i = 1
        while i <= length(G)
			G[i], G[end] = G[end], G[i]
			temp = polyrem(G[end], G[1:end-1])
			if temp != G[end]
				changed = true
				G[end] = temp
			end

			# Remove redundant polynomial
			if iszero(G[end])
				changed = true
				pop!(G)
				continue
			end
			G[end] = G[end] / leading_coefficient(G[end])
			removed = false
			for j in 1:length(G)-1
				if divides(leading_term(G[j]), leading_term(G[end]))
					changed = true
					removed = true
					pop!(G)
					break
				end
			end
			removed && continue
			
			# Swap back if we are keeping it
			G[i], G[end] = G[end], G[i]
			i += 1
		end

		
		!changed && break
	end

	sort!(G, by=p -> leading_term(p))
	return G
end

