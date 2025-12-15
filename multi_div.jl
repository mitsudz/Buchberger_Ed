"""
Computes a remainder after multivariate polynomial division of `f` by the divisors.

If needed, we could also keep track of the quotients (useful for the 
Extended Buchberger's Algorithm - analogous to Extended Euclidean Algorithm)
"""
function polyrem(f::T, divisors::Vector{T})::T where {T <: AbstractPolynomial}
    r = 0*f
    while f != 0
        lt_f = leading_term(f)
        found_divisor = false

        # Divide leading term by a divisor and update f
        for divisor in divisors
            lt_div = leading_term(divisor)
            if divides(lt_div, lt_f)
                found_divisor = true
                f -= divisor * div_multiple(lt_f, lt_div)
                break
            end
        end

        # Add to remainder if no divisor found
        if found_divisor == false
            r += lt_f
            f -= lt_f
        end
    end

    return r
end

