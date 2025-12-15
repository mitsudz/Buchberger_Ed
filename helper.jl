###############################################################
# Helper functions for manipulating multivariate polynomials. # 
###############################################################

import Base: lcm
using MultivariatePolynomials

"""
Computes the least common multiple of the leading terms of two polynomials.
"""
function lcm(f::T, g::T)::T where {T <: AbstractPolynomial}
    lm_f = leading_monomial(f)
    lm_g = leading_monomial(g)
    lcm(lm_f, lm_g)
end

"""
Computes the S-polynomial of `f` and `g`
spoly(f, g) = (lcm/lt(f))*f - (lcm/lt(g))*g
"""
function spoly(f::T, g::T)::T where {T <: AbstractPolynomial}
    lcm_fg = lcm(f, g)

    lt_f = leading_term(f)
    lt_g = leading_term(g)

    return div_multiple(lcm_fg, lt_f) * f - div_multiple(lcm_fg, lt_g) * g
end

