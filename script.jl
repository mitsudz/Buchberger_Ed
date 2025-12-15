##########################################
# A basic script to test this repository #
##########################################

include("helper.jl")
include("multi_div.jl")
include("reduce_grobner.jl")
include("buchberger.jl")

using DynamicPolynomials

@polyvar x y z

f = (3//1)x^2 + 2x*y
g = (2//1)x - y
h = (1//1)x^3*y^2 + z^2 + 5z^5*x^3
@show f
@show g
@show h 

S = spoly(f, g)
println("S-polynomial of f and g: ", S, "\n")

println("Dividing f by y^2 and 2x - y")
println(polyrem(f, [(1//1)y^2, (2//1)x-y]), "\n")

grobner_basis = buchberger([f, g, h])
println("Grobner basis of f, g, h:")
println(grobner_basis, "\n")

red_grob = reduce_gb(grobner_basis)
println("Reduced grobner basis of f, g, h:")
println(red_grob)
;
