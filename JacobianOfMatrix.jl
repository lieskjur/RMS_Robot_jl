using ForwardDiff
using Tullio
using BenchmarkTools
f(x) = [ 	x[1]^2 x[1]*x[2]                                                                                  
       		x[2]*x[1] x[2]^2	]

fake(u,x) = [	2*u[1]*x[1] + u[2]*x[2]		u[2]*x[1]
				u[1]*x[2]					u[1]*x[1]+2*u[2]*x[2]	]

J(x) = @benchmark ForwardDiff.jacobian(f,x)
x = [3,1]
SH = (2,2,2)
dJdx = reshape(J(x),SH)

u = [2,4]
dJudx = zeros(2,2)
@tullio dJudx[i,k] = dJdx[i,j,k]*u[i]


m = 3
J = 2
include("dwdx.jl")

dwdx(rand(4),rand(2))
