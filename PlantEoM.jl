# Robot Creation
include("RobotModel.jl"); using .RobotModel
R = Robot(1,1)

# Gardening Tools
Plant(x,u) = R.ODE(x,u)
LinPlant(z,w) = [ z[3:4] ; w ]
Trowel(z,w) = R.M(z)*w+R.C(z)

# Eliptical trajectory
z_d(t) = [t,t,1,1]
w_d(t) = z_d(t)[3:4]

#=x_d(t) = [t,t,1,1]
w = xx_d(1)[3:4]=#

# System
using DifferentialEquations

function f(x,p,t)
	u = Trowel(z_d(t),w_d(t))
	dx = Plant(x,u)
end

prob = ODEProblem(f,z_d(0),[0,1])
solve(prob)

