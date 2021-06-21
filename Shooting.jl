using Plots
using LinearAlgebra
include("PontryaginMinEnergy.jl"); using .PontryaginMinEnergy

## Pohybove rovnice robota
M(m,I,q,dq) = diagm([ I+m*q[2]^2 , m ])
C(m,I,q,dq) = [	2*m*q[2]*dq[2]*q[1]
				m*q[2]*dq[1]^2		]
J = I(2)

RobotEoM = PontryaginMinEnergy.EquationsOfMotion((q,dq)->M(1,1,q,dq),(q,dq)->C(1,1,q,dq),J,2)

## Uloha
# Okrajove podminky
tspan = (0.0,1.5);
x_0 = [0,0,0,0];
x_f = [1,1,0,0];
#=tspan = (0.5,3.);
x_0 = [0,0,1,0];
x_f = [3,1,0,0];=#

# Kriterium optimality
R = I(2)

## Vypocet energeticky optimalni trajektorie
sol_opt = PontryaginMinEnergy.ShootingMethod.OptimalTrajectory(RobotEoM,tspan,x_0,x_f,R)

## Prubeh stavovych velicin
fig1 = plot(legend=:outertopright,dpi=200,title="Energeticky optimalní trajektorie");
for i in 1:4 plot!(fig1,sol_opt,ls=:solid, vars=(0,i),lc=i,label="x"*string(i)) end
#savefig(fig1,"FIGs/Shooting_trajectory")
display(fig1)
