using LinearAlgebra
using BoundaryValueDiffEq
using DifferentialEquations
using ForwardDiff
using Plots
include("RobotDyn.jl"); using .RobotDyn
include("ParabolicTrajectory.jl"); using .ParabolicTrajectory
pt = ParabolicTrajectory
include("RegulaceNQR.jl"); using .RegulaceNQR

## Nelinearni system
robot = Robot(1,1) # (m,I) => (:M,:C,:J)
n = 2; r = 2*n

## Uloha
# Okrajove podminky
#=tspan = (0.0,1.5);
x_0 = [0,0,0,0];
x_f = [1,1,0,0];=#
tspan = (0.5,3.);
x_0 = [0,0,1,0];
x_f = [3,1,0,0];

# Kriterium opt. trajektorie
Q_opt = zeros(r,r);
R_opt = I(n);

## Popis systemu
# idealni system
ref(x,u) = [ x[3:4] ; robot.M(x) \ ( u - robot.C(x) ) ] # bez treni

# pozorovatel : dx = f(x) + g(x)*u
f(x) = [ x[3:4] ; -robot.M(x) \ robot.C(x) ]
g(x) = [ zeros(2,2) ; inv( robot.M(x) ) ]

# nemodelovane tlumeni 
B(x) = I(n)*x[n+1:r] # M*dx + C(x) + B(x) = u 
b(x) = [ zeros(n) ; - robot.M(x) \ B(x) ] # dx = f(x)+b(x)+g(x)*u 

# potryaginuv princip maxima
u_opt(x,p) = 0.5 * R_opt \ [ zeros(n,n) ; inv(robot.M(x)) ]' * p
w(x,u) = robot.M(x) \ ( u - robot.C(x) )
function dpdt(x,p,u)
	dwdx(x) = ForwardDiff.jacobian((x)->w(x,u),x)
	2*Q_opt*x - [ zeros(n,n) I(n) ; dwdx(x) ]'p
end

## Urceni optimalni trajektorie
# ODE system pro optimalizaci trajektorie
function S!(dy,y,ρ,t)
	x = y[1:r]
	p = y[r+1:2*r]
	dy[1:2*r] = [	ref(x,u_opt(x,p))
					dpdt(x,p,u_opt(x,p))	]
end

# optimalni trajektorie linearizovaneho systemu
λ,α,β = pt.opt_λ(n,R_opt,tspan.-tspan[1],x_0,x_f)
y_λ(t) = [	pt.z(n,R_opt,λ,t)
			pt.p(α,β,t)		]

# Nastrel dle parabolickeho rizeni lin sys. (pouze ilustracni)
prob_λ = ODEProblem(S!,y_λ(0),tspan)
sol_λ = solve(prob_λ);

fig1 = plot(legend=:outertopright, vars=(0,1:n),dpi=200,title="Nástřel vs. Pontryagin");
for i in 1:4 plot!(fig1,sol_λ,ls=:dot, vars=(0,i),lc=i,label="x"*string(i)*" - 0.") end

# funkce residui (podminka splneni okrajovych podm.)
function bc!(residual,sol,ρ,t)
	residual[1:r] = sol(tspan[1])[1:r] - x_0
	residual[r+1:2*r] = sol(tspan[2])[1:r] - x_f
end

# Reseni okrajove ulohy (optimalni rizeni)
bvp_opt = BVProblem(S!,bc!,y_λ(tspan[1]),tspan)
sol_opt = solve(bvp_opt, Shooting(Tsit5()));

for i in 1:4 plot!(fig1,sol_opt,ls=:dash, vars=(0,i),lc=i,label="x"*string(i)*" - opt") end
#savefig(fig1,"FIGs/Shooting_NQR_1_LinVsPontr")

## Simulace systemu
# ODE system
function S(R,Q,y,ρ,t)
	ξ = y[1:r]; # stavove veliciny rostlinky
	x = y[r+1:2*r]; # stavove veliciny pozorovatele
	p = y[2*r+1:3*r]; # konjugovane promenne
	# Idealni model
	μ_opt = u_opt(x,p)
	dx = ref(x,μ_opt)
	dp = dpdt(x,p,μ_opt)
	# Regulator
	K = NQR(R,Q,f,g,ξ)
	Δξ = ξ-x
	u = μ_opt - K*Δξ
	# Rostlinka
	dξ = f(ξ)+b(ξ)+g(ξ)*u
	#
	return [dξ;dx;dp]
end

# Pocatecni podminky simulace
p_0 = sol_opt.u[1][5:8]
y_0 = [x_0;x_0;p_0]

# NQR regulace
j = 0 #>>
k = 2
R_NQR = I(n)*10.0^j
Q_NQR = diagm([1.5,0.7,0.7,0.1])*10.0^k
#Q_NQR = diagm([1.5,1.5,0.5,0.5])*10.0^k # pro prikladovou trajektorii

prob_reg = ODEProblem((y,ρ,t)->S(R_NQR,Q_NQR,y,ρ,t),y_0,tspan)
sol_reg = solve(prob_reg);

fig2 = plot(legend=:outertopright,dpi=200,
	title="Teoretická trajektorie vs. Tlumení + NQR:\nRᵢᵢ = 1e"*string(j)*", Qᵢᵢ ≈ 1e"*string(k)
);
for i in 1:4 plot!(fig2,sol_reg,ls=:dash, vars=(0,r+i),lc=i,label="x"*string(i)*" - opt") end
for i in 1:4 plot!(fig2,sol_reg,ls=:solid, vars=(0,i),lc=i,label="x"*string(i)*" - reg") end
#savefig(fig2,"FIGs/Shooting_NQR_2_OptVsNQR_1e"*string(j)*"_1e"*string(k))
display(fig2) #<