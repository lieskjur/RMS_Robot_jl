module PontryaginMinEnergy

	struct EquationsOfMotion
		M;C;J
		n::Int
		w::Function
		f::Function
	end

	function EquationsOfMotion(M::Function,C::Function,J::AbstractArray,n::Int)
		EquationsOfMotion(
			(x)->M(x[1:n],x[n+1:2*n]),
			(x)->C(x[1:n],x[n+1:2*n]),
			J,
			n,
			(x,u)->M(x[1:n],x[n+1:2*n]) \ ( J*u - C(x[1:n],x[n+1:2*n]) ),
			(x,u)->[ x[n+1:2*n] ; M(x[1:n],x[n+1:2*n]) \ ( J*u - C(x[1:n],x[n+1:2*n]) ) ]
			)
	end
	function EquationsOfMotion(M::Function,C::Function,J::Function,n::Int)
		EquationsOfMotion(
			(x)->M(x[1:n],x[n+1:2*n]),
			(x)->C(x[1:n],x[n+1:2*n]),
			(x)->J(x[1:n],x[n+1:2*n]),
			n,
			(x,u)->M(x[1:n],x[n+1:2*n]) \ ( J(x[1:n],x[n+1:2*n])*u - C(x[1:n],x[n+1:2*n]) ),
			(x,u)->[ x[n+1:2*n] ; M(x[1:n],x[n+1:2*n]) \ ( J(x[1:n],x[n+1:2*n])*u - C(x[1:n],x[n+1:2*n]) ) ]
			)
	end

	module LinearizationMethod
		using LinearAlgebra

		κ(ν,t) = [	1/6*I(ν)*t^3	1/2*I(ν)*t^2	I(ν)*t	I(ν)
					1/2*I(ν)*t^2	I(ν)*t 			I(ν)	zeros(ν,ν) ]

		# Synteza trajektorie
		function opt_λ(ν,R,T,z_0,z_f)
			K = [	κ(ν,T[1])
					κ(ν,T[2])	]
			Ω = kron(I(2*ν),R)
			ζ = [z_0;z_f]
			λ = 2*K\Ω*ζ
			return λ
		end

		# Prubeh stavovych velicin
		z(ν,R,λ,t) = 0.5*kron(I(ν),R)\(κ(ν,t)*λ)

		# prubeh konjugovanych promennych
		p(ν,λ,t) = [ -λ[1:ν] ; λ[1:ν]*t + λ[ν+1:2*ν] ]
	end

	module ShootingMethod
		using ForwardDiff
		using BoundaryValueDiffEq
		using DifferentialEquations
		using LinearAlgebra
		using ..PontryaginMinEnergy: EquationsOfMotion
		using ..LinearizationMethod

		# casova derivace konjugovanych promennych
		function dpdt(EoM::EquationsOfMotion,x::AbstractArray,p::AbstractArray,u::AbstractArray)
			dwdx(x) = ForwardDiff.jacobian((x)->EoM.w(x,u),x)
			return -[ zeros(EoM.n,EoM.n) I(EoM.n) ; dwdx(x) ]'p
		end

		# optimalni akcni zasah
		dfdu(M::Function,J::AbstractArray,n,x) = [ zeros(n,n) ; M(x)\J ]
		dfdu(M::Function,J::Function,n,x) = [ zeros(n,n) ; M(x)\J(x) ]

		function u_opt(EoM::EquationsOfMotion,n::Int,R::AbstractArray,x::AbstractArray,p::AbstractArray)
			0.5 * R \ dfdu(EoM.M,EoM.J,n,x)' * p
		end

		# Vypocet optimalni trajektorie
		function OptimalTrajectory(
			EoM::EquationsOfMotion,
			tspan::Tuple{Float64,Float64},
			x_0::AbstractArray,x_f::AbstractArray,
			R::AbstractArray
			)

			@assert length(x_0) == length(x_f)
			n = EoM.n; r = n*2


			# optimalni trajektorie linearizovaneho systemu
			λ = LinearizationMethod.opt_λ(n,R,tspan.-tspan[1],x_0,x_f)
			y_λ(t) = [	LinearizationMethod.z(n,R,λ,t)
						LinearizationMethod.p(n,λ,t)	]

			# ODE system pro optimalizaci trajektorie
				function S!(dy,y,ρ,t)
					x = y[1:r]
					p = y[r+1:2*r]
					u = u_opt(EoM,n,R,x,p)
					dy[1:2*r+1] = [	EoM.f(x,u)
									dpdt(EoM,x,p,u)
									u'R*u		]
				end

			# funkce residui (podminka splneni okrajovych podm.)
			function bc!(residual,sol,ρ,t)
				residual[1:r] = sol(tspan[1])[1:r] - x_0
				residual[r+1:2*r] = sol(tspan[2])[1:r] - x_f
			end

			# Reseni okrajove ulohy (optimalni rizeni)
			bvp_opt = BVProblem(S!,bc!,append!(y_λ(tspan[1]),0),tspan)
			sol_opt = solve(bvp_opt, Shooting(Tsit5()));

			return sol_opt
		end
	end
end