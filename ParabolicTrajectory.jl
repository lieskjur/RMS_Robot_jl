module ParabolicTrajectory
	#export p,opt_λ
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
		return λ,[ λ[ν*k+1:ν*(k+1)] for k = 0:2*ν-1 ]... # λ,α,β,γ,δ
	end

	# Prubeh stavovych velicin
	z(ν,R,λ,t) = 0.5*kron(I(ν),R)\(κ(ν,t)*λ)

	# prubeh konjugovanych promennych
	p(α,β,t) = [ -α ; α*t + β ]
end