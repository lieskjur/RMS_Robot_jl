module RegulaceNQR
	export NQR
	using MatrixEquations

	ϕ(f,x,i) = f([zeros(i);x[i+1:end]])

	function Rozklad(f,x)
		ξ = x; ξ[ξ.==0] .= floatmin(Float64)
		lx = length(x)
		F = [ ϕ(f,x,i) for i in 0:lx ]
		A = reduce(hcat,[ (F[i]-F[i+1])/ξ[i] for i in 1:lx ])
	end

	function NQR(R,Q,f,g,x)
		A = Rozklad(f,x)
		B = g(x)
		P,CLSEIG = arec(A,B,R,Q)
		K = R\(B'P)
	end
end