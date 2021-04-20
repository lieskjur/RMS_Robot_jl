module RobotDyn
	export Robot
	using LinearAlgebra

	struct Robot
		m; I;
		M; C;
		Robot(m,I) = new(
			m,I,
			(x)->M(m,I,x),
			(x)->C(m,I,x),
			)
	end

	M(m,I,x) = diagm([ I+m*x[2]^2 , m ])
	C(m,I,x) = [	2*m*x[2]*x[4]*x[1]
					m*x[2]*x[3]^2		]
end