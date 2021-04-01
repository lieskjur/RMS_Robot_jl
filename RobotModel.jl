module RobotModel
	export Robot
	using LinearAlgebra

	struct Robot
		m; J;
		M; C;
		ODE;
		Robot(m,J) = new( m,J,
			(x)->M(m,J,x),
			(x)->C(m,J,x),
			(x,u)->ODE(m,J,x,u)	)
	end

	M(m,J,x) = diagm([ J+m*x[1]^2 , m ])
	C(m,J,x) = [	2*m*x[2]*x[4]*x[1]
					m*x[2]*x[3]^2		]
	ODE(m,J,x,u) = [ x[3:4] ; M(m,J,x)\(u-C(m,J,x)) ]
end