using Pkg
Pkg.activate(".")

using Catalyst, DifferentialEquations, Plots 

# Initial Catalyst model 
ErkModel = @reaction_network begin
    (k₁,k₋₁),  M + MKK   <--> M_MKK
    k₂,        M_MKK    --> Mp + MKK
    (k₃,k₋₃),  Mp + MKK  <--> Mp_MKK
    k₄,        Mp_MKK   --> Mpp + MKK 
    (h₁,h₋₁),  Mpp + MKP <--> Mpp_MKP
    h₂,        Mpp_MKP  --> Mp + MKP 
    (h₃,h₋₃),  Mp + MKP  <--> Mp_MKP
    h₄,        Mp_MKP   --> M + MKP
end k₁ k₂ k₃ k₄ h₁ h₂ h₃ h₄ k₋₁ k₋₃ h₋₁ h₋₃

# Create Graph 
Graph(ErkModel)

# Default parameter values
p = [:k₁ => 0.00166, :k₂ => 0.0001, :k₃ => 0.1,:k₄ => 0.00166, 
			:h₁ => 0.02, :h₂ => 0.001,  :h₃ => 0.02,  :h₄ => 0.02, 			
			:k₋₁ => 0.0001, :k₋₃ => 0.1,  :h₋₁ => 0.02, :h₋₃ => 0.02]

# Inital conditions for states            
u0 = [:M => 201, :MKK => 100, :M_MKK => 10, :Mp_MKK => 20, :MKP => 2, 
		:Mpp => 23, :Mp => 10, :Mpp_MKP => 10, :Mp_MKP => 10]

# Simulation time
tspan = (0.0, 4.0)

# ODEProblem convertion 
ode_prob = ODEProblem(ErkModel, u0, tspan, p)
sol  = solve(ode_prob)
p1 = plot(sol, legend = :outerright, grid = "off")

# SDEProblem convertion 
sde_prob = SDEProblem(ErkModel, u0, tspan, p)
sol2  = solve(sde_prob,  LambaEM(), tstops = range(0.0, step = 4e-3, length = 1001))
p2 = plot(sol2, legend = :outerright, grid = "off")

# JumpProblem
d_prob = DiscreteProblem(ErkModel, u0, tspan, p)
j_prob = JumpProblem(ErkModel, d_prob, Direct())
sol3 = solve(j_prob, SSAStepper())
p3 = plot(sol3, legend = :outerright, grid = "off")