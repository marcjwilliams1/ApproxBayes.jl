using DifferentialEquations
using Distributions
srand(1)
println("Test ABC SMC with ODE model")

f = @ode_def LV begin
  dx = a*x - b*x*y
  dy = b*x*y - y
end a=>1.0 b=>1.0

inputa = 0.8
inputb = 2.2
g = LV(a = inputa, b = inputb)
x0 = [1.0; 0.5]
tspan = (0.0, 15.0)
prob = ODEProblem(g, x0, tspan)
sol = solve(prob)

#generate target data
times = 1.0:1.0:15.0
x = map(x -> x[1], sol(times).u) .+ rand(Normal(0.0, 0.25^2), length(times))
y = map(x -> x[2], sol(times).u) .+ rand(Normal(0.0, 0.25^2), length(times))
targetdata = [x, y]

function simLV(params, constants, targetdata)
  a = params[1]
  b = params[2]
  x0 = [1.0; 0.5]
  tspan = (0.0, 15.0)
  times = 1.0:1.0:15.0
  h = LV(a = a, b = b)
  prob = ODEProblem(h, x0, tspan)
  sol = solve(prob)
  x1 = map(x -> x[1], sol(times).u)
  y1 = map(x -> x[2], sol(times).u)
  d = sum((x1 .- targetdata[1]).^2 + (y1 .- targetdata[2]).^2)
  return d, sol
end

setup = ABCSMC(simLV,
  2,
  0.1,
  Prior([Uniform(0.0, 5.0), Uniform(0.0, 5.0)]),
  maxiterations = 10^6,
  convergence = 0.05,
  nparticles = 500
  )
@time ressmc = runabc(setup, targetdata, verbose = false, progress = false);

println("\t Check parameters are inferred correctly for Lotka Volterra model")
@test isapprox(mean(ressmc.parameters, weights(ressmc.weights), 1)[1], inputa, rtol = 0.05)
@test isapprox(mean(ressmc.parameters, weights(ressmc.weights), 1)[2], inputb, rtol = 0.05)
