using DifferentialEquations
using Distributions
using Distances
using Plots
using ApproxBayes

f = @ode_def LV begin
  dx = a*x - b*x*y
  dy = b*x*y - y
end a=>1.0 b=>1.0

x0 = [1.0; 0.5]
tspan = (0.0, 15.0)
prob = ODEProblem(f, x0, tspan)
sol = solve(prob)
plot(sol)

#generate target data
times = 1.0:1.0:15.0
x = map(x -> x[1], sol(times)) .+ rand(Normal(0.0, 0.5^2), length(times))
y = map(x -> x[1], sol(times)) .+ rand(Normal(0.0, 0.5^2), length(times))
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
  #println(params)

  x1 = map(x -> x[1], sol(times))
  y1 = map(x -> x[1], sol(times))

  d = sum((x1 .- targetdata[1]).^2 + (y1 .- targetdata[2]).^2)
  return d, sol
end


setup = ABCSMC(simLV,
  2,
  0.1,
  Prior([Uniform(0.0, 5.0), Uniform(0.0, 5.0)]),
  maxiterations = 10^6,
  convergence = 0.01,
  nparticles = 1000
  )
@time ressmc = runabc(setup, targetdata, verbose = true, progress = true);
show(ressmc)
