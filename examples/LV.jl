using DifferentialEquations
using Distributions
using Distances
using Gadfly
using ApproxBayes
using DataFrames

function getsolution(sol, times)
  #function to get species abundance from ODE solution
  if VERSION < v"0.6-"
    x1 = map(x -> x[1], sol(times))
    y1 = map(x -> x[2], sol(times))
  else
    x1 = map(x -> x[1], sol(times).u)
    y1 = map(x -> x[2], sol(times).u) #Need .u syntax as of 0.6
  end
  return x1, y1
end

# define Lotka Voltera ODE where a and b can be modified
f = @ode_def LV begin
  dx = a*x - b*x*y
  dy = b*x*y - y
end a=>0.8 b=>2.0

x0 = [1.0; 0.5]
tspan = (0.0, 15.0)
prob = ODEProblem(f, x0, tspan)
sol = solve(prob)

#generate target data by sampling 15 points and then adding Gaussian noise
times = 1.0:2.0:15.0
x, y = getsolution(sol, times)
x += rand(Normal(0.0, 0.3^2), length(times))
y += rand(Normal(0.0, 0.3^2), length(times))
targetdata = [x, y]

#simulations function for ABC. return distance (sum of squared distances) and solution
function simLV(params, constants, targetdata)
  a = params[1]
  b = params[2]
  x0 = [1.0; 0.5]
  tspan = (0.0, 15.0)
  times = 1.0:2.0:15.0
  h = LV(a = a, b = b)
  prob = ODEProblem(h, x0, tspan)
  sol = solve(prob)
  x1, y1 = getsolution(sol, times)
  d = sum((x1 .- targetdata[1]).^2 + (y1 .- targetdata[2]).^2)
  return d, sol
end

#define ABC setup type
setup = ABCSMC(simLV,
  2,
  0.1,
  Prior([Uniform(0.0, 5.0), Uniform(0.0, 5.0)]),
  maxiterations = 10^6,
  convergence = 0.01,
  nparticles = 1000
  )
  #run ABC SMC algorithm
@time ressmc = runabc(setup, targetdata, verbose = true, progress = true);

#show results
show(ressmc)

#plot posterior parameters
plotparameterposterior(ressmc)

#plots results using Gadfly
#Create data frame for target data
DFt = DataFrame(x = targetdata[1], y = targetdata[2], time = times)
d1 = stack(DFt, [:x, :y])

#get solution for
res = mean(ressmc.parameters, weights(ressmc.weights), 1)
x0 = [1.0; 0.5]
tspan = (0.0, 15.0)
h = LV(a = res[1], b = res[2])
prob = ODEProblem(h, x0, tspan)
sol = solve(prob)
times2 = 1.0:0.001:15.0
x2, y2 = getsolution(sol, times2)

DF1 = DataFrame(x = x2, y = y2, time = times2)
d2 = stack(DF1, [:x, :y])

l1 = layer(d1,x = :time, y = :value, color = :variable, Geom.point)
l2 = layer(d2,x = :time, y = :value, color = :variable, Geom.line)
Gadfly.plot(l1, l2)
