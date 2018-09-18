# Test implementation of parallel algorithms using multithreading

using DifferentialEquations
using Base.Threads
using Distributions
using Distances
using Plots
gr()

function getsolution(sol, times)
    x1 = map(x -> sol(x)[1], times)
    y1 = map(x -> sol(x)[2], times)
  return x1, y1
end

# define Lotka Voltera ODE where a and b can be modified
f = @ode_def LV begin
  dx = a*x - b*x*y
  dy = b*x*y - y
end a b

x0 = [1.0; 0.5]
tspan = (0.0, 15.0)
p = [2.0, 0.8]
prob = ODEProblem(f, x0, tspan, p)
sol = solve(prob)
Plots.plot(sol)

#generate target data by sampling 15 points and then adding Gaussian noise
times = 1.0:2.0:15.0
x, y = getsolution(sol, times)
x .+= rand(Normal(0.0, 1), length(times))
y .+= rand(Normal(0.0, 1), length(times))
targetdata = [x, y]


#simulations function for ABC. return distance (sum of squared distances) and solution
function simLV(params, constants, targetdata)
  a = params[1]
  b = params[2]
  x0 = [1.0; 0.5]
  tspan = (0.0, 15.0)
  times = 1.0:2.0:15.0
  prob = ODEProblem(f, x0, tspan, [a, b])
  sol = solve(prob)
  x1, y1 = getsolution(sol, times)
  d = sum((x1 .- targetdata[1]).^2 + (y1 .- targetdata[2]).^2)
  return d, sol
end

function simSerial(N,tol)
  iCnt = 0
  for ii=1:N
    (d,sol) = simLV(p+0.05*randn(size(p)), 1, targetdata)
    if d < tol
      iCnt += 1;
    end
  end
  return iCnt
end

function simParallel(N,tol)
  dVec = Array{Float64}(undef,N)
  iCnt = Atomic{Int64}(0)
  @threads for ii=1:N
    (d,sol) = simLV(p+0.05*randn(size(p)), 1, targetdata)
    dVec[ii] = d
    if d < tol
      atomic_add!(iCnt, 1)
    end
    if iCnt[] > 1000
      break
    end
  end
  return dVec
end
