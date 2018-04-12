using DifferentialEquations
using Distributions
using Distances
using Gadfly
using ApproxBayes
using DataFrames
using Plots
using DiffEqBayes

f1 = @ode_def_nohes LotkaVolterraTest1 begin
 dx = a*x - x*y
 dy = -3*y + x*y
end a

p = [1.5]
u0 = [1.0,1.0]
tspan = (0.0,10.0)
prob1 = ODEProblem(f1,u0,tspan,p)

σ = 0.01                         # noise, fixed for now
t = collect(linspace(1,10,10))   # observation times
sol = solve(prob1,Tsit5())

randomized = VectorOfArray([(sol(t[i]) + σ * randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)

bayesian_result_stan = stan_inference(prob1,t,data,priors;num_samples=300,
                                num_warmup=500,likelihood=Normal,
                                vars =(StanODEData(),InverseGamma(3,2)))

bayesian_result_turing = turing_inference(prob1,Tsit5(),t,data,[Normal(1.5, 1)];num_samples=500)

bayesian_result_hmc = dynamichmc_inference(prob1, data, [Normal(1.5, 1)], t, [bridge(ℝ, ℝ⁺, )])

function createabcfunction(prob, t, distancefunction, slver)
    function simfunc(params, constants, targetdata)
        prob.p = params
        sol = solve(prob, slver)
        randomized = VectorOfArray([sol(t[i]) for i in 1:length(t)])
        data = convert(Array, randomized)
        distancefunction(targetdata, data), data
    end
end

function abc_inference(prob::DEProblem, data, priors, t; ϵ=0.001, distancefunction = euclidean, slvr = Tsit5(), ABCalgorithm = ABCSMC, verbose = false, progress = false, kwargs...)

    setup = ABCalgorithm(createABCfunction(prob, t, distancefunction, slvr),
      length(priors),
      ϵ,
      Prior(priors);
      kwargs...
      )

    abcresult = runabc(setup, data, progress = progress)

    return abcresult
end

result_abc = abc_inference(prob1, data, [Normal(5.5, 1)], t)




f1 = @ode_def_nohes LotkaVolterraTest4 begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
u0 = [1.0,1.0]
tspan = (0.0,10.0)
p = [1.5,1.0,3.0,1.0]
prob1 = ODEProblem(f1,u0,tspan,p)
sol = solve(prob1,Tsit5())
t = collect(linspace(1,10,10))
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)
priors = [Truncated(Normal(1.5,0.01),0,2),Truncated(Normal(1.0,0.01),0,1.5),
          Truncated(Normal(3.0,0.01),0,4),Truncated(Normal(1.0,0.01),0,2)]

priors = [Uniform(0.0, 5.0), Uniform(0.0, 5.0), Uniform(0.0, 5.0), Uniform(0.0, 5.0)]

result_abc = abc_inference(prob1, data, priors, t, verbose = true, progress = true)

bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=500,epsilon = 0.001)
