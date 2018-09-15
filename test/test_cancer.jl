function probabilityM1(s₁::BigInt, t₁::BigFloat, n::BigInt)
    pM1 = factorial(s₁) ./ (exp(t₁) * (n + 1) ^(s₁ + 1))
    pM2 = (factorial(n) * factorial(s₁)) ./ factorial(n + s₁ + 1)
    p = pM1 ./(pM1 + pM2)
    p = Float64(p)
    return p
end

s1calc(x) = sum(x)
t1calc(x::Array{BigInt, 1}) = Float64(sum(log.(map(a -> factorial(a), x))))
t1calc(x::Array{Int64, 1}) = sum(log.(map(a -> factorial(a), x)))

function poissonsimulator(params, cst, targetdata)

    λ = params[1]
    simdata = rand(Poisson(λ), 100)
    s₁ = s1calc(simdata)
    if sum(simdata.>19) > 0
        t₁ = t1calc(map(x -> BigInt(x), simdata))
    else
        t₁ = t1calc(simdata)
    end
    euclidean(targetdata, [s₁, t₁]), 1, rand([true, false])
    #euclidean(targetdata[2], [s₁, t₁][2]), 1
end

function geometricsimulator(params, cst, targetdata)

    μ = params[1]
    simdata = rand(Geometric(μ), 100)
    s₁ = s1calc(simdata)
    if sum(simdata.>19) > 0
        t₁ = t1calc(map(x -> BigInt(x), simdata))
    else
        t₁ = t1calc(simdata)
    end
    euclidean(targetdata, [s₁, t₁]), 1, rand([true, false])
    #euclidean(targetdata[2], [s₁, t₁][2]), 1
end

function generatedata()
    λ = 0.5
    n = 100
    data = rand(Poisson(λ), n)
    s₁ = s1calc(data)
    t₁ = t1calc(data)
    pM1 = probabilityM1(BigInt(s₁), BigFloat(t₁), BigInt(n))

  return [s₁, t₁], pM1
end

################################################################
println("Test Bayes factors are calculated correctly (within 5% of true value)")
Random.seed!(1234)

cst = [[i] for i in 1:2]

ABCsetupsmc = ABCSMCModel([poissonsimulator, geometricsimulator], [1, 1], 1.0,
[Prior([Exponential(1)]), Prior([Uniform(0.0, 1.0)])];
nparticles = 500, maxiterations = 10^7)

ABCsetuprej = ABCRejectionModel([poissonsimulator, geometricsimulator], [1, 1], 1.0,
[Prior([Exponential(1)]), Prior([Uniform(0.0, 1.0)])];
nparticles = 500, maxiterations = 10^7)

td, pM1 = generatedata()
println("\t Checking ABC SMC model selection with algorithm that returns true or false for each simulations")
@time abcressmc = ApproxBayes.runabcCancer(ABCsetupsmc, td);
@test isapprox(pM1, abcressmc.modelprob[1], rtol = 0.05)






function getnormal2(params, constants, targetdata)

  simdata = rand(Normal(params...), 1000)
  ApproxBayes.ksdist(simdata, targetdata), 1, rand([true, false])
end

println("Test parameters of normal distribution are inferred correctly (mean within 5% of true value)")

Random.seed!(1)
p1 = 2.0
p2 = 0.4
targetdata = rand(Normal(p1, p2), 1000)

#setup ABC alogrithm specifications for Rejection algorithm
setup = ABCRejection(getnormal2,
  2,
  0.1,
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]);
  maxiterations = 10^6,
  )

#now run SMC algorithm
setup = ABCSMC(getnormal2,
  2,
  0.05,
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]),
  )
ressmc = ApproxBayes.runabcCancer(setup, targetdata, verbose = false);

println("\t Check ABC SMC with true/false returned from simulations correctly infers parameters")
# test that mean value of posterior is within 10% of true value
@test isapprox(mean(ressmc.parameters, weights(ressmc.weights), dims = 1)[1], p1, rtol = 0.05)
@test isapprox(mean(ressmc.parameters, weights(ressmc.weights), dims = 1)[2], p2, rtol = 0.05)
