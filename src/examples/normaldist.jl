srand(1)

targetdata = rand(Normal(2, 0.4), 100)

function getnormal(params, constants, targetdata)

  simdata = rand(Normal(params...), 100)

  ksdist(simdata, targetdata)

end

#@time particles, parameters = runabc(ABCRejection(getnormal, 2, 0.1, 100, 1.0, 1000000, PriorUniform([0 20; 0 2.0]), targetdata);
