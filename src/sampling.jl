
function getproposal(p::Prior, nparams)

  newparams = zeros(Float64, nparams)

  for i in 1:nparams
    newparams[i] = rand(Uniform(p.p[i,:]...))
  end

  return newparams
end
