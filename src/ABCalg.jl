
function runabc(ABCsetup::ABCRejection, pr::Prior)

  #initalize array of particles
  particles = Array(Particlerejection, ABCsetup.nparticles)

  #set particle indicator to 1
  i = 1

  while i < ABCsetup.nparticles

    #get new proposal parameters
    newparams = getproposal(p::Prior, nparams)

    #simulate with new parameters
    dist = ABCsetup.sim_func(newparams, ABCsetup.constants)

    #if simulated data is less than target tolerance accept particle
    if dist < ABCsetup.ϵ
      particles[i] = ParticleRejection(newparams)
    end


  end

  return particles

end
