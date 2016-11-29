
function runabc(ABCsetup::ABCRejection, targetdata)

  #initalize array of particles
  particles = Array(ParticleRejection, ABCsetup.nparticles)

  #parameters = zeros(Float64, ABCsetup.nparticles, ABCsetup.nparams)


  i = 1 #set particle indicator to 1
  its = 0 #keep track of number of iterations

  while (i < (ABCsetup.nparticles + 1)) & (its < ABCsetup.maxiterations)

    its += 1

    #get new proposal parameters
    newparams = getproposal(ABCsetup.prior, ABCsetup.nparams)

    #simulate with new parameters
    dist = ABCsetup.sim_func(newparams, ABCsetup.constants, targetdata)

    #if simulated data is less than target tolerance accept particle
    if dist < ABCsetup.ϵ
      particles[i] = ParticleRejection(newparams)
      #parameters[i,:] = newparams
      i +=1
    end


  end

  return ABCrejectionresults(particles, its, ABCsetup)

end
