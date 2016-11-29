
function runabc(ABCsetup::ABCRejection, targetdata)

  #initalize array of particles
  particles = Array(ParticleRejection, ABCsetup.nparticles)

  i = 1 #set particle indicator to 1
  its = 0 #keep track of number of iterations
  distvec = zeros(Float64, ABCsetup.nparticles) #store distances in an array

  while (i < (ABCsetup.nparticles + 1)) & (its < ABCsetup.maxiterations)

    its += 1

    #get new proposal parameters
    newparams = getproposal(ABCsetup.prior, ABCsetup.nparams)

    #simulate with new parameters
    dist = ABCsetup.sim_func(newparams, ABCsetup.constants, targetdata)

    #if simulated data is less than target tolerance accept particle
    if dist < ABCsetup.ϵ
      particles[i] = ParticleRejection(newparams)
      distvec[i] = dist
      i +=1
    end


  end

  i > ABCsetup.nparticles || error("Only accepted $(i-1) particles with ϵ <  $(ABCsetup.ϵ). \nDecrease ϵ or increase maxiterations ")

  return ABCrejectionresults(particles, its, ABCsetup, distvec)

end
