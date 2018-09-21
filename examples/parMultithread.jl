# Test implementation of parallel algorithms using multithreading

using Base.Threads
include("parProblem.jl")

function simSerial(N,tol)
  dVec = Array{Float64}(undef,N)
  iCnt = 0
  for ii=1:N
    (d,sol) = simLV(p+0.05*randn(size(p)), 1, targetdata)
    if d < tol
      iCnt += 1;
    end
    if iCnt > 1000
      break
    end
  end
  return dVec
end

function simParallel(N,tol)
  dVec = Array{Float64}(undef,N)
  iCnt = Atomic{Int64}(1)
  @threads for ii=1:N
    d,sol = simLV(p+0.05*randn(size(p)), 1, targetdata)
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

@time simSerial(100000,10.0);
@time dV = simParallel(100000,10.0);
