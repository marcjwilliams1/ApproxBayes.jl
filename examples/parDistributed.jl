# Test implementation of parallel algorithms using multithreading

using Distributed
using SharedArrays
addprocs(2)
@everywhere include("/home/rowan/Documents/github_repos/ApproxBayes/examples/parProblem.jl")

function simParallel(N,tol)
  dVec = SharedArray{Float64}(N)
  iCnt = SharedArray{Int64}(1)
  iCnt[1] = 0
  @distributed for ii=1:N
    (d,sol) = simLV(p+0.05*randn(size(p)), 1, targetdata)
    dVec[ii] = d
    # if d < tol
      iCnt[1] += 1
    # end
  end
  return dVec, iCnt
end

@time dV,ii = simParallel(10000,15.0);
