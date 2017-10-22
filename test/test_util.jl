#random numbers from the same distribution should have small ks dist
@test ApproxBayes.ksdist(rand(10^5), rand(10^5)) < 0.01
