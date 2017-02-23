Pro testjem_avsigclip
  pattern = dist(100,200)
  pattern = rebin(pattern, 100, 200, 10)
  seed=10
  err = randomn(seed,100,200,10)*sqrt(pattern)
  array = err+pattern
  weight = 1./pattern
  pull = 1.0
  maxiter = 10
  ave = jem_avsigclip(array, weight, pull, pull, maxiter)
  help, ave
end
