;; a is pdf of some random variable
;; p is the confidence interval probability fraction
;; function finds the smallest interval containing fraction p of the values of the pdf
Function ConfidenceInterval, a, p

na = n_elements(a)
j = floor(p*na)

sorted = a[sort(a)]
usemin = min( (sorted[j:*] - sorted[0:na-1-j]), minhere)

return, sorted[[minhere, minhere+j]]
end
