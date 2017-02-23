;; simple function to make sure trailing / is present in directory strings.

Function directoryify, dir
  outdir = dir
  if strmid( dir, strlen(dir)-1, 1 ) ne '/' then $
     outdir +='/'
  return, outdir
end
