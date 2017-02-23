Function get_jem_lib, verbose=verbose
  common get_jem_lib_blk, shlib

  build_lib = n_elements(shlib) eq 0
  if (not build_lib) then build_lib = not file_test(shlib, /read)
  if (build_lib) then begin
     call_ex_dir = getenv('JEM_IDL')+'call_external/'
     source = [ 'drizzle1d' ]
     export_rtns = [ 'drizzle1d_natural', 'drizzle1d' ]
     make_dll, source, 'jem_lib', export_rtns, $
               input_dir=call_ex_dir, dll_path=shlib, $
               verbose=verbose, show_all_output=verbose
  endif
  return, shlib
end
