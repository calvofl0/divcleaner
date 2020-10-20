function div, Bb1, Bb2, Bb3, cell_sizes=h
 dim1 = size(Bb1, /dimensions)
 dim2 = size(Bb2, /dimensions)
 dim3 = size(Bb3, /dimensions)
 dim  = lonarr(3)
 dim(0) = dim1(0)-1
 dim(1) = dim1(1)
 dim(2) = dim1(2)
 D = dblarr(dim)

 if n_elements(h) eq 0 then begin
  h = ptr_new()
  no_h=1
 endif else begin
  h=double(h)
  no_h=0
 endelse

 s = call_external('divcleaner_idl.so', 'divergence', D, $
                   double(Bb1), double(Bb2), double(Bb3), $
                   long(dim), h, $
                   /unload, value=[0,0,0,0,0,no_h])
 return, D
end

function fixTopBottomBoundary, Bb1, Bb2, Bb3, cell_sizes=h
 dim1 = size(Bb1, /dimensions)
 dim2 = size(Bb2, /dimensions)
 dim3 = size(Bb3, /dimensions)
 dim  = lonarr(3)
 dim(0) = dim1(0)-1
 dim(1) = dim1(1)
 dim(2) = dim1(2)
 Bb3_out = dblarr(dim3)

 if n_elements(h) eq 0 then begin
  h = ptr_new()
  no_h=1
 endif else begin
  h = double(h)
  no_h=0
 endelse

 s = call_external('divcleaner_idl.so', 'fixTopBottomBoundary', Bb3_out, $
                   double(Bb1), double(Bb2), double(Bb3), $
                   long(dim), h, $
                   /unload, value=[0,0,0,0,0,no_h])
 return, Bb3_out
end

pro phiToBbc, Bb1, Bb2, Bb3, phi, cell_sizes=h
 dim = size(phi, /dimensions)
 Bb1 = dblarr([dim[0]+1, dim[1], dim[2]])
 Bb2 = dblarr([dim[0], dim[1]+1, dim[2]])
 Bb3 = dblarr([dim[0], dim[1], dim[2]+1])

 if n_elements(h) eq 0 then begin
  h = ptr_new()
  no_h=1
 endif else begin
  h = double(h)
  no_h=0
 endelse

 s = call_external('divcleaner_idl.so', 'phiToBbc', Bb1, Bb2, Bb3, $
                   double(phi), long(dim), h, $
                   /unload, value=[0,0,0,0,0,no_h])
end

function poisson, x0, bc, thr, maxiter, b=b, boundary=d, cell_sizes=h
 dim = size(x0, /dimensions)
 x_out = dblarr(dim)

 if n_elements(b) eq 0 then begin
  b = ptr_new()
  no_b=1
 endif else begin
  b = double(b)
  no_b=0
 endelse

 if n_elements(d) eq 0 then begin
  d = ptr_new()
  no_d=1
 endif else begin
  d = double(d)
  no_d=0
 endelse

 if n_elements(h) eq 0 then begin
  h = ptr_new()
  no_h=1
 endif else begin
  h = double(h)
  no_h=0
 endelse

 s = call_external('divcleaner_idl.so', 'poisson', x_out, $
                   double(x0), long(dim), long(bc), $
                   double(thr), long(maxiter), $
                   b, d, h, /unload, value=[0,0,0,0,0,0,no_b,no_d,no_h])
 return, x_out
end

pro killdiv, Bb1_out, Bb2_out, Bb3_out, Bb1, Bb2, Bb3, $
             cell_sizes=h, boundary_conditions=bc, boundary=d, $
             initial_guess=x0, threshold=thr, maxiter=maxiter
 dim1 = size(Bb1, /dimensions)
 dim2 = size(Bb2, /dimensions)
 dim3 = size(Bb3, /dimensions)
 dim  = lonarr(3)
 dim(0) = dim1(0)-1
 dim(1) = dim1(1)
 dim(2) = dim1(2)
 Bb1_out = dblarr(dim1)
 Bb2_out = dblarr(dim2)
 Bb3_out = dblarr(dim3)

 if n_elements(h) eq 0 then begin
  h = ptr_new()
  no_h=1
 endif else begin
  h = double(h)
  no_h=0
 endelse

 if n_elements(bc) eq 0 then begin
  bc = ptr_new()
  no_bc=1
 endif else begin
  bc = long(bc)
  no_bc=0
 endelse

 if n_elements(d) eq 0 then begin
  d = ptr_new()
  no_d=1
 endif else begin
  d = double(d)
  no_d=0
 endelse

 if n_elements(x0) eq 0 then begin
  x0 = ptr_new()
  no_x0=1
 endif else begin
  x0 = double(x0)
  no_x0=0
 endelse

 if n_elements(thr) eq 0 then begin
  thr = ptr_new()
  no_thr=1
 endif else begin
  thr = double(thr)
  no_thr=0
 endelse

 if n_elements(maxiter) eq 0 then begin
  maxiter = ptr_new()
  no_maxiter=1
 endif else begin
  maxiter = long(maxiter)
  no_maxiter=0
 endelse

 s = call_external('divcleaner_idl.so', 'killdiv', Bb1_out, Bb2_out, Bb3_out, $
                   double(Bb1), double(Bb2), double(Bb3), $
                   long(dim), h, bc, d, x0, thr, maxiter, $
                   /unload, value=[0,0,0,0,0,0,0,no_h,no_bc,no_d,no_x0,no_thr,no_maxiter])
end
