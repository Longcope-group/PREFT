pro straighten_tube, tube, cf_smooth=cf_smooth

; if cf_smooth > 0 then smooth the velocity over the section of tube 
;  a distance cf_smooth from conduction fronts.

if( not keyword_set( cf_smooth ) ) then cf_smooth=0.0

calc_tv, tube
calc_len, tube

; find midpoint
nh = tube.n/2
x0 = tube.x[*,nh]
l0 = tube.l[nh]

; now create a stright tube
tube.x[0,*] = x0[0] + ( tube.l - l0 )
tube.x[1,*] = x0[1]
tube.x[2,*] = x0[2]

tv = tube.tv_e
vpar = total( tube.v*tv, 1 )

;  implement smoothing
if( cf_smooth gt 1.0d-3 ) then begin
  print, 'smoothing velocity over ', cf_smooth, ' Mm from conduction fronts'
  n = n_elements( tube.l )
  nh = n/2
  l0 = intersect( tube.l[0:nh], tube.t[0:nh], 0.1 );  use 100,000 for CF location
  dl = tube.l - l0[0]
  ii = where( tube.l lt l0[0]+cf_smooth, nii )
  if( nii gt 1 ) then begin
    vpar[ii] = 0.0; smooth( vpar[ii], 10, /edge_trunc )
  endif
  l0 = intersect( tube.l[nh:*], tube.t[nh:*], 0.1 );  use 100,000 for CF location
  dl = tube.l - l0[0] 
  ii = where( tube.l gt l0[0]-cf_smooth, nii )
  if( nii gt 1 ) then begin
    vpar[ii] = 0.0; smooth( vpar[ii], 10, /edge_trunc )
  endif
endif

; keep only parallel velocity

tube.v[0,*] = vpar
tube.v[1,*] = 0.0
tube.v[2,*] = 0.0

calc_tv, tube
calc_len, tube

return
end