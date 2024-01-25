; beginning with a straight horizontal tube, introduce a bend 
; about the mid-point.

pro bend_tube, tube, dth, zero=zero, bend_frac=bend_frac

nmp = tube.n/2
if( keyword_set( bend_frac ) ) then begin
  calc_len, tube
  lbend = bend_frac*max( tube.l )
  mm = min( abs( tube.l - lbend ), nmp )
endif
print, nmp

xmp = tube.x[*,nmp]
xo = tube.x[0,*] - xmp[0]

if( keyword_set( zero ) ) then xmp = [ 0.0, 0.0, 0.0 ]

tube.x[0,0:nmp] = xmp[0] + xo[0:nmp]*cos(0.5*dth/!radeg)
tube.x[2,0:nmp] = xmp[2] + xo[0:nmp]*sin(0.5*dth/!radeg)
tube.x[0,nmp:*] = xmp[0] + xo[nmp:*]*cos(0.5*dth/!radeg)
tube.x[2,nmp:*] = xmp[2] - xo[nmp:*]*sin(0.5*dth/!radeg)

round_rx_bend, tube, nmp
calc_tv, tube
set_rho, tube;  maintain same profile even after lengths have changed

return
end