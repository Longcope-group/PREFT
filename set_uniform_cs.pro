; the initialization program also contains evaluation routine
;

pro field_at_points, tube
;  the field strength and its gradients at n points
common uniform_cs, b0

tube.b = replicate( b0, tube.n )
; tube.db = fltarr( 3, tube.n )

return
end

; -----------------------------

pro init_field_line, tube, len=len, dtheta=dtheta

if( n_elements( dtheta ) lt 1 ) then dtheta = 90.0;  in degrees
if( not keyword_set( len ) ) then len = 10.0

n = tube.n
tube.time = 0.0d0
; set_consts, tube

dl = len/float(n)

tube.x[0,*] = dl*cos( dtheta*0.5/!radeg )*( findgen(n)- 0.5*n )
tube.x[1,*] = 0.0
tube.x[2,*] = dl*sin( dtheta*0.5/!radeg )*( findgen(n)- 0.5*n )
tube.x[2,*] = abs( tube.x[2,*] )
; rb = 10*dl
; tube.x[2,*] = sqrt( tube.x[2,*]^2 + rb^2 )

mz = max( tube.x[2,*] )
tube.x[2,*] = mz - tube.x[2,*]

round_rx_bend, tube, n/2, nb=10

tube.v = 0*tube.x

calc_tv, tube
field_at_points, tube

return
end

; --------------------------------------------------------------

pro set_uniform_cs, b0_
common uniform_cs, b0

b0 = b0_

return
end