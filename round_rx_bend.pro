pro round_rx_bend_old, tube, ib, nb=nb
;  circular (or helical) bend between points ib-nb and ib+nb

if( n_elements( nb ) lt 1 ) then nb=4;  9 points in circle

calc_tv, tube

x0 = reform( tube.x[*,ib-nb], 3 )
x1 = reform( tube.x[*,ib+nb], 3 )
t0 = reform( tube.x[*,ib-nb]-tube.x[*,ib-nb-1], 3 )
t1 = reform( tube.x[*,ib+nb+1]-tube.x[*,ib+nb], 3 )
t0 = t0/sqrt( total( t0^2 ) )
t1 = t1/sqrt( total( t1^2 ) )

dth_tot = acos( total( t0*t1 ) ); angle between tangents
if( dth_tot lt 0.03 ) then return;   less than 2 degree deflection

; center of circle
alpha = 0.5*total( ( x1 - x0 )*t0 )/( 1.0 - total( t0*t1 ) )
xc = 0.5*( x1 + x0 ) + alpha*( t1 - t0 )
rad = sqrt( total( (x0-xc)^2 ) );  radius of circle

d0 = ( x0 - xc )/rad;  unit vector to point 0

dth = dth_tot/float(2*nb)
for i=ib-nb+1, ib+nb-1 do begin
  th = dth*( i-ib+nb )
  tube.x[*,i] = xc + rad*cos(th)*d0 + rad*sin(th)*t0
endfor

return
end

; -----------------------------

pro round_rx_bend, tube, ib, nb=nb
;  cubic spline between points ib-nb and ib+nb

if( n_elements( nb ) lt 1 ) then nb=4;  9 points 
calc_tv, tube

x0 = reform( tube.x[*,ib-nb], 3 )
x1 = reform( tube.x[*,ib+nb], 3 )
dx0 = reform( tube.x[*,ib-nb]-tube.x[*,ib-nb-1], 3 )
dx1 = reform( tube.x[*,ib+nb+1]-tube.x[*,ib+nb], 3 )

s = findgen(2*nb+1)/float(2*nb);  normalized coordinate: 0 <= s <= 1

tube.x[0,(ib-nb):(ib+nb)] = s^2*(3.-2*s)*x1[0] + (1.-s)^2*(1.+2*s)*x0[0] $
                          + s^2*(s-1.)*(2*nb)*dx1[0] + s*(1.-s)^2*(2*nb)*dx0[0]
tube.x[2,(ib-nb):(ib+nb)] = s^2*(3.-2*s)*x1[2] + (1.-s)^2*(1.+2*s)*x0[2] $
                          + s^2*(s-1.)*(2*nb)*dx1[2] + s*(1.-s)^2*(2*nb)*dx0[2]

return
end
