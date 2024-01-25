; a model of a flaring loop from 2011-09-13T23:15

forward_function tube_erg

l_fin = 60.0;  final length in Mm
b0 = 100.0;  field strength
derg = 1.0d11;  erg/cm^2 released into entire loop
dth = 90.0;  reconnection angle in degrees
tmax=4.0;  maximum temperature
tmin=0.01;  chromospheric temperature

; now compute properties

l_rx = 1.0d-8*2*!pi*derg/(b0*sin(0.25*dth/!radeg)^2)^2;  length of retracted section
d_len = 2*l_rx*sin(0.25*dth/!radeg)^2;   change due to retaction

l0 = l_fin + d_len

init_lam_rlf, tmin=1.0d6*(1.01*tmin)
preft, /verb
set_uniform_cs, b0
tube = rtv_tube_equilib( l0, tmin=tmin, chr=2.0, tmax=tmax, n=200 )
bend_tube, tube, dth; , /zero
set_heat_for_equilib, tube
print, 'n=', tube.n
print, 'min(dl) = ', min(tube.dl)

va = b0/sqrt(4*!pi*tube.rho[tube.n/2])

t_rx = 0.5*l_rx/va;  time to permit reconnection

; report on run
print, 'tube length = ', l0
print, '  reconn. length = ', l_rx
print, 'reconnection time = ', t_rx

window, 0, xs=600, ys=500

plot, tube.x[0,*], tube.x[2,*], /iso, yr=[-30,2], yst=1, tit=string(tube.time)

; now the run

dt = 0.5
nst = ceil( t_rx/dt )
tarr = [ tube ]

erg_i = tube_erg( tarr[0] )

window, 1, xs=600, ys=500
window, 2, xs=600, ys=900 
for i=1, nst do begin
  adv_tube, tube, dt, max=100000L
  tarr = [ tarr, tube ]
  wset, 1
  plot, tube.x[0,*], tube.x[2,*], /iso, yr=[-30,2], yst=1, tit=string(tube.time)
  wset, 2
  tube_plot, tube, fp=7.0, tr=[7.0d3, 50.0d6 ]
  erg = tube_erg( tube )
  print, ' erg: ', 1.0d8*b0*( erg.kin_par + erg.therm - erg_i.therm )
endfor

; now straightn and continue run
straighten_tube, tube
dt = 1.0
nst = 30

;window, 1, xs=600, ys=900
;window, 2, xs=600, ys=500 
for i=1, nst do begin
  adv_tube, tube, dt, max=500000L
  tarr = [ tarr, tube ]
  wset, 2
  tube_plot, tube, fp=7.0, vfp=0.25, tr=[7.0d3, 50.0d6 ]
  print, ' time:  ', tube.time
endfor

n = n_elements( tarr )

erg_f = tube_erg( tarr[n-1] )

terg = b0*( erg_f.kin + erg_f.therm - erg_i.therm )
print, terg*1.0d8

save, file='bend_run.sav', tarr

end
