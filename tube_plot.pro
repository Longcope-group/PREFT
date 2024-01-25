pro tube_plot, tube, fp=fp, vmax=vmax, vfp=vfp, denmax=denmax, tr=tr, $
  lhalf=lhalf, init=init, pr=pr

if( not keyword_set( fp ) ) then fp=5.0
if( not keyword_set( denmax ) ) then denmax=1.0d13
if( not keyword_set( pr ) ) then pr=[0.1, 1000.0 ]
if( not keyword_set( tr ) ) then tr=[0.01, 100.0 ]*1.0d6

tit = 't = ' + string( tube.time ) + ' [s]'

calc_len, tube
if( n_elements( init ) ge 1 ) then calc_len, init

lh = 0.5*max(tube.l)
if( not keyword_set( lhalf ) ) then lhalf=lh

ypos=[0.07,0.29, 0.51, 0.73, 0.95 ]
n_e = tube.rho*tube.epamu/1.67d-8
vpar = total( tube.v*tube.tv_e, 1 )
if( not keyword_set( vmax ) ) then vmax = max( vpar )
if( not keyword_set( vfp ) ) then vfp = 0.25

plot_io, tube.l, n_e, xr=[0.0,fp], xst=1, $
  pos=[0.1, ypos[0], 0.4, ypos[1]], ytit='n!de!n [cm!u-3!n]', $
  xtit='x [ Mm ]', yr=[1.0d8,denmax]
if( n_elements( init ) ge 1 ) then begin
  n_e0 = init.rho*init.epamu/1.67d-8
  oplot, init.l, n_e0, lines=2
endif
plot_io, tube.l, n_e, xr=[0.0,lhalf], xst=1, yr=[1.0d8,denmax], $
  pos=[0.45, ypos[0], 0.9, ypos[1]], /noerase, xtit='x [ Mm ]', $
  ytickf=''
if( n_elements( init ) ge 1 ) then oplot, init.l, n_e0, lines=2
oplot, lh*[1,1], [ 1.0d8, denmax ], lines=1
; axis, /yax, yst=1, ytickf=''

plot_io, tube.l, 1.0d6*tube.t, xr=[0.0,fp], xst=1, /noerase, $
  pos=[0.1, ypos[1], 0.4, ypos[2]], xtickf='', $
  ytickf='', yr=tr
if( n_elements( init ) ge 1 ) then oplot, init.l, init.t*1.0d6, lines=2
plot_io, tube.l, 1.0d6*tube.t, xr=[0.0, lhalf], xst=1, $
  pos=[0.45, ypos[1], 0.9, ypos[2]], yst=8, /noerase, $
  xtickf='', ytickf='', yr=tr
if( n_elements( init ) ge 1 ) then oplot, init.l, init.t*1.0d6, lines=2
oplot, lh*[1,1], tr, lines=1

axis, /yax, yst=1, ytit='T [K]'

plot, tube.l, vpar, xr=[0.0,fp], xst=1, /noerase, $
  pos=[0.1, ypos[2], 0.4, ypos[3]], xtickf='', $
  yr=vfp*[-1.0,1.0], ytit='v [ Mm/s ]', yst=1
oplot, [0.0,fp], [0,0], lines=1
plot, tube.l, vpar, xr=[0.0,lhalf], xst=1, $
  pos=[0.45, ypos[2], 0.9, ypos[3]], /noerase, $
  xtickf='', yr=vmax*[-1.,1.], yst=1
oplot, [0.0,lhalf], [0,0], lines=1
oplot, lh*[1,1], vmax*[-1,1], lines=1

plot_io, tube.l, tube.p, xr=[0.0,fp], xst=1, /noerase, $
  pos=[0.1, ypos[3], 0.4, ypos[4]], xtickf='', $
  ytickf='', yr=pr, yst=1
if( n_elements( init ) ge 1 ) then oplot, init.l, init.p, lines=2
plot_io, tube.l, tube.p, xr=[0.0,lhalf], xst=1, /noerase, $
  pos=[0.45, ypos[3], 0.9, ypos[4]], yst=9, xtickf='', $
  ytickf='', yr=pr, tit=tit
if( n_elements( init ) ge 1 ) then oplot, init.l, init.p, lines=2
oplot, lh*[1,1], pr, lines=1

axis, /yax, yst=1, ytit='p [erg/cm!u3!n]'

return
end
