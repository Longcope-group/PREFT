; Version 2.5
;    2013-08-15, 2014-01-22, 2014-05-15, 2014-07-11, 2014-12-23, 2015-02-12
; Version 1.0
;    Dana Longcope 2012-12-18
;
; This file contains code necessary to advance TFT equations: 
;   IDL> adv_tube, tube, dt
; calls are made to the following externally defined procedures:
;   pro field_at_points, tube
;   
; ---------------------
;
;  units:  
;          length        [ Mm = 1.0e8 cm ]
;          time          [ sec ]
;          energy        [ microflare = uf = 1.0e24 erg ]
;          temperature   [ MK = 1.0e6 K ]
;
;          velocity [ Mm/sec = 1.0e8 cm/sec ]
;          pressure [ uf/Mm^3 = erg/cm^3 ]
;          B        [ G ]
;          mass     [ HT = uf s^2 / Mm^2 = 1.0e8 gm ]
;          mass density [ HT/Mm^2 = 1.0e-16 gm/cm^3 ]
;
;  grid:
;
;
;   x, v, k, rho_e  
; edges        0     1     2     3     n-3   n-2   n-1      
;              +--o--+--o--+--o--+-->>--+--o--+--o--+
; cells           0     1     2            n-3  n-2
;   t, dm, p, 
;   rho, dl, tv
;   
;
;    shifting:   cell -> edge:  rho_e = 0.5*( shift( rho, 1 ) + rho )
;                edge -> cell:  v_c = 0.5*( shift( v, -1 ) + v )
;                               dx = shift( tube.x, 0, -1 ) - tube.x
;
; ---------------------------------

function ion_consts
; set the physical constants in the tube
;  consider only H & He => 0.7% error
;  assume complete ionization of both H & He

a_he = 10.0d0^( -1.075 );    solar He abundance = 10.925
mph = 1.0 + 4*a_he;          mass per H nucleus [ amu ]
pph = 2.0 + 3*a_he;          particles per H nucleus
eph = 1.0 + 2*a_he;          electrons per H nucleus

mpp = mph/pph
epamu = eph/mph

; constants for viscosity & conductivity
mu0 = 0.12;  mu = mu0*(t)^(2.5)
kap0 = 0.1;  kap = map0*(t)^(2.5)

res = { mpp:mpp, epamu:epamu, mu0:mu0, kap0:kap0 }

return, res
end

; -----------------------------

function make_tube, n

if( n_params() lt 1 ) then n=25

;  primative quantities
x = dblarr( 3, n );  position [ Mm ]
v = fltarr( 3, n );  velocity [ Mm/sec ]
t = dblarr( n );  temperature [ 1.0e6 K ]
dm = dblarr( n ); differential mass

time = 0.0d0;   [ sec ]
gam = 5.0/3.0;
; mpp = 0.5;     mass per particle [ amu ]  
; epamu = 1.0;   number of electrons per amo

;  derived quantities
rho = fltarr( n );  mass density [ 1.0e-16 gm/cm^2 ]
dl = dblarr( n ); differential length 
dl_e = dblarr( n ); differential length @ edges
l = dblarr( n ); total length
tv = dblarr( 3, n ); the tangent vector
tv_e = dblarr( 3, n ); the tangent vector @ edges
k = dblarr( 3, n ); the curvature vector
p = fltarr( n );  pressure [ erg/cm^3 ]
b = fltarr( n );  the field strength @ edges
db = fltarr( 3, n ); the gradient of field strength (@ edges) 
gpar = fltarr( n );  parallel component of grav. acceleration [ Mm/s/s ]
mu = fltarr( n );  dynamic viscosity  [ gm/cm/sec ]
kap = fltarr( n );  thermal conductivity  [ 1.0e10 erg/cm/K/sec ]
dvn = fltarr( n );  tv.(dv/dl)  @ centers
heat = fltarr( n );  heating rate integrated over a cell [ 1.0e8 erg/sec/cm^2 ]
rad_loss = fltarr( n );  radiative losses integrated over a cell [ 1.0e8 erg/sec/cm^2 ]
;ne = fltarr ( n  ); number density ; ashfield 3_25

const_str = ion_consts()

tube = { n:n, time:time, gam:gam, x:x, v:v, l:l, dl:dl, dl_e:dl_e, $
         dm:dm, tv:tv, tv_e:tv_e, k:k, rho:rho, t:t, p:p, b:b, db:db, gpar:gpar, $
         dvn:dvn, mu:mu, kap:kap, heat:heat, rad_loss:rad_loss, net_erg_loss:0.0, $
         mpp:const_str.mpp, epamu:const_str.epamu, mu0:const_str.mu0, kap0:const_str.kap0, $
         inv_hflf:1.0 };, ne:ne}

return, tube
end

; -----------------------------

pro calc_tv, tube
;  geometric factors: dl, dl_e, tv, tv_e, k

; tangent vector: at centers
dx = shift( tube.x, 0, -1 ) - tube.x
dx[*,tube.n-1] = dx[*,tube.n-2]

tube.dl = sqrt( total( dx^2, 1 ) )
tube.tv = dx/( [1,1,1] # tube.dl )

;  tangent vector at edges: defined to be perp. to k at edges
tube.tv_e = 0.5*( tube.tv + shift( tube.tv, 0, 1 ) )
tube.tv_e[*,0] = tube.tv[*,0]
tm = sqrt( total( tube.tv_e^2, 1 ) )
tube.tv_e = tube.tv_e/( [1,1,1] # tm )

;  distance between centers: on edges
tube.dl_e = 0.5*( tube.dl + shift( tube.dl, 1 ) )
tube.dl_e[0] = tube.dl_e[1]

; curvature vector: on edges
dt = tube.tv - shift( tube.tv, 0, 1 )
dt[*,0] = dt[*,1]

tube.k = dt/( [1,1,1] # tube.dl_e )

return
end

; -----------------------------

pro calc_prho, tube
;  the density & then pressure from primative variables
;   t and dm.  
;  ************ calc_tv and field_at_points both must be called first

r = 0.01*(1.38/1.67/tube.mpp);        [ erg / cm^3 / MK / 1.0e-16 gm ]

tube.rho = tube.dm*tube.b/tube.dl;     mass density [ 1.0e-16 gm/cm^3 ]
tube.p = r*tube.rho*tube.t;            pressure     [ ergs ]
;tube.ne = tube.epamu*tube.rho/1.67d-8

return
end

; -----------------------------

pro set_rho, tube
;  tube.rho is set.  compute tube.dm to enforce this

calc_tv, tube;   in case this is not current
field_at_points, tube
tube.dm = tube.dl*tube.rho/tube.b

return
end 

; -----------------------------

pro calc_len, tube

tube.l[0] = 0.0
for i=1, tube.n-1 do tube.l[i] = tube.l[i-1] + tube.dl[i]

return
end

; -----------------------------

pro calc_rad_loss, tube

n_e = tube.epamu*tube.rho/1.67d-8;                 [ cm^{-3} ]
lam = lam_rlf( 1.0d6*tube.t );                     [ erg cm^3/sec ]
n_eTOP = n_e*n_e; < 1.0d22 ; set uppper threshold of ne to be 10d11, fix 'rad instabiliteis' - whaiv
tube.rad_loss = n_eTOP*lam*tube.dl;               [ 1.0e8 erg/cm^2/sec ]

return
end

; -----------------------------

pro set_kappa, tube, inv_hflf
; impose the free-streaming limit on kappa

; temp at edges
t_e = 0.5*( tube.t + shift( tube.t, 1 ) )
t_e[0] = t_e[1]

te32 = t_e*sqrt(t_e);  t_e^(1.5)  built for speed
kap_sp = tube.kap0*t_e*te32;       spitzer conductivity [ 1.0e16 erg/cm/sec/MK ]

if( inv_hflf gt 0.01 ) then begin
  grad_t = ( tube.t - shift( tube.t, 1 ) )/tube.dl_e;        [ 1.0e-8 MK/cm ]
  grad_t[0] = 0.0;
  grad_t[tube.n-1] = 0.0

  n_ee = 0.5*( tube.rho + shift( tube.rho, 1 ) )*tube.epamu/1.67d-8;              [ cm^{-3} ]
  ;   kb^(1.5) * m_e^(-0.5) = 5.37e-11 erg cm / s / K^(1.5)
  ;                         = 0.0537  erg cm / s / (MK)^(1.5)
  f_fs = 1.5d-8*0.0537*n_ee*te32;                    [ 1.0d8 erg/cm^2/s ]
  f_fs[0] = f_fs[1]

  fctr = inv_hflf*kap_sp*grad_t/f_fs
  rat = sqrt( 1.0 + fctr*fctr )
endif else rat = 1.0

tube.kap = kap_sp/rat;                       set to limit

return
end


; -----------------------------


pro set_min_visc, tube, fctr
;  need mu > rho*c_s*dl = \sqrt( gamma p * rho )* dl

mu_min = fctr*tube.dl*sqrt( tube.gam*tube.p*tube.rho )
tube.mu = tube.mu > mu_min

return
end

; -----------------------------

pro calc_dv, tube, dv

; prepare standard tube variables
calc_tv, tube
field_at_points, tube
calc_prho, tube

; density at edges
; rho_e = 0.5*( tube.rho + shift( tube.rho, 1 ) )
mb = tube.dm*tube.b
rho_e = 0.5*( mb + shift( mb, 1 ) )/tube.dl_e;     mass density [ 1.0e-16 gm/cm^3 ]
rho_e[0] = rho_e[1]

; pressure gradient
dp = ( tube.p - shift( tube.p, 1 ) )/tube.dl_e

; viscous stress
tube.mu = tube.mu0*tube.t*tube.t*sqrt(tube.t);   dynamic viscosity [ gm/cm/sec ]
set_min_visc, tube, 0.5;    changed from 0.05 - whaiv       impose a minimum for resolution

dvdl = shift( tube.v, 0, -1 ) - tube.v
dvdl[*,tube.n-1] = dvdl[*,tube.n-2]
tube.dvn = total( dvdl*tube.tv, 1 )/tube.dl

sig = tube.mu*tube.dvn
dsig = ( sig - shift( sig, 1 ) )/tube.dl_e

tension = (tube.b*tube.b-4*!pi*tube.p) > (0.05*tube.b*tube.b);  an ad-hoc limit to prevent fire hose
dv = 0.25*tube.k*([1,1,1]#tension)/!pi - tube.tv_e*([1,1,1]#(dp-dsig))

; ----- other viscosity term [added 7/23/13]
sig_e = 0.5*(  sig + shift( sig, 1 ) )
dv = dv + tube.k*([1,1,1]#sig_e)

dv = dv/([1,1,1]#rho_e)

; gravity
dv = dv + ([1,1,1]#tube.gpar)*tube.tv_e

; apply BC
dv[*,0] = 0.0
dv[*,tube.n-1] = 0.0

return
end

; -----------------------------

function calc_dfc_direct, tube, fc=fc
; direct computation.  used to check matrix computation
;  fc is conductive flux [ 1.0e8 erg/cm^2/s ]

; temp at edges
t_e = 0.5*( tube.t + shift( tube.t, 1 ) )
t_e[0] = t_e[1]

tube.kap = 0.1*(t_e)^(2.5);    conductivity [ 1.0e10 erg/cm/sec/K ]

fc = tube.kap*( tube.t - shift( tube.t, 1 ) )/tube.dl_e
fc[0] = 0.0;           no heat flux across ends
fc[tube.n-1] = 0.0
dfc = ( shift( fc, -1 ) - fc )/tube.dl;  div( Fc ): in cell

dfc = (tube.gam-1.0)*tube.t*dfc/tube.p

; apply BC
dfc[0] = 0.0
dfc[tube.n-2] = 0.0
dfc[tube.n-1] = 0.0

return, dfc
end

; -----------------------------

pro calc_dtemp, tube, dtemp
; return DT/dt
; all changes except thermal conduction
; *** variables are prepared from previous call to calc_dv.  including tube.dvn

; conductive flux: at edges cell -- for explicit advance
; dfc = calc_dfc_direct( tube )

; viscous heating: in cells
vheat = tube.mu*tube.dvn^2

; adiabatic contribution
dlnBdt = total( tube.v*tube.db, 1 )/tube.b;  at edges
dlnTdt_ad = 0.5*( shift( dlnBdt, -1 ) + dlnBdt ) - tube.dvn

calc_rad_loss, tube
dtemp = (tube.gam-1.0)*tube.t*( dlnTdt_ad + ( vheat + ( tube.heat - tube.rad_loss )/tube.dl )/tube.p )
; :::::: for adiabatic test:
; dtemp = (tube.gam-1.0)*tube.t*dlnTdt_ad

; ################### attempt to mimic ionization of H & He
; dtemp = dtemp*ionization_comp( tube )

; apply BC
dtemp[0] = 0.0
dtemp[tube.n-2] = 0.0
dtemp[tube.n-1] = 0.0

return
end

; -----------------------------

function dfc_mat, tube
; compute the tridiagonal matrix which used to compute div( Fc )
;   dfc = mat[0,*]*shift( tube.t, 1 ) + mat[1,*]*tube.t + mat[2,*]*shift( tube.t, -1 )
;

; temp at edges
t_e = 0.5*( tube.t + shift( tube.t, 1 ) )
t_e[0] = t_e[1]

; tube.kap = tube.kap0*(t_e)^(2.5);    conductivity [ 1.0e10 erg/cm/sec/K ]
set_kappa, tube, tube.inv_hflf

mat = fltarr( 3, tube.n )

; fc = tube.kap*( tube.t - shift( tube.t, 1 ) )/tube.dl_e
; dfc = ( shift( fc, -1 ) - fc )/tube.dl;  div( Fc ): in cell
; dfc = (tube.gam-1.0)*tube.t*dfc/tube.p

; coefficient of shift( tube.t, 1 ) : sub-diagonal
mat[0,*] = tube.kap/tube.dl_e/tube.dl

; coefficient of shift( tube.t, -1 ) : super-diagonal
mat[2,*] = shift( tube.kap/tube.dl_e, -1 )/tube.dl

; coefficient of tube.t = shift( shift( tube.t, 1 ), -1 )
mat[1,*] = -tube.kap/tube.dl_e/tube.dl - shift( tube.kap/tube.dl_e, -1 )/tube.dl

mat = (tube.gam-1.0)*mat*( [1,1,1] # ( tube.t/tube.p ) )

; ################### attempt to mimic ionization of H & He
; mat = mat*( [1,1,1] #ionization_comp( tube ) )

; BC
mat[*,0] = 0.0
mat[*,tube.n-2] = 0.0
mat[*,tube.n-1] = 0.0

return, mat
end

; -----------------------------

pro temp_adv_impl, tube, dt

mat = -dt*dfc_mat( tube )
mat[1,*] = mat[1,*] + 1.0; add identity

tnp1 = trisol( reform( mat[0,*] ), reform( mat[1,*] ), reform( mat[2,*] ), tube.t )
tube.t = tnp1

return
end


; -----------------------------

pro one_step, tube, dt, isotherm=isotherm

safe = 0.51

; advance positions
xt = tube.x
vt = tube.v

calc_dv, tube, dv
tube.x = tube.x + safe*dt*tube.v
tube.v = vt + safe*dt*dv

calc_dv, tube, dv
tube.x = xt + dt*tube.v
tube.v = vt + dt*dv

; advance temperature
if( not keyword_set( isotherm ) ) then begin
; full explicit time step for everything but conduction
  tt = tube.t
  calc_dtemp, tube, dtemp
  tube.t = tube.t + safe*dt*dtemp

  calc_dtemp, tube, dtemp
  tube.t = tt + dt*dtemp;  tube.t now contains result of full time step.
; full implicit time step for conduction
  temp_adv_impl, tube, dt
endif

tube.net_erg_loss = tube.net_erg_loss + dt*total( tube.rad_loss - tube.heat )
tube.time = tube.time + dt

return
end

; -----------------------------

pro one_step_expl, tube, dt, isotherm=isotherm
; explicit time-stepping

safe = 0.51

; advance positions
xt = tube.x
vt = tube.v

calc_dv, tube, dv
tube.x = tube.x + safe*dt*tube.v
tube.v = vt + safe*dt*dv

calc_dv, tube, dv
tube.x = xt + dt*tube.v
tube.v = vt + dt*dv

; advance temperature
if( not keyword_set( isotherm ) ) then begin
  tt = tube.t
  calc_dtemp, tube, dtemp
  dfc = calc_dfc_direct( tube )
  tube.t = tube.t + safe*dt*( dtemp + dfc )

  calc_dtemp, tube, dtemp
  dfc = calc_dfc_direct( tube )
  tube.t = tt + dt*( dtemp + dfc )
endif

tube.time = tube.time + dt

return
end

; -----------------------------

function cfl_step, tube, include
; include[i] is array of 1 & 0 to include/exclude terms:
;   i=0: advection
;     1: alfven
;     2: sound
;     3: viscosity
;     4: thermal conduction
;     5: radiative losses

va = tube.b/sqrt(4*!pi*tube.rho)
cs = sqrt(tube.gam*tube.p/tube.rho)

f_adv = max( abs( tube.dvn ) );                      i=0
f_cs = max(cs/tube.dl);			          1
f_va = max(va/tube.dl);			          2
f_visc = max( tube.mu/tube.dl^2/tube.rho );            3
f_kap = max( tube.kap*tube.t/tube.dl^2/tube.p );       4
f_rlf = max( tube.rad_loss/tube.dl/tube.p )

f_tot = total( [ f_adv, f_cs, f_va, f_visc, f_kap, f_rlf ]*include )

return, 1.0/f_tot
end

; -----------------------------

pro adv_tube, tube, dtime, dtmax=dtmax, maxstep=maxstep, safe=safe, isotherm=isotherm, cnt=cnt

if( not keyword_set( dtmax ) ) then dtmax=0.01
if( not keyword_set( maxstep ) ) then maxstep=10000L
if( not keyword_set( safe ) ) then safe=0.1
if( not keyword_set( isotherm ) ) then isotherm=0

dtmin=1.0d-15
t=dtmin
cnt = 0L 

; take small step to prepare for clf_step
one_step, tube, dtmin, isoth=isotherm

repeat begin
  if( ( cnt mod 5 ) eq 0 ) then $;  only reset time step every 5 steps, to save time
  dtc = cfl_step( tube, [1,1,1,1,0,1] )
  dtf = dtime - t
  dt = min( [ dtmax, dtc*safe, dtf ] ) > dtmin
  one_step, tube, dt, isoth=isotherm
  t = t + dt
  cnt = cnt + 1
endrep until( ( t ge dtime ) or ( cnt gt maxstep ) )

return
end

; -----------------------------

function tube_erg, tube
; energy in units of 1.0e8 erg/Mx

; kinetic
; rho_e = 0.5*( tube.rho + shift( tube.rho, 1 ) )
mb = tube.dm*tube.b
rho_e = 0.5*( mb + shift( mb, 1 ) )/tube.dl_e;     mass density [ 1.0e-16 gm/cm^3 ]
rho_e[0] = rho_e[1]
b_e = 0.5*( tube.b + shift( tube.b, 1 ) )
b_e[0] = b_e[1]
vm2 = total( tube.v^2, 1 )
kin = total( 0.5*rho_e*vm2*tube.dl_e/tube.b )

vpl = total( tube.v*tube.tv_e, 1 )
kin_par = total( 0.5*rho_e*vpl*vpl*tube.dl_e/tube.b )

therm = total( tube.p*tube.dl/tube.b )/(tube.gam-1.0)
mag = total( tube.b*tube.dl )/(4*!pi)

tot = kin + therm + mag

; gravitational potential energy
phi = total( tube.gpar*tube.dl_e, /cum )
grav = total( tube.dm*phi )

; energy loss from radiation - whaiv
rad = total(tube.rad_loss/tube.b)


erg = { tot:tot, kin:kin, kin_par:kin_par, therm:therm, mag:mag, grav:grav, rad:rad }

return, erg
end

; -----------------------------

pro write_tube, tube, append=append, cl=cl

if( not keyword_set( append ) ) then openw, 1, 'preft.ftb'

printf, 1, tube.n, tube.time, form=' (i8, d15.9)'
for i=0, tube.n-1 do printf, 1, tube.x[*,i], tube.v[*,i], tube.t[i], tube.dm[i], $
  form='(8e14.6)'

if( keyword_set( cl ) ) then close, 1

return
end

; -----------------------------

pro preft, verbos=verbose

if( keyword_set( verbose ) ) then begin
  print, systime( 0 )
  openr, 1, 'preft.pro'
  line = 'abc'
  readf, 1, line
  print, line
  readf, 1, line
  print, line
  close, 1
endif

return
END
