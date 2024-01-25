; set the parameters in a tube to be an RTV, constant pressure,
; equilibrium.  feet of the loop are at tmin located a distance chr from
; each end.  The apex, midway between these
; in physical length, is at tmax.

function rtv_tube_equilib, len, chr=chr, tmax=tmax, tmin=tmin, n=n, $
  debug=debug

if( not keyword_set(n) ) then n=400;  number of points in tube interior
if( not keyword_set(tmin) ) then tmin=0.03;  30,000K chromosphere
if( not keyword_set(tmax) ) then tmax=2.0
if( not keyword_set(chr) ) then chr=2.0;  length of chromopshere in Mm

nh = n/2

cs = ion_consts()
chi = cs.mpp*cs.epamu/(1.38d-10);         [ MK / erg ]
kap0 = 1.0d-6;                            [ erg / cm / s / K ]

; the temperature grid
fctr = 0.9
nb = 4000
nc = 240
; middle range
tb = (tmin/fctr)*exp(findgen(nb)*alog(fctr*fctr*tmax/tmin)/float(nb-1) )
na = ceil( 2*(tb[0]-tmin)/(tb[1]-tb[0]) )
ta = tmin + (tb[0]-tmin)*( findgen(na)/float(na-1) )^2
nc = ceil( 2*(tmax-tb[nb-1])/(tb[nb-1]-tb[nb-2]) )
tc = tmax + (tb[nb-1]-tmax)*( findgen(nc)/float(nc-1) )^2
tc = reverse( tc )
t = [ ta, tb[1:(nb-2)], tc ]
nt = n_elements( t )
dlt = 2*( t - shift(t,1) )/( t + shift(t,1) )
dlt[0] = dlt[1]

;  the computations
lam = lam_rlf( t*1.0d6 );        [ erg cm^3 / s ]
f = 0.1*total( t^(1.5)*lam*dlt/kap0, /cum );   [ 1.0e10 cm^4 K^5 ]
rt = ( f - f[nt-1]*( t^(3.5)-tmin^(3.5) )/( tmax^(3.5)-tmin^(3.5) ) ) > 1.0d-22
xi = 100*0.707107*t^(2.5)/sqrt(rt);            [ 1.0e8 cm^(-2) ]
ig = 0.5*( xi*dlt + shift( xi*dlt, 1 ) )
ig[0] = 0.5*xi[0]*dlt[0]
col_e = total( ig, /cum );                     [ 1.0e8 cm^(-2) ]
; col_e = shift( col_e, 1 )
; col_e[0] = 0.0

if( keyword_set( debug ) ) then begin
  window, 0
  plot_oi, t, dlt
  window, 1
  plot_oi, t, col_e
endif

p = (2.0/chi/len)*total( t*xi*dlt );           [ erg / cm^3 ]
p = 1.02*p;  help smooth out ends

; n_e0 = chi*p/tmin
; n_ea = chi*p/tmax
; print, n_e0, n_ea

; create the grid
col_norm = col_e/max(col_e);               runs from 0 to 1
grid_dm = (dindgen(nh)+3.0)^(0.5);         column increases
; grid_dm = replicate( 1.0, nh )
grid_col = [ 0.0, total( grid_dm[0:(nh-2)], /cum ) ]
tgh = interpol( t, col_norm, grid_col/grid_col[nh-1] )
;  dl ~ dm/rho ~ dm*T/p
dlh = tgh*grid_dm
l_tot = total( dlh )
dlh = 0.5*len*dlh/l_tot

; arrays for the central portion
dl_c = [ dlh, reverse( dlh[0:(nh-2)] ) ]
t_c = [ tgh, reverse( tgh[0:(nh-2)] ) ]

; chromosphere
chr_exp = 0.5
nchr = ceil( 3.0^(chr_exp)*( (chr_exp+1.0)*chr/dlh[0] )^(1.0/(chr_exp+1.0)) );      
dl_chr = dlh[0]*( ( 3.0 + findgen(nchr) )/3.0 )^(chr_exp)


dl = [ reverse( dl_chr ), dl_c, dl_chr ]
l = [ 0.0, total( dl, /cum ) ];  the full length array

t_chr = t_c[0]; the temperature of the chromosphere; MK
t_arr = [ replicate( t_chr, nchr ), t_c, replicate( t_chr, nchr ) ]

z_chr = l[0:(nchr-1)] - l[nchr-1];  depth from top of chromosphere
gp = 2.74d-4;   g in Mm/s/s
r = 0.01*(1.38/1.67/cs.mpp)
h_ch = t_chr*r/gp
p_chr = p*exp( -z_chr/h_ch );   exponential atmosphere


; place result into a tube
n_tot = n_elements(l)
tube = make_tube( n_tot )


tube.x[0,*] = l
tube.p = [ p_chr, replicate( p, n_elements( dl_c ) ), reverse( p_chr ) ]

tube.t = [ t_arr, t_arr[0] ]
tube.p[n_tot-1] = tube.p[n_tot-2]
tube.gpar[0:(nchr-1)] = -gp
tube.gpar[(n_tot-nchr):(n_tot-1)] = gp

tube.rho = tube.p/r/tube.t
set_rho, tube
calc_len, tube
set_heat_for_equilib, tube

return, tube
end
