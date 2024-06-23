%load_water

mat_water.set('family', 'water');
mat_water.propertyGroup('def').set('thermalexpansioncoefficient', '');
mat_water.propertyGroup('def').set('bulkviscosity', '');
mat_water.propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(T)' '0' '0' '0' 'alpha_p(T)' '0' '0' '0' 'alpha_p(T)'});
mat_water.propertyGroup('def').set('bulkviscosity', 'muB(T)');
mat_water.propertyGroup('def').descr('thermalexpansioncoefficient_symmetry', '');
mat_water.propertyGroup('def').descr('bulkviscosity_symmetry', '');
mat_water.propertyGroup('def').set('dynamicviscosity', 'eta(T)');
mat_water.propertyGroup('def').descr('dynamicviscosity_symmetry', '');
mat_water.propertyGroup('def').set('ratioofspecificheat', 'gamma_w(T)');
mat_water.propertyGroup('def').descr('ratioofspecificheat_symmetry', '');
mat_water.propertyGroup('def').set('electricconductivity', {'5.5e-6[S/m]' '0' '0' '0' '5.5e-6[S/m]' '0' '0' '0' '5.5e-6[S/m]'});
mat_water.propertyGroup('def').descr('electricconductivity_symmetry', '');
mat_water.propertyGroup('def').set('heatcapacity', 'Cp(T)');
mat_water.propertyGroup('def').descr('heatcapacity_symmetry', '');
mat_water.propertyGroup('def').set('density', 'rho(T)');
mat_water.propertyGroup('def').descr('density_symmetry', '');
mat_water.propertyGroup('def').set('thermalconductivity', {'k(T)' '0' '0' '0' 'k(T)' '0' '0' '0' 'k(T)'});
mat_water.propertyGroup('def').descr('thermalconductivity_symmetry', '');
mat_water.propertyGroup('def').set('soundspeed', 'cs(T)');
mat_water.propertyGroup('def').descr('soundspeed_symmetry', '');
mat_water.propertyGroup('def').addInput('temperature');
mat_water.propertyGroup('def').func.create('eta', 'Piecewise');
mat_water.propertyGroup('def').func('eta').set('arg', 'T');
mat_water.propertyGroup('def').func('eta').set('pieces', {'273.15' '413.15' '1.3799566804-0.021224019151*T^1+1.3604562827E-4*T^2-4.6454090319E-7*T^3+8.9042735735E-10*T^4-9.0790692686E-13*T^5+3.8457331488E-16*T^6'; '413.15' '553.75' '0.00401235783-2.10746715E-5*T^1+3.85772275E-8*T^2-2.39730284E-11*T^3'});
mat_water.propertyGroup('def').func('eta').set('argunit', 'K');
mat_water.propertyGroup('def').func('eta').set('fununit', 'Pa*s');
mat_water.propertyGroup('def').func.create('Cp', 'Piecewise');
mat_water.propertyGroup('def').func('Cp').set('arg', 'T');
mat_water.propertyGroup('def').func('Cp').set('pieces', {'273.15' '553.75' '12010.1471-80.4072879*T^1+0.309866854*T^2-5.38186884E-4*T^3+3.62536437E-7*T^4'});
mat_water.propertyGroup('def').func('Cp').set('argunit', 'K');
mat_water.propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
mat_water.propertyGroup('def').func.create('rho', 'Piecewise');
mat_water.propertyGroup('def').func('rho').set('arg', 'T');
mat_water.propertyGroup('def').func('rho').set('smooth', 'contd1');
mat_water.propertyGroup('def').func('rho').set('pieces', {'273.15' '293.15' '0.000063092789034*T^3-0.060367639882855*T^2+18.9229382407066*T-950.704055329848'; '293.15' '373.15' '0.000010335053319*T^3-0.013395065634452*T^2+4.969288832655160*T+432.257114008512'});
mat_water.propertyGroup('def').func('rho').set('argunit', 'K');
mat_water.propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
mat_water.propertyGroup('def').func.create('k', 'Piecewise');
mat_water.propertyGroup('def').func('k').set('arg', 'T');
mat_water.propertyGroup('def').func('k').set('pieces', {'273.15' '553.75' '-0.869083936+0.00894880345*T^1-1.58366345E-5*T^2+7.97543259E-9*T^3'});
mat_water.propertyGroup('def').func('k').set('argunit', 'K');
mat_water.propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
mat_water.propertyGroup('def').func.create('cs', 'Interpolation');
mat_water.propertyGroup('def').func('cs').set('table', {'273' '1403';  ...
'278' '1427';  ...
'283' '1447';  ...
'293' '1481';  ...
'303' '1507';  ...
'313' '1526';  ...
'323' '1541';  ...
'333' '1552';  ...
'343' '1555';  ...
'353' '1555';  ...
'363' '1550';  ...
'373' '1543'});
mat_water.propertyGroup('def').func('cs').set('interp', 'piecewisecubic');
mat_water.propertyGroup('def').func('cs').set('argunit', 'K');
mat_water.propertyGroup('def').func('cs').set('fununit', 'm/s');
mat_water.propertyGroup('def').func('cs').set('allowrand', 'off');
mat_water.propertyGroup('def').func.create('an1', 'Analytic');
mat_water.propertyGroup('def').func('an1').set('funcname', 'alpha_p');
mat_water.propertyGroup('def').func('an1').set('expr', '-1/rho(T)*d(rho(T),T)');
mat_water.propertyGroup('def').func('an1').set('args', {'T'});
mat_water.propertyGroup('def').func('an1').set('argunit', 'K');
mat_water.propertyGroup('def').func('an1').set('fununit', '1/K');
mat_water.propertyGroup('def').func('an1').set('plotargs', {'T' '273.15' '373.15'});
mat_water.propertyGroup('def').func.create('an2', 'Analytic');
mat_water.propertyGroup('def').func('an2').set('funcname', 'gamma_w');
mat_water.propertyGroup('def').func('an2').set('expr', '1+(T/Cp(T))*(alpha_p(T)*cs(T))^2');
mat_water.propertyGroup('def').func('an2').set('args', {'T'});
mat_water.propertyGroup('def').func('an2').set('argunit', 'K');
mat_water.propertyGroup('def').func('an2').set('fununit', '1');
mat_water.propertyGroup('def').func('an2').set('plotargs', {'T' '273.15' '373.15'});
mat_water.propertyGroup('def').func.create('an3', 'Analytic');
mat_water.propertyGroup('def').func('an3').set('funcname', 'muB');
mat_water.propertyGroup('def').func('an3').set('expr', '2.79*eta(T)');
mat_water.propertyGroup('def').func('an3').set('args', {'T'});
mat_water.propertyGroup('def').func('an3').set('argunit', 'K');
mat_water.propertyGroup('def').func('an3').set('fununit', 'Pa*s');
mat_water.propertyGroup('def').func('an3').set('plotargs', {'T' '273.15' '553.75'});