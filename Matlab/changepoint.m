close all
clear

x = load('obitosBRmar2015jul2020.txt');
ydiff = load('ydiff', '-mat'); 
ydiff1 = struct2array(ydiff);
ydiff1 = ydiff1';
wt = dwt(ydiff1,'haar',5);

[pts_Opt, kopt, t_est] = wvarchg(ydiff1,2,1);
sprintf('The estimated change points are %d\n',pts_Opt)

decomp = load('LA8decomp', '-mat')

a4 = wrcoef('a',decomp.coefs,decomp.longs,'sym8',4);
d4 = wrcoef('d',decomp.coefs,decomp.longs,'sym8',4);
d3 = wrcoef('d',decomp.coefs,decomp.longs,'sym8',3);
d2 = wrcoef('d',decomp.coefs,decomp.longs,'sym8',2);
d1 = wrcoef('d',decomp.coefs,decomp.longs,'sym8',1);

figure
subplot(611),plot(ydiff1);
subplot(612),plot(d1);
subplot(613),plot(d2);
subplot(614),plot(d3);
subplot(615),plot(d2);
subplot(616),plot(a4);