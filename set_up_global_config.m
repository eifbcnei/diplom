n = 10;
dlambda = 0.5;
num_tgt = 2;
ang1 = -5;
ang2 = -ang1;
theta1 = deg2rad(ang1);
theta2 = deg2rad(ang2);
v = 5/n;
v1=v;v2=v;
sigma = 1;
l = 30;
angs = deg2rad(-30:0.01:30);
total = 100;
specific_plots = 0;
rng default;

fouriervalid = 1;
logger = '';
yplot = 'b-';