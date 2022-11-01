n = 10;
dlambda = 0.5;
num_tgt = 2;
theta1 = deg2rad(-10);
theta2 = deg2rad(15);
v = 10;
v1=v;v2=v;
sigma = 1;
l = 5000;
rng default;
malt = zeros(n);

for i = 1:l
    a1 = (v1 / 2)^0.5 * (normrnd(0,sigma)+normrnd(0,sigma)*1j);
    a2 = (v2 / 2)^0.5 * (normrnd(0,sigma)+normrnd(0,sigma)*1j);
    s1 = exp(1j*2*3.14*dlambda *((1:n)-1)' * sin(theta1));
    s2 = exp(1j*2*3.14*dlambda *((1:n)-1)' * sin(theta2));
    z = (sigma / 2)^0.5 * (normrnd(0, 1, [n, 1]) + normrnd(0, 1, [n, 1])*1j);
    x = a1 * s1 + a2 * s2 + z;
    m = sigma^2 * (eye(n) + v1 * (s1 * s1') + v2 * (s2 * s2'));
    e = eig(m);
    malt = malt + x * x';
    if (i == 1)
       figure
       scatter(e, zeros(size(e)), 'ro')
       hold on;
       grid on
    end
end
malt = malt / l;
ealt = eig(malt);

scatter((ealt), zeros(size(ealt)), 'bx')
hold off