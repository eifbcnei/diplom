malt = zeros(n);
ltheta1 = zeros([total, 1]);
ltheta2 = zeros([total, 1]);
threshold = 40;
y = zeros([length(angs), 1]);

for q = 1:total
    disp(q);
    for i = 1:l
        a1 = (v1 / 2)^0.5 * (normrnd(0,sigma)+normrnd(0,sigma)*1j);
        a2 = (v2 / 2)^0.5 * (normrnd(0,sigma)+normrnd(0,sigma)*1j);
        s1 = exp(1j*2*3.14*dlambda *((1:n)-1)' * sin(theta1));
        s2 = exp(1j*2*3.14*dlambda *((1:n)-1)' * sin(theta2));
        z = sigma*(1 / 2)^0.5 * (normrnd(0, 1, [n, 1]) + normrnd(0, 1, [n, 1])*1j);
        x = a1 * s1 + a2 * s2 + z;
        m = sigma^2 * (eye(n) + v1 * (s1 * s1') + v2 * (s2 * s2'));
        e = eig(m);
        malt = malt + x * x';
        if (i == 1 && q==1)
           figure
           scatter(e, zeros(size(e)), 'ro')
           hold on;
           grid on
        end

    end
    malt = malt / l;
    ealt = eig(malt);
    if (q==1)
        scatter((ealt), zeros(size(ealt)), 'bx')
        legend('lambda', 'lambda^');
        hold off
    end

    ytest = zeros([length(angs), 1]);
    for k = 1:length(angs)
        angi = angs(k);
        sphi = exp(1j*2*3.14*dlambda *((1:n)-1)' * sin(angi));
        wn = (1 / (sphi' * inv(malt) * sphi)) * (inv(malt) * sphi);
        y(k) =  y(k) + wn' * malt * wn;
        ytest(k) = wn' * malt * wn;
    end

    [~,rtheta]=findpeaks(real(ytest), rad2deg(angs));

    if (valid == 1)
        ltheta1(q) = rtheta(1);
        ltheta2(q) = rtheta(2);
    end
end
y = y / total;
figure;
plot(rad2deg(angs), y ./ max(y), 'bx-');
title('Capon')
grid on;

caponnorm = y ./ max(y);

meantheta1 = mean(transpose(ltheta1));
meantheta2 = mean(transpose(ltheta2));

sigmasqrtheta1 = 1/(total-1) * sum((ltheta1-meantheta1).^2);
sigmasqrtheta2 = 1/(total-1) * sum((ltheta2-meantheta2).^2);
