figure;
plot(rad2deg(angs), fouriernorm, 'b-');
hold on;
plot(rad2deg(angs), caponnorm, 'g-');
plot(rad2deg(angs), heatnoisenorm, 'r-');
legend('Фурье', 'Кейпон', 'Тепловой шум');
grid on;