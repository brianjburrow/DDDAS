subplot(2,2,1)
scatter(oDam(1,:), oCAPS, 0.1, 'black')
xlabel("Ply 1 Damage")
ylabel("Max Load Factor")
hold on
subplot(2,2,2)
scatter(oDam(2,:), oCAPS,  0.1,  'black')
xlabel("Ply 2 Damage")
ylabel("Max Load Factor")

hold on 
subplot(2,2,3)
scatter(oDam(3,:), oCAPS, 0.1,  'black')

xlabel("Ply 3 Damage")
ylabel("Max Load Factor")
hold on
subplot(2,2,4)
scatter(oDam(4,:), oCAPS, 0.1,'black')
xlabel("Ply 4 Damage")
ylabel("Max Load Factor")

figure(2)
scatter(vecnorm(oDam), oCAPS, 2.0, 'black')
xlabel('Damage Magnitude')
ylabel('Max Load Factor')