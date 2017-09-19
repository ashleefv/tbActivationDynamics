sol = WDetModel('Wigginton-model-parameters.mat');

figure(1)
plot(sol.x, sol.y(1,:), 'r-', sol.x, sol.y(2,:), 'g-', sol.x, sol.y(3,:), 'b-');
legend('M_R', 'M_I', 'M_A');

figure(2)
plot(sol.x, sol.y(4,:), 'r-', sol.x, sol.y(5,:), 'g-', sol.x, sol.y(6,:), 'b-');
legend('T_0', 'T_1', 'T_2');

figure(3)
plot(sol.x, sol.y(7,:), 'r-', sol.x, sol.y(8,:), 'g-', sol.x, sol.y(9,:), 'b-', sol.x, sol.y(10,:), 'k-');
legend('I_y', 'I_12', 'I_10', 'I_4');

figure(4)
plot(sol.x, sol.y(11,:), 'r-', sol.x, sol.y(12,:), 'g-');
legend('B_E', 'B_I');