clear all 

% Parameters
beta = 0.99;
kappa = 0.13;
rho = 0.8;

% Coefficients from numerical solution (approximated)
a = 0.701502; % 1.43991; 
b = 0.262662; % -1.70648; 
c = 0.885;
d = 0.115;  % Not used here since eta_t = 0

% Simulation settings
T = 21;                     % Number of periods
y = zeros(T, 1);            % Output
delta_M = zeros(T, 1);      % Money growth
eps = zeros(T, 1);          % Monetary shocks
eta = zeros(T, 1);          % Productivity shocks

% Apply one-unit shock at t = 0 -- or 1 standard deviation for eps
eps(1) = 1;
%eps(1) = 0.00066;

% Simulate model
for t = 1:T
    if t == 1
        delta_M(t) = eps(t);
        y(t) = c * eps(t);
    else
        delta_M(t) = rho * delta_M(t-1) + eps(t);
        y(t) = a * y(t-1) + b * delta_M(t-1) + c * eps(t);
    end
end

% Plot impulse response
figure;
plot(0:T-1, y, 'o-', 'LineWidth', 2);
title('True Impulse Response of Output to a 1-Monetary Shock','FontSize',9);
xlabel('Time after shock (t)');
ylabel('Output (y_t)');
grid on;

% Save figure
saveas(gcf, 'out/2_IRF_theoretical_y.png');
