clear all 

% Parameters from previous solutions
a = 0; %0.701502; % 1.43991; 
b = 0.49232; %0.262662; % -1.70648; 
c = 0.6154; %0.885;
d = 0.115;  % Not used here since eta_t = 0
rho = 0.8;
T = 600; % I will burn the first 100
N_sim = 500;

rng(456,'twister');
% rng(123456,'twister');

% Standard deviations (set so that Var(ΔM) = Var(η))
sigma_eta = 0.007;
sigma_eps = 0.00066;

H = 20;  % horizon

% shock to be simulated in Jorda's beta_h's
shock_eps = 0.00066;

% Storage for LP with control (y_{t-1})
IRF_LP_control = zeros(H+1, N_sim);

for s = 1:N_sim
    % Simulate from true model
    eta = sigma_eta * randn(T + H + 1, 1);
    % eta = zeros(T+H+1,1);
    eps = sigma_eps * randn(T + H + 1, 1);
    deltaM = zeros(T + H + 1, 1);
    y = zeros(T + H + 1, 1);

    for t = 2:T + H + 1
        deltaM(t) = rho * deltaM(t-1) + eps(t);
        y(t) = a * y(t-1) + b * deltaM(t-1) + c * eps(t) + d * eta(t);
    end

    % LP estimation at each horizon h
    for h = 0:(H+1)
        Y_h = y((1+h):(T+h));                % y_{t+h}
        eps_t = eps(2:T+1);                    % ε_t
        y_lag = y(1:T);                      % y_{t-1}
        X = [ones(T,1), eps_t, y_lag];       % Add constant + controls
        b_lp = pinv(X)*Y_h;
        IRF_LP_control(h+1, s) = b_lp(2)*shock_eps;    % β_h on ε_t
    end
end

% Median IRF
median_LP_control = median(IRF_LP_control, 2);

% --- Compute true impulse response ---
horizon = 21; 
true_irf = zeros(horizon,1);
deltaM = zeros(horizon,1);
eps_irf = zeros(horizon,1); eps_irf(1) = shock_eps; % for Jorda shock

for t = 1:horizon
    if t == 1
        deltaM(t) = eps_irf(t);
        true_irf(t) = c * eps_irf(t);
    else
        deltaM(t) = rho * deltaM(t-1) + eps_irf(t);
        true_irf(t) = a * true_irf(t-1) + b * deltaM(t-1) + c * eps_irf(t);
    end
end


% Plot comparison with true IRF
figure; hold on;
plot(0:H, true_irf(1:H+1), 'k--', 'LineWidth', 3);
plot(0:H, median_LP_control(2:end), 'g-', 'LineWidth', 2);
legend('True IRF', 'LP with y_{t-1} control', 'Location', 'NorthEast');
xlabel('Horizon (h)');
ylabel('Impulse Response of y_{t+h}');
title('Jordà LP with Control vs. True IRF');
grid on;

% Save figure
saveas(gcf, 'out/5_IRF_Jorda_control_y.png');
