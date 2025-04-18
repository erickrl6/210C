clear all

% Parameters from previous solutions
a = 0; %0.701502; % 1.43991; 
b = 0.615; %0.262662; % -1.70648; 
c = 0.492; %0.885;
d = 0.115;  % Not used here since eta_t = 0
rho = 0.8;
T = 600; % I will burn the first 100
N_sim = 500;

rng(456,'twister');
%rng(123456,'twister');

% Standard deviations (set so that Var(ΔM) = Var(η))
sigma_eta = 0.007;
sigma_eps = 0.00066;

% Storage for impulse responses with ε_t and 6 lags of ΔM
horizon = 21;
IRF_dM6 = zeros(horizon, N_sim);

for s = 1:N_sim
    % Generate shocks
    eta = sigma_eta * randn(T+6, 1);  % +6 for lag space
    % eta = zeros(T+6,1);
    eps = sigma_eps * randn(T+6, 1);

    % Simulate ΔM, y
    deltaM = zeros(T+6,1);
    y = zeros(T+6,1);
    for t = 2:T+6
        deltaM(t) = rho * deltaM(t-1) + eps(t);
        y(t) = a * y(t-1) + b * deltaM(t-1) + c * eps(t) + d * eta(t);
    end

    % Build regression dataset
    Y = y(107:end);
    Y_lag = y(106:end-1);
    E = eps(107:end);  % contemporaneous ε_t

    % ΔM_{t-1} to ΔM_{t-6}
    deltaM_lags = zeros(T-100, 6);
    for j = 1:6
        deltaM_lags(:, j) = deltaM(107-j:end-j);  % deltaM(t-j)
    end


    % Stack regressor matrix
    X = [ones(T-100,1), Y_lag, deltaM_lags, E];
    b_est = pinv(X)*Y;

    % Compute impulse response
    y_irf = zeros(horizon, 1);
    eps_irf = zeros(horizon, 1); eps_irf(1) = 1;%0.00066;
    deltaM_irf = zeros(horizon,1); deltaM_irf(1) = eps_irf(1); % initial shock just contemporary
    y_irf(1) = b_est(1) + b_est(end)*eps_irf(1); % initial response: alpha + gamma*shock

    for t = 2:horizon % get deltaM series given an initial contempraneous shock of 1 sd
        deltaM_irf(t) = rho * deltaM_irf(t-1) + eps_irf(t);
    end

    for t = 2:horizon
        y_irf(t) = b_est(1) + b_est(2) * y_irf(t-1) + b_est(end) * eps_irf(t);  % lag of y and contemp monetary shock
        for j = 1:min(6, t-1)
            y_irf(t) = y_irf(t) + b_est(2+j) * deltaM_irf(t-j);
        end
    end
    IRF_dM6(:, s) = y_irf;
end

% Median IRF
median_dM6 = median(IRF_dM6, 2);

% --- Compute true impulse response ---
true_irf = zeros(horizon,1);
deltaM = zeros(horizon,1);
eps_irf = zeros(horizon,1); eps_irf(1) = 1;%0.00066;

for t = 1:horizon
    if t == 1
        deltaM(t) = eps_irf(t);
        true_irf(t) = c * eps_irf(t);
    else
        deltaM(t) = rho * deltaM(t-1) + eps_irf(t);
        true_irf(t) = a * true_irf(t-1) + b * deltaM(t-1) + c * eps_irf(t);
    end
end

% Plot comparison
figure; hold on;
plot(0:horizon-1, true_irf, 'k--', 'LineWidth', 3);
plot(0:horizon-1, median_dM6, '-', 'LineWidth', 2);
legend('True IRF', 'Model: y_{t-1} + ε_t + 6 lags of  ΔM', 'Location', 'NorthEast');
xlabel('Time (t)');
ylabel('Impulse Response of y_t');
title('Misspecified Model: y_{t-1} + ε_t + 6 lags of ΔM');
grid on;

% Save figure
saveas(gcf, 'out/4_IRF_monetary_y.png');
