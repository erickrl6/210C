clear all 

% Parameters from previous solutions
a = 0; %0.701502; % 1.43991; 
b = 0.49232; %0.262662; % -1.70648; 
c = 0.6154; %0.885;
d = 0.115;  % Not used here since eta_t = 0
rho = 0.8;
T = 600; % I will burn the first 100
N_sim = 500;

%rng(12345,'twister');  % for reproducibility - seed
%rng(123456,'twister'); 
rng(456,'twister');

% Standard deviations (set so that Var(ΔM) = Var(η))
sigma_eta = 0.007;
sigma_eps = 0.00066;

% IR storage: [horizon x simulations]
horizon = 21;
IRF_1lag = zeros(horizon, N_sim);
IRF_4lag = zeros(horizon, N_sim);
IRF_12lag = zeros(horizon, N_sim);

% Loop over simulations
for s = 1:N_sim
    %rng(456+s,'twister');
    % Generate shocks
    eta = sigma_eta * randn(T+1, 1); % productivity shocks
    %eta = zeros(T+1,1);
    eps = sigma_eps * randn(T+1, 1);

    % Simulate ΔM, y
    deltaM = zeros(T+1,1);
    y = zeros(T+1,1);

    for t = 2:T+1
        deltaM(t) = rho * deltaM(t-1) + eps(t);
        y(t) = a * y(t-1) + b * deltaM(t-1) + c * eps(t) + d * eta(t);
    end

    % Trim initial obs + burn in
    Y = y(102:end);
    E = eps(102:end);
    
    %% Misspecified model: 1 lag of y
    X1 = [ones(T-100,1), y(101:end-1), E];
    b1 = pinv(X1)*Y;
    % Impulse response: simulate y_t after ε_0 = 1 standard dev
    y1 = zeros(horizon,1);
    eps_irf = zeros(horizon,1); eps_irf(1) = 0.00066; % 0.00066;
    y1(1) = b1(1) + b1(end)*eps_irf(1); % initial response: alpha + gamma*shock
    for t = 2:horizon
        y1(t) = b1(1) + b1(2) * y1(t-1) + b1(3) * eps_irf(t);
    end
    IRF_1lag(:, s) = y1;

    %% Misspecified model: 4 lags of y
    Y4lags = lagmatrix(y(102:end), 1:4);
    X4 = [ones(T-104,1), Y4lags(5:end,:), E(5:end)];
    Y4 = Y(5:end);
    b4 = pinv(X4)*Y4;

    y4 = zeros(horizon,1);
    y4(1) = b4(1) + b4(end)*eps_irf(1); % initial response: alpha + gamma*shock

    for t = 2:horizon
        lag_sum = 0;
        % Use up to 4 lags (or fewer if t-1 < 4)
        for j = 1:min(4, t-1)
            lag_sum = lag_sum + b4(1+j)*y4(t-j);
        end
        y4(t) = b4(1) + lag_sum + b4(end) * eps_irf(t);
    end

    IRF_4lag(:, s) = y4;

    %% Misspecified model: 12 lags of y
    Y12lags = lagmatrix(y(102:end), 1:12);
    X12 = [ones(T-112,1), Y12lags(13:end,:), E(13:end)];
    Y12 = Y(13:end);
    b12 = pinv(X12) * Y12;

    y12 = zeros(horizon,1);
    y12(1) = b12(1) + b12(end)*eps_irf(1); % initial response: alpha + gamma*shock

    for t = 2:horizon
        lag_sum = 0;
        % Use up to 12 lags (or fewer if t-1 < 12)
        for j = 1:min(12, t-1)
            lag_sum = lag_sum + b12(1+j)*y12(t-j);
        end
        y12(t) = b12(1) + lag_sum + b12(end) * eps_irf(t);
    end

    IRF_12lag(:, s) = y12;
end

%% Plot median IRFs
median_1 = median(IRF_1lag, 2);
median_4 = median(IRF_4lag, 2);
median_12 = median(IRF_12lag, 2);

% --- Compute true impulse response ---
true_irf = zeros(horizon,1);
deltaM = zeros(horizon,1);
eps_irf = zeros(horizon,1); eps_irf(1) = 0.00066; %0.00066;

for t = 1:horizon
    if t == 1
        deltaM(t) = eps_irf(t);
        true_irf(t) = c * eps_irf(t);
    else
        deltaM(t) = rho * deltaM(t-1) + eps_irf(t);
        true_irf(t) = a * true_irf(t-1) + b * deltaM(t-1) + c * eps_irf(t);
    end
end

figure; hold on;
plot(0:horizon-1, true_irf, 'k--', 'LineWidth', 2);
plot(0:horizon-1, median_1, 'r', 'LineWidth', 2);
plot(0:horizon-1, median_4, 'g', 'LineWidth', 2);
plot(0:horizon-1, median_12, 'b', 'LineWidth', 2);
legend('True IRF','1 lag', '4 lags', '12 lags','Location','north');
xlabel('Time (t)');
ylabel('Impulse Response of y_t');
title('Median Impulse Response from Misspecified Models','FontSize',9);
grid on;

% Save figure
saveas(gcf, 'out/3_IRF_misspec_y.png');