close all
clear
clc

LUT_U_to_F = 2.25;
x_pulzi_na_meter = 22755;
fi_pulzi_na_stopinjo = 45.5;

%----------- Vzbujanje -------------
% Define input force (example: sinusoidal input)
Tsimulation = 15;
time = 0:0.01:Tsimulation;  % Simulation time vector (s)

% Pulse
Astep = -3;        % Amplitude of the pulse (U)
Tstep = 1.2;        % Duration of the pulse (s)
U = zeros(size(time)); % Initialize force vector
U(time <= Tstep) = Astep; % Set force to Astep for the duration of Tstep

% Convert input voltage to force s
F = LUT_U_to_F .* U;

%----------- Meritve -------------
meritve = load('meritve 6_12/plus_A3_T1p2_t15.mat');
% Access specific variables
t_meritve = meritve.out.simout.time;
fi_meritve = meritve.out.simout.data(:,1);
x_meritve = meritve.out.simout.data(:,2);

u_meritve = meritve.u;
t_vzb_meritve = meritve.t;

% Convert x output from pulses to meters or degrees
x_meritve = x_meritve / x_pulzi_na_meter;
fi_meritve = fi_meritve / fi_pulzi_na_stopinjo;
% Convert angle to rad
fi_meritve = fi_meritve * pi/180;

% Find start x position and substract from meritve
x_offset = x_meritve(1);
x_meritve = x_meritve - x_offset;

%----------- Parametri modela -------------

% Define parameters
m1 = 4;            % Mass of the cart (kg)
m2 = 0.36;         % Mass of the pendulum (kg)
l = 0.451;         % Length of the pendulum (m)
f1 = 10;           % Damping constant for the cart (kg/s)
f2 = 0.00145;      % Damping constant for the pendulum (kg·m^2/s)
J = 0.08433;       % Pendulum inertia (kg·m^2)
g = 9.81;          % Gravitational acceleration (m/s^2)

% Initialize state variables
x = zeros(size(time));       % Position of the cart (m)
x_dot = zeros(size(time));   % Velocity of the cart (m/s)
x_ddot = zeros(size(time));  % Acceleration of the cart (m/s^2)
phi = zeros(size(time));     % Pendulum angle (rad)
phi_dot = zeros(size(time)); % Angular velocity of the pendulum (rad/s)
phi_ddot = zeros(size(time));% Angular acceleration of the pendulum (rad/s^2)

% Precompute constants
c1 = 1 / (m1 + m2);
c2 = 1 / (J + (m2^2 * l^2) / (m1 + m2));

%----------- Izracun modela -------------
% Simulation loop
for i = 1:length(time)-1
    % Compute cart acceleration
    x_ddot(i) = c1 * (F(i) - f1 * x_dot(i) + (m2 * l / J) * (f2 * phi_dot(i) + m2 * l * g * phi(i)));

    % Compute pendulum angular acceleration
    phi_ddot(i) = -c2 * (f2 * phi_dot(i) + m2 * l * g * phi(i) - (m2 * l / (m1 + m2)) * F(i) + (m2 * l * f1 / (m1 + m2)) * x_dot(i));

    % Update state variables using Euler integration
    x_dot(i+1) = x_dot(i) + x_ddot(i) * 0.01;
    x(i+1) = x(i) + x_dot(i) * 0.01;
    phi_dot(i+1) = phi_dot(i) + phi_ddot(i) * 0.01;
    phi(i+1) = phi(i) + phi_dot(i) * 0.01;
end


%----------- Plotting -------------

% Create a new figure
figure;

% Plot Input Force F
subplot(2,1,1);
plot(time, F, 'b', 'LineWidth', 1.5); % Plot input force
title('Input Force (F)');
xlabel('Time (s)');
ylabel('Force (N)');
grid on;

% Plot Input Voltage U and Measured Voltage u_meritve
subplot(2,1,2);
plot(time, U, 'r', 'LineWidth', 1.5); % Plot calculated voltage U
hold on;
plot(t_vzb_meritve, u_meritve, 'k--', 'LineWidth', 1.5); % Plot measured voltage u_meritve
title('Input Voltage (U) and Measured Voltage (u_{meritve})');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Calculated Voltage (U)', 'Measured Voltage (u_{meritve})');
grid on;


% Plot Cart Position
figure;
subplot(1,2,1);
plot(time, x, 'r', 'LineWidth', 1.5); % Plot calculated cart position
hold on; % Hold on to add measured data
plot(t_meritve, x_meritve, 'k--', 'LineWidth', 1.5); % Plot measured cart position
title('Cart Position (x)');
xlabel('Time (s)'); 
ylabel('Position (m)');
legend('Calculated Position', 'Measured Position');
grid on;

% Plot Pendulum Angle
subplot(1,2,2);
plot(time, phi, 'g', 'LineWidth', 1.5); % Plot calculated pendulum angle
hold on; % Hold on to add measured data
plot(t_meritve, fi_meritve, 'k--', 'LineWidth', 1.5); % Plot measured pendulum angle
title('Pendulum Angle (phi)');
xlabel('Time (s)'); 
ylabel('Angle (rad)');
legend('Calculated Angle', 'Measured Angle');
grid on;