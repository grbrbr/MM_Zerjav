
close all
clc


% Cas stopnice
T_stop = 0.6;
% Amplituda
A = 0.1;

% Vzbujanje 
t = 0:0.01:15;
u = rectangularPulse(0,T_stop,t) * A;
ts = timeseries(u,t);

simIn = Simulink.SimulationInput("model_func");
simIn = simIn.setVariable("F_vzbujanje", ts);
out = sim(simIn);

