
encoder_reset

% Cas stopnice
T_stop = 0.6;
% Amplituda
A = 0;

% Vzbujanje 
t = 0:0.01:15;
u = rectangularPulse(0,T_stop,t) * A;
ts = timeseries(u,t);

simIn = Simulink.SimulationInput("moedl1");
simIn = simIn.setVariable("vzbujanje", ts);
out = sim(simIn);


simIn1 = Simulink.SimulationInput("model");
out1 = sim(simIn);