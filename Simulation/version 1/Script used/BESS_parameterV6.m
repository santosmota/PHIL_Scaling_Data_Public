%Detailed BESS 2MW parameters
%Prepared by Joseph Banda, PhD at NTNU
clc;
%Simulation time steps
%Tsample = 1/24000;
%Tplant = Tsample*0.01;
Tplant = 10e-6;
Tdelay = 1e-6;

%Time information from Template
Ts_sim = 10E-6;
Ts_phy = Ts_sim;
TsFPGA = 10E-9;
k_Ts = 1; %ratio between simulation time step and control time step
Ts_c = k_Ts * Ts_phy; % control time step [s]
FPGA_clock_freq = 1E8; % Update rate for FPGA. ( used by FPGA signal generator, filters, integrators etc.)

%% run script for configuring electrical simulation model
ABB_BESS_50_initial;
