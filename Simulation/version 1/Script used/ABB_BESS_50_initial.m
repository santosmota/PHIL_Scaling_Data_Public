%%%% NTNU - BEGIN %%%%

%% BESS-STATCOM System
systemBase = {};
systemBase.Sn = 5e6; % VA
systemBase.Vn = 13800;  % V rms phase to phase
systemBase.fn = 50; % Hz
systemBase.Vn_phase = systemBase.Vn/sqrt(3); % V
systemBase.In   = systemBase.Sn/(systemBase.Vn*sqrt(3)) ; % A
systemBase.Zn   = systemBase.Vn*systemBase.Vn/systemBase.Sn ; % Ohm
systemBase.Rn   = systemBase.Zn; % Ohm
systemBase.Wn   = 2*pi()*systemBase.fn ; % rad/s
systemBase.Ln   = systemBase.Zn/systemBase.Wn ; % H
systemBase.Cn   = 1/(systemBase.Zn*systemBase.Wn) ; % F

%% Grid parameters
grid.rss = 0.02;    % Maximum steady-state frequency deviation [pu]
grid.rtr = 0.05;    % Maximum transient frequency deviation [pu]

%% ESS - Grid converter specification
ess = {};
% Sampling time, system frequency, switching frequency
ess.Fn = systemBase.fn;
ess.Fs_control = 1/Ts_phy;          % [Hz] Sampling frequency controller
ess.Fsw = 3e3;                      % [Hz] Switching frequency PWM
ess.Ts_control = 1/ess.Fs_control;  % [s] Sampling time controller A/D converter
ess.Tsw = 1/ess.Fsw;                % [s] Switching period

% Converter Specifications
ess.Pn = 3e6;                       % [W] Active power
ess.Sn = 5e6;                       % [VA] Apparent power
ess.Un = 690;                       % [Vrms] Line voltage (at LV side of trafo)
ess.Uanp = ess.Un*sqrt(2/3);        % [Vpeak] Phase to ground voltage (at LV side of trafo)
ess.Uabp = ess.Un*sqrt(2);          % [Vpeak] Line voltage (at LV side of trafo)
ess.Iap = ess.Sn/ess.Un*sqrt(2/3);  % [Apeak] Line current (at LV side of trafo)
ess.Udc = 1100;                     % [Vdc] Rated dc voltage
ess.Idc = ess.Sn/ess.Udc;           % [Adc] Rated dc current
%ess.Cdc = ess.Idc/(2*ess.Fn*ess.Udc); % [F] DC link capacitor % JB: Not sure, check with Daniel Mota
ess.Cdc = 20e-3;     %reduced capacitance - Yet to check the design
ess.LCL.l2 = 0.08;                  % [pu] Short circuit reactance of trafo
ess.LCL.r2 = 0.005;                 % [pu] Short circuit resistance of trafo
ess.LCL.r1 = 0.01;                  % [pu] LCL main reactance resistance
ess.dampDC = 0.8; %sqrt(2)/2;       % [-] damping Udc controller 

%% ESS - Grid converter design
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% ESS - Grid converter');
[designOk, ess.LCL] = GC_LCLdesign(ess);
if designOk == 1 
    [ess.CurrContSingle, ess.DCVContSingle] = GC_PItuning(ess,0);
    [ess.CurrContDual, ess.DCVContDual] = GC_PItuning(ess,1);
end

% ESS - Grid converter controllers
ess.filtertype = 1;                 % None = 0; Notch = 1; DSC = 2

%Current 
ess.CurrCont.PImax = 1.5;           % [pu] max value of PI output
ess.CurrCont.PImin = -1.5;          % [pu] min value of PI output
ess.CurrCont.x = ess.LCL.l1;        % [pu] reactance for dq decoupling
ess.CurrCont.f_LPF_noise = ess.LCL.fres;    % [Hz] low pass cutout filter
ess.CurrCont.zeta = sqrt(2)/2;      % [pu] damping coefficient for notch filters

%DC link voltage
ess.DCVCont.PImax = 1.5;            % [pu] Max value of PI output
ess.DCVCont.PImin = -1.5;           % [pu] Min value of PI output
ess.DCVCont.Start = 0.1;

%AC voltage controller
ess.ACVCont.PImax = 1.5;            % [pu] Max value of PI output
ess.ACVCont.PImin = -1.5;           % [pu] Min value of PI output
ess.ACVCont.kg = 1;                 % [pu/pu] Loop gain kg(kp + 1/sTi)
ess.ACVCont.kp = 0.00;              % [pu/pu] Proportional gain kg(kp + 1/sTi)
ess.ACVCont.Ti = 0.15;               % [s] Integral time
ess.ACVCont.kiq_droop = 0.00;

%%% The end
