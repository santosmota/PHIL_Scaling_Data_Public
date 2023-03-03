% Magnitude invariant PU normalisation - Iterative algorithm to arrive at operating point of 60 kVA scaled down converter


clc;
clear all;
%Laboratory Scaled down Converter component ratings

%LCL filter 
Lr_sd_lab = 500e-6; %Main reactor
Lt_sd_lab = 2.38e-4; % Transformer leakage inductance viewed from LV side
Rt_sd_lab = 4.94e-2; % Transformer leakage resistance viewed from LV side
Cac_sd_lab = 5e-5;  %Shunt branch capacitance
Rac_sd_lab = 0.33577; % Shunt branch capacitance resistance
%DC link
Cdc_sd_lab = 1.40e-02; % DC link Capacitance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Full Scale BESS Ratings and PU calculations
Sb_fs = 5e6; % BESS rating
Vac_fs = 690 ; %LV Voltage rating L-L
Iac_fs = Sb_fs/(sqrt(3)*Vac_fs); % Converter current LV side
%DC side
Vdc_fs = 1000; 
Pf_fs = 1;
Idc_fs = (Sb_fs*Pf_fs)/(Vdc_fs);
%Transformer
Transfo_ratio = 13800/690;
Xt_fs_pu = 0.08;
Rt_fs_pu = 0.005;
Zs_fs_base =  Vac_fs/(sqrt(3)*Iac_fs);
freq_fs_base = 60;
Ls_fs_base = Zs_fs_base/(2*pi()*freq_fs_base);
Cs_fs_base = 1/((2*pi()*freq_fs_base*Zs_fs_base));
Xt_fs = Xt_fs_pu*Zs_fs_base;
Lt_fs = Xt_fs/(2*pi()*freq_fs_base);
Lt_fs_pu = Lt_fs/Ls_fs_base;
Rt_fs = Rt_fs_pu*Zs_fs_base;
Fsw_fs = 3000; 
%LCL filter
Lr_fs = 1.06e-04;
Cac_fs = 1.53e-03;
Rac_fs = 3.51e-02;
Lr_fs_pu = Lr_fs/Ls_fs_base;
Cac_fs_pu = Cac_fs/Cs_fs_base;
Rac_fs_pu = Rac_fs/Zs_fs_base;
LCL_fs_resonance = sqrt((Lr_fs+Lt_fs)/(Lr_fs*Lt_fs*Cac_fs))/(2*pi());
del_Ilr_fs = Vdc_fs/(8*Lr_fs*Fsw_fs);
del_Ilr_fs_pu = del_Ilr_fs/(sqrt(2)*Iac_fs);
%DC Capacitance
Cdc_fs = 1.85e-2;
Cdc_fs_pu = Cdc_fs/Cs_fs_base;
H_fs = (Cdc_fs*Vdc_fs^2)/(2*Sb_fs);

%Magnitude Invariant Method Calculations for scaled down converter
for(k = 50:1:400)
 Vac_sd(k) = k; 
for(p = 5:1:72)
 Iac_sd(p) = p;  
 freq_sd_base = 60; %Might try to change to 50Hz
 Fsw_sd = 3000; 
 Vdc_fs_to_Vac_fs = Vdc_fs/Vac_fs;
 Vdc_sd(k) = Vdc_fs_to_Vac_fs*Vac_sd(k);
 Sb_sd(k,p) = sqrt(3)*Vac_sd(k)*Iac_sd(p);
 Zs_sd_base(k,p) = Vac_sd(k)/(sqrt(3)*Iac_sd(p));
 Ls_sd_base(k,p) = Zs_sd_base(k,p)/(2*pi()*freq_sd_base);
 Cs_sd_base(k,p) = 1/((2*pi()*freq_sd_base*Zs_sd_base(k,p)));
 %LCL filter - Main reactor
 Lr_sd_calcu(k,p) = (Lr_fs/Ls_fs_base)*Ls_sd_base(k,p);
 Lr_sd_pu(k,p) = Lr_sd_lab/Ls_sd_base(k,p);
 del_Ilr_sd(k,p) = Vdc_sd(k)/(8*Lr_sd_lab*Fsw_sd);
 del_Ilr_sd_pu(k,p) = del_Ilr_sd(k,p)/(sqrt(2)*Iac_sd(p));
 %Transformer
  Lt_sd_calcu(k,p) = (Lt_fs/Ls_fs_base)*Ls_sd_base(k,p);
  Lt_sd_pu(k,p) = Lt_sd_lab/Ls_sd_base(k,p);
  Rt_sd_cal(k,p) = (Rt_fs/Zs_fs_base)*Zs_sd_base(k,p);
  Rt_sd_pu(k,p) = Rt_sd_lab/Zs_sd_base(k,p);
 %LCL filter - Shunt branch
  Cac_sd_calcu(k,p) = (Cac_fs/Cs_fs_base)*Cs_sd_base(k,p);
  Cac_sd_pu(k,p) = Cac_sd_lab/Cs_sd_base(k,p);
  LCL_sd_resonance(k,p) = sqrt((Lr_sd_lab+Lt_sd_lab)/(Lr_sd_lab*Lt_sd_lab*Cac_sd_lab))/(2*pi()); 
  Rac_sd_calcu(k,p) = (Rac_fs/Zs_fs_base)*Zs_sd_base(k,p);
  Rac_sd_pu(k,p) = Rac_sd_lab/Zs_sd_base(k,p);
  %DC link capacitance
  Cdc_sd_calcu(k,p) = Cdc_fs*(Vdc_fs/Vdc_sd(k))^2*(Sb_sd(k,p)/Sb_fs);
  Cdc_sd_pu(k,p) = Cdc_sd_lab/Cs_sd_base(k,p);
  H_sd(k,p) = (Cdc_sd_lab*Vdc_sd(k)^2)/(2*Sb_sd(k,p));
 %Pu mismatch errors over Full scale converter
 %loss less components
  err_Lr(k,p) = abs((Lr_sd_pu(k,p) - Lr_fs_pu)/Lr_fs_pu)*100;
  err_Lt(k,p) = abs((Lt_sd_pu(k,p) - Lt_fs_pu)/Lt_fs_pu)*100;
  err_Cac(k,p) = abs((Cac_sd_pu(k,p) - Cac_fs_pu)/Cac_fs_pu)*100;
  err_delIlr(k,p) = abs((del_Ilr_sd_pu(k,p)-del_Ilr_fs_pu)/del_Ilr_fs_pu)*100;
  err_H(k,p) = abs((H_sd(k,p)-H_fs)/H_fs)*100;

 %lossy components
  err_Rt(k,p) = abs((Rt_sd_pu(k,p) - Rt_fs_pu)/Rt_fs_pu)*100;
  err_Rac(k,p)  = abs((Rac_sd_pu(k,p) - Rac_fs_pu)/Rac_fs_pu)*100;
end
end
limit = 0.001;
%Finding minimal values
err_Lr(err_Lr<limit) = NaN;
err_Lt(err_Lt<limit) = NaN;
err_Cac(err_Cac<limit) = NaN;
err_delIlr(err_delIlr<limit) = NaN;
err_H(err_H<limit) = NaN;
err_Rt(err_Rt<limit) = NaN;
err_Rac(err_Rac<limit) = NaN;

%At this point focussing to minimise the error for only four points Lr, Lt,Cac and ripple current
minValue = min(err_Lr(:));
[row, column] = find(err_Lr == minValue);
X = ['minimum value of error for Lr_pu in %'];
disp(X)
minValue
Voltage = row
Current = column

minValue = min(err_Lt(:));
[row, column] = find(err_Lt == minValue);
X = ['minimum value of error for Lt_pu in %'];
disp(X)
minValue
Voltage = row
Current = column

minValue = min(err_Cac(:));
[row, column] = find(err_Cac == minValue);
X = ['minimum value of error for shunt cap Cac_pu in %'];
disp(X)
minValue
Voltage = row
Current = column

minValue = min(err_delIlr(:));
[row, column] = find(err_delIlr == minValue);
X = ['minimum value of error for ripple current_pu in %'];
disp(X)
minValue
Voltage = row
Current = column

minValue = min(err_H(:));
[row, column] = find(err_H == minValue);
X = ['minimum value of error for H_pu in %'];
disp(X)
minValue
Voltage = row
Current = column






%plots
surf(err_Lr)
hold on
%surf(err_Lt)
hold on
surf(err_Cac)
hold on
surf(err_delIlr)

legend('errLr','errLt','errCac','errdelIlr')
%End of script






































%Backup
%Vac_sd = 150 ; % Volts - this is a range that needs to be iterated
% Iac_sd = 10 ; % Amps - this is a range that needs to be iterated
% freq_sd_base = 60; 
% Fsw_sd = 2500; 
% Vdc_fs_to_Vac_fs = Vdc_fs/Vac_fs;
% Vdc_sd = Vdc_fs_to_Vac_fs*Vac_sd;
% Sb_sd = sqrt(3)*Vac_sd*Iac_sd;
% Zs_sd_base = Vac_sd/(sqrt(3)*Iac_sd);
% Ls_sd_base = Zs_sd_base/(2*pi()*freq_sd_base);
% Cs_sd_base = 1/((2*pi()*freq_sd_base*Zs_sd_base));
% %LCL filter - Main reactor
% Lr_sd_calcu = (Lr_fs/Ls_fs_base)*Ls_sd_base;
% Lr_sd_pu = Lr_sd_lab/Ls_sd_base;
% del_Ilr_sd = Vdc_sd/(8*Lr_sd_lab*Fsw_sd);
% del_Ilr_sd_pu = del_Ilr_sd/(sqrt(2)*Iac_sd);
% %Transformer
% Lt_sd_calcu = (Lt_fs/Ls_fs_base)*Ls_sd_base;
% Lt_sd_pu = Lt_sd_lab/Ls_sd_base;
% Rt_sd_cal = (Rt_fs/Zs_fs_base)*Zs_sd_base;
% Rt_sd_pu = Rt_sd_lab/Zs_sd_base;
% %LCL filter - Shunt branch
% Cac_sd_calcu = (Cac_fs/Cs_fs_base)*Cs_sd_base;
% Cac_sd_pu = Cac_sd_lab/Cs_sd_base;
% LCL_sd_resonance = sqrt((Lr_sd_lab+Lt_sd_lab)/(Lr_sd_lab*Lt_sd_lab*Cac_sd_lab))/(2*pi());
% Rac_sd_calcu = (Rac_fs/Zs_fs_base)*Zs_sd_base;
% Rac_sd_pu = Rac_sd_lab/Zs_sd_base;
% % DC link capacitance
% Cdc_sd_calcu = Cdc_fs*(Vdc_fs/Vdc_sd)^2*(Sb_sd/Sb_fs);
% Cdc_sd_pu = Cdc_sd_lab/Cs_sd_base;
% H_sd = (Cdc_sd_lab*Vdc_sd^2)/(2*Sb_sd);
% 
% %Pu mismatch errors over Full scale converter
% %loss less components
% err_Lr = ((Lr_sd_pu - Lr_fs_pu)/Lr_fs_pu)*100
% err_Lt = ((Lt_sd_pu - Lt_fs_pu)/Lt_fs_pu)*100
% err_Cac = ((Cac_sd_pu - Cac_fs_pu)/Cac_fs_pu)*100
% err_delIlr = ((del_Ilr_sd_pu-del_Ilr_fs_pu)/del_Ilr_fs_pu)*100
% err_H = ((H_sd-H_fs)/H_fs)*100
% %lossy components
% err_Rt = ((Rt_sd_pu - Rt_fs_pu)/Rt_fs_pu)*100
% err_Rac = ((Rac_sd_pu - Rac_fs_pu)/Rac_fs_pu)*100




































