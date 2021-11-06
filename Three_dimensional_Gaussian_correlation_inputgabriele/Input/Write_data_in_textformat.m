clear all;
clc;
close all;
cd('/home/prateek/DATA_CD/codes/final_PL/Wall_pressure_calculation_final_code/Three_dimensional_Gaussian_correlation_inputgabriele/Input');
load('PL_input.mat');
Data = [ypos/delta99 velgrad ufluc (1.4*4/sqrt(pi).*Lambda)/delta99]  ; % 1.4*4/np.sqrt(np.pi) taken by gabriele
dlmwrite('input_gabriele.txt', Data,'precision',16); % sixteen digit precision corresponds to double precision. Since 
% fortran code is double precision save the data in double precision mode. 

