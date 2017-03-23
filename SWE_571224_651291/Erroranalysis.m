function Erroranalysis()
clear all; close all; clc;
x=-1:0.01:1;
Amp_error=sqrt((1-(x.^2)/2+(x.^4)/24).^2+(x-(x.^3)/6).^2);
Phase_error=atan((x-(x.^3)/6)./(1-(x.^2)/2+(x.^4)/24))-x;
hold on 
figure(1) 
plot(x,Amp_error); 
xlabel('Lamda*Delta t'); 
ylabel('Amplitude Error'); 
grid on 
title('Amplitude Error Plot'); 
figure(2) 
plot(x,Phase_error); 
xlabel('Lamda*Delta t'); 
ylabel('Phase Error'); 
grid on 
title('Phase Error Plot'); 
