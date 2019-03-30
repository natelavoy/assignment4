%% Assignment 4
% Nathan Lavoy
% 100995612
% Submitted: March 30, 2019
clc; close all;  clear all;  %initialization of the matlab environment

global G C b; %define global variables

G = zeros(6,6); 
C = zeros(6,6); 
b = zeros(6,1); 

vol(1,0,1);
res(1,2,1);
cap(1,2,0.25);
res(2,0,2);
ind(2,3,0.2);
res(3,4,10);
ccvs2(5,0,4,0,100);
res(5,6,0.1);
res(6,0,1000);
disp('The G matix is:')
disp(G)
disp('The C matix is:')
disp(C)

steps = 50;
Vin = linspace(-10,10, 50);
for n=1:steps
    b(7) = Vin(n);
    X = G\b;
    V3(n) = X(3);
    VO(n) = X(6);
    gain(n) = VO(n)/Vin(n);
end
figure(1);
plot(Vin,V3, Vin, VO);
title('DC Sweep')
xlabel('Input Voltage (V)')
ylabel('Ouput Voltage (V)')
legend(['V_3';'V_O'])
%% Frequency Response
F = logspace(0, 9, 5000);
OutputNode = 6;  
b(7) = 1;
for n=1:length(F)
    w = 2*pi*F(n);
    s = 1i*F(n);
    A = G + s*C;   

    X = A\b;
    Vout(n) = abs(X(OutputNode));
    gain(n) = 20*log(abs(X(OutputNode)));
end

figure(2);
semilogx(F, Vout);
xlabel('Frequency (Hz)');
ylabel('V_O (V)');
title('Frequncy Response');
figure(3);
semilogx(F, gain);
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
title('Gain Response');
%% Part D - Monte Carlo
Cvals = 0.25 + 0.05.*randn(1,1000);
for n = 1:1000
    s = 1i*pi;
    C(1,1) = Cvals(n);
    C(2,2) = Cvals(n);
    C(1,2) = Cvals(n)*-1;
    C(2,1) = Cvals(n)*-1;
    
    A = G + s*C;   
    X = A\b;
    gainMC(n) = 20*log10(abs(X(OutputNode)));
end
figure(4)
histogram(gainMC);

%% Transient Analysis
Inputs;
Fs = 1000;
T = 1/Fs;
L = 1000;
t = (0:L-1)*T;
[stepOut,stepinF,stepOutF] = transient(G,C,b,stepin,t,T);
[sinOut,sininF,sinOutF] = transient(G,C,b,sinin,t,T);
[gaussOut,gaussinF,gaussOutF] = transient(G,C,b,gaussin,t,T);

% FFT of step
figure(5);
subplot(3,1,1)
plot(t,stepin,t,stepOut);
title('Step Function Time Domain Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n','V_O');
disp('Note: the infinite fourier transform is not shown due to the low sampling rate')
subplot(3,1,2)
plot(stepinF) 
title('Single-Sided Amplitude Spectrum of Step Function Input')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([450 550])
subplot(3,1,3)
plot(stepOutF,'red') 
title('Single-Sided Amplitude Spectrum of Step Function Output')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([450 550])

% FFT of sine 
figure(6);
subplot(3,1,1)
plot(t,sinin,t,sinOut);
title('Sine Function Time Domain Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n','V_O');
subplot(3,1,2)
plot(sininF) 
title('Single-Sided Amplitude Spectrum of Sine Function Input')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([300 700])
subplot(3,1,3)
plot(sinOutF,'red') 
title('Single-Sided Amplitude Spectrum of Sine Function Output')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([300 700])

% FFT of Gauss
figure(7);
subplot(3,1,1)
plot(tg,gaussin,t,gaussOut);
title('Gaussian Pulse Time Domain Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n','V_O');
subplot(3,1,2)
plot(gaussinF) 
title('Single-Sided Amplitude Spectrum of Gaussian Pulse Input')
xlabel('Frequency (f)')
ylabel('|P(f)|')
xlim([300 700])
subplot(3,1,3)
plot(gaussOutF) 
title('Single-Sided Amplitude Spectrum of Gaussian Pulse Output')
xlabel('f (Hz)')
ylabel('|P(f)|')
xlim([300 700])
%% Noise
assignment4_part3
%% Varying Bandwidth
VaryingBW
%% VaryingStepSize
varyingStepSize
%% Non-linear Gain
% If the response given is non-linear in the form V = aI3 + bI3^2 + cI3^3
% than the addition of a F(X) vector must be added. This will be solved by
% a Newton-Raphson method which would be implemented to solve the F(X)
% equation. It would then be used in the linear decomposition or
% trapezoidal rule. 
%% Tapezoidal Rule
function [Vout,ffts1, ffts2] = transient(G,C,b,input,time,h)
    xold = zeros(length(b),1);
    Vout = zeros(numel(time),1);
    bnew = zeros(length(b),1);
    
    for n=2:numel(time)
        bmat = b;
        bmat(7) = input(n-1);
        bnew(7) = input(n);
        A = (2*C/h + G);
        xnew = A\((2*C/h - G)*xold + bmat + bnew);
        Vout(n) = xnew(6);
        xold = xnew;
    end
    ffts1 = fftshift(abs(fft(input)));
    ffts2 = fftshift(abs(fft(Vout)));
end