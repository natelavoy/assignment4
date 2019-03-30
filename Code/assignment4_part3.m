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
cap(3,4,0.00001);
cur(3,4,1);
ccvs2(5,0,4,0,100);
res(5,6,0.1);
res(6,0,1000);
disp('The new C Matrix including nominal 1 values for the current source is:')
disp(C)
%% Frequency Response
F = logspace(0, 9, 5000);
OutputNode = 6; 
noise = randn(1,5000);
noise = noise./max(abs(noise)).*0.01;
for n=1:length(F)
    b(3) = noise(n);
    b(4) = -noise(n);
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
figure(4);
plot(F,noise)
title('In')
xlabel('Time (s)');
ylabel('Current (A)');

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

function [Vout,ffts1, ffts2] = transient(G,C,b,input,time,h)
    xold = zeros(length(b),1);
    Vout = zeros(numel(time),1);
    bnew = zeros(length(b),1);
    
    for n=2:numel(time)
        noise = randn();
        noise = noise./max(abs(noise)).*0.01;
        b(3) = noise;
        b(4) = -noise;   
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