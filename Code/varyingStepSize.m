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
%% Transient Analysis
Inputs;
Fs = 1000;
T = 1/Fs;
L = 1000;
t = (0:L-1)*T;
[stepOut,stepinF,stepOutF] = transient(G,C,b,stepin,t,T);
[sinOut,sininF,sinOutF] = transient(G,C,b,sinin,t,T);
[gaussOut,gaussinF,gaussOutF] = transient(G,C,b,gaussin,t,T);

[stepOut2,stepinF2,stepOutF2] = transient(G,C,b,stepin,t,T/4);
[sinOut2,sininF2,sinOutF2] = transient(G,C,b,sinin,t,T/4);
[gaussOut2,gaussinF2,gaussOutF2] = transient(G,C,b,gaussin,t,T/4);

% FFT of step
figure(5);
subplot(2,1,1)
plot(t,stepin,t,stepOut);
title('Step Function Time Domain Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n','V_O');
subplot(2,1,2)
plot(t,stepin,t,stepOut2);
title('Step Function Time Domain Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n','V_O');

% FFT of sine 
figure(6);
subplot(2,1,1)
plot(t,sinin,t,sinOut);
title('Sine Function Time Domain Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n','V_O');
subplot(2,1,2)
plot(t,sinin,t,sinOut2);
title('Sine Function Time Domain Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n','V_O');


% FFT of Gauss
figure(7);
subplot(2,1,1)
plot(tg,gaussin,t,gaussOut);
title('Gaussian Pulse Time Domain Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n','V_O');
subplot(2,1,2)
plot(tg,gaussin,t,gaussOut2);
title('Gaussian Pulse Time Domain Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_i_n','V_O');

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
%%
% The results show how the larger step size results in a less accurate
% result with a slower response to change in the circuit.