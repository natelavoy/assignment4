clear; clc;
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
Cval = [0.00001,0.001,0.01];
for i = 1:3
    F = logspace(0, 3, 5000);
    OutputNode = 6; 
    noise = randn(1,5000);
    noise = noise./max(abs(noise)).*0.01;
    C(3,3) = Cval(i);
    C(4,4) = Cval(i);
    C(3,4) = -Cval(i);
    C(4,3) = -Cval(i);
    for n=1:length(F)
        b(3) = noise(n);
        b(4) = -noise(n);
        w = 2*pi*F(n);
        s = 1i*F(n);
        A = G + s*C;   

        X = A\b;
        Vout(i,n) = abs(X(OutputNode));
        gain(i,n) = 20*log(abs(X(OutputNode)));
    end
end
figure(1);
subplot(2,1,1)
semilogx(F, Vout(1,:), F, Vout(2,:), F, Vout(3,:));
xlabel('Frequency (Hz)');
ylabel('V_O (V)');
title('Frequncy Response');
legend('C = 10e-6','C = 1e-3','C = 0.01');
subplot(2,1,2)
semilogx(F, gain(1,:), F, gain(2,:), F, gain(3,:));
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
title('Gain Response');
legend('C = 10e-6','C = 1e-3','C = 0.01');
%%
% The figure shows that increasing the value of the capacitor will increase
% the bandwidth of the signal.