%% Step Function
Fs = 1000;
T = 1/Fs;
L = 1000;
t = (0:L-1)*T;
stepin = zeros(1,1000);
for i = 1:1000
    if i <= 30
        stepin(i) = 0;
    else
        stepin(i) = 1;
    end
end
% figure(1)
% plot(t,stepin)
% title('Step Function in Time Domain')
% xlabel('Time (t)')
% xlim([0,1.01]);
% ylabel('V(t)')
% ylim([-0.01,1.01])

Y = fft(stepin);
P2 = abs(Y/L);
P1_stepin = P2(1:L/2+1);
P1_stepin(2:end-1) = 2*P1_stepin(2:end-1);
f = Fs*(0:(L/2))/L;
% figure(2)
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of Step Function')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
%% Sin Input
Fs = 1000;
T = 1/Fs;
L = 1000;
t = (0:L-1)*T;
sinin = sin(2*pi*(1/0.03)*t);
% figure(3)
% plot(t,sinin)
% title('Sine Function in Time Domain')
% xlabel('Time (t)')
% ylabel('V(t)')

Y = fft(sinin);
P2 = abs(Y/L);
P1_sinein = P2(1:L/2+1);
P1_sinein(2:end-1) = 2*P1_sinein(2:end-1);
f = Fs*(0:(L/2))/L;
% figure(4)
% plot(f,P1_sinein) 
% title('Single-Sided Amplitude Spectrum of Sine Function')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

%% Gaussian Pulse
Fs = 1000;           % Sampling frequency
tg = 0:1/Fs:1;  % Time vector 
L = length(tg);      % Signal length
dis = [0:1e-3:1];
gaussin = normpdf(dis,0.53,0.03);
gaussin = gaussin./max(gaussin);
% figure(5);
% plot(tg,gaussin)
% title('Gaussian Pulse in Time Domain')
% xlabel('Time (t)')
% ylabel('X(t)')

n = 2^nextpow2(L);
Y = fft(gaussin,n);
f = Fs*(0:(n/2))/n;
P = abs(Y/n);
% figure(6);
% plot(f,P(1:n/2+1)) 
% title('Gaussian Pulse in Frequency Domain')
% xlabel('Frequency (f)')
% ylabel('|P(f)|')