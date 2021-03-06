clc
clear all
close all
% dPhi(t)/dt = w(t) = w0*t 
% Phi(t) = 0.5*w0*t^2 + C

%parameters
T=15e-3;%duration of signal 
fs=2.5*100e3;%f sample
Ts=1/fs;%T sample
f0=15e3;%Low freq
f1=13e3;%High freq
t=0:Ts:T; %time
Fpass = 2e3;
Fstop = 5.5e3;
HarmonyAmplitudedb = [10,20,7,12,8,9];
HarmonyAmplitude = 10.^(HarmonyAmplitudedb/20);
HarmonyAmplitude = HarmonyAmplitude/sum(HarmonyAmplitude);
descimationfactor=5;
NumOfFiles = 1e3 + 1;
Results = {};
% clip_val = 0.3;

Amplitude_Variance = 5; %dB
freq_Variance = 50; %Hz

for j =1:NumOfFiles
    %Chirp Generation
    x=0;
    HarmonyAmplitudedb_temp = HarmonyAmplitudedb + sqrt(Amplitude_Variance)*randn(1,6);
    HarmonyAmplitude = 10.^(HarmonyAmplitudedb_temp/20);
    HarmonyAmplitude = HarmonyAmplitude/sum(HarmonyAmplitude);
    for i =1:6
        z=(t-mean(t)) / max(t-mean(t));
        f0_temp = f0 + sqrt(freq_Variance) * randn(1);
        f1_temp = f1 + sqrt(freq_Variance) * randn(1);
        phi = 0.5*(f0_temp+f1)*t - ((f0_temp-f1)/2) * ((2*t/T -1).^4 * (T/8)); %we assume w is 11-13 kHz with most oftern in 12kHz. so we did 3thOrder interpulation
        x = x + HarmonyAmplitude(i)*sin(2*i*pi*phi);%Ith-harmony
    end
    
    %Filter
    LPF = RealFilter(fs,Fpass,Fstop);
    
    %Noise adding
    x=x + filter(LPF.Numerator,1,wgn(1,length(x),-20)) + wgn(1,length(x),-30); % Total Low Freq Noise + Noise    
%     x(x>clip_val) = clip_val;
%     x(x<(-1*clip_val)) = -1*clip_val;
    Results{j,1} = x;
    %Results{j,2} = x_descimated;
    %Results{j,2} = T;
end

csvwrite(['ResultsSimulation',num2str(NumOfFiles),'.csv'],Results)

%save('ResultsSimulation','Results'); %matfile

%FFT Plot
N=length(x);%Nfft
f=(-N/2 : N/2 -1)*fs/N;
Y=fft(x);
Y=fftshift(Y);
figure;plot(f(f>0),db(abs(Y(f>0))));xlabel('frequancy');ylabel('db');title('FFT of signal');

%Spectogram
figure;stft(x,fs,'Window',kaiser(256,5),'OverlapLength',220,'FFTLength',512);ylim([0,fs/(2e3)]);
title("Bat signal smapled in 250Khz Clipped")

%Time Domain
figure;plot(t,x);xlabel("T[ms]");
title("Bat signal smapled in 250Khz - Time Domain")

