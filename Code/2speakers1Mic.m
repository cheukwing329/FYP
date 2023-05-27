close all; clear; clc;

bit_tx = [1,0,1,0];
A_c = 1;
f_c = 500;

bitpersec =1;
f_speaker = 44100;
Numdelay = 25;

t = linspace(0,1,f_speaker);
playRec = audioPlayerRecorder('Device', 'STUDIO-CAPTURE', 'PlayerChannelMapping', [1,2], 'RecorderChannelMapping', [1]);
s_A = zeros(ceil(10/bitpersec)*44100,1); 
s_B = zeros(ceil(10/bitpersec)*44100,1); 

s_0 = 0;
s_A_1 = sin(2*pi*f_c*t);
s_A_1 = s_A_1(1:f_speaker/bitpersec)';
s_B_1 = sin(2*pi*f_c*t);
s_B_1 = s_B_1(1:f_speaker/bitpersec)';

for i = 1:1:length(bit_tx) %1(Longer)
    if bit_tx(i) == 1
        s_A((i-1)*f_speaker/bitpersec+1:i*f_speaker/bitpersec,1) = s_A_1;
    else 
        s_A((i-1)*f_speaker/bitpersec+1:i*f_speaker/bitpersec,1) = s_0;
    end 
end

for i = 1:1:length(bit_tx) %2
    if bit_tx(i) == 1
        s_B((i-1)*f_speaker/bitpersec+1+Numdelay:i*f_speaker/bitpersec+Numdelay,1) = s_B_1;
    else 
        s_B((i-1)*f_speaker/bitpersec+1+Numdelay:i*f_speaker/bitpersec+Numdelay,1) = s_0;
    end 
end

for i = 1:f_speaker:length(s_A)
    A = s_A(i:i+f_speaker-1,1);
    B = zeros(f_speaker,1);
    audioToDevice = [B,A];
    input = playRec(audioToDevice);
    s_receive_A(i:i+f_speaker-1,1)= input(:); 
end

figure();
plot(s_receive_A);
title('A only');

for i = 1:f_speaker:length(s_A)
    A = zeros(f_speaker,1);
    B = s_B(i:i+f_speaker-1,1);
    %%output(:,2) = s(i:i+f_speaker-1,1);
    audioToDevice = [B,A];
    input = playRec(audioToDevice);
    s_receive_B(i:i+f_speaker-1,1)= input(:); 
end

figure();
plot(s_receive_B);
title('B only');

for i = 1:f_speaker:length(s_A)
    A = s_A(i:i+f_speaker-1,1);
    B = s_B(i:i+f_speaker-1,1);
    %%output(:,2) = s(i:i+f_speaker-1,1);
    audioToDevice = [B,A];
    input = playRec(audioToDevice);
    s_receive_AB(i:i+f_speaker-1,1)= input(:); 
end

figure();
plot(s_receive_AB);
title('2Speakers');