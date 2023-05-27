close all; clear; clc;

%% Read Image "Lenna"
im = im2double(imread('Lenna.png'));
im = rgb2gray(im);
im = imresize(im,[64 64]);
figure();
imshow(im);
title('Orginal Image');

%% Setup
Q = 4;
A_c = 1;
f_c = 500;
bitpersec = 70;

f_speaker = 44100;
t = linspace(0,1,f_speaker);

playRec = audioPlayerRecorder('Device', 'STUDIO-CAPTURE', 'PlayerChannelMapping', [1], 'RecorderChannelMapping', [1]);

error=0;
error_percent=0;


s_0 = -cos(2*pi*f_c*t);
s_0 = s_0(1:f_speaker/bitpersec)';
s_1 = cos(2*pi*f_c*t);
s_1 = s_1(1:f_speaker/bitpersec)';
%% Image --> bit:
bit_im = im(:);
bit_im = bit_im./(max(abs(bit_im)));
bit_im = round(bit_im.*(Q-1));
bit_tx = de2bi(bit_im,'left-msb');
bit_tx = bit_tx';
bit_tx = bit_tx(:);
totalbit = length(bit_tx);
bit_tx = [1;1;bit_tx];
totalbit_tx = length(bit_tx);

s = zeros((ceil(totalbit_tx/bitpersec)+2)*f_speaker,1); 
s_receive = zeros((ceil(totalbit_tx/bitpersec)+2)*f_speaker,1);
s_rrx = zeros((ceil(totalbit_tx/bitpersec)+2)*f_speaker,1);

%% PSK Modulation
for i = 1:1:totalbit_tx
    if bit_tx(i) == 1
        s((i-1)*f_speaker/bitpersec+1:i*f_speaker/bitpersec,1)=s_1;
    else 
        s((i-1)*f_speaker/bitpersec+1:i*f_speaker/bitpersec,1)=s_0;
    end 
end

%% Tx & Rx
for i = 1:f_speaker:length(s)
    output = s(i:i+f_speaker-1,1);
    input = playRec(output);
    s_receive(i:i+f_speaker-1,1)= input(:); 
end

%% Correlation
s_receive=s_receive(f_speaker:end);
[corr,lags]=xcorr(s_1,s_receive(1:3*f_speaker/bitpersec));
delay = find(corr==max(corr));
delay = abs(lags(delay));
s_rrx = s_receive((delay+1):(ceil(totalbit_tx/bitpersec)*f_speaker+delay));

%% Normalize
s_rrx(:) = s_rrx./max(s_rrx((3/10)*f_speaker/bitpersec:(7/10)*f_speaker/bitpersec));

%% Bit Detection
for i=1:totalbit_tx
    x_1 = s_rrx((i-1)*f_speaker/bitpersec+1:(i)*f_speaker/bitpersec).*s_1;
    x_0 = s_rrx((i-1)*f_speaker/bitpersec+1:(i)*f_speaker/bitpersec).*s_0;
    r_1 = sum(x_1);
    r_0 = sum(x_0);
    if (r_1-r_0)>=0
        bit_rx(i) = 1;
    else
        bit_rx(i)=0;
    end
end

figure();
plot(bit_tx,'o','Color','g');hold on
plot(bit_rx,'o','Color','r');hold off
xlim([0.5 totalbit+0.5]);
ylim([-0.5 1.5]);
title('Red=Received,Green=Original');
for i = 1:totalbit_tx
    if bit_tx(i) ~= bit_rx(i)
        error = error+1;
    end
end

error_percent = error/totalbit_tx;

%% bit --> Image
bit_ref = bit_tx(3:end);
bit_ref = reshape(bit_ref, log2(Q),[])';
im_ref = bi2de(bit_ref,'left-msb');
im_ref = im_ref ./((Q-1));
im_ref = reshape(im_ref,64,64);

bit_r = bit_rx(3:end);
bit_r = reshape(bit_r, log2(Q),[])';
im_rx = bi2de(bit_r,'left-msb');
im_rx = im_rx ./((Q-1));
im_rx = reshape(im_rx,64,64);

figure();
imshow(im_ref);
title('Reference Image');

figure();
imshow(im_rx);
title('Received Image');