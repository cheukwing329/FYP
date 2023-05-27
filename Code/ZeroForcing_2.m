close all; clear; clc;

%% Setup
A_c = 1;
f_c = 1000;
f_speaker = 44100;
N = 10;                      %number of bit
t = linspace(0,N,N*f_speaker);

%% H matrix
att = 0.8;                  %attenuation
dis = [1,1.118;1.118,1];    %distance
det_t = dis/343;            %time delay
ang_diff = 2*pi*f_c*det_t;  %phase diff
H = att*exp(-1j*ang_diff);
H_inv = inv(H);

%% Message
% m_1 = randi([0 1],1,N) + 1j*randi([0 1],1,N);
% m_2 = randi([0 1],1,N) + 1j*randi([0 1],1,N);
m_1 = [1+1j,1+1j,1+1j,0+0j,0+0j,0+0j,0+1j,1+1j,0+0j,0+0j]; % 1010 and 0101
m_2 = [0+0j,0+1j,1+0j,1+0j,1+0j,0+0,0+1j,1+1j,1+1j,1+0j]; % 1100 and 0011

m_1t=zeros(1,N*f_speaker);
m_2t=zeros(1,N*f_speaker);

for i=1:N
    if real(m_1(i))==1
        m_1t(1,(i-1)*f_speaker+1 :i*f_speaker)=m_1t(1,(i-1)*f_speaker+1 :i*f_speaker) + ones(1,f_speaker);
    elseif real(m_1(i))==0
        m_1t(1,(i-1)*f_speaker+1 :i*f_speaker)=m_1t(1,(i-1)*f_speaker+1 :i*f_speaker) + zeros(1,f_speaker);
    end
    if imag(m_1(i))==1
        m_1t(1,(i-1)*f_speaker+1 :i*f_speaker)=m_1t(1,(i-1)*f_speaker+1 :i*f_speaker) + j*ones(1,f_speaker);
    elseif imag(m_1(i))==0
        m_1t(1,(i-1)*f_speaker+1 :i*f_speaker)=m_1t(1,(i-1)*f_speaker+1 :i*f_speaker) + j*zeros(1,f_speaker);
    end
    if real(m_2(i))==1
        m_2t(1,(i-1)*f_speaker+1 :i*f_speaker)=m_2t(1,(i-1)*f_speaker+1 :i*f_speaker) + ones(1,f_speaker);
    elseif real(m_2(i))==0
        m_2t(1,(i-1)*f_speaker+1 :i*f_speaker)=m_2t(1,(i-1)*f_speaker+1 :i*f_speaker) + zeros(1,f_speaker);
    end
    if imag(m_2(i))==1
        m_2t(1,(i-1)*f_speaker+1 :i*f_speaker)=m_2t(1,(i-1)*f_speaker+1 :i*f_speaker) + j*ones(1,f_speaker);
    elseif imag(m_2(i))==0
        m_2t(1,(i-1)*f_speaker+1 :i*f_speaker)=m_2t(1,(i-1)*f_speaker+1 :i*f_speaker) + j*zeros(1,f_speaker);
    end
end

M = [m_1t;m_2t];

figure();
subplot(2,1,1);
plot(t,real(m_1t));
% ylim([-0.1 1.1]);
title("Real m_1(t)");
subplot(2,1,2);
plot(t,imag(m_1t));
% ylim([-0.1 1.1]);
title("Imag m_1(t)");

figure();
subplot(2,1,1);
plot(t,real(m_2t));
% ylim([-0.1 1.1]);
title("Real m_2(t)");
subplot(2,1,2);
plot(t,imag(m_2t));
% ylim([-0.1 1.1]);
title("Imag m_2(t)");

%% Transmit signal
S = H_inv*M;
s_1 = S(1,:);
s_2 = S(2,:);

x_1=real(s_1.*exp(j*2*pi*f_c*t));
x_2=real(s_2.*exp(j*2*pi*f_c*t));

figure();
subplot(2,1,1);
plot(t,x_1);
title("x_1");
subplot(2,1,2);
plot(t,x_2);
title("x_2");

%% Receive
det_t_op = round(1/343*44100);           %opposite
det_t_dia = round(1.118/343*44100);      %diagonally

x_1_1 = [zeros(1,det_t_op),x_1(1:end-det_t_op)];
x_1_2 = [zeros(1,det_t_dia),x_1(1:end-det_t_dia)];
x_2_1 = [zeros(1,det_t_dia),x_2(1:end-det_t_dia)];
x_2_2 = [zeros(1,det_t_op),x_2(1:end-det_t_op)];

y_1 = att.*x_1_1 + att.*x_2_1;
y_2 = att.*x_1_2 + att.*x_2_2;

figure();
subplot(2,1,1);
plot(t,y_1);
title("y_1");
subplot(2,1,2);
plot(t,y_2);
title("y_2");

%% Demodulation

y_1_dm1 = y_1.*(2.*cos(2*pi*f_c*t));
y_1_dm2 = y_1.*(-2.*sin(2*pi*f_c*t));
y_2_dm1 = y_2.*(2.*cos(2*pi*f_c*t));
y_2_dm2 = y_2.*(-2.*sin(2*pi*f_c*t));
for i=1:N
        y_1_real(i) = round(mean(y_1_dm1((i-1)*f_speaker+1:i*f_speaker)));
        y_1_imag(i) = round(mean(y_1_dm2((i-1)*f_speaker+1:i*f_speaker)));
        y_2_real(i) = round(mean(y_2_dm1((i-1)*f_speaker+1:i*f_speaker)));
        y_2_imag(i) = round(mean(y_2_dm2((i-1)*f_speaker+1:i*f_speaker)));
end


%% Check
figure();
subplot(2,1,1);
plot(real(m_1),'o','Color','r');hold on
plot(y_1_real,'o','Color','g');hold off
xlim([0.5 N+0.5]);
ylim([-0.5 1.5]);
title('Real m_1 checking');
subplot(2,1,2);
plot(imag(m_1),'o','Color','r');hold on
plot(y_1_imag,'o','Color','g');hold off
xlim([0.5 N+0.5]);
ylim([-0.5 1.5]);
title('Imag m_1 checking');

figure();
subplot(2,1,1);
plot(real(m_2),'o','Color','r');hold on
plot(y_2_real,'o','Color','g');hold off
xlim([0.5 N+0.5]);
ylim([-0.5 1.5]);
title('Real m_2 checking');
subplot(2,1,2);
plot(imag(m_2),'o','Color','r');hold on
plot(y_2_imag,'o','Color','g');hold off
xlim([0.5 N+0.5]);
ylim([-0.5 1.5]);
title('Imag m_2 checking');

