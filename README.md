# FYP: Real-Time Wireless Acoustic Communication System
Produced by Cheuk Wing Sin

# Brief introduction
This final year project is to build an acoustic communication system. 

Microphone is used as system input and speaker is used as system output. As well as audio interface is used to connect microphone speaker with the computer.

The system use the sound that human can hear as a signal to transfer image(Lenna.png) from output tot input.

It is a lab-based project. In this project, MATLAB is used in programming for the purpose of component control and data analysis. 

There are two phases in the project: 

## Phase 1: SISO 
In Phase 1, only 1 microphone and 1 speaker are used to simulate a single input single output, SISO system. 

There are three basic communication schemes are used, which include Amplitude Shift Keying (ASK), Frequency Shift Keying (FSK) and Phase Shift Keying (PSK).
The corresponding code are SendPic_ASK.m, SendPic_FSK.m, SendPic_PSK.m in the "Code" folder.

## Phase 2: MIMO
In Phase 2, 2 microphones and 1 speaker are used first to simulate the signal interference. 
The corresponding code is 2speakers1Mic.m in the "Code" folder.

Then, 2 microphones and 2 speakers are used to simulate the zero-forcing beamforming in multiple input multiple output, MIMO system.
The corresponding code is ZeroForcing_2.m in the "Code" folder.

# For more details
Please refer to "FYP_Report_SINCheukWing.pdf"
