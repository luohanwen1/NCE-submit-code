% Copyright (c) 2023, Hanwen Luo, Mengfan Cheng in  National Engineering Research Center for Next Generation Internet Access System
% for Optoelectronics and School of Optical and Electronic Information, 
% Huazhong University of Science and Technology, Wuhan 430074, China
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation 
%    and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its contributors may
%    be used to endorse or promote products derived from this software without 
%    specific prior written permission.
%
% 4. In case results are published that rely on this source code, please cite
% our paper published in <Nature Communications Engineering> entitled 
% "Device-compatible Ultra-high-order Quantum Noise Stream Cipher Based on Delta-Sigma Modulator and Optical Chaos" [1]. 

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
% NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
% OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
% WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%
% [1] https://www.nature.com/articles/...
%%
close all
clear 
%% parameter
D_M=6;                      %Plaintext QAM signal   4=log(16QAM)   6=log(64QAM)
encryptedLevel = 20;        %Ciphertext QAM signal
C_N=encryptedLevel-D_M/2;   %Encryption added bits
L = 16-1;
M = 8;
osr = L/M; %oversampling ratio
bitM = 2;  %DSM quantization bit
%% Tx data
load PRBS_15
PRBS = PRBS';
PRBS = [PRBS ; 0];
tx_data = repmat(PRBS,[D_M 1]);
tx_data = reshape(tx_data,D_M,length(tx_data)/D_M).';
%% Y00 mapping         cipher_I=data_I xor SI  +   BI    cipher_Q=data_Q xor SQ  +   BQ 
load PRBS_25
encry_data=PRBS_25(1:length(PRBS)*(D_M/2+C_N)*2);
encry_data1=reshape(encry_data,(D_M/2+C_N)*2,length(encry_data)/(D_M/2+C_N)/2);
encry_data1=encry_data1';
CI=zeros(length(PRBS),D_M/2+C_N);
CQ=CI;
SI=tx_data(:,1:D_M/2);
RI=encry_data1(:,1:D_M/2);
BI=encry_data1(:,D_M/2+1:D_M/2+C_N);

SQ=tx_data(:,D_M/2+1:end);
RQ=encry_data1(:,D_M/2+C_N+1:D_M/2+D_M/2+C_N);
BQ=encry_data1(:,D_M/2+D_M/2+C_N+1:end);

CI(:,1:D_M/2)=xor(SI,RI);
CI(:,1+D_M/2:end)=BI;

CQ(:,1:D_M/2)=xor(SQ,RQ);
CQ(:,1+D_M/2:end)=BQ;
CI = (bi2de(CI,'left-msb'))*2-2^(D_M/2+C_N)+1;
CQ = (bi2de(CQ,'left-msb'))*2-2^(D_M/2+C_N)+1;
%% QAM modulate
y_Tx=CI.'+1j*CQ.';
scatterplot(y_Tx);
title('QNSC encryption')
%% RRC QNSC
sps = 2; % Number of samples per symbol (oversampling factor)
filtlen = 20; % Filter length in symbols
rolloff = 0.1; % Filter rolloff factor
rrcFilter = rcosdesign(rolloff,filtlen,sps);
Tx = upfirdn(y_Tx,rrcFilter,sps,1);
%% QNSCTX before
tx_I=real(Tx);
tx_Q=imag(Tx);
%% DC 
tx_I = tx_I - mean(tx_I);
tx_Q = tx_Q - mean(tx_Q);
%% DSM
m=3;
tx_Idsm = DSMfunction(tx_I,0.28,L,M,2^(bitM),20,m);
tx_Qdsm = DSMfunction(tx_Q,0.28,L,M,2^(bitM),20,m);
%% scrambler
N = 8;
scrambler = comm.Scrambler(N,'1 + z^-2 + z^-3 + z^-5 + z^-7', ...
    [3 3 1 1 3 3 1]);
scrDataI = scrambler(tx_Idsm.');
scrDataI = scrDataI-4;
scrDataQ = scrambler(tx_Qdsm.');
scrDataQ = scrDataQ-4;
tx_dsm = scrDataI.' +1i*scrDataQ.';
tx_I = real(tx_dsm);
tx_Q = imag(tx_dsm);
save ('testdata.mat',"tx_Q","tx_I");
%%
function v =DSMfunction(sig,band,L,M,bitM,inputM,k)
osr = L/M;      
y_Tx1 = sig-mean(sig); 
datain = y_Tx1/(2^inputM);
datatemp = resample(datain,L,1);
dataOS = resample(datatemp,1,M);
order = 4;     
osr2 = osr/2/band;  
nlev = bitM;
f0 = 0;         
opt = 1;        
Hinf = 3.5; 
form= 'CRFF';
ntf = synthesizeNTF(order,osr2,opt,Hinf,f0);		% Optimized zero placement
[a,g,b,c] = realizeNTF(ntf,form);
ABCD = stuffABCD(a,g,b,c,form);
l=1+0.1*k*bitM;
dataOStemp = dataOS*l;
v = simulateDSM(dataOStemp,ABCD,nlev);
end

