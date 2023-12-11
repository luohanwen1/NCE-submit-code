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
D_M=6;              
encryptedLevel = 20;       
C_N=encryptedLevel-D_M/2;  
L = 16-1;
M = 8;
osr = L/M; 
bitM = 2;  

%% Tx data
load PRBS_15
PRBS = PRBS';
PRBS = [PRBS ; 0];
tx_data = repmat(PRBS,[D_M 1]);
tx_data = reshape(tx_data,D_M,length(tx_data)/D_M).';
%% Y00 mapping         cipher_I=data_I xor SI  +   BI    cipher_Q=data_Q xor SQ  +   BQ 
load PRBS_25
% PRBS_25 = circshift(PRBS_25,10);  %decrypt-Eve
encry_data=PRBS_25(1:length(PRBS)*(D_M/2+C_N)*2);
encry_data1=reshape(encry_data,(D_M/2+C_N)*2,length(encry_data)/(D_M/2+C_N)/2);
encry_data1=encry_data1';
CCI=zeros(length(PRBS),D_M/2+C_N);
CCQ=CCI;
SI=tx_data(:,1:D_M/2);
RI=encry_data1(:,1:D_M/2);
BI=encry_data1(:,D_M/2+1:D_M/2+C_N);

SQ=tx_data(:,D_M/2+1:end);
RQ=encry_data1(:,D_M/2+C_N+1:D_M/2+D_M/2+C_N);
BQ=encry_data1(:,D_M/2+D_M/2+C_N+1:end);

CCI(:,1:D_M/2)=xor(SI,RI);
CCI(:,1+D_M/2:end)=BI;

CCQ(:,1:D_M/2)=xor(SQ,RQ);
CCQ(:,1+D_M/2:end)=BQ;
CI = (bi2de(CCI,'left-msb'))*2-2^(D_M/2+C_N)+1;
CQ = (bi2de(CCQ,'left-msb'))*2-2^(D_M/2+C_N)+1;
y_Tx=CI+1j*CQ; 
%% FEC RS（255，239）
load testdata.mat;
rx_dsm = tx_I+1i*tx_Q;
rx_dsm = QAMFECDecode(rx_dsm,16,255,239,add);
rx_Idsm = real(rx_dsm);
rx_Qdsm = imag(rx_dsm);
%% Descrambler
N = 8;
descrambler = comm.Descrambler(N,'1 + z^-2 + z^-3 + z^-5 + z^-7', ...
    [3 3 1 1 3 3 1]);
descrDataI = descrambler(rx_Idsm.');
descrDataQ = descrambler(rx_Qdsm.');
descrDataI = descrDataI-4;
descrDataQ = descrDataQ-4;
m=3;
l=1+0.1*m*2^(bitM);
rx_Idsm = descrDataI.'/l;
rx_Qdsm = descrDataQ.'/l;
%% Filter
fback = 0.28;
d1 = designfilt('lowpassfir', ...
    'PassbandFrequency',0.99*fback,'StopbandFrequency',1*fback, ...
    'PassbandRipple',5,'StopbandAttenuation',80, ...
    'DesignMethod','kaiserwin','SampleRate',osr);
dataIOUT = filtfilt(d1,rx_Idsm);
%I
datatemp = resample(dataIOUT,M,1);
dataIDS = resample(datatemp,1,L);
rxI_filter_downsample = dataIDS*(2^encryptedLevel);
%Q
dataQOUT = filtfilt(d1,rx_Qdsm);
dataQtemp = resample(dataQOUT,M,1);
dataQDS = resample(dataQtemp,1,L);
rxQ_filter_downsample = dataQDS*(2^encryptedLevel);
rxI = rxI_filter_downsample;
rxQ = rxQ_filter_downsample;
y_Rx=rxI+1j*rxQ;
%% QNSC RRC filter
sps = 2; % Number of samples per symbol (oversampling factor)
filtlen = 20; % Filter length in symbols
rolloff = 0.1; % Filter rolloff factor
rrcFilter = rcosdesign(rolloff,filtlen,sps);
RX = upfirdn(y_Rx,rrcFilter,1,sps);
RX = RX(filtlen + 1:end - filtlen); % Account for delay
y_Rx1 = RX;
scatterplot(y_Rx1); 
%%  decrypt
E_data=((bi2de(BI,'left-msb'))*2-2^C_N+1)+1j*((bi2de(BQ,'left-msb'))*2-2^C_N+1);
D_data=y_Rx1-E_data.';
D_dataI=real(D_data);
D_dataQ=imag(D_data);
D_data=D_dataI+1j*D_dataQ;
hold on
scatterplot(D_data,1,0,'r.',gcf); 
hold off
%% BER
D_dataI=D_dataI/2^C_N;
D_dataQ=D_dataQ/2^C_N;
D_data=D_dataI+1j*D_dataQ;
D_data = D_data/std(D_data)*sqrt(42);
D_dataI = real(D_data);
D_dataQ = imag(D_data);
load  vik.mat
figure
b1=binscatter(D_dataI(1:length(D_dataI)/2),D_dataQ(1:length(D_dataI)/2),150);
colormap(vik);
figureMagic([-12 12],4,1, [-12 12],4,1 );
axis([-8,8,-8,8])
canvasX=4;     
canvasY=4;    
canvasL=10;      
canvasH= 8;    
set(gcf,'unit','centimeters','position',[canvasX canvasY canvasL canvasH])
set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold','linewidth',2);
colorbar off
grid off
figure
b2=binscatter(D_dataI(length(D_dataI)/2+1:end),D_dataQ(length(D_dataI)/2+1:end),150);
colormap(vik);
figureMagic([-12 12],4,1, [-12 12],4,1 );
axis([-8,8,-8,8])
set(gcf,'unit','centimeters','position',[canvasX canvasY canvasL canvasH])
set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold','linewidth',2);
colorbar off
grid off
S_data = qamdemod(D_data,2^D_M,'bin');
S_data_bit=de2bi(S_data,'left-msb');
SI_rx=S_data_bit(:,1:D_M/2);
SQ_rx=~S_data_bit(:,D_M/2+1:end);
SI_d=double(xor(SI_rx,RI));
SQ_d=double(xor(SQ_rx,RQ));
count1=abs(SI-SI_d);
count1=reshape(count1,1,(length(S_data)*D_M/2));
count2=abs(SQ-SQ_d);
count2=reshape(count2,1,(length(S_data)*D_M/2));
BER=sum(count1+count2);
BER_ratio=BER/(length(count1)+length(count2));
fprintf('Error count: %d/%d, uncoded BER=%.2e, -log10(BER)=%0.3f\n',BER,length(count1)+length(count2),BER_ratio,-log10(BER_ratio));

function rxSig = QAMFECDecode(sig,M,N,K,addbit)
demodSig = qamdemod(sig,M,'bin','OutputType','bit');% Demodulated the noisy signal.
rsDecoder = comm.RSDecoder(N,K,'BitInput',true);
rxData = rsDecoder(demodSig);% Decode the data.
rxData = rxData(1:end-addbit);
rxData1 = reshape(rxData,log2(M),length(rxData)/log2(M));
rxData1 = rxData1.';
rxData2 = bi2de(rxData1,'left-msb');
rxSig =  qammod(rxData2,M,'bin');
rxSig = rxSig.';
end

