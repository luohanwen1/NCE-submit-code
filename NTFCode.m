%% NTF design
clear; close all;
order = 4;      %The order of DSM 
osr = 3.3482;   %oversampling rate
f0 = 0;         % 0 means baseband signal
opt = 1;        %A flag used to request optimized NTF zeros.
Hinf = 3.5;     %The maximum out-of-band gain of the NTF
form= 'CRFF';   %The form of DSM

% Synthesize NTF
ntf = synthesizeNTF(order,osr,opt,Hinf,f0);		% Optimized zero placement
[a,g,b,c] = realizeNTF(ntf,form);   % Convert the NTF into a set of coefficients
% Here we get the parameters values of NTF

%% Next we plot the zero-poles plot and NTF curve
% plot the zero-poles plot
figure;
plotPZ(ntf,'b',8);   
set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold','linewidth',2);
canvasX=4;      
canvasY=4;      
canvasL=10;       
canvasH= 8;    
set(gcf,'unit','centimeters','position',[canvasX canvasY canvasL canvasH])

% plot the NTF curve
figure;
f = [linspace(0,0.75/osr,100) linspace(0.75/osr,0.5,100)];
z = exp(2i*pi*f);
magH = dbv(evalTF(ntf,z));
% magH = evalTF(ntf,z);
plot(f,magH,'linewidth',2);
figureMagic([0 0.5],0.1,1, [-50 20],10,1 );
xlabel('Normalized frequency');
ylabel('Magnitude(dB)');
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','NTF Magnitude Response');
set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold','linewidth',2);
set(gcf,'unit','centimeters','position',[canvasX canvasY canvasL canvasH])
