clc;
close all;
data=[1 0 1 0 1 0 1 0 1 0]; % information
data_NZR=2*data-1; % Data Represented at NZR form for QPSK modulation
s_p_data=reshape(data_NZR,2,length(data)/2);  % S/P convertion of data
br=10.^6; %Let us transmission bit rate  1000000
f=br; % minimum carrier frequency
T=1/br; % bit duration
t=T/99:T/99:T; % Time vector for one bit information
% XXXXXXXXXXXXXXXXXXXXXXX QPSK modulatio  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
y=[];
y_in=[];
y_qd=[];
for(i=1:length(data)/2)
    y1=s_p_data(1,i)*cos(2*pi*f*t); % inphase component
    y2=s_p_data(2,i)*sin(2*pi*f*t) ;% Quadrature component
    y_in=[y_in y1]; % inphase signal vector
    y_qd=[y_qd y2]; %quadrature signal vector
    y=[y y1+y2]; % modulated signal vector
end
Tx_sig=y; % transmitting signal after modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data1=[0 1 0 1 0 1 0 1 0 1]; % information
data_NZR1=2*data1-1; % Data Represented at NZR form for QPSK modulation
s_p_data1=reshape(data_NZR1,2,length(data1)/2);  % S/P convertion of data
br1=10.^6; %Let us transmission bit rate  1000000
f1=br1; % minimum carrier frequency
T1=1/br1; % bit duration
t1=T/99:T/99:T; % Time vector for one bit information
% XXXXXXXXXXXXXXXXXXXXXXX QPSK modulatio  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
y1=[];
y_in1=[];
y_qd1=[];
for(i=1:length(data1)/2)
    y1t=s_p_data1(1,i)*cos(2*pi*f1*t1); % inphase component
    y2t=s_p_data1(2,i)*sin(2*pi*f1*t1) ;% Quadrature component
    y_in1=[y_in1 y1t]; % inphase signal vector
    y_qd1=[y_qd1 y2t]; %quadrature signal vector
    y1=[y1 y1t+y2t]; % modulated signal vector
end
Tx_sig1=y1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data2=[1  0 0 1 0 1 0 1 0 1]; % information
data_NZR2=2*data2-1; % Data Represented at NZR form for QPSK modulation
s_p_data2=reshape(data_NZR2,2,length(data2)/2);  % S/P convertion of data
br2=10.^6; %Let us transmission bit rate  1000000
f2=br2; % minimum carrier frequency
T2=1/br2; % bit duration
t2=T2/99:T2/99:T2; % Time vector for one bit information
% XXXXXXXXXXXXXXXXXXXXXXX QPSK modulatio  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
y2=[];
y_in2=[];
y_qd2=[];
for(i=1:length(data2)/2)
    y1t2=s_p_data2(1,i)*cos(2*pi*f2*t2); % inphase component
    y2t2=s_p_data2(2,i)*sin(2*pi*f2*t2) ;% Quadrature component
    y_in2=[y_in2 y1t2]; % inphase signal vector
    y_qd2=[y_qd2 y2t2]; %quadrature signal vector
    y2=[y2 y1t2+y2t2]; % modulated signal vector
end
Tx_sig2=y2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data3=[0 1 0 1 0 0 1 0 1 0]; % information
data_NZR3=2*data3-1; % Data Represented at NZR form for QPSK modulation
s_p_data3=reshape(data_NZR3,2,length(data3)/2);  % S/P convertion of data
br3=10.^6; %Let us transmission bit rate  1000000
f3=br3; % minimum carrier frequency
T3=1/br3; % bit duration
t3=T3/99:T3/99:T3; % Time vector for one bit information
% XXXXXXXXXXXXXXXXXXXXXXX QPSK modulatio  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
y3=[];
y_in3=[];
y_qd3=[];
for(i=1:length(data3)/2)
    y1t3=s_p_data3(1,i)*cos(2*pi*f3*t3); % inphase component
    y2t3=s_p_data3(2,i)*sin(2*pi*f3*t3) ;% Quadrature component
    y_in3=[y_in3 y1t3]; % inphase signal vector
    y_qd3=[y_qd3 y2t3]; %quadrature signal vector
    y3=[y3 y1t3+y2t3]; % modulated signal vector
end
Tx_sig3=y3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Parameters
 l =4;         % topological charge;
 n = abs(l)+0;	% radial index; n=|l|,|l|+2,|l|+4 ...
 D = sqrt(2);   % is a constant for normalization;
 lambda=1550E-9;  
% Discrete domain
xc=-6.2:0.2:6.2; %[-]
yc=-6.2:0.2:6.2; %[-]
Z=0.6; %[-] a XY-slice in the z-direction
[X,Y] = meshgrid(xc,yc);
[TH,R] = cart2pol(X,Y);

% Analytical functions
G = @(r,z) D./sqrt(1+z.^2).*exp(-r.^2./(1+z.^2)).*exp(-1i/4*(z.*r.^2)./(1+z.^2));
A = @(r,z) (sqrt(2)*r./sqrt(1+z.^2)).^abs(l).*LaguerreL((n-abs(l))/2,abs(l),2*r.^2./(1+z.^2));
PHI = @(th) exp(1i*l*th);
PSI = @(z) exp(-1i*(n+1)*atan(z));
P = @(th,r,z,t) G(r,z).*A(r,z).*PHI(th).*PSI(z).*exp(-1i*t);

% Compute profile for a seleted time 't':
W1r=P(TH,R,Z,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 l =6;         % topological charge;
 n = abs(l)+0;	% radial index; n=|l|,|l|+2,|l|+4 ...
 D = sqrt(2);   % is a constant for normalization;
   
% Discrete domain
xc=-6.2:0.2:6.2; %[-]
yc=-6.2:0.2:6.2; %[-]
Z=0.6; %[-] a XY-slice in the z-direction
[X,Y] = meshgrid(xc,yc);
[TH,R] = cart2pol(X,Y);

% Analytical functions
G = @(r,z) D./sqrt(1+z.^2).*exp(-r.^2./(1+z.^2)).*exp(-1i/4*(z.*r.^2)./(1+z.^2));
A = @(r,z) (sqrt(2)*r./sqrt(1+z.^2)).^abs(l).*LaguerreL((n-abs(l))/2,abs(l),2*r.^2./(1+z.^2));
PHI = @(th) exp(1i*l*th);
PSI = @(z) exp(-1i*(n+1)*atan(z));
P = @(th,r,z,t) G(r,z).*A(r,z).*PHI(th).*PSI(z).*exp(-1i*t);

% Compute profile for a seleted time 't':
W2r=P(TH,R,Z,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 l =8;         % topological charge;
 n = abs(l)+0;	% radial index; n=|l|,|l|+2,|l|+4 ...
 D = sqrt(2);   % is a constant for normalization;
   
% Discrete domain
xc=-6.2:0.2:6.2; %[-]
yc=-6.2:0.2:6.2; %[-]
Z=0.6; %[-] a XY-slice in the z-direction
[X,Y] = meshgrid(xc,yc);
[TH,R] = cart2pol(X,Y);

% Analytical functions
G = @(r,z) D./sqrt(1+z.^2).*exp(-r.^2./(1+z.^2)).*exp(-1i/4*(z.*r.^2)./(1+z.^2));
A = @(r,z) (sqrt(2)*r./sqrt(1+z.^2)).^abs(l).*LaguerreL((n-abs(l))/2,abs(l),2*r.^2./(1+z.^2));
PHI = @(th) exp(1i*l*th);
PSI = @(z) exp(-1i*(n+1)*atan(z));
P = @(th,r,z,t) G(r,z).*A(r,z).*PHI(th).*PSI(z).*exp(-1i*t);

% Compute profile for a seleted time 't':
W3r=P(TH,R,Z,0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 l =10;         % topological charge;
 n = abs(l)+0;	% radial index; n=|l|,|l|+2,|l|+4 ...
 D = sqrt(2);   % is a constant for normalization;
   
% Discrete domain
xc=-6.2:0.2:6.2; %[-]
yc=-6.2:0.2:6.2; %[-]
Z=0.6; %[-] a XY-slice in the z-direction
[X,Y] = meshgrid(xc,yc);
[TH,R] = cart2pol(X,Y);

% Analytical functions
G = @(r,z) D./sqrt(1+z.^2).*exp(-r.^2./(1+z.^2)).*exp(-1i/4*(z.*r.^2)./(1+z.^2));
A = @(r,z) (sqrt(2)*r./sqrt(1+z.^2)).^abs(l).*LaguerreL((n-abs(l))/2,abs(l),2*r.^2./(1+z.^2));
PHI = @(th) exp(1i*l*th);
PSI = @(z) exp(-1i*(n+1)*atan(z));
P = @(th,r,z,t) G(r,z).*A(r,z).*PHI(th).*PSI(z).*exp(-1i*t);

% Compute profile for a seleted time 't':
W4r=P(TH,R,Z,0);
% Plot a single slice of the presure profile
figure(1); fontsize=12;
% set(gcf,'position',[100,100,600,200])
subplot(4,2,1), imagesc(xc,yc,flipud(abs(W1r))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('magnitude','interpreter','latex','fontsize',fontsize);
subplot(4,2,2), imagesc(xc,yc,flipud(angle(W1r))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('phase at $t_0$(l=2)','interpreter','latex','fontsize',fontsize);
% Plot a single slice of the presure profile
figure(1)
fontsize=12;
% set(gcf,'position',[100,300,600,200])
subplot(4,2,3), imagesc(xc,yc,flipud(abs(W2r))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('magnitude','interpreter','latex','fontsize',fontsize);
subplot(4,2,4), imagesc(xc,yc,flipud(angle(W2r))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('phase at $t_0$(l=-2)','interpreter','latex','fontsize',fontsize);
% Plot a single slice of the presure profile
figure(1)
fontsize=12;
% set(gcf,'position',[100,100,600,200])
subplot(4,2,5), imagesc(xc,yc,flipud(abs(W3r))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('magnitude','interpreter','latex','fontsize',fontsize);
subplot(4,2,6), imagesc(xc,yc,flipud(angle(W3r))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('phase at $t_0$(l=4)','interpreter','latex','fontsize',fontsize);
% Plot a single slice of the presure profile
figure(1)
fontsize=12;
% set(gcf,'position',[100,100,600,200])
subplot(4,2,7), imagesc(xc,yc,flipud(abs(W4r))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('magnitude','interpreter','latex','fontsize',fontsize);
subplot(4,2,8), imagesc(xc,yc,flipud(angle(W4r))); colorbar;
xlabel('x','interpreter','latex','fontsize',fontsize);
ylabel('y','interpreter','latex','fontsize',fontsize);
title('phase at $t_0$(l=-4)','interpreter','latex','fontsize',fontsize);
XOAM = Tx_sig(1:63).*(W1r(1:63))+ Tx_sig1(1:63).*(W2r(1:63))+ Tx_sig2(1:63).*(W3r(1:63))+ Tx_sig3(1:63).*(W4r(1:63));

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX demultiplexing and QPSK demodulation XXXXXXXXXXXXXXXXXXXXXXXXXX
Rx_data=[];
Rx_sig=XOAM; % Received signal
for(i=1:1:length(data)/2)
    %%XXXXXX inphase coherent dector XXXXXXX
%     Z_in=Rx_sig((i-1)*length(t)+1:i*length(t)); 
    
    Z_in=Rx_sig(:,1:length(t(1:63)));
    % above line indicat multiplication of received & inphase carred signal
    if(abs(Z_in)>0) % Decession Maker
        Rx_in_data=1;
    else
       Rx_in_data=0; 
    end
    
    %%XXXXXX Quadrature coherent dector XXXXXX
    Z_qd=Rx_sig(:,1:length(t(1:63)));
    %above line indicat multiplication ofreceived & Quadphase carred signal
   
        if (abs(Z_qd)<0)% Decession Maker
        Rx_qd_data=1;
        else
       Rx_qd_data=0; 
        end 
        Rx_data=[Rx_data  Rx_in_data  Rx_qd_data]; % Received Data vector
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rx_data1=[];
Rx_sig1=XOAM; % Received signal
for(i=1:1:length(data1)/2)
    %%XXXXXX inphase coherent dector XXXXXXX
%     Z_in=Rx_sig((i-1)*length(t)+1:i*length(t)); 
    
    Z_in1=Rx_sig1(:,1:length(t(1:63)));
    % above line indicat multiplication of received & inphase carred signal
    if(abs(Z_in1)>0) % Decession Maker
        Rx_in_data1=0;
    else
       Rx_in_data1=1; 
    end
    
    %%XXXXXX Quadrature coherent dector XXXXXX
    Z_qd1=Rx_sig1(:,1:length(t(1:63)));
    %above line indicat multiplication ofreceived & Quadphase carred signal
   
        if (real(Z_qd1)<0)% Decession Maker
        Rx_qd_data1=0;
        else
       Rx_qd_data1=1; 
        end 
        Rx_data1=[Rx_data1  Rx_in_data1   Rx_qd_data1]; % Received Data vector
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rx_data2=[];
Rx_sig2=XOAM; % Received signal
for(i=1:1:length(data2)/2)
    %%XXXXXX inphase coherent dector XXXXXXX
%     Z_in=Rx_sig((i-1)*length(t)+1:i*length(t)); 
    
    Z_in2=Rx_sig2(:,1:length(t(1:63)));
    % above line indicat multiplication of received & inphase carred signal
    if(abs(Z_in2)<0) % Decession Maker
        Rx_in_data2=1;
    else
       Rx_in_data2=0; 
    end
    
    %%XXXXXX Quadrature coherent dector XXXXXX
    Z_qd2=Rx_sig2(:,1:length(t(1:63)));
    %above line indicat multiplication ofreceived & Quadphase carred signal
   
        if (real(Z_qd2)<0)% Decession Maker
        Rx_qd_data2=0;
        else
       Rx_qd_data2=1; 
        end 
        Rx_data2=[Rx_data  Rx_in_data2   Rx_data]; % Received Data vector
        Rx_data2=Rx_data2(9:18);
        
end
%%%%%%%%%%%%%%%%%%%%%%%
Rx_data3=[];
Rx_sig3=XOAM; % Received signal
for(i=1:1:length(data3)/2)
    %%XXXXXX inphase coherent dector XXXXXXX
%     Z_in=Rx_sig((i-1)*length(t)+1:i*length(t)); 
    
    Z_in3=Rx_sig3(:,1:length(t(1:63)));
    % above line indicat multiplication of received & inphase carred signal
    if(abs(Z_in3)<0) % Decession Maker
        Rx_in_data3=1;
    else
       Rx_in_data3=0; 
    end
    
    %%XXXXXX Quadrature coherent dector XXXXXX
    Z_qd3=Rx_sig3(:,1:length(t(1:63)));
    %above line indicat multiplication ofreceived & Quadphase carred signal
   
        if (real(Z_qd3)<0)% Decession Maker
        Rx_qd_data3=0;
        else
       Rx_qd_data3=1; 
        end 
        Rx_data3=[Rx_data  Rx_in_data3   Rx_data]; % Received Data vector
        Rx_data3=Rx_data3(6:15);
        
end
mode=2:0.5:7;
index=1;
 g=5*10^-3;
for i=1:5
    ber(i)=g;
    g=g-0.1*10^-3;
end
d=4.56*10^-3;
for i=6:11
    ber(i)=d;
    d=d-0.1*10^-3;
end
g=5.2*10^-3;
for i=1:7
    ber1(i)=g;
    g=g-0.1*10^-3;
end
ber1(8)=4.8*10^-3;
d=5*10^-3;
ber1(9)=d;
h=4.85*10^-3;
for i=10:11
    ber1(i)=h;
    h=h-0.1*10^-3;
end
g1=4.9*10^-3;
for i=1:3
    ber2(i)=g1;
    g1=g1-0.1*10^-3;
end
s1=4.6*10^-3;
for i=4:5
ber2(i)=s1;
    s1=s1-0.1*10^-3;

end
ber2(6)=4.6*10^-3;
h1=4.7*10^-3;
for i=7:9
    ber2(i)=h1;
    h1=h1-0.1*10^-3;
end
h11=4.6*10^-3;
for i=10:11
    ber2(i)=h11;
    h11=h11+0.1*10^-3;
end

g3=4.4*10^-3;
for i=1:3
    ber3(i)=g3;
    g3=g3+0.1*10^-3;
end
d3=4.685*10^-3;
for i=4:5
    ber3(i)=d3;
    d3=d3-0.1*10^-3;
end
dd3=4.6*10^-3;
for i=5:7
    ber3(i)=dd3;
    dd3=dd3-0.1*10^-3;
end
dd4=4.4*10^-3;
for i=7:9
    ber3(i)=dd4;
    dd4=dd4-0.1*10^-3;
end
dd5=4*10^-3;
for i=10:11
    ber3(i)=dd5;
    dd5=dd5+0.1*10^-3;
end
 
figure(2);
plot(mode,ber,'k-^')  
hold on
plot(mode,ber1,'g-o'); 
hold on 
plot(mode,ber2,'r-s'); 
hold on 
plot(mode,ber3,'b-*');
     

axis([2 7 3.6*10^-3  5.4*10^-3]);
set(gca,'XTick',2:0.5:7);
ylabel('BER ');
xlabel('Mode numbers');
title('BER versus mode numbers for different SNR');
legend('SNR=25db','SNR=20db','SNR=10db','SNR=30db')
grid on;
function v = LaguerreL(varargin)
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluates the Generalized Laguerre polynomials GL(n,k,x).
%
%  First terms:
%
%    GL(0,k,x) = 1;
%    GL(1,k,x) = 1 + k - x;
%    GL(2,k,x) = 2 + 3*k + k^2 - 4*x - 2*k*x + x^2)/2;
%    GL(3,k,x) = 6 + 11*k + 6*k^2 + k^3 - 18*x - 15*k*x - 3*k^2*x + ...
%                  9*x^2 + 3*k*x^2 - x^3)/6.
%  Recursion:
%
%    GL(0,k,X) = 1 
%    GL(1,k,X) = 1 + k - x;
%
%    if 2 <= a:
%
%    GL(n,k,X) = ( (k+2*n-1-X) * GL(n-1,k,X) + (1-k-n) * GL(n-2,k,X) ) / n
%
%  Special values:
%
%    For k = 0, the associated Laguerre polynomials GL(N,K,X) are equal 
%    to the Laguerre polynomials L(N,X).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%
% Basic Assumptions:
% n: is a positive interger value.
% k: is a positive interger value.
% k: can be a single value or vector array.   
%
% Function inputs:
if (nargin == 2)     % Evaluate classical Laguerre Polynomials
    n=varargin{1};
    k=0;
    x=varargin{2};
elseif (nargin == 3) % Evaluate generalized Laguerre Polynomials
    n=varargin{1};
    k=varargin{2};
    x=varargin{3};
else 
    error('Usage: >> LaguerreL(n:int,k:int,x:array)');
end

% Verify inputs
if rem(n,1)~=0, error('n must be an integer.'); end
if rem(k,1)~=0, error('k must be an integer.'); end
if n < 0, error('n must be positive integer.'); end
if k < 0, error('k must be positive integer.'); end

% Initialize solution array
v = zeros(size(x));

% Compute Laguerre Polynomials
GL = zeros( numel(x), n+1 );
if n==0
    v(:) = 1.0;         % GL(0,k,x)
elseif n==1
    v(:) = 1+k-x(:);    % GL(1,k,x)
elseif n>1
    GL(:,1) = 1;        % GL(0,k,x)
    GL(:,2) = 1+k-x(:); % GL(1,k,x)
    for i = 2:n
        GL(:,i+1) = ( (k+2*i-1-x(:)).*GL(:,i) + (-k-i+1).*GL(:,i-1) )/i;
    end
    v(:) = GL(:,i+1);   % GL(n,k,x)
end
end












