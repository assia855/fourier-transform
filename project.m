%SQUARE SIGNAL:
clear all;close all;clc 

t = linspace(0,1*2*pi)';%define the time

x = square(t);%the build function of square wave

%hold on;

syms t n

T=2;%the period

n=1:20%the number of harmonic

w0=2*pi/T;%

 
%the coeffiecient of the equation
a0 = (1/T).*int(1,t,0,2);
an = (2/T)*int(1*cos(n*w0.*t),t,0,1);
bn = (2/T)*int(1*sin(n*w0.*t),t,0,1);

M=[];%the matrix

t = linspace(0,5*2*pi)';

plot(t,x,'.-'),title('square signal')  %plotting the scare signal

figure;
stem(bn);title('the harmonics');

t2 = linspace(0,pi,512)';

for n=1:20

    S=bn(n).*sin(n*w0.*t2);
    M(n,:)=S';

end

figure

plot(M')

%3D representation:

figure

Rcarre=cumsum(M);

mesh(Rcarre);

% The reconstructed signal:

figure;

plot(Rcarre(20,:)),title("The reconsctructed signal");

Fs=1000;
N=4;

%frequency axis:
f = (0:n-1)*(Fs/n-1);
%Module calculation:anlysis fft of the fundamental signal
x_fft=abs(fft(x));
v=max(x_fft);x_fft = abs(fft(x))./v;

%the fft of polluted signal:
xn_fft=abs(fft(Rcarre(20,:)));
vv=max(x_fft);xn_fft=abs(fft(Rcarre(20,:)))./vv;

%to plot fft:
noSamples = length(t2);
f = 0:Fs/noSamples:Fs - Fs/noSamples;

noSamples = length(t);
f1 = 0:Fs/noSamples:Fs - Fs/noSamples;

figure;
subplot(2,1,1);
plot(f,xn_fft),title("spectrum of the reconstructued signal")

subplot(2,1,2);
plot(f1,x_fft),title("spectrum of the fundamental signal")
 
%%
%SWATOOTH SIGNAL:
clear all;close all;clc;

t = linspace(0,5*2*pi)';
x = sawtooth(t);% the sawtooth build function 

hold on

syms t n

T=2;%perid

n=1:20;%number of harmonic

w0=2*pi/T;

A=2; %amplitude

%Coefficient: 
a0 = A/2;
an =0;
bn = -A ./(n*pi);

z=[];%Matrix

t = linspace(0,5*2*pi)';
plot(t,x,'.-'),title('sawtooth signal')  

figure;
stem(bn),title('The harmonics');

t2 = linspace(0,pi,512)';

for n=1:20

    S=bn(n).*sin(n*w0.*t2);
    z(n,:)=S';

end

figure

plot(z');

%3D signal:

figure
sawtooth=cumsum(z);
mesh(sawtooth),title("3D")

%reconstructed signal:

figure
plot(sawtooth(20,:));title("reconstructed signal")

%FFT:
Fs=1000;
N=4;
%frequency axis:
f = (0:n-1)*(Fs/n-1);
%Module calculation:anlysis fft of thefundamental signal
x_fft=abs(fft(x));
v=max(x_fft);x_fft = abs(fft(x))./v;

%let's the fft of polluted signal:
xn_fft=abs(fft(sawtooth(20,:)));
vv=max(x_fft);xn_fft=abs(fft(sawtooth(20,:)))./vv;

%to plot fft:
noSamples = length(t2);
f = 0:Fs/noSamples:Fs - Fs/noSamples;

noSamples = length(t);
f1 = 0:Fs/noSamples:Fs - Fs/noSamples;

figure;
subplot(2,1,1);
plot(f,xn_fft),title("spectrum of the reconstructued signal")
subplot(2,1,2);
plot(f1,x_fft),title("spectrum of the fundamental signal")


%%
%TRIANGLE SIGNAL:

clear all; close all; clc 

%t=0:10;
t = linspace(0,5*2*pi)';
T=2;
w0=2*pi/T;
x=abs(sawtooth(t));

plot(t,x),title("Triangular signal")
H=15;
n=1:15;
b=(8./n.^2*pi.^2).*sin(n*pi./2);%all even harmonics egual zero
figure;stem(b);title("The harmonics")

t2 = linspace(0,pi,512);
a0=0;
figure;
hold on;
for n=1:H;
    harmonic=b(n)*sin(n*w0*t2);
    fnum=a0+harmonic;
    pause(0.01);
    Z(n,:)= fnum';
   
end
plot( Z');title("the accumulated signal");


%3D

figure;

Rcarre=cumsum(Z );

mesh(Rcarre);title("3D ")
figure;plot(Rcarre(n,:));title("the reconstracted signal");

%FFT:
Fs=1000;
N=4;
%frequency axis:
f0 = (0:n-1)*(Fs/n-1);
%Module calculation:anlysis fft of the fundamental signal
x_fft=abs(fft(x));
v=max(x_fft);x_fft = abs(fft(x))./v;

%let's the fft of polluted signal:
xn_fft=abs(fft(Rcarre(n,:)));
vv=max(x_fft);xn_fft=abs(fft(Rcarre(n,:)))./vv;

%to plot fft:
noSamples = length(t2);
f2 = 0:Fs/noSamples:Fs - Fs/noSamples;

noSamples = length(t);
f1 = 0:Fs/noSamples:Fs - Fs/noSamples;

figure;
subplot(2,1,1);
plot(f2,xn_fft),title("spectrum of the reconstructued signal")

subplot(2,1,2);
plot(f1,x_fft),title("spectrum of the fundamental signal")
%% 
%full rectied signal :
clear all; close all; clc;

 t = 0:1e-5:1e-3;%time
 y = abs(0.2 * sin(2*pi*2e3*t));%equation
 plot(t,y),title("full wave rectifier ")
 
 H=25;%harmonics i choose 25
 A=2;%amplitude
 n=1:H

T=2;%the period
w0=2*pi/T;
%the coefficient
a0 = 2.*A./pi;
k=pi.*(4*n.^2-1);
an = -4.*A./k;
bn = 0;

figure;stem(an);title("The Harmonics")

%t2 = linspace(0,pi,512);
t2=[0:0.01:6];
z=[];
figure;
hold on;
for n=1:H;
    
    harmonic=a0+an(n)*cos(n*w0*t2);
    Z(n,:)= bn+harmonic';
   
end
plot( Z');title("The accumulated signal");

%3D

figure;

Rcarre=cumsum(Z);
mesh(Rcarre);title("3D ")

figure;
plot(Rcarre(n,:));title("Reconstructed signal")

%FFT:
Fs=1000;
N=4;
%frequency axis:
f0 = (0:n-1)*(Fs/n-1);
%Module calculation:anlysis fft of the fundamental signal
x_fft=abs(fft(y));
v=max(x_fft);x_fft = abs(fft(y))./v;

%let's the fft of polluted signal:
xn_fft=abs(fft(Rcarre(n,:)));
vv=max(x_fft);xn_fft=abs(fft(Rcarre(n,:)))./vv;

%to plot fft:
noSamples = length(t2);
f2 = 0:Fs/noSamples:Fs - Fs/noSamples;

noSamples = length(t);
f1 = 0:Fs/noSamples:Fs - Fs/noSamples;

figure;
subplot(2,1,1);
plot(f2,xn_fft),title("spectrum of the reconstructued signal")

subplot(2,1,2);
plot(f1,x_fft),title("spectrum of the fundamental signal")

%%
%Half wave rectification:
clear all; close all; clc;

l=linspace(0,10,100);

sig=sin(2*pi*50*l);


for t=1:100

if sin(2*pi*50*l(t))<=0

sig(t)=0;

else

sig(t) = sin(2*pi*50*l(t));

end

end

plot(sig);title("The half rectified signal")





