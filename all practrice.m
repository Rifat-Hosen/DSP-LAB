%unit Step
n = -5:1:10;
u = (n>=0);
subplot(4,4,1);
stem(n,u);
ylabel("u(t)");
xlabel("t");
xlim([-6 12]);
ylim([-.2 1.2]);
title("Unit Step");

%unit Impulse ut=1, t=0 & ut=0 for else

ui = (n==0)
subplot(4,4,2);
stem(n,ui);
ylabel("ui(t)");
xlabel("t");
xlim([-6 12]);
ylim([-.2 1.2]);
title("Unit Impulse");


%Ramp signal is denoted by r(t), and it is defined as r(t) = t when t >= 0 else 0

r = (n>=0).*n;
subplot(4,4,3);
stem(n,r);
ylabel("Amplitude");
xlabel("ramp signal");
xlim([-6 12]);
ylim([-1 12]);
title("Ramp Signal");


%expotential 
x= 0:1:10
a =0.8;
x2 = a.^x;
subplot(4,4,4);
stem(x,x2);
ylabel("x2(n)");
xlabel("n");
title("Expotential");

%Sinusoidal sequence:-
t=0:0.01:pi;
y=sin(2*pi*t);
subplot(4,4,5);
plot(t,y);
ylabel('Amplitude');
xlabel('e');
title('Sinusoidal Sequence');

% Cosine Sequence:-
t=0:0.01:pi;
y=cos(2*pi*t);
subplot(4,4,6);
plot(t,y);
ylabel('Amplitude');
xlabel('f');
title('Cosine Sequence');



%*****************************************************************************

%Allisins effect........................................
%x = 5*sin(2*pi*1000*t+pi/2);

a = 2;
f = 1000;
t = linspace(0,.01,1000)
x = a*sin(2*pi*f*t);
subplot(4,4,7)
plot(t,x)
xlim([-.0001 .011])

%fs = f frequency..................
fs = 800;
T = 1/fs;
nmin = ceil(0/T);
nmax = floor(.01/T);
n = nmin:nmax;
x4 = a*sin(2*pi*f*n*T);
subplot(4,4,8)
plot(t,x)
hold on
plot(n*T,x4,'o')
xlim([-.0001 .011])
hold off

%2*f frequency...............
fs = 3000;
T = 1/fs;
nmin = ceil(0/T);
nmax = floor(.01/T);
n = nmin:nmax;
x3 = a*sin(2*pi*f*n*T);
subplot(4,4,9)
plot(t,x)
hold on
plot(n*T,x3,'o')
xlim([-.0001 .011])
hold off

%uppper frequency.........
fs = 8000;
T = 1/fs;
nmin = ceil(0/T);
nmax = floor(.01/T);
n = nmin:nmax;
x2 = a*sin(2*pi*f*n*T);
subplot(4,4,10)
plot(t,x)
hold on
plot(n*T,x2,'o')

xlim([-.0001 .011])
hold off




%TEST............
% n=0:1:30; fs=8000; f=f/fs; y=1*sin(2*pi*f*n); stem(n,y);







