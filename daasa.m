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


—-------------------------------------------------------------------------------------------


%convulation
%x(n)*h(n) = y(n)
%(h(n) = impulse
clc;
clear all;
x=[1 2 3 4 0];
h=[0 0 1 5 2];
lenx=length(x);
lenh=length(h);
X=[x , zeros(1,lenx)];
H=[h , zeros(1,lenh)];
for i=1:lenx+lenh-1
    Y(i)=0;
    for j=1:lenx
        if(i-j+1>0)
            Y(i)=Y(i)+X(j)*H(i-j+1);
        else;
        end
    end
end

subplot(3,1,1);
stem(x);
title('x(n)');
subplot(3,1,2);
stem(h);
title('h(n)');
subplot(3,1,3);
stem(Y);
title('Y(n)');

%...........................................
%corelation.......

clc;
clear all;
x=[1 2 3 4 0];
h=[0 0 1 5 2];
%flip h..........
h= flip(h);
lenx=length(x);
lenh=length(h);
X=[x , zeros(1,lenx)];
H=[h , zeros(1,lenh)];
for i=1:lenx+lenh-1
    Y(i)=0;
    for j=1:lenx
        if(i-j+1>0)
            Y(i)=Y(i)+X(j)*H(i-j+1);
        else;
        end
    end
end

subplot(3,1,1);
stem(x);
title('x(n)');
subplot(3,1,2);
stem(h);
title('h(n)');
subplot(3,1,3);
stem(Y);
title('Y(n)');


%autocorelation..................

clc;
clear all;
x=[1 2 3 4 0];
h=x;
%h = x..........
h= flip(h);
lenx=length(x);
lenh=length(h);
X=[x , zeros(1,lenx)];
H=[h , zeros(1,lenh)];
for i=1:lenx+lenh-1
    Y(i)=0;
    for j=1:lenx
        if(i-j+1>0)
            Y(i)=Y(i)+X(j)*H(i-j+1);
        else;
        end
    end
end

subplot(3,1,1);
stem(x);
title('x(n)');
subplot(3,1,2);
stem(h);
title('h(n)');
subplot(3,1,3);
stem(Y);
title('Y(n)');



_________________________________________________________

%DFT...............


clc;
clear al;
x=[0.3535 0.3535 0.6464 1.0607 0.3535 -1.0607 -1.3535 -0.3535]
N=length(x);
X=zeros(N,1);
for m=1:N
   for n=1:N
       X(m)=X(m) + x(n)*exp((-2j*pi*(n-1)*(m-1))/N);
   end
end
h=0:N-1;
disp(X);
subplot(3,1,1);
plot(h,abs(x));
title('X(n)');
subplot(3,1,2);
plot(h,abs(X));
title('Amplitude Spectrum');
subplot(3,1,3);
plot(h,angle(X)*180/pi)
title('Phase Spectrum');

