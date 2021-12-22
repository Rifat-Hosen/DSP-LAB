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