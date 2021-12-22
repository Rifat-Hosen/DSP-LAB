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



