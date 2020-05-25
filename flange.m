function [y]=flange(fs, v, x, r)
%[y] = flange(fs, v, x, r)
%

%Coded by: Stephen G. McGovern, date: 08.03.03
md= ceil(v*fs);
n=1:length(x)+md;
v=round(v*fs);
z=zeros(md,1);
m=max(abs(x));
x=[z;x;z];
rr=2*pi/round(fs*r);
b=round((v/2)*(1-cos(rr.*n)));
x=x(n+md)+x(n+md-b);
m=m/max(abs(x));
y=m*x;