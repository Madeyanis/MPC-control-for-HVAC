function [A, B] = StateJcbFcnBc(x ,u, ts)
mpump = 50; 
Volume = 84;
DensiteAir = 1.204;
e = 0.355;
k = 2781;
s1 = (6*4) * 2;
s2 = (6*3.5) * 2;
s3 = (4*3.5) * 2;
s = s1 + s2 + s3;
R = e * k/s ;
M = Volume * DensiteAir;
c = 1256;
sfe = 1.5;
Q = 100;
Az  = (-(1/M)*mpump)-1/(M*c*R);
theta = 0.9;
a = 0.9;
beta = 0.9;
sigma = 0.8;

A = [Az (mpump/M) 0;
       theta a 0;
       0 0 1];
   
B = [0 0 0 (1/M*c*R) (sfe/M*c) (Q/M*c);
       beta 0 0 sigma 0 0;
       0 ts 0 0 0 0];
  

end

