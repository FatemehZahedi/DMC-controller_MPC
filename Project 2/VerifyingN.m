%%................................................................ ........
%..........................verifying parameter.............................
% verify N
% N=40
clear
 clc
[n1,d1,n2,d2]=Inputsys(1);
Gs1 = tf(n1,d1);
Ts=0.1;
Gd1 = c2d(Gs1,Ts,'zoh');
[num1,den1]=tfdata(Gd1,'v');
Gs2 = tf(n2,d2);
Gd2 = c2d(Gs2,Ts,'zoh');
[num2,den2]=tfdata(Gd2,'v');
sys_info = stepinfo(Gd1);
ts1 = sys_info.SettlingTime;
tr1=sys_info.RiseTime; 
sys_info = stepinfo(Gd2);
ts2 = sys_info.SettlingTime;
tr2=sys_info.RiseTime; 
%...................................................
t=1:Ts:10;
[g1,t1] = step(Gd1,t);
[g2,t2] = step(Gd2,t);
P1=floor(tr1/Ts);
P2=floor(tr2/Ts);
N1=floor( ts1/Ts);
N2=floor( ts2/Ts);
P=P2;
N=40;
M=P;
%.....................Toeplitz Matrix.................................
b1 = zeros(1,P); b1(1,1)= g1(2);
a1 = g1(2:P+1);
G1 = toeplitz(a1,b1);
G1(:,M) = G1(:,M:P)*ones(P-M+1,1);
G1 = G1(:,1:M);
%........................................................
b2 = zeros(1,P); b2(1,1)= g2(2);
a2 = g2(2:P+1);
G2 = toeplitz(a2,b2);
G2(:,M) = G2(:,M:P)*ones(P-M+1,1);
G2 = G2(:,1:M);
G=[G1 G2];
%....................Hankel Matrix.....................................
c1= [g1(3:P+2)];
r1 = [(g1(P+2:N+1))' zeros(1,P-1)];
G1_ = hankel(c1,r1);
%..........................................................
c2 = [g2(3:P+2)];
r2 = [(g2(P+2:N+1))' zeros(1,P-1)];
G2_ = hankel(c2,r2);
G_=[G1_ G2_];
%....................................................................................
%%....................... Designing........................................
gamma =1;
gain_DC=(num1(1)+num1(2)+num1(3))/(den1(1)+den1(2)+den1(3));
gain_DC2=(num2(1)+num2(2)+num2(3))/(den2(1)+den2(2)+den2(3));
Q = eye(P);
R1 =((1.2)^2)*gamma*gain_DC^2*eye(M);
R2=gamma*gain_DC2^2*eye(M);
R=[R1 zeros(M); zeros(M) R2];
Kdmc=(G'*Q*G+R)\(G'*Q);
alpha=0.5;
x01=0.0882;
x02=441.2;
%..........................................................................................
U1_ = zeros(P,length(t));
U2_ = zeros(P,length(t));
dU1_=zeros(N-1,length(t));
dU2_=zeros(N-1,length(t));
dU_=[dU1_;dU2_];
% U_=[U1_ ; U2_];
Ud1=zeros(N+P-1,length(t));
Ud2=zeros(N+P-1,length(t));
d=zeros(1,length(t));
%y1=0; %linear
y1=441.2;
u_1=[];
u_2=[];
ym=[];
y=0;
Y_d=zeros(P,length(t));
Y_past=zeros(P,length(t));
Y_m=zeros(P,length(t));
D=zeros(P,length(t));
E=zeros(P,length(t));
U1=zeros(M,length(t));
U2=zeros(M,length(t));
% U=[U1;U2];
dU1=zeros(M,length(t));
dU2=zeros(M,length(t));
dU=[dU1;dU2];
%..................step...........................
r =ones(length(t),1);
for i=1:length(t)-1 
    
for j=1:P
  Y_d(j,i+1)=(alpha^j)*y+(1-(alpha)^j)*r(i+1); % Programmed
end 
    
Y_past(:,i+1)=G_*dU_(:,i+1)+g1(N+1)*U1_(:,i+1)+g2(N+1)*U2_(:,i+1);
D(:,i+1)=d(i+1)*ones(P,1);
    
E(:,i+1)=Y_d(:,i+1)-Y_past(:,i+1)-D(:,i+1);
dU(:,i+1)=Kdmc*E(:,i+1);
dU1(:,i+1)=dU(1:M,i+1);
dU2(:,i+1)=dU(M+1:2*M,i+1);
U1(1,i+1)=dU1(1,i+1)+U1(1,i);
U2(1,i+1)=dU2(1,i+1)+U2(1,i);
dU(:,i+1)=[dU1(:,i+1);dU2(:,i+1)];
    
Y_m(:,i+1)=G*dU(:,i+1)+Y_past(:,i+1);

Ud1(2:N+P-1,i+2)=Ud1(1:N+P-2,i+1);
Ud1(1,i+2)=U1(1,i+1);
U1_(:,i+2)=Ud1(N:N+P-1,i+2);

Ud2(2:N+P-1,i+2)=Ud2(1:N+P-2,i+1);
Ud2(1,i+2)=U2(1,i+1);
U2_(:,i+2)=Ud2(N:N+P-1,i+2);


dU1_(2:N-1,i+2) = dU1_(1:N-2,i+1);
dU1_(1,i+2)=dU1(1,i+1);
dU2_(2:N-1,i+2) = dU2_(1:N-2,i+1);
dU2_(1,i+2)=dU2(1,i+1);


dU_(:,i+2)=[dU1_(:,i+2);dU2_(:,i+2)];

u1=U1(1,i+1);
u2=U2(1,i+1);
sim('Model')
%d(i+1)=yl(end)-Y_m(1,i); %linear
d(i+1)=y(end)-Y_m(1,i);
%y=yl(end); % linear
y=y(end);%+dist(i,1);    % nonlinear
%y1=[y1;yl(end)];  % linear
y1=[y1; y+441.2];
ym=[ym; Y_m(1,i)];
u_1=[u_1; u1];
u_2=[u_2; u2];
x01=x1(end);
x02=x2(end);
end
figure(1);
plot(y1,'b');
grid on
title('Response of the nonlinear system');
xlabel('sample');
figure(2);
plot(u_1,'b');
grid on
xlabel('sample');
title('Control law for input 1 without bias');
figure(3);
plot(u_2,'b');
grid on
xlabel('sample');
title('Control law for input 2 without bias');
%................................................
% N=17
clear
 clc
[n1,d1,n2,d2]=Inputsys(1);
Gs1 = tf(n1,d1);
Ts=0.1;
Gd1 = c2d(Gs1,Ts,'zoh');
[num1,den1]=tfdata(Gd1,'v');
Gs2 = tf(n2,d2);
Gd2 = c2d(Gs2,Ts,'zoh');
[num2,den2]=tfdata(Gd2,'v');
sys_info = stepinfo(Gd1);
ts1 = sys_info.SettlingTime;
tr1=sys_info.RiseTime; 
sys_info = stepinfo(Gd2);
ts2 = sys_info.SettlingTime;
tr2=sys_info.RiseTime; 
%...................................................
t=1:Ts:10;
[g1,t1] = step(Gd1,t);
[g2,t2] = step(Gd2,t);
P1=floor(tr1/Ts);
P2=floor(tr2/Ts);
N1=floor( ts1/Ts);
N2=floor( ts2/Ts);
P=P2;
N=N1;
M=P;
%.....................Toeplitz Matrix.................................
b1 = zeros(1,P); b1(1,1)= g1(2);
a1 = g1(2:P+1);
G1 = toeplitz(a1,b1);
G1(:,M) = G1(:,M:P)*ones(P-M+1,1);
G1 = G1(:,1:M);
%........................................................
b2 = zeros(1,P); b2(1,1)= g2(2);
a2 = g2(2:P+1);
G2 = toeplitz(a2,b2);
G2(:,M) = G2(:,M:P)*ones(P-M+1,1);
G2 = G2(:,1:M);
G=[G1 G2];
%....................Hankel Matrix.....................................
c1= [g1(3:P+2)];
r1 = [(g1(P+2:N+1))' zeros(1,P-1)];
G1_ = hankel(c1,r1);
%..........................................................
c2 = [g2(3:P+2)];
r2 = [(g2(P+2:N+1))' zeros(1,P-1)];
G2_ = hankel(c2,r2);
G_=[G1_ G2_];
gamma =1;
gain_DC=(num1(1)+num1(2)+num1(3))/(den1(1)+den1(2)+den1(3));
gain_DC2=(num2(1)+num2(2)+num2(3))/(den2(1)+den2(2)+den2(3));
Q = eye(P);
R1 =((1.2)^2)*gamma*gain_DC^2*eye(M);
R2=gamma*gain_DC2^2*eye(M);
R=[R1 zeros(M); zeros(M) R2];
Kdmc=(G'*Q*G+R)\(G'*Q);
alpha=0.5;
x01=0.0882;
x02=441.2;
%..........................................................................................
U1_ = zeros(P,length(t));
U2_ = zeros(P,length(t));
dU1_=zeros(N-1,length(t));
dU2_=zeros(N-1,length(t));
dU_=[dU1_;dU2_];
% U_=[U1_ ; U2_];
Ud1=zeros(N+P-1,length(t));
Ud2=zeros(N+P-1,length(t));
d=zeros(1,length(t));
%y1=0; %linear
y1=441.2;
u_1=[];
u_2=[];
ym=[];
y=0;
Y_d=zeros(P,length(t));
Y_past=zeros(P,length(t));
Y_m=zeros(P,length(t));
D=zeros(P,length(t));
E=zeros(P,length(t));
U1=zeros(M,length(t));
U2=zeros(M,length(t));
% U=[U1;U2];
dU1=zeros(M,length(t));
dU2=zeros(M,length(t));
dU=[dU1;dU2];
%..................step...........................
r =ones(length(t),1);
for i=1:length(t)-1 
    
for j=1:P
  Y_d(j,i+1)=(alpha^j)*y+(1-(alpha)^j)*r(i+1); % Programmed
end 
    
Y_past(:,i+1)=G_*dU_(:,i+1)+g1(N+1)*U1_(:,i+1)+g2(N+1)*U2_(:,i+1);
D(:,i+1)=d(i+1)*ones(P,1);
    
E(:,i+1)=Y_d(:,i+1)-Y_past(:,i+1)-D(:,i+1);
dU(:,i+1)=Kdmc*E(:,i+1);
dU1(:,i+1)=dU(1:M,i+1);
dU2(:,i+1)=dU(M+1:2*M,i+1);
U1(1,i+1)=dU1(1,i+1)+U1(1,i);
U2(1,i+1)=dU2(1,i+1)+U2(1,i);
dU(:,i+1)=[dU1(:,i+1);dU2(:,i+1)];
    
Y_m(:,i+1)=G*dU(:,i+1)+Y_past(:,i+1);

Ud1(2:N+P-1,i+2)=Ud1(1:N+P-2,i+1);
Ud1(1,i+2)=U1(1,i+1);
U1_(:,i+2)=Ud1(N:N+P-1,i+2);

Ud2(2:N+P-1,i+2)=Ud2(1:N+P-2,i+1);
Ud2(1,i+2)=U2(1,i+1);
U2_(:,i+2)=Ud2(N:N+P-1,i+2);


dU1_(2:N-1,i+2) = dU1_(1:N-2,i+1);
dU1_(1,i+2)=dU1(1,i+1);
dU2_(2:N-1,i+2) = dU2_(1:N-2,i+1);
dU2_(1,i+2)=dU2(1,i+1);


dU_(:,i+2)=[dU1_(:,i+2);dU2_(:,i+2)];

u1=U1(1,i+1);
u2=U2(1,i+1);
sim('Model')
%d(i+1)=yl(end)-Y_m(1,i); %linear
d(i+1)=y(end)-Y_m(1,i);
%y=yl(end); % linear
y=y(end);%+dist(i,1);    % nonlinear
%y1=[y1;yl(end)];  % linear
y1=[y1; y+441.2];
ym=[ym; Y_m(1,i)];
u_1=[u_1; u1];
u_2=[u_2; u2];
x01=x1(end);
x02=x2(end);
end
figure(1);
hold on
plot(y1,'m');
figure(2);
hold on
plot(u_1,'m');
figure(3);
hold on
plot(u_2,'m');
%.....................................................................
% N=5
clear
 clc
[n1,d1,n2,d2]=Inputsys(1);
Gs1 = tf(n1,d1);
Ts=0.1;
Gd1 = c2d(Gs1,Ts,'zoh');
[num1,den1]=tfdata(Gd1,'v');
Gs2 = tf(n2,d2);
Gd2 = c2d(Gs2,Ts,'zoh');
[num2,den2]=tfdata(Gd2,'v');
sys_info = stepinfo(Gd1);
ts1 = sys_info.SettlingTime;
tr1=sys_info.RiseTime; 
sys_info = stepinfo(Gd2);
ts2 = sys_info.SettlingTime;
tr2=sys_info.RiseTime; 
%...................................................
t=1:Ts:10;
[g1,t1] = step(Gd1,t);
[g2,t2] = step(Gd2,t);
P1=floor(tr1/Ts);
P2=floor(tr2/Ts);
N1=floor( ts1/Ts);
N2=floor( ts2/Ts);
P=P2;
N=5;
M=P;
%.....................Toeplitz Matrix.................................
b1 = zeros(1,P); b1(1,1)= g1(2);
a1 = g1(2:P+1);
G1 = toeplitz(a1,b1);
G1(:,M) = G1(:,M:P)*ones(P-M+1,1);
G1 = G1(:,1:M);
%........................................................
b2 = zeros(1,P); b2(1,1)= g2(2);
a2 = g2(2:P+1);
G2 = toeplitz(a2,b2);
G2(:,M) = G2(:,M:P)*ones(P-M+1,1);
G2 = G2(:,1:M);
G=[G1 G2];
%....................Hankel Matrix.....................................
c1= [g1(3:P+2)];
r1 = [(g1(P+2:N+1))' zeros(1,P-1)];
G1_ = hankel(c1,r1);
%..........................................................
c2 = [g2(3:P+2)];
r2 = [(g2(P+2:N+1))' zeros(1,P-1)];
G2_ = hankel(c2,r2);
G_=[G1_ G2_];
gamma =1;
gain_DC=(num1(1)+num1(2)+num1(3))/(den1(1)+den1(2)+den1(3));
gain_DC2=(num2(1)+num2(2)+num2(3))/(den2(1)+den2(2)+den2(3));
Q = eye(P);
R1 =((1.2)^2)*gamma*gain_DC^2*eye(M);
R2=gamma*gain_DC2^2*eye(M);
R=[R1 zeros(M); zeros(M) R2];
Kdmc=(G'*Q*G+R)\(G'*Q);
alpha=0.5;
x01=0.0882;
x02=441.2;
%..........................................................................................
U1_ = zeros(P,length(t));
U2_ = zeros(P,length(t));
dU1_=zeros(N-1,length(t));
dU2_=zeros(N-1,length(t));
dU_=[dU1_;dU2_];
% U_=[U1_ ; U2_];
Ud1=zeros(N+P-1,length(t));
Ud2=zeros(N+P-1,length(t));
d=zeros(1,length(t));
%y1=0; %linear
y1=441.2;
u_1=[];
u_2=[];
ym=[];
y=0;
Y_d=zeros(P,length(t));
Y_past=zeros(P,length(t));
Y_m=zeros(P,length(t));
D=zeros(P,length(t));
E=zeros(P,length(t));
U1=zeros(M,length(t));
U2=zeros(M,length(t));
% U=[U1;U2];
dU1=zeros(M,length(t));
dU2=zeros(M,length(t));
dU=[dU1;dU2];
%..................step...........................
r =ones(length(t),1);
for i=1:length(t)-1 
    
for j=1:P
  Y_d(j,i+1)=(alpha^j)*y+(1-(alpha)^j)*r(i+1); % Programmed
end 
    
Y_past(:,i+1)=G_*dU_(:,i+1)+g1(N+1)*U1_(:,i+1)+g2(N+1)*U2_(:,i+1);
D(:,i+1)=d(i+1)*ones(P,1);
    
E(:,i+1)=Y_d(:,i+1)-Y_past(:,i+1)-D(:,i+1);
dU(:,i+1)=Kdmc*E(:,i+1);
dU1(:,i+1)=dU(1:M,i+1);
dU2(:,i+1)=dU(M+1:2*M,i+1);
U1(1,i+1)=dU1(1,i+1)+U1(1,i);
U2(1,i+1)=dU2(1,i+1)+U2(1,i);
dU(:,i+1)=[dU1(:,i+1);dU2(:,i+1)];
    
Y_m(:,i+1)=G*dU(:,i+1)+Y_past(:,i+1);

Ud1(2:N+P-1,i+2)=Ud1(1:N+P-2,i+1);
Ud1(1,i+2)=U1(1,i+1);
U1_(:,i+2)=Ud1(N:N+P-1,i+2);

Ud2(2:N+P-1,i+2)=Ud2(1:N+P-2,i+1);
Ud2(1,i+2)=U2(1,i+1);
U2_(:,i+2)=Ud2(N:N+P-1,i+2);


dU1_(2:N-1,i+2) = dU1_(1:N-2,i+1);
dU1_(1,i+2)=dU1(1,i+1);
dU2_(2:N-1,i+2) = dU2_(1:N-2,i+1);
dU2_(1,i+2)=dU2(1,i+1);


dU_(:,i+2)=[dU1_(:,i+2);dU2_(:,i+2)];

u1=U1(1,i+1);
u2=U2(1,i+1);
sim('Model')
%d(i+1)=yl(end)-Y_m(1,i); %linear
d(i+1)=y(end)-Y_m(1,i);
%y=yl(end); % linear
y=y(end);%+dist(i,1);    % nonlinear
%y1=[y1;yl(end)];  % linear
y1=[y1; y+441.2];
ym=[ym; Y_m(1,i)];
u_1=[u_1; u1];
u_2=[u_2; u2];
x01=x1(end);
x02=x2(end);
end
figure(1);
hold on
plot(y1,'c');
hold on
plot(r+441.22,'r')
legend('N=40','N=17','N=5','r');
figure(2);
hold on
plot(u_1,'c');
legend('N=40','N=17','N=5');
figure(3);
hold on
plot(u_2,'c');
legend('N=40','N=17','N=5');
%............................................................
