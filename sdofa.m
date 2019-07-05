
clear;
m=1;
M=m;
g=9.8;
K=10^2;
k=K;
C=0;
u=0.3;
h=1*10^-3;
T=1;
n=floor(T/h)+1;
t=zeros(n,1);
r=10;
%广义a方法参数
rho=0.8;
alpha=(2*rho-1)/(rho+1);
beta=1/(1+rho)^2;
gamma=(3-rho)/(2*(1+rho));
delta=rho/(rho+1); 
%TR方法
% alpha=0;
% beta=0.25;
% gamma=0.5;
% delta=0;
for i=1:n
    t(i,1)=(i-1)*h;
end
qd=zeros(n,1);
qv=zeros(n,1);
qa=zeros(n,1);
qf=zeros(n,1);
%初始条件
qv(1,1)=1.0;
qf(1,1)=-u*m*g;
qa(1,1)=-u*g;
R_C=(1-alpha)/(gamma*h).*M+(1-delta)*beta/gamma*h.*K+(1-delta).*C;
iR_C=inv(R_C);
for i=2:n%i-1步到i步
    %迭代初始值
    qd(i,1)=qd(i-1,1);
    qv(i,1)=qv(i-1,1);
    qa(i,1)=qa(i-1,1);
    qf(i,1)=qf(i-1,1);
    
    %迭代过程中不变的量
    W_T=1.0;
    R_1=delta.*(-K*qd(i-1,1)-C*qv(i-1,1)+W_T*qf(i-1,1))-alpha.*M*qa(i-1,1);
    R_2=1/(gamma*h).*qv(i-1,1)+(1/gamma-1).*qa(i-1,1);
    R_3=R_1+(1-alpha)*M*R_2-(1-delta)*K*(qd(i-1,1)+...
        (1-beta/gamma)*h.*qv(i-1,1)+(0.5-beta/gamma)*h*h.*qa(i-1,1));
    error=1.0;
    while (error>10^(-15))
        qv(i,1)=(1-delta).*iR_C*W_T*qf(i,1)+iR_C*R_3;
        g_T=W_T'*qv(i,1);
        Projc=projc(qf(i,1)-r.*g_T,m*g,u);
        dProjc=dprojc(qf(i,1)-r.*g_T,m*g,u);
        H=qf(i,1)-Projc;
        D=r*(1-delta)*dProjc*W_T'*iR_C*W_T+1-dProjc;
        qf(i,1)=qf(i,1)-inv(D)*H;
        error=abs(H);              
    end
    qv(i,1)=(1-delta).*iR_C*W_T*qf(i,1)+iR_C*R_3;
    qd(i,1)=qd(i-1,1)+beta*h./gamma*qv(i,1)+...
        (1-beta./gamma).*h.*qv(i-1,1)+(1/2-beta/gamma)*h*h*qa(i-1,1);
    qa(i,1)=1/(gamma*h).*(qv(i,1)-qv(i-1,1))-(1/gamma-1).*qa(i-1,1);
%     disp(t(i,1));
end

%理论解
ed=zeros(n,1);
ev=zeros(n,1);
ea=zeros(n,1);
ef=zeros(n,1);
for i=1:n
    A1=1/10;B1=u*m*g/k;C1=-B1;
    tn1=atan(A1/B1)/10;
    if(t(i,1)<tn1)
        ed(i,1)=A1*sin(10*t(i,1))+B1*cos(10*t(i,1))+C1;
        ev(i,1)=10*A1*cos(10*t(i,1))-10*B1*sin(10*t(i,1));
        ea(i,1)=-100*A1*sin(10*t(i,1))-100*B1*cos(10*t(i,1));
        ef(i,1)=-u*m*g;
        x1=A1*sin(10*tn1)+B1*cos(10*tn1)+C1;
        C2=u*m*g/k;
        A2=x1-C2;
        tn2=pi/10+tn1;
    else if(t(i,1)<tn2)
        ed(i,1)=A2*cos(10*t(i,1)-10*tn1)+C2;
        ev(i,1)=-10*A2*sin(10*t(i,1)-10*tn1);
        ea(i,1)=-100*A2*cos(10*t(i,1)-10*tn1);
        ef(i,1)=u*m*g;
        else
        ed(i,1)=A2*cos(10*tn2-10*tn1)+C2;
        ev(i,1)=0.0;
        ea(i,1)=0.0;
        ef(i,1)=0.0;    
        end
    end  
end
% for i=1:n
%     if (t(i,1)<atan(450/147)/100)
%     ed(i,1)=9/10000*sin(100*t(i,1))+147/500000*cos(100*t(i,1))-147/500000;
%     ev(i,1)=9/100*cos(100*t(i,1))-147/5000*sin(100*t(i,1));
%     ea(i,1)=-9*sin(100*t(i,1))-147/50*cos(100*t(i,1));
%     ef(i,1)=-u*m*g;
%     x1=9/10000*sin(100*(atan(450/147)/100))+147/500000*cos(100*(atan(450/147)/100))-147/500000-(147/500000);
%     else
%         if (t(i,1)<(atan(450/147)+pi)/100)
% %             ed(i,1)=147/500000+x1*cos(100*(t(i,1)-atan(450/147)/100));
% %             ev(i,1)=-100*x1*sin(100*(t(i,1)-atan(450/147)/100));
% %             ea(i,1)=-10000*x1*cos(100*(t(i,1)-atan(450/147)/100));
%             ed(i,1)=147/500000+cos(100*(t(i,1)-atan(450/147)/100))*(24733551/68933500000);
%             ev(i,1)=-sin(100*(t(i,1)-atan(450/147)/100))*(24733551/689335000);
%             ea(i,1)=-cos(100*(t(i,1)-atan(450/147)/100))*(24733551/6893350);
% %             ed(i,1)=147/500000+cos(100*(t(i,1)-atan(450/147)/100))*(1119/3125000);
% %             ev(i,1)=-sin(100*(t(i,1)-atan(450/147)/100))*(1119/31250);
% %             ea(i,1)=-cos(100*(t(i,1)-atan(450/147)/100))*(11190/3125);
%             ef(i,1)=u*m*g;
%         else
% %             ed(i,1)=147/500000-x1;
%             ed(i,1)=147/500000-24733551/68933500000;
% %             ed(i,1)=147/500000-1119/3125000;
%             ev(i,1)=0;
%             ea(i,1)=0;
%             ef(i,1)=k*ed(i,1);
%         end
%     end
% end
%计算误差
errd1=0;errd2=0;
errv1=0;errv2=0;
erra1=0;erra2=0;
for i=1:floor(n/10):n
    errd1=errd1+abs(qd(i,1)-ed(i,1));
    errd2=errd2+abs(ed(i,1));
    errv1=errv1+abs(qv(i,1)-ev(i,1));
    errv2=errv2+abs(ev(i,1));
    erra1=erra1+abs(qa(i,1)-ea(i,1));
    erra2=erra2+abs(ea(i,1));
end
errd=errd1/errd2
errv=errv1/errv2
erra=erra1/erra2
% ----------------------画图----------------------
figure;
grid on
hold on 
plot(t,ev(:,1),'k-','linewidth',2); % 接触点1相对切向速度      
xlabel('$t$ [s]'); ylabel('$\dot{g}_T$ [m/s]');

figure;
grid on
hold on 
plot(t,qd(:,1),'k-','linewidth',2); % 接触点1位移     
xlabel('$t$ [s]'); ylabel('$x$ [m/s]');

figure;
grid on
hold on 
plot(t,qf(:,1),'k-','linewidth',2); % 接触点1摩擦力      
xlabel('$t$ [s]'); ylabel('$F_f$ [N]');

figure;
grid on
hold on 
plot(t,qa(:,1),'k-','linewidth',2); % 接触点1加速度     
xlabel('$t$ [s]'); ylabel('$\ddot{g}_T$ [m/s]');
