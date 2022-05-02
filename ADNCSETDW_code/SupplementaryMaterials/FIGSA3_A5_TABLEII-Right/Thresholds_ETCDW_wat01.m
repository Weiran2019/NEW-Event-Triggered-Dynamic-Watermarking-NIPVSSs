%clear all; clc;
%%导入数据
C=[1,0,0,0;0,1,0,0];
B=[0.0001;0.0002;0.0100;0.0300];
len_vec=zeros(1,6);
trirate_vec=zeros(1,6);
%%第1组数据
Data = importdata('TVKF_wat01_del00001_ETCDW_nor1.txt');
len=length(Data.data(:,1));
len_vec(1)=len;
pos=Data.data(:,2);
ang=Data.data(:,3);
gamma=Data.data(:,4);
Im=Data.data(:,5);
res1=Data.data(:,6);
res2=Data.data(:,7);
Psi11=Data.data(:,8);
Psi12=Data.data(:,9);
Psi21=Data.data(:,10);
Psi22=Data.data(:,11);
wat=Data.data(:,12);
%%计算触发率
count=0;
for i=1:len
    if(gamma(i)==1)
        count=count+1;
    end
end
trirate=count/len*100;
trirate_vec(1)=trirate;
%%计算检测式
bR=zeros(2,2);
bD=zeros(2,1);
ind1=zeros(1,len);
ind2=zeros(1,len);
for i=2:len
    bR=bR+([res1(i)*res1(i),res1(i)*res2(i);res2(i)*res1(i),res2(i)*res2(i)]-[Psi11(i),Psi12(i);Psi21(i),Psi22(i)]);
    bD=bD+(wat(i-1)*[res1(i);res2(i)]+(1-gamma(i))*C*B*0.0001);
    ind1(i)=norm(bR/(i-1));
    ind2(i)=norm(bD/(i-1));
end

Fonts=18;
k=1:len;
% fig1=figure(1);
% plot(k,pos(k),'b','LineWidth',2);hold on;
% 
% fig2=figure(2);
% plot(k,ang(k),'b','LineWidth',2);hold on;

fig3=figure(3);
set(fig3, 'Position', [0 0 1000 350]); 
set(fig3, 'PaperSize', [29.7000 21.0000]); 
set(fig3,'PaperPosition',[1 1 21 7]);
semilogy(k,ind1(k),'b','LineWidth',2);hold on;
le11=xlabel('$i$');
le12=ylabel('$\frac{1}{i}{\tilde \mathcal{R}_{i}}$');
set(gca,'FontSize',Fonts);
set(le11,'Interpreter','latex');
set(le12,'Interpreter','latex');

fig4=figure(4);
set(fig4, 'Position', [0 0 1000 350]); 
set(fig4, 'PaperSize', [29.7000 21.0000]); 
set(fig4,'PaperPosition',[1 1 21 7]);
semilogy(k,ind2(k),'b','LineWidth',2);hold on;
le11=xlabel('$i$');
le12=ylabel('$\frac{1}{i}{\mathcal{D}_{i}}$');
set(gca,'FontSize',Fonts);
set(le11,'Interpreter','latex');
set(le12,'Interpreter','latex');

%%第2组数据
Data = importdata('TVKF_wat01_del00001_ETCDW_nor2.txt');
len=length(Data.data(:,1));
len_vec(2)=len;
pos=Data.data(:,2);
ang=Data.data(:,3);
gamma=Data.data(:,4);
Im=Data.data(:,5);
res1=Data.data(:,6);
res2=Data.data(:,7);
Psi11=Data.data(:,8);
Psi12=Data.data(:,9);
Psi21=Data.data(:,10);
Psi22=Data.data(:,11);
wat=Data.data(:,12);
%%计算触发率
count=0;
for i=1:len
    if(gamma(i)==1)
        count=count+1;
    end
end
trirate=count/len*100;
trirate_vec(2)=trirate;
%%计算检测式
bR=zeros(2,2);
bD=zeros(2,1);
ind1=zeros(1,len);
ind2=zeros(1,len);
for i=2:len
    bR=bR+([res1(i)*res1(i),res1(i)*res2(i);res2(i)*res1(i),res2(i)*res2(i)]-[Psi11(i),Psi12(i);Psi21(i),Psi22(i)]);
    bD=bD+(wat(i-1)*[res1(i);res2(i)]+(1-gamma(i))*C*B*0.0001);
    ind1(i)=norm(bR/(i-1));
    ind2(i)=norm(bD/(i-1));
end

k=1:len;
% figure(1);
% plot(k,pos(k),'r','LineWidth',2);hold on;
% 
% figure(2);
% plot(k,ang(k),'r','LineWidth',2);hold on;

figure(3);
semilogy(k,ind1(k),'r','LineWidth',2);hold on;

figure(4);
semilogy(k,ind2(k),'r','LineWidth',2);hold on;

%%第3组数据
Data = importdata('TVKF_wat01_del00001_ETCDW_nor3.txt');
len=length(Data.data(:,1));
len_vec(3)=len;
pos=Data.data(:,2);
ang=Data.data(:,3);
gamma=Data.data(:,4);
Im=Data.data(:,5);
res1=Data.data(:,6);
res2=Data.data(:,7);
Psi11=Data.data(:,8);
Psi12=Data.data(:,9);
Psi21=Data.data(:,10);
Psi22=Data.data(:,11);
wat=Data.data(:,12);
%%计算触发率
count=0;
for i=1:len
    if(gamma(i)==1)
        count=count+1;
    end
end
trirate=count/len*100;
trirate_vec(3)=trirate;
%%计算检测式
bR=zeros(2,2);
bD=zeros(2,1);
ind1=zeros(1,len);
ind2=zeros(1,len);
for i=2:len
    bR=bR+([res1(i)*res1(i),res1(i)*res2(i);res2(i)*res1(i),res2(i)*res2(i)]-[Psi11(i),Psi12(i);Psi21(i),Psi22(i)]);
    bD=bD+(wat(i-1)*[res1(i);res2(i)]+(1-gamma(i))*C*B*0.0001);
    ind1(i)=norm(bR/(i-1));
    ind2(i)=norm(bD/(i-1));
end

k=1:len;
% figure(1);
% plot(k,pos(k),'c','LineWidth',2);hold on;
% 
% figure(2);
% plot(k,ang(k),'c','LineWidth',2);hold on;

figure(3);
semilogy(k,ind1(k),'c','LineWidth',2);hold on;

figure(4);
semilogy(k,ind2(k),'c','LineWidth',2);hold on;

%%第4组数据
Data = importdata('TVKF_wat01_del00001_ETCDW_nor4.txt');
len=length(Data.data(:,1));
len_vec(4)=len;
pos=Data.data(:,2);
ang=Data.data(:,3);
gamma=Data.data(:,4);
Im=Data.data(:,5);
res1=Data.data(:,6);
res2=Data.data(:,7);
Psi11=Data.data(:,8);
Psi12=Data.data(:,9);
Psi21=Data.data(:,10);
Psi22=Data.data(:,11);
wat=Data.data(:,12);
%%计算触发率
count=0;
for i=1:len
    if(gamma(i)==1)
        count=count+1;
    end
end
trirate=count/len*100;
trirate_vec(4)=trirate;
%%计算检测式
bR=zeros(2,2);
bD=zeros(2,1);
ind1=zeros(1,len);
ind2=zeros(1,len);
for i=2:len
    bR=bR+([res1(i)*res1(i),res1(i)*res2(i);res2(i)*res1(i),res2(i)*res2(i)]-[Psi11(i),Psi12(i);Psi21(i),Psi22(i)]);
    bD=bD+(wat(i-1)*[res1(i);res2(i)]+(1-gamma(i))*C*B*0.0001);
    ind1(i)=norm(bR/(i-1));
    ind2(i)=norm(bD/(i-1));
end

k=1:len;
fig1=figure(1);
set(fig1, 'Position', [0 0 1000 350]); 
set(fig1, 'PaperSize', [29.7000 21.0000]); 
set(fig1,'PaperPosition',[1 1 21 7]);
plot(pos(1:500),'m','LineWidth',2);
le11=xlabel('$k$');
le12=ylabel('$\alpha$ (m)');
set(gca,'FontSize',Fonts);
set(le11,'Interpreter','latex');
set(le12,'Interpreter','latex');

fig2=figure(2);
set(fig2, 'Position', [0 0 1000 350]); 
set(fig2, 'PaperSize', [29.7000 21.0000]); 
set(fig2,'PaperPosition',[1 1 21 7]);
plot(ang(1:500),'m','LineWidth',2);
le11=xlabel('$k$');
le12=ylabel('$\theta$ (rad)');
set(gca,'FontSize',Fonts);
set(le11,'Interpreter','latex');
set(le12,'Interpreter','latex');

fig5=figure(5);
set(fig5, 'Position', [0 0 1000 350]); 
set(fig5, 'PaperSize', [29.7000 21.0000]); 
set(fig5,'PaperPosition',[1 1 21 7]);
stairs(gamma(1:500),'m','LineWidth',2);
le11=xlabel('$k$');
le12=ylabel('$\gamma(k)$');
ylim([-1,2]);
set(gca,'FontSize',Fonts);
set(le11,'Interpreter','latex');
set(le12,'Interpreter','latex');

figure(3);
semilogy(k,ind1(k),'m','LineWidth',2);hold on;

figure(4);
semilogy(k,ind2(k),'m','LineWidth',2);hold on;

%%第5组数据
Data = importdata('TVKF_wat01_del00001_ETCDW_nor5.txt');
len=length(Data.data(:,1));
len_vec(5)=len;
pos=Data.data(:,2);
ang=Data.data(:,3);
gamma=Data.data(:,4);
Im=Data.data(:,5);
res1=Data.data(:,6);
res2=Data.data(:,7);
Psi11=Data.data(:,8);
Psi12=Data.data(:,9);
Psi21=Data.data(:,10);
Psi22=Data.data(:,11);
wat=Data.data(:,12);
%%计算触发率
count=0;
for i=1:len
    if(gamma(i)==1)
        count=count+1;
    end
end
trirate=count/len*100;
trirate_vec(5)=trirate;
%%计算检测式
bR=zeros(2,2);
bD=zeros(2,1);
ind1=zeros(1,len);
ind2=zeros(1,len);
for i=2:len
    bR=bR+([res1(i)*res1(i),res1(i)*res2(i);res2(i)*res1(i),res2(i)*res2(i)]-[Psi11(i),Psi12(i);Psi21(i),Psi22(i)]);
    bD=bD+(wat(i-1)*[res1(i);res2(i)]+(1-gamma(i))*C*B*0.0001);
    ind1(i)=norm(bR/(i-1));
    ind2(i)=norm(bD/(i-1));
end

k=1:len;
% figure(1);
% plot(k,pos(k),'y','LineWidth',2);hold on;
% 
% figure(2);
% plot(k,ang(k),'y','LineWidth',2);hold on;

figure(3);
semilogy(k,ind1(k),'y','LineWidth',2);hold on;

figure(4);
semilogy(k,ind2(k),'y','LineWidth',2);hold on;

%%第6组数据
Data = importdata('TVKF_wat01_del00001_ETCDW_nor6.txt');
len=length(Data.data(:,1));
len_vec(6)=len;
pos=Data.data(:,2);
ang=Data.data(:,3);
gamma=Data.data(:,4);
Im=Data.data(:,5);
res1=Data.data(:,6);
res2=Data.data(:,7);
Psi11=Data.data(:,8);
Psi12=Data.data(:,9);
Psi21=Data.data(:,10);
Psi22=Data.data(:,11);
wat=Data.data(:,12);
%%计算触发率
count=0;
for i=1:len
    if(gamma(i)==1)
        count=count+1;
    end
end
trirate=count/len*100;
trirate_vec(6)=trirate;
%%计算检测式
bR=zeros(2,2);
bD=zeros(2,1);
ind1=zeros(1,len);
ind2=zeros(1,len);
for i=2:len
    bR=bR+([res1(i)*res1(i),res1(i)*res2(i);res2(i)*res1(i),res2(i)*res2(i)]-[Psi11(i),Psi12(i);Psi21(i),Psi22(i)]);
    bD=bD+(wat(i-1)*[res1(i);res2(i)]+(1-gamma(i))*C*B*0.0001);
    ind1(i)=norm(bR/(i-1));
    ind2(i)=norm(bD/(i-1));
end

k=1:len;
% figure(1);
% plot(k,pos(k),'g','LineWidth',2);hold off;
% 
% figure(2);
% plot(k,ang(k),'g','LineWidth',2);hold off;

figure(3);
semilogy(k,ind1(k),'g','LineWidth',2);hold on;

figure(4);
semilogy(k,ind2(k),'g','LineWidth',2);hold on;

%%检测阈值
max_len=max(len_vec);
thd1=zeros(1,max_len);
thd2=zeros(1,max_len);
kap1=1e-6;
rho1=1;
kap2=1e-5;
rho2=1;
for k=2:max_len
   thd1(k)=sqrt((1+rho1)*kap1*log((k-1))/(k-1))+1e-3;
   thd2(k)=sqrt((1+rho2)*kap2*log((k-1))/(k-1));
end

k=1:max_len;

figure(3);
semilogy(k,thd1(k),'--k','LineWidth',2);hold off;

figure(4);
semilogy(k,thd2(k),'--k','LineWidth',2);hold off;

saveas(fig3,'Thd_rr_ETCDW_wat01','pdf');
saveas(fig4,'Thd_rd_ETCDW_wat01','pdf');
saveas(fig1,'ETCDW01_nor4_x1','pdf');
saveas(fig2,'ETCDW01_nor4_x2','pdf');
saveas(fig5,'ETCDW01_nor4_ga','pdf');

sum=0;
for i=1:6
   trirate_vec(i)
   sum=sum+trirate_vec(i);
end
sum/6

RR=1;
NN=7000;
y=[pos,ang]';
J=0;
for i=1:NN
   J=J+(y(:,i))'*y(:,i);
end
J/NN