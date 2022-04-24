close all;
%% meth1
% ==============需要读取.mat文件================
  load('data/G_meth1_800.mat');%加载保存的迭代梯度信息x_k_store{}
  load('data/X_meth1_800.mat');%加载保存的迭代解信息gradient()
  load('data/C_meth1_800.mat');%加载保存的邻接矩阵信息C_store
  lamuda1=5*10^(-4);
  lamuda2=5*10^(-4);
  load('cov.mat');
  load('L_cov.mat');
  A(435756:581009,:)=[];
L(435756:581009)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
  C=C_store;
% 参数设置
Maxgen=size(x_k_store,2);%迭代次数
agent_num=size(C);%智能体个数

% 训练集
for k=1:Maxgen
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  % 目标函数
  fi=sum(1./(1+exp(L'.*A*x_k_1)),1);
  obj_m1(k)=fi/size(A,1)+lamuda1*norm(x_k_1,1)+lamuda2*norm(x_k_1,2)^2; 
  %  XLX
  XLX=0;
  for i=1:agent_num
    temp_X=zeros(size(A,2),1);
    for j=1:agent_num
        temp_X=temp_X+C(i,j)*(x_k(:,i)-x_k(:,j));
    end
    XLX=XLX+x_k(:,i)'*temp_X;
  end
  w_m1(k)=XLX;
  % 梯度值
   %取所有智能体中第一个
   g_k=gradient_sto{k}+lamuda1*sign(x_k_1);
   zz_m1(k)=norm(g_k(:,1));
end
% 测试集
  load('cov.mat');
  load('L_cov.mat');
  A(1:435755,:)=[];
L(1:435755)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
for k=1:Maxgen
     % 测试集正确率
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  result=A*x_k_1;
  result(result>=0)=1;
  result(result<0)=-1;
  yy_m1(k)=sum((result==L'))/size(L,2);
  
end

%% meth2
% ==============需要读取.mat文件================
  load('data/G_meth2_800.mat');%加载保存的迭代梯度信息x_k_store{}
  load('data/X_meth2_800.mat');%加载保存的迭代解信息gradient()
  load('data/C_meth2_800.mat');%加载保存的邻接矩阵信息C_store
  lamuda1=5*10^(-4);
  lamuda2=5*10^(-4);
  load('cov.mat');
  load('L_cov.mat');
  A(435756:581009,:)=[];
L(435756:581009)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
  C=C_store;
% 参数设置
Maxgen=size(x_k_store,2);%迭代次数
agent_num=size(C);%智能体个数

% 训练集
for k=1:Maxgen
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  % 目标函数
  fi=sum(1./(1+exp(L'.*A*x_k_1)),1);
  obj_m2(k)=fi/size(A,1)+lamuda1*norm(x_k_1,1)+lamuda2*norm(x_k_1,2)^2; 
  %  XLX
  XLX=0;
  for i=1:agent_num
    temp_X=zeros(size(A,2),1);
    for j=1:agent_num
        temp_X=temp_X+C(i,j)*(x_k(:,i)-x_k(:,j));
    end
    XLX=XLX+x_k(:,i)'*temp_X;
  end
  w_m2(k)=XLX;
  % 梯度值
   %取所有智能体中第一个
   g_k=gradient_sto{k}+lamuda1*sign(x_k_1);
   zz_m2(k)=norm(g_k(:,1));
end
% 测试集
  load('cov.mat');
  load('L_cov.mat');
  A(1:435755,:)=[];
L(1:435755)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
for k=1:Maxgen
     % 测试集正确率
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  result=A*x_k_1;
  result(result>=0)=1;
  result(result<0)=-1;

  yy_m2(k)=sum((result==L'))/size(L,2);
  
end

%% meth1_sw
% ==============需要读取.mat文件================
  load('data/G_meth1_sw_800.mat');%加载保存的迭代梯度信息x_k_store{}
  load('data/X_meth1_sw_800.mat');%加载保存的迭代解信息gradient()
  load('data/C_meth1_sw_800.mat');%加载保存的邻接矩阵信息C_store
  lamuda1=5*10^(-4);
  lamuda2=5*10^(-4);
  load('cov.mat');
  load('L_cov.mat');
  A(435756:581009,:)=[];
L(435756:581009)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
  C=C_store{1};
% 参数设置
Maxgen=size(x_k_store,2);%迭代次数
agent_num=size(C);%智能体个数

% 训练集
for k=1:Maxgen
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  % 目标函数
  fi=sum(1./(1+exp(L'.*A*x_k_1)),1);
  obj_m1s(k)=fi/size(A,1)+lamuda1*norm(x_k_1,1)+lamuda2*norm(x_k_1,2)^2; 
  %  XLX
  XLX=0;
  for i=1:agent_num
    temp_X=zeros(size(A,2),1);
    for j=1:agent_num
        temp_X=temp_X+C(i,j)*(x_k(:,i)-x_k(:,j));
    end
    XLX=XLX+x_k(:,i)'*temp_X;
  end
  w_m1s(k)=XLX;
  % 梯度值
   %取所有智能体中第一个
   g_k=gradient_sto{k}+lamuda1*sign(x_k_1);
   zz_m1s(k)=norm(g_k(:,1));
end
% 测试集
  load('cov.mat');
  load('L_cov.mat');
  A(1:435755,:)=[];
L(1:435755)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
for k=1:Maxgen
     % 测试集正确率
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  result=A*x_k_1;
  result(result>=0)=1;
  result(result<0)=-1;

  yy_m1s(k)=sum((result==L'))/size(L,2);
end

%% meth3
% ==============需要读取.mat文件================
  load('data/G_meth3_sw_800.mat');%加载保存的迭代梯度信息x_k_store{}
  load('data/X_meth3_sw_800.mat');%加载保存的迭代解信息gradient()
%   load('data/C_meth3_sw_800.mat');%加载保存的邻接矩阵信息C_store
  lamuda1=5*10^(-4);
  lamuda2=5*10^(-4);
  load('cov.mat');
  load('L_cov.mat');
  A(435756:581009,:)=[];
L(435756:581009)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
%   C=C_store{1};
% 参数设置
Maxgen=size(x_k_store,2);%迭代次数
agent_num=size(C);%智能体个数

% 训练集
for k=1:Maxgen
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  % 目标函数
  fi=sum(1./(1+exp(L'.*A*x_k_1)),1);
  obj_m3(k)=fi/size(A,1)+lamuda1*norm(x_k_1,1)+lamuda2*norm(x_k_1,2)^2; 
  %  XLX
  XLX=0;
  for i=1:agent_num
    temp_X=zeros(size(A,2),1);
    for j=1:agent_num
        temp_X=temp_X+C(i,j)*(x_k(:,i)-x_k(:,j));
    end
    XLX=XLX+x_k(:,i)'*temp_X;
  end
  w_m3(k)=XLX;
  % 梯度值
   %取所有智能体中第一个
   g_k=gradient_sto{k}+lamuda1*sign(x_k_1);
   zz_m3(k)=norm(g_k(:,1));
end
% 测试集
  load('cov.mat');
  load('L_cov.mat');
  A(1:435755,:)=[];
L(1:435755)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
for k=1:Maxgen
     % 测试集正确率
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  result=A*x_k_1;
  result(result>=0)=1;
  result(result<0)=-1;
  yy_m3(k)=sum((result==L'))/size(L,2);
 
end

%% meth4
% ==============需要读取.mat文件================
  load('data/G_meth4_sw_800.mat');%加载保存的迭代梯度信息x_k_store{}
  load('data/X_meth4_sw_800.mat');%加载保存的迭代解信息gradient()
%   load('data/C_meth3_sw_800.mat');%加载保存的邻接矩阵信息C_store
  lamuda1=5*10^(-4);
  lamuda2=5*10^(-4);
  load('cov.mat');
  load('L_cov.mat');
  A(435756:581009,:)=[];
L(435756:581009)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
%   C=C_store{1};
% 参数设置
Maxgen=size(x_k_store,2);%迭代次数
agent_num=size(C);%智能体个数

% 训练集
for k=1:Maxgen
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  % 目标函数
  fi=sum(1./(1+exp(L'.*A*x_k_1)),1);
  obj_m4(k)=fi/size(A,1)+lamuda1*norm(x_k_1,1)+lamuda2*norm(x_k_1,2)^2; 
  %  XLX
  XLX=0;
  for i=1:agent_num
    temp_X=zeros(size(A,2),1);
    for j=1:agent_num
        temp_X=temp_X+C(i,j)*(x_k(:,i)-x_k(:,j));
    end
    XLX=XLX+x_k(:,i)'*temp_X;
  end
  w_m4(k)=XLX;
  % 梯度值
   %取所有智能体中第一个
   g_k=gradient_sto{k}+lamuda1*sign(x_k_1);
   zz_m4(k)=norm(g_k(:,1));
end
% 测试集
  load('cov.mat');
  load('L_cov.mat');
  A(1:435755,:)=[];
L(1:435755)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
for k=1:Maxgen
     % 测试集正确率
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  result=A*x_k_1;
  result(result>=0)=1;
  result(result<0)=-1;
  yy_m4(k)=sum((result==L'))/size(L,2);
 
end

%% meth5
% ==============需要读取.mat文件================
  load('data/G_meth5_sw_800.mat');%加载保存的迭代梯度信息x_k_store{}
  load('data/X_meth5_sw_800.mat');%加载保存的迭代解信息gradient()
%   load('data/C_meth3_sw_800.mat');%加载保存的邻接矩阵信息C_store
  lamuda1=5*10^(-4);
  lamuda2=5*10^(-4);
  load('cov.mat');
  load('L_cov.mat');
  A(435756:581009,:)=[];
L(435756:581009)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
%   C=C_store{1};
% 参数设置
Maxgen=size(x_k_store,2);%迭代次数
agent_num=size(C);%智能体个数

% 训练集
for k=1:Maxgen
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  % 目标函数
  fi=sum(1./(1+exp(L'.*A*x_k_1)),1);
  obj_m5(k)=fi/size(A,1)+lamuda1*norm(x_k_1,1)+lamuda2*norm(x_k_1,2)^2; 
  %  XLX
  XLX=0;
  for i=1:agent_num
    temp_X=zeros(size(A,2),1);
    for j=1:agent_num
        temp_X=temp_X+C(i,j)*(x_k(:,i)-x_k(:,j));
    end
    XLX=XLX+x_k(:,i)'*temp_X;
  end
  w_m5(k)=XLX;
  % 梯度值
   %取所有智能体中第一个
   g_k=gradient_sto{k}+lamuda1*sign(x_k_1);
   zz_m5(k)=norm(g_k(:,1));
end
% 测试集
  load('cov.mat');
  load('L_cov.mat');
  A(1:435755,:)=[];
L(1:435755)=[];
  L=double(L);
  L(L==1)=-1;
  L(L==2)=1;
for k=1:Maxgen
     % 测试集正确率
  x_k=x_k_store{k};
  x_k_1=x_k(:,1);%取所有智能体中第一个
  result=A*x_k_1;
  result(result>=0)=1;
  result(result<0)=-1;
  yy_m5(k)=sum((result==L'))/size(L,2);
 
end
figure(1);
plot(1:Maxgen,w_m1,'r-.','linewidth',1),hold on;
plot(1:Maxgen,w_m2,'b--','linewidth',1),hold on;
plot(1:Maxgen,w_m1s,'r','linewidth',1),hold on;
plot(1:Maxgen,w_m3,'k:','linewidth',1.5);
legend('meth1','meth2');
ylabel('$$\|D(\bar{x})\|$$','Interpreter','latex')
xlabel('iterations');
%% 绘图  
figure(3);
% subplot(2,2,1)
% plot(1:Maxgen,w_m1,'r-.','linewidth',1),hold on;
% plot(1:Maxgen,w_m2,'b--','linewidth',1),hold on;
plot(1:Maxgen,w_m1s,'r','linewidth',1),hold on;
plot(1:Maxgen,w_m3,'k:','linewidth',1.5), hold on;
plot(1:Maxgen,w_m4,'b-.','linewidth',1), hold on;
plot(1:Maxgen,w_m5,'g--','linewidth',1), hold on;
legend('meth1-s','meth2-s','meth3-s','meth4-s');
ylabel('$$\|D(\bar{x})\|$$','Interpreter','latex')
xlabel('iterations')
figure(2);
% subplot(2,2,1)
% plot(1:Maxgen,w_m1,'r-.','linewidth',1),hold on;
% plot(1:Maxgen,w_m2,'b--','linewidth',1),hold on;
% plot(1:Maxgen,w_m1s,'r','linewidth',1),hold on;
% plot(1:Maxgen,w_m3,'k:','linewidth',1.5);
% legend('meth1','meth2','meth1-s','meth3-s');
% ylabel('$$\|D(\bar{x})\|$$','Interpreter','latex')
% xlabel('iterations')
% title('(a)')

subplot(1,2,1)
% plot(1:Maxgen,zz_m1,'r-.','linewidth',1),hold on;
% plot(1:Maxgen,zz_m2,'b--','linewidth',1),hold on;
plot(1:Maxgen,zz_m1s,'r','linewidth',1),hold on;
plot(1:Maxgen,zz_m3,'k:','linewidth',1), hold on;
plot(1:Maxgen,zz_m4,'b-.','linewidth',1),hold on;
plot(1:Maxgen,zz_m5,'g--','linewidth',1),hold on;
%ylim([0,0.22])
legend('meth1-s','meth2-s','meth3-s','meth4-s');
ylabel('$$\|\nabla f(x)\|$$','Interpreter','latex')
xlabel('iterations')


subplot(1,2,2)
% plot(1:Maxgen,obj_m1,'r-.','linewidth',1),hold on;
% plot(1:Maxgen,obj_m2,'b--','linewidth',1),hold on;
plot(1:Maxgen,obj_m1s,'r','linewidth',1),hold on;
plot(1:Maxgen,obj_m3,'k:','linewidth',1), hold on;
plot(1:Maxgen,obj_m4,'b-.','linewidth',1), hold on;
plot(1:Maxgen,obj_m5,'g--','linewidth',1), hold on;
legend('meth1-s','meth2-s','meth3-s','meth4-s');
ylabel('$$f(x)$$','Interpreter','latex')
xlabel('iterations')
axes('Position',[0.67,0.34,0.22,0.18]);   % set the small figure     

plot(300:10:650,obj_m1s(300:10:650),'r','linewidth',1), hold on;      
plot(300:10:650,obj_m4(300:10:650),'b-.','linewidth',1), hold on;% plot the local small figure                                                                                                               
xlim([300,650]);                         % set the axes range 
% subplot(2,2,4)
% plot(1:Maxgen,yy_m1,'r-.','linewidth',1),hold on;
% plot(1:Maxgen,yy_m2,'b--','linewidth',1),hold on;
% plot(1:Maxgen,yy_m1s,'r','linewidth',1),hold on;
% plot(1:Maxgen,yy_m3,'k:','linewidth',1.5);
% legend('meth1','meth2','meth1-s','meth3-s');
% ylabel('Test accuracy')
% xlabel('iterations')
% title('(d)')
