clear all
%% 数据加载
load('cov.mat');%载入A：A(581010x54):581010个数据
load('L_cov.mat');%载入L：A(1x581010):581010个结果
load('data/C_meth1_smote_sw2_800');%载入C_store
A=double(A);
L=double(L);
L(L==1)=-1;
L(L==2)=1;
A(435756:581009,:)=[];
L(435756:581009)=[];
%% 参数设置
agent_num=10;% agent个数
Maxgen=800;% 迭代次数

load('data/C_meth1_smote_sw2_800.mat');%载入C_store
C=C_store;
x_k_i_last=zeros(54,agent_num);% 第k次迭代，每个智能体x的值
q=zeros(54,agent_num);% Q阵
lamuda1=0.5*10^(-5);
lamuda2=0.5*10^(-5);
global v;% V阵
%% 数据预处理
%根据智能体个数裁剪数据，每个智能体十分之一的数据
for i=1:agent_num
    L_cut(i,:)=L((i-1)*floor(size(A,1)/agent_num)+1:i*floor(size(A,1)/agent_num));
    A_cut(:,:,i)=A((i-1)*floor(size(A,1)/agent_num)+1:i*floor(size(A,1)/agent_num),:); 
end
%% GUROBI设置
% x_k_sdp=sdpvar(123,1,'full');
% ops = sdpsettings('verbose',0,'solver','GUROBI');
%%
local_n=floor(size(A,1)/agent_num);

%% 算法主体
for k=1:Maxgen
    k   
    alpha=0.6;
    if rem(k,5)==0
        Ck=C{1,5};
    else
        Ck=C{1,rem(k,5)};
    end
    % 读取上次迭代的数据
    
    if k==1
        x_k_last=zeros(54,agent_num);
        xx_k_1=x_k_last;
        y_k_last=zeros(54,agent_num);
    elseif k==2
      x_k_last0=xx_k_1;
      x_k_last=x_k_store{k-1};
      y_k_last=y_k_store{k-1};
    else
      x_k_last0=x_k_store{k-2};
      x_k_last=x_k_store{k-1};
      y_k_last=y_k_store{k-1};
     end
        
     if k==1
        y=zeros(54,agent_num);
        for i=1:agent_num
            x_k_i_last=x_k_last(:,i);
            mid=L_cut(i,:)'.*A_cut(:,:,i); 
            gradient_k0(:,i)=zeros(54,1);
            gradient=-mid.*exp(mid*x_k_i_last)./(1+exp(mid*x_k_i_last)).^2;
            gradient=sum(gradient,1)/floor(size(A,1)/agent_num)+2*lamuda2*x_k_i_last';
            gradient_k(:,i)=gradient';
            y(:,i)=y_k_last(:,i)+gradient_k(:,i)+lamuda1*sign(x_k_i_last);
        end
        y_k_store{k}=y;
     else
             for i=1:agent_num 
               x_k_i_last0=x_k_last0(:,i);
               x_k_i_last=x_k_last(:,i);

               %%--求梯度
               mid=L_cut(i,:)'.*A_cut(:,:,i); 
                gradient0=-mid.*exp(mid*x_k_i_last0)./(1+exp(mid*x_k_i_last0)).^2;
                gradient0=sum(gradient0,1)/floor(size(A,1)/agent_num)+2*lamuda2*x_k_i_last0';
                gradient_k0(:,i)=gradient0';

                %-----求梯度--------
        %        mid=L_cut(i,:)'.*A_cut(:,:,i); 
                gradient=-mid.*exp(mid*x_k_i_last)./(1+exp(mid*x_k_i_last)).^2;
                gradient=sum(gradient,1)/floor(size(A,1)/agent_num)+2*lamuda2*x_k_i_last';
                gradient_k(:,i)=gradient';
                clear mid;

                %%----y 更新
                y(:,i)=y_k_last(:,i)+gradient_k(:,i)+lamuda1*sign(x_k_i_last)-gradient_k0(:,i)-lamuda1*sign(x_k_i_last0);
            end
            y_k_store{k}=y;
     end
%     y_sum=zeros(123,agent_num);
    xx=zeros(54,agent_num);
    for i=1:agent_num
        x_k_i_last=x_k_last(:,i);
        y_sum=zeros(54,1);
        for j=1:agent_num
            y_sum=y_sum+Ck(i,j)*y(:,j);
        end
        xx(:,i)=x_k_i_last-alpha*y_sum;
    end
    for i=1:agent_num
        %%----求合并信息
        x_consens=zeros(54,1);
        for j=1:agent_num
            x_consens=x_consens+Ck(i,j)*xx(:,j);
        end
        %---------求q----------
        x_k_i_new(:,i)=x_consens;%-alpha*gradient'-alpha*lamuda1*sign(x_k_i_last);%q(123x10)
    end
    gradient_sto{k}=gradient_k;
     x_k_store{k}=x_k_i_new;%统一更新旧的待优化参数值
     
end

