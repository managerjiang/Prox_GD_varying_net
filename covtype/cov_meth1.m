clear all
%% 数据加载
load('cov.mat');%载入A：A(581012x54):581012个数据
load('L_cov.mat');%载入L：A(1x581012):581012个结果
A=double(A);
L=double(L);
L(L==1)=-1;
L(L==2)=1;
A(435756:581009,:)=[];
L(435756:581009)=[];
%% 参数设置
agent_num=10;% agent个数
Maxgen=800;% 迭代次数
% C=doubly_stochastic(agent_num);% 生成随机的邻接矩阵（行和列和都是1）
% C_store=C
load('data/C_meth1_800.mat');%载入C_store
C=C_store
x_k_i_last=zeros(54,agent_num);% 第k次迭代，每个智能体x的值
q=zeros(54,agent_num);% Q阵
alpha=0.9;% 步长
lamuda1=5*10^(-5);
lamuda2=5*10^(-5);
global v;% V阵
%% 数据预处理
%根据智能体个数裁剪数据，每个智能体十分之一的数据
for i=1:agent_num
    L_cut(i,:)=L((i-1)*floor(size(A,1)/agent_num)+1:i*floor(size(A,1)/agent_num));
    A_cut(:,:,i)=A((i-1)*floor(size(A,1)/agent_num)+1:i*floor(size(A,1)/agent_num),:); 
end
%% GUROBI设置
x_k_sdp=sdpvar(54,1,'full');
ops = sdpsettings('verbose',0,'solver','GUROBI');
%% 算法主体
for k=1:Maxgen
    k   
    % 读取上次迭代的数据
    if k==1
      x_k_last=zeros(54,agent_num);
    else
      x_k_last=x_k_store{k-1};
    end
    % 开始循环算法
    for i=1:agent_num       
        x_k_i_last=x_k_last(:,i);
        %-----求梯度--------
        mid=L_cut(i,:)'.*A_cut(:,:,i); 
        gradient=-mid.*exp(mid*x_k_i_last)./(1+exp(mid*x_k_i_last)).^2;
        gradient=sum(gradient,1)/floor(size(A,1)/agent_num)+2*lamuda2*x_k_i_last';
        gradient_k(:,i)=gradient';
        clear mid;
        %---------求q----------
        q(:,i)=x_k_i_last-alpha*gradient';%q(54x10)
    end
    gradient_sto{k}=gradient_k;
     %所有个体异步更新q值，再更新V和x
    for i=1:agent_num
         %---------求v-----------
          mid=0;%v(54x1)
        for j=1:agent_num
           mid=mid+C(i,j)*q(:,j);
        end
        v=mid;
        clear mid;
        %----------方法三：yalmip+gurobi求解--------------
     f=lamuda1*norm(x_k_sdp,1)+norm(x_k_sdp-v,2)^2/(2*alpha);
     reuslt = optimize([],f,ops);
     x_k_i_new(:,i)=value(x_k_sdp);
    end
     x_k_store{k}=x_k_i_new;%统一更新旧的待优化参数值
     
end

