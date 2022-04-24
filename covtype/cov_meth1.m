clear all
%% ���ݼ���
load('cov.mat');%����A��A(581012x54):581012������
load('L_cov.mat');%����L��A(1x581012):581012�����
A=double(A);
L=double(L);
L(L==1)=-1;
L(L==2)=1;
A(435756:581009,:)=[];
L(435756:581009)=[];
%% ��������
agent_num=10;% agent����
Maxgen=800;% ��������
% C=doubly_stochastic(agent_num);% ����������ڽӾ����к��кͶ���1��
% C_store=C
load('data/C_meth1_800.mat');%����C_store
C=C_store
x_k_i_last=zeros(54,agent_num);% ��k�ε�����ÿ��������x��ֵ
q=zeros(54,agent_num);% Q��
alpha=0.9;% ����
lamuda1=5*10^(-5);
lamuda2=5*10^(-5);
global v;% V��
%% ����Ԥ����
%��������������ü����ݣ�ÿ��������ʮ��֮һ������
for i=1:agent_num
    L_cut(i,:)=L((i-1)*floor(size(A,1)/agent_num)+1:i*floor(size(A,1)/agent_num));
    A_cut(:,:,i)=A((i-1)*floor(size(A,1)/agent_num)+1:i*floor(size(A,1)/agent_num),:); 
end
%% GUROBI����
x_k_sdp=sdpvar(54,1,'full');
ops = sdpsettings('verbose',0,'solver','GUROBI');
%% �㷨����
for k=1:Maxgen
    k   
    % ��ȡ�ϴε���������
    if k==1
      x_k_last=zeros(54,agent_num);
    else
      x_k_last=x_k_store{k-1};
    end
    % ��ʼѭ���㷨
    for i=1:agent_num       
        x_k_i_last=x_k_last(:,i);
        %-----���ݶ�--------
        mid=L_cut(i,:)'.*A_cut(:,:,i); 
        gradient=-mid.*exp(mid*x_k_i_last)./(1+exp(mid*x_k_i_last)).^2;
        gradient=sum(gradient,1)/floor(size(A,1)/agent_num)+2*lamuda2*x_k_i_last';
        gradient_k(:,i)=gradient';
        clear mid;
        %---------��q----------
        q(:,i)=x_k_i_last-alpha*gradient';%q(54x10)
    end
    gradient_sto{k}=gradient_k;
     %���и����첽����qֵ���ٸ���V��x
    for i=1:agent_num
         %---------��v-----------
          mid=0;%v(54x1)
        for j=1:agent_num
           mid=mid+C(i,j)*q(:,j);
        end
        v=mid;
        clear mid;
        %----------��������yalmip+gurobi���--------------
     f=lamuda1*norm(x_k_sdp,1)+norm(x_k_sdp-v,2)^2/(2*alpha);
     reuslt = optimize([],f,ops);
     x_k_i_new(:,i)=value(x_k_sdp);
    end
     x_k_store{k}=x_k_i_new;%ͳһ���¾ɵĴ��Ż�����ֵ
     
end

