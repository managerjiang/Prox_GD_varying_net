clear all
%% ���ݼ���
load('cov.mat');%����A��A(581010x54):581010������
load('L_cov.mat');%����L��A(1x581010):581010�����
load('data/C_meth1_smote_sw2_800');%����C_store
A=double(A);
L=double(L);
L(L==1)=-1;
L(L==2)=1;
A(435756:581009,:)=[];
L(435756:581009)=[];
%% ��������
agent_num=10;% agent����
Maxgen=800;% ��������
C=C_store%Ĭ�϶Ա��㷨�뷽��һ���ڽӾ�����һ����
% C{1}=doubly_stochastic(10);% ����5��������ڽӾ����к��кͶ���1��
% C{2}=doubly_stochastic(10);
% C{3}=doubly_stochastic(10);
% C{4}=doubly_stochastic(10);
% C{5}=doubly_stochastic(10);
% C_store=C;
tau_i=5;% ����
lamuda1=0.5*10^(-5);
lamuda2=0.5*10^(-5);
global v;% V��
%% ����Ԥ����
%��������������ü����ݣ�ÿ��������ʮ��֮һ������
for i=1:agent_num
    L_cut(i,:)=L((i-1)*floor(size(A,1)/agent_num)+1:i*floor(size(A,1)/agent_num));
    A_cut(:,:,i)=A((i-1)*floor(size(A,1)/agent_num)+1:i*floor(size(A,1)/agent_num),:); 
end
%% GUROBI����
% x_i_sdp=sdpvar(54,1,'full');
% ops = sdpsettings('verbose',0,'solver','GUROBI');
%% ��ʼ��
x=zeros(54,agent_num);% x��
for i=1:agent_num
 %-----���ݶ�--------
 mid=L_cut(i,:)'.*A_cut(:,:,i);
 gradient_i_new=-mid.*exp(mid*x(:,i))./(1+exp(mid*x(:,i))).^2;
 gradient_i_new=sum(gradient_i_new,1)/floor(size(A,1)/agent_num)+2*lamuda2*x(:,i)';
 clear mid;
 gradient(:,i)=gradient_i_new';
end
y=gradient;% y��
pi=agent_num*y-gradient;% pai��
z=zeros(54,agent_num);% z��

%% �㷨����
for k=1:Maxgen
    k   
    % ��ʼѭ���㷨
    for i=1:agent_num  
         %----------��������yalmip+gurobi���--------------
    % f=lamuda1*norm(x_i_sdp,1)+pi(:,i)'*(x_i_sdp-x(:,i));
     mid=L_cut(i,:)'.*A_cut(:,:,i);
     mid=-1./(1+exp(mid*x(:,i)));
       f=@(x_i_sdp)lamuda1*norm(x_i_sdp,1)+pi(:,i)'*(x_i_sdp-x(:,i))+sum(mid,1)/floor(size(A,1)/agent_num)+gradient(:,i)'*(x_i_sdp-x(:,i))+tau_i/2*norm(x_i_sdp-x(:,i),2)^2;
     clear mid;   
%      optimize([],f,ops);
     [xi,fval]=fminunc(f,x(:,i));
     xx(:,i)=xi;
        z(:,i)=x(:,i)+1/k^0.6*(xx(:,i)-x(:,i));
    end
     
     C_k=lamda(C,k);%lameda(10x10)
    for i=1:agent_num     
        x_i_new=0;
        y_i_new=0;
        for j=1:agent_num
           x_i_new=x_i_new+C_k(i,j)*z(:,j);
           y_i_new=y_i_new+C_k(i,j)*y(:,j);
        end
        x(:,i)=x_i_new;
        %-----���ݶ�--------
        mid=L_cut(i,:)'.*A_cut(:,:,i); 
        gradient_i_new=-mid.*exp(mid*x_i_new)./(1+exp(mid*x_i_new)).^2;
        gradient_i_new=sum(gradient_i_new,1)/floor(size(A,1)/agent_num)+2*lamuda2*x_i_new';
        clear mid;   
        y_i_new=y_i_new+gradient_i_new'-gradient(:,i);
        y(:,i)=y_i_new;
        pi(:,i)=agent_num*y_i_new-gradient_i_new';
        gradient(:,i)=gradient_i_new';
    end
    x_k_store{k}=x;
    gradient_sto{k}=gradient;
end   


