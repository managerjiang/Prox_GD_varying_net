function [A1,L1]=combine_data(A,L)
%输入：A(32561x123)/L(1x32561)
%输出：混合好的样本
num_posi=sum(L==1,2);
num_nega=sum(L==-1,2);
A_posi=A;
A_posi(L==-1,:)=[];
A_nega=A;
A_nega(L==1,:)=[];
A1=zeros(2*num_posi,size(A,2));
L1=zeros(1,2*num_posi);
choose=randperm(2*num_posi,num_posi);
L1(choose)=1;
L1(L1==0)=-1;
A1(choose,:)=A_posi;
A1(L1==-1,:)=A_nega(randperm(num_nega,num_posi),:);
end