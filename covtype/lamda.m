function y=lamda(C,k)
%C�Ǹ�����5��ͨ�����˽ṹ
%k�ǵ�ǰ��������
%y�ǽ��
C1=C{1};
C2=C{2};
C3=C{3};
C4=C{4};
C5=C{5};
num=fix(k/5);
yu=rem(k,5);
y=(C1*C2*C3*C4*C5)^num;
if yu==0
    y=y;
elseif yu==1
    y=y*C1;
elseif yu==2
    y=y*C1*C2;
elseif yu==3
    y=y*C1*C2*C3;
elseif yu==4
    y=y*C1*C2*C3*C4;
end
end