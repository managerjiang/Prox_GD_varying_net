%未知量a(1,2),a(2,3),a(3,4),a(4,5),a(5,6),a(6,7),a(6,8),a(7,8)
%n矩阵维度
function a=matrix_general(n)
    for num=1:100
        %n=10;
        a=rand(n,n)*(2/n);%随机生成的0-1的数据
        %对角为0
        for i=1:n
            a(i,i)=0;
        end
        a;
        %手动选择一个范围，使得该范围的数据设置为零，达到随机赋值零的目的
        %保证节点联通的效果，用眼看吧
        %调整稀疏度x(0-10)
        x=6;
        for i=1:12
            for j=i:n
                if ((a(i,j)>0)&&(a(i,j)<x/(5*n)))
                    a(i,j)=0;
                end
            end
        end
        %对称化
        for i=1:n
            for j=i:n
                a(j,i)=a(i,j);
            end
        end

%         %计算a(1,2),a(2,3),a(3,4),a(4,5),a(5,6),
        for i=1:n-3
            a(i,i+1)=1;
            for j=1:n
                if j ~= i+1
                    a(i,i+1)=a(i,i+1)-a(i,j);
                    a(i+1,i)=a(i,i+1);
                end
            end
        end
     %   需要计算的未知量还剩3个
%         a(n-2,n-1);a(n-2,n);a(n-1,n)
%         a(n-2,n-1)+a(n-2,n)=suman_2;
%         a(n-2,n-1)+a(n-1,n)=suman_1;
%         a(n-2,n)+a(n-1,n)=suman_0;
        suman_2=1;
        suman_1=1;
        suman_0=1;
        for i=1:n-3
            suman_2=suman_2-a(n-2,i);
            suman_1=suman_1-a(n-1,i);
            suman_0=suman_0-a(n,i);
        end
        a(n-2,n-1)=(suman_2+suman_1-suman_0)/2;
        a(n-2,n)=(suman_2-suman_1+suman_0)/2;
        a(n-1,n)=(-suman_2+suman_1+suman_0)/2;

        a(n,n-2)=a(n-2,n);
        a(n-1,n-2)=a(n-2,n-1);
        a(n,n-1)=a(n-1,n);

        suma=zeros(n,1);
        for i=1:n
            for j=1:n
              suma(i)=suma(i)+a(i,j);
            end 
        end

        mina=min(a(:));
        if (mina>=0)
            a
            suma
            break;
        end
    end  
end