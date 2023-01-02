clc
clear all
vpa(10);
A=zeros(31,59);
A(27,23:37)=1;
a0=-10
a1=4
a2=1
a3=4
a4=1
a2new=(a2+a4)*a2/(4.5*a4+a2);
a4new=4.5*(a2+a4)*a2/(4.5*a4+a2);
for q=1:10000
    A2=A;
for i=2:30
    for j=2:58
        if i==27 
            
            if ~ismember([[16:22],[38:44]],j)
                continue
        
            else
                A(i,j)=(A(i,j+1)*a1+A(i-1,j)*a2new+A(i,j-1)*a3+A(i+1,j)*a4new)/(-1*a0);

            end
         end
        A(i,j)=(A(i,j+1)*a1+A(i-1,j)*a2+A(i,j-1)*a3+A(i+1,j)*a4)/(-1*a0);
    end
end
%max((A-A2)./A));
 if max(max(abs((A-A2)./A)))<00000.1
     disp("HALT");
     break;
 end
end
% [X, Y] = meshgrid(linspace(0,30,31), linspace(0,29,59));
% Z = griddata(,A(:,2),X,Y);
% Z(isnan(Z)) = 0;
% imshow(Z)
for l=1:1
    col=31-l;
B=zeros(col,59-l);
for i=1:col
    for j=1:(59-l)
       u0=A(i+1,j);
       u1=A(i+1,j+1);
       u2=A(i,j+1);
       u3=A(i,j);
       B(i,j)=2/3*( (u0-u1)^2 +(u3-u2)^2 +(u0-u1)*(u3-u2))+1/6*((u0-u3)^2+(u1-u2)^2 +(u0-u3)*(u1-u2));
    end
end
%A=B;
end
sum=0;
e0=8.8543*(10^(-12));
for i=1:30
    for j=1:58
        if(i<27)
         sum=sum+e0*B(i,j);   
        else
         sum=sum+4.5*e0*B(i,j);   
        end
    end
end
A_0=zeros(31,59);
A_0(27,23:37)=1;

for q=1:10000
    A_2=A_0;
for i=2:30
    for j=2:58
         if i==27 && ~any(ismember([[16:22],[38:44]],j))
                continue
         else
        A_0(i,j)=(A_0(i,j+1)*a1+A_0(i-1,j)*a2+A_0(i,j-1)*a3+A_0(i+1,j)*a4)/(-1*a0);
         end
    end
end
%max((A-A2)./A));
 if max(max(abs((A_0-A_2)./A_0)))<00000.1
     disp("HALT");
     break;
 end
end
% [X, Y] = meshgrid(linspace(0,30,31), linspace(0,29,59));
% Z = griddata(,A(:,2),X,Y);
% Z(isnan(Z)) = 0;
% imshow(Z)
B_0=zeros(30,58);
for i=1:30
    for j=1:58
       u0=A_0(i+1,j);
       u1=A_0(i+1,j+1);
       u2=A_0(i,j+1);
       u3=A_0(i,j);
       B_0(i,j)=2/3*( (u0-u1)^2 +(u3-u2)^2 +(u0-u1)*(u3-u2))+1/6*((u0-u3)^2+(u1-u2)^2 +(u0-u3)*(u1-u2));
    end
end
sum_0=0;
e0=8.8543*(10^(-12));
for i=1:30
    for j=1:58
         sum_0=sum_0+e0*B_0(i,j);   
    end
end
vp_a=299792458;
z0=1/(vp_a*sqrt(sum*sum_0));
vp=vp_a*sqrt(sum_0/sum)
e_f=sum/sum_0;
%A
E0_x=zeros(31,59);
E0_y=zeros(31,59);
for i=2:30
    for j=2:58
        E0_x(i,j)=A(i,j-1)-A(i,j+1);
        E0_y(i,j)=(A(i+1,j)-A(i-1,j))/2;
    end
end
E0_x_0=zeros(31,59);
E0_y_0=zeros(31,59);

for i=2:30
    for j=2:58
        E0_x_0(i,j)=A_0(i,j-1)-A_0(i,j+1);
        E0_y_0(i,j)=(A_0(i+1,j)-A_0(i-1,j))/2;
    end
end
test=E0_y-E0_y_0;
find(test~= 0)
subplot(2,2,1)

heatmap(A,'Colormap',   jet,'Title','CPWG with Substrate');
subplot(2,2,2)
timeout = 10
heatmap(A_0,'Colormap',   jet,'Title','CPWG with No Substrate');
subplot(2,2,3)
title('CPWG with Substrate');hold on;
quiver(E0_x(end:-1:1,:),E0_y(end:-1:1,:));hold off;
subplot(2,2,4)
title('CPWG with No Substrate');hold on;
quiver(E0_x_0(end:-1:1,:),E0_y_0(end:-1:1,:));hold off;
