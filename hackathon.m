clc
clear all
close all
vpa(10);
ylen=30;
xlen=29;
maxiter=1000000;
error=0.00000001;

subslen=4;
airlen=26;
groundlen=7;
gaplen=4;
linelen=7;

h1=0.5;
h3=0.5;
h2=1;
h4=1;
eps0=4.5;
e0=8.8543*(10^(-12));
a0 = (-2) * ( ( 1/ (h1*h3) ) + ( 1 / (h2 * h4) ) );
a1 = ( 2 / ( h1 * (h1 + h3 )));
a2 = ( 2 / ( h2 * (h2 + h4 )));
a3 = ( 2 / ( h3 * (h1 + h3 )));
a4 = ( 2 / ( h4 * (h2 + h4 )));
a2new=(a2+a4)*a2/(4.5*a4+a2);
a4new=eps0*(a2+a4)*a2/(4.5*a4+a2);
% ai constants
col=xlen/h1 +1;
row=ylen/h2 +1;
%grid dots per y and x axis
A=zeros(row,col);
A(row-4,23:37)=1;
for q=1:maxiter
    A3=A;
for i=2:(row-1)
    for j=2:(col-1)
        if i==row-subslen/h2 
            
            if ~any(ismember([[(2+groundlen/h1):((groundlen+gaplen)/h1)],[(((groundlen+gaplen+linelen)/h1)+2):(col-groundlen/h1 -1)]],j ))
                continue;
        
            else
                A(i,j)=(A(i,j+1)*a1+A(i-1,j)*a2new+A(i,j-1)*a3+A(i+1,j)*a4new)/(-1*a0);

            end
        
        else
        A(i,j)=(A(i,j+1)*a1+A(i-1,j)*a2+A(i,j-1)*a3+A(i+1,j)*a4)/(-1*a0);
        end
    end
end
 if max(max(abs((A-A3 )./A)))<error
     fprintf("\n approached the acceptable solution for the Substrate Case %d th iteration ",q);
     disp("HALT! The potential distribition of Substrate case is ready");
     break;
 end
end

B=zeros(row-1,col-1);
for i=1:(row-1)
    for j=1:(col-1)
       u0=A(i+1,j);
       u1=A(i+1,j+1);
       u2=A(i,j+1);
       u3=A(i,j);
       B(i,j)=(h2/(3*h1))*( (u0-u1)^2 +(u3-u2)^2 +(u0-u1)*(u3-u2))+(h1/(3*h2))*((u0-u3)^2+(u1-u2)^2 +(u0-u3)*(u1-u2));
    end
end
sum=0;
for i=1:row-1
    for j=1:col-1
        if i<row-subslen/h2
         sum=sum+e0*B(i,j);   
        else
         sum=sum+eps0*e0*B(i,j);   
        end
    end
end
A_0=zeros(row,col);
A_0(row-subslen/h2,(1+(groundlen+gaplen)/h1):(1+(groundlen+gaplen+linelen)/h1))=1;

for q=1:maxiter
    A_2=A_0;
for i=2:row-1
    for j=2:col-1
         if i==row-subslen/h2 && ~any(ismember([[16:22],[38:44]],j))
                continue
         else
        A_0(i,j)=(A_0(i,j+1)*a1+A_0(i-1,j)*a2+A_0(i,j-1)*a3+A_0(i+1,j)*a4)/(-1*a0);
         end
    end
end

 if max(max(abs((A_0-A_2)./A_0)))<error
     fprintf("\n approached the acceptable solution for the no Substrate Case %d th iteration ",q);
     disp("HALT! The potential distribition of no Substrate case is ready");
     break;
 end
end

%applying finite differences method to find the value for each mesh of the
%field
B_0=zeros(row-1,col-1);
for i=1:row-1
    for j=1:col-1
       u0=A_0(i+1,j);
       u1=A_0(i+1,j+1);
       u2=A_0(i,j+1);
       u3=A_0(i,j);
       B_0(i,j)=(h2/(3*h1))*( (u0-u1)^2 +(u3-u2)^2 +(u0-u1)*(u3-u2))+(h1/(3*h2))*((u0-u3)^2+(u1-u2)^2 +(u0-u3)*(u1-u2));
        
       %computing the integral for each grid area
    end
end
sum_0=0;

for i=1:row-1
    for j=1:col-1
         sum_0=sum_0+e0*B_0(i,j); 
         %finding the capacitance by doing the sigma function for
         %capacitance
    end
end

vp_a=299792458;
z0=1/(vp_a*sqrt(sum*sum_0));%Characteristic line impedance 
vp=vp_a*sqrt(sum_0/sum);%phase velocity for given substrate problem  
e_f=sum/sum_0;%Effective relative permittivity of the problem

E0_x=zeros(row,col);
E0_y=zeros(row,col);

for i=2:row-1
    for j=2:col-1
        E0_x(i,j)=A(i,j-1)-A(i,j+1);
        E0_y(i,j)=(A(i+1,j)-A(i-1,j))/2;
    end
end
E0_x_0=zeros(row,col);
E0_y_0=zeros(row,col);

for i=2:row-1
    for j=2:col-1
        E0_x_0(i,j)=A_0(i,j-1)-A_0(i,j+1);
        E0_y_0(i,j)=(A_0(i+1,j)-A_0(i-1,j))/2;
        %calculation electrical field of no substance CPWG
    end
end
%ploting Potential Dist. and Vector form of electrical field for CPWG with
%substance and no substance

subplot(2,2,1)
heatmap(A,'Colormap',   jet,'Title','Color map of Potential distribution,CPWG with Substrate');

subplot(2,2,2)
heatmap(A_0,'Colormap',   jet,'Title','Color map of Potential distribution,CPWG with no Substrate');
subplot(2,2,3)
title('Electrical field vector form ,CPWG with Substrate')
hold on;
quiver(E0_x(end:-1:1,:),E0_y(end:-1:1,:));hold off;
subplot(2,2,4)
title('Electrical field vector form ,CPWG with no Substrate');hold on;
quiver(E0_x_0(end:-1:1,:),E0_y_0(end:-1:1,:));hold off;
