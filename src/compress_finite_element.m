function [x1,y1,ss,Range,time]=compress_finite_element(x,y,num,E,Po,density,percentage,Maxit,Tol,Thick,force)
%描述：基于有限元原理对轨迹数据进行压缩
%输入：
%     x:            double列向量       x坐标
%     y:            double列向量       y坐标
%     num:          int               坐标总数
%     E:            double            弹性模量 elastic modulus  20
%     Po:           double            泊松比   Poisson's ratio  0.2
%     density:      double            质量密度  Mass density    2.7e10
%     percentage:   double            分位数   95
%     Maxit:        int               最大迭代次数 100
%     Tol:          double            误差阈值 1e-6
%     Thick:        int               松弛因子   1
%     force:        double            力单位  1
%输出：
%     x1:           double列向量       压缩后x坐标
%     y1:           double列向量       压缩后y坐标
%     ss:           int               坐标个数
%     Range:        二维列向量          原坐标与压缩后坐标的映射关系
%     time:         double列向量       每个坐标的时间戳
%例子：y=func(1,'type1');

%data=[x' y'];
data=[x y];
%time0=double(time0');
time0 = [1:1:size(x)]';
time0=double(time0);
A=data;
[r,c]=size(data);
No_nel=3;                 % number of nodes per element
No_dof=2;                 % number of dofs per node
Prop(1)=E;                % elastic modulus  20 
Prop(2)=Po;               % Poisson's ratio  0.2
Prop(3)=density;             % Mass density
t=Thick;                     % Thickness
%----------------------------------------------------------
% Input data for nodal coordinate values
% gcoord(i,j) where i-> node no. and j-> x or y
%----------------------------------------------------------
% Input data for nodal connectivity for each element 
% nodes(i,j) where i-> element no. and j-> connected nodes
%--------------------------------------------------------
[m1,n1]=max(data(:,1));
[m2,n2]=min(data(:,1));
[m3,n3]=max(data(:,2));
[m4,n4]=min(data(:,2));
Range=[m2 m1;m4 m3];
cx = mean(x);
cy = mean(y);
step1=(m1-m2)/29;
step2=(m3-m4)/29;
for i=1:30
    data(i+r,1)=(i-1)*step1+m2;
    data(i+r,2)=m4;
    data(i+r+30,1)=(i-1)*step1+m2;
    data(i+r+30,2)=m3;
    data(i+r+60,1)=m2;
    data(i+r+60,2)=(i-1)*step2+m4;
    data(i+r+90,1)=m1;
    data(i+r+90,2)=(i-1)*step2+m4;    
end
data(r+121:r+1020,1)=m2+(m1-m2)*rand([900,1]);
data(r+121:r+1020,2)=m4+(m3-m4)*rand([900,1]);
%data(r+121:r+9120,1)=m2+(m1-m2)*rand([9000,1]);
%data(r+121:r+9120,2)=m4+(m3-m4)*rand([9000,1]);
[Node_number,cc]=size(data);
[gcoord(:,2),index]=sort(data(:,2));
gcoord(:,1)=data(index(:,1),1);  
nodes = delaunay(gcoord(:,1),gcoord(:,2));
[row,col]=size(nodes);
Element_number=row;
%------------------------------------------------------
% Inout data for boundary conditions
%------------------------------------------------------
ed(1:Node_number,1:2)=1;            %elment_displacement 
constraint=[];
% force vector
f=zeros(Node_number*No_dof,1);
for i=1:Node_number
    f(2*i-1) = (cx - gcoord(i,1))*10000*force;
    f(2*i) = (cy - gcoord(i,2))*10000*force;
end
% for i=1:Node_number
%     if(gcoord(i,1)==m2)
%         f(2*i-1)=force;
%     else
%         f(2*i-1)=-force;
%     end
%     if(gcoord(i,2)==m4)
%         f(2*i)=force;
%     else
%         f(2*i)=-force;
%     end
% end
for loopi=1:length(constraint);
    ed(constraint(loopi,1),constraint(loopi,2))=0;
end
dof=0;
for loopi=1:Node_number
     for loopj=1:2
        if ed(loopi,loopj)~=0
            dof=dof+1;
            ed(loopi,loopj)=dof;
        end
    end
end
%--------------------------------------------------------
%Initilization of matrices 
%--------------------------------------------------------
k=sparse(dof,dof);         % system stiffness matrix
disp=zeros(dof,1);         % system displacement vector
e2s(1:6)=0;    % index of transform the elment displament number to structural
%---------------------------
Number_con=-1;                % apply boundary conditions
 for loopi=1:length(constraint)   
     Number_con= Number_con+1;
     f((constraint(loopi,1)-1)*No_dof+constraint(loopi,2)-Number_con)=[];
 end
%--------------------------------------------------------
% Compute system stiffness and mass matrices
%-------------------------------------------------------
 for loopi=1:Element_number            %  loop for the total number of elements
     for zi=1:3
        e2s((zi-1)*2+1)=ed(nodes(loopi,zi),1);
        e2s((zi-1)*2+2)=ed(nodes(loopi,zi),2);
     end
   for loopj=1:3
       xycoord(loopj,1)=gcoord(nodes(loopi,loopj),1);
       xycoord(loopj,2)=gcoord(nodes(loopi,loopj),2);
   end
  iopt=1;                               % plane stress analysis
  Opt_mass=1;                           % Consistent mass matrix
  %[dk,dm]=LineartriElement1(Prop,No_nel,No_dof,xycoord,t,iopt,Opt_mass);
  xycoord=double(xycoord);
  dk=LineartriElement1(Prop,No_nel,No_dof,xycoord,t,iopt);
    for jx=1:6
        for jy=1:6
            if(e2s(jx)*e2s(jy)~=0)
                k(e2s(jx),e2s(jy))=k(e2s(jx),e2s(jy))+dk(jx,jy);
            end
        end
    end
 end
 %--------------------------------------
 % compute system displacement
 %--------------------------------------
 %solve the matrix equation
 %disp=k\f;
 disp=gmres(k,f,Maxit,Tol,5);
 %disp=gmres(k,f);
 %------------------------------------------------------
 % Construct the taotal displacement
 %------------------------------------------------------
 %  print fem solutions
 %---------------------------------------------------
   for i=1:Node_number*No_dof
      if rem(i,2)~=0
          gc((i+1)/2,1)=disp(i);
      else
          gc(i/2,2)=disp(i);
      end
  end
  [v,index1]=sort(index);
  gcc(:,1)=gc(index1(:,1),1);
  gcc(:,2)=gc(index1(:,1),2);
  di(:,1)=gcoord(index1(:,1),1);
  di(:,2)=gcoord(index1(:,1),2);
 for i=1:r
    xy(i,1)=di(i,1)+gcc(i,1);
    xy(i,2)=di(i,2)+gcc(i,2);    
 end
dx=diff(xy,1,1);
dx1=sqrt(dx(:,1).^2+dx(:,2).^2);
inmax1=find(diff(sign(diff(A(:,1))))==-2)+1;
inmin1=find(diff(sign(diff(A(:,1))))==2)+1;
inmax2=find(diff(sign(diff(A(:,2))))==-2)+1;
inmin2=find(diff(sign(diff(A(:,2))))==2)+1;
local=[inmax1;inmin1;inmax2;inmin2];
local=sort(local);
local=unique(local);
ss=1;total=0;
x1(ss)=di(1,1);
y1(ss)=di(1,2);
time(1)=0;
dx2=(dx1-min(dx1))/(max(dx1)-min(dx1));
for i=1:length(dx2)-1
    if dx2(i,1)==0
        total=total+1;
    else
        total=total+1;
    if mod(i,10)~=0        
    if dx2(i,1)<double(prctile(dx2,double(percentage)))
        for j=1:double(length(local))
            if (i+1)==local(j,1)
                ss=ss+1;
                x1(ss)=di(i+1,1);
                y1(ss)=di(i+1,2);
                if total==1
                    time(ss)=time0(i+1,1)-time0(i,1);
                    total=0;
                else
                    time(ss)=time0(i+1,1)-time0(i-total+1,1);
                    total=0;
                end
            end
        end
    else
        ss=ss+1;
        x1(ss)=di(i+1,1);
        y1(ss)=di(i+1,2);
        if total==1
            time(ss)=time0(i+1,1)-time0(i,1);
            total=0;
        else 
            time(ss)=time0(i+1,1)-time0(i-total+1,1);
            total=0;
        end              
    end
    else
        ss=ss+1;
        x1(ss)=di(i+1,1);
        y1(ss)=di(i+1,2);
        if total==1
            time(ss)=time0(i+1,1)-time0(i,1);
            total=0;
        else
            time(ss)=time0(i+1,1)-time0(i-total+1,1);
            total=0;
        end
    end
    end
end
ss=ss+1;
total=total+1;
x1(ss)=di(length(dx1)+1,1);
y1(ss)=di(length(dx1)+1,2);
if total==1
    time(ss)=time0(length(dx1)+1,1)-time0(length(dx1),1);
    total=0;
else
    time(ss)=time0(length(dx1)+1,1)-time0(length(dx1)-total+1,1);
    total=0;
end