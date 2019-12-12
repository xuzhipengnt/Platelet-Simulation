clear;
a=10.5; %Ellipsoid - a
b=10.5; %Ellipsoid - b
c=3.5;  %Ellipsoid - c
c2=c/a*(a-0.5);
c3=c/a*(a-1);
deltalength=1;
fd=@(p) p(:,1).^2/(a*a)+p(:,2).^2/(b*b)+p(:,3).^2/(c*c)-1;  % Generate Mesh Data --Layer 1
fd2=@(p) p(:,1).^2/((a-0.5)*(a-0.5))+p(:,2).^2/((b-0.5)*(b-0.5))+p(:,3).^2/(c2*c2)-1;  % Generate Mesh Data --Layer 2
fd3=@(p) p(:,1).^2/((a-1)*(a-1))+p(:,2).^2/((b-1)*(b-1))+p(:,3).^2/(c3*c3)-1;  % Generate Mesh Data --Layer 3
[p1,t]=distmeshsurface(fd,@huniform,deltalength,[-a-2,-b-2,-c-2; a+2,b+2,c+2]);  % p -- Position Data  t -- Triangle Data
[p2,t]=distmeshsurface(fd2,@huniform,deltalength-0.1,[-a-2,-b-2,-c-2; a+2,b+2,c+2]);  % p -- Position Data  t -- Triangle Data
%[p3,t]=distmeshsurface(fd3,@huniform,deltalength-0.2,[-a-2,-b-2,-c-2; a+2,b+2,c+2]);  % p -- Position Data  t -- Triangle Data
p=[p1;p2];
np=size(p,1)
na=0;
at=2;
bt=0;
ant=0;
%%%%%Boundary Size%%%%%
xlo=-a*10; xhi=a*10;
ylo=-b*10; yhi=b*10;
zlo=-c*10; zhi=c*19;
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Generate DPD Particles%%%%%%%%%%%%
dpd_n=100000;
dpd_x=rand(dpd_n,1)*(xhi-xlo)+xlo;
dpd_y=rand(dpd_n,1)*(yhi-ylo)+ylo;
dpd_z=rand(dpd_n,1)*(zhi-zlo)+zlo;
dpd_p=[dpd_x dpd_y dpd_z];
%%%%%Remove the particles overlap with platelet%%%%%%%%
rm_list=(dpd_p(:,1).^2/((a+0.5)*(a+0.5))+dpd_p(:,2).^2/((b+0.5)*(b+0.5))+dpd_p(:,3).^2/((c+0.5)*(c+0.5))-1<0);
dpd_p(rm_list,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('topo_data');
load('topo_data');
dpd_new_size=size(dpd_p,1);
fid=fopen('data.txt','w');
fprintf(fid,'%s\n','   ');
fprintf(fid,'%d %s\n',dpd_new_size+np,'atoms');
fprintf(fid,'%s\n','# 1 WH');
%fprintf(fid,'%d %s\n',na,'angles');
fprintf(fid,'%d %s\n',at,'atom types');
%fprintf(fid,'%d %s\n',ant,'angle types');
fprintf(fid,'%d %d %s\n',xlo,xhi,'xlo xhi');
fprintf(fid,'%d %d %s\n',ylo,yhi,'ylo yhi');
fprintf(fid,'%d %d %s\n',zlo,zhi,'zlo zhi');
fprintf(fid,'\n');
fprintf(fid,'%s\n','Atoms');
fprintf(fid,'\n');
nk=0;
for k=1:np
   nk=nk+1;
   fprintf(fid,'%d %d %d %d %d %d %d\n',nk,1,1,0,p(k,1),p(k,2),p(k,3)); 
end

for k=1:dpd_new_size
   nk=nk+1;
   fprintf(fid,'%d %d %d %d %d %d %d\n',nk,2,2,0,dpd_p(k,1),dpd_p(k,2),dpd_p(k,3)); 
end
%
%for k=1:(nb/2);
%  nbk=nbk+1;
%  fprintf(fid,'%d %d %d %d\n',nbk,1,t(k,1),t(k,2));
%  nbk=nbk+1;
%  fprintf(fid,'%d %d %d %d\n',nbk,1,t(k,2),t(k,3));
%end
%fprintf(fid,'%s\n','Angles');
fclose(fid);