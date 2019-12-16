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
pl2=p;
pl3=p;
pl2(:,3)=pl2(:,3)+9;
pl3(:,3)=pl3(:,3)+18;
np=size(p,1);
np2=size(pl2,1);
np3=size(pl3,1);
na=0;
at=4;
bt=0;
ant=0;
%%%%%Boundary Size%%%%%
xlo=-a*5; xhi=a*5;
ylo=-b*5; yhi=b*5;
zlo=-c*40; zhi=c*40;
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Generate DPD Particles%%%%%%%%%%%%
dpd_n=200000;
dpd_x=rand(dpd_n,1)*(xhi-xlo)+xlo;
dpd_y=rand(dpd_n,1)*(yhi-ylo)+ylo;
dpd_z=rand(dpd_n,1)*(zhi-zlo)+zlo;
dpd_p=[dpd_x dpd_y dpd_z];
%%%%%Remove the particles overlap with platelet%%%%%%%%
rm_list=(dpd_p(:,1).^2/((a+0.5)*(a+0.5))+dpd_p(:,2).^2/((b+0.5)*(b+0.5))+dpd_p(:,3).^2/((c+0.5)*(c+0.5))-1<0);
dpd_p(rm_list,:)=[];
rm_list=(dpd_p(:,1).^2/((a+0.5)*(a+0.5))+dpd_p(:,2).^2/((b+0.5)*(b+0.5))+(dpd_p(:,3)-9).^2/((c+0.5)*(c+0.5))-1<0);
dpd_p(rm_list,:)=[];
rm_list=(dpd_p(:,1).^2/((a+0.5)*(a+0.5))+dpd_p(:,2).^2/((b+0.5)*(b+0.5))+(dpd_p(:,3)-18).^2/((c+0.5)*(c+0.5))-1<0);
dpd_p(rm_list,:)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('topo_data');
load('topo_data');
dpd_new_size=size(dpd_p,1);
fid=fopen('data.txt','w');
fprintf(fid,'%s\n','   ');
fprintf(fid,'%d %s\n',dpd_new_size+np+np2+np3,'atoms');
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
for k=1:np2
   nk=nk+1;
   fprintf(fid,'%d %d %d %d %d %d %d\n',nk,2,2,0,pl2(k,1),pl2(k,2),pl2(k,3)); 
end
for k=1:np3
   nk=nk+1;
   fprintf(fid,'%d %d %d %d %d %d %d\n',nk,3,3,0,pl3(k,1),pl3(k,2),pl3(k,3)); 
end
for k=1:dpd_new_size
   nk=nk+1;
   %if dpd_p(k,2)>(yhi-10) 
  %  fprintf(fid,'%d %d %d %d %d %d %d\n',nk,5,5,0,dpd_p(k,1),dpd_p(k,2),dpd_p(k,3)); 
  % elseif dpd_p(k,2)<(ylo+10) 
  %  fprintf(fid,'%d %d %d %d %d %d %d\n',nk,6,6,0,dpd_p(k,1),dpd_p(k,2),dpd_p(k,3)); 
  % else
         fprintf(fid,'%d %d %d %d %d %d %d\n',nk,4,4,0,dpd_p(k,1),dpd_p(k,2),dpd_p(k,3));  
  % end
 
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