function [ q,C,L,sol,vars,tempos,emin,emind,soll ] = lagrange( f,var,c )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[r,~]=size(var);
sym('lamda','real');
for i=1:r
    vars(i,1)=sym(var(i,1),'real');
end
[r1,~]=size(c);
for i=1:r1
    s=horzcat('lamda',num2str(i));
    tvars(i,1)=sym(s,'real');
end
C=0;
for i=1:r1
    C=C+(tvars(i,1)*c(i,1));
end
L=f+C;
for i=1:r
    q(i,1)=diff(L,var(i,1));
end
q=cat(1,q,c);
sol=solve(q);
temp=fieldnames(sol);




[r22,~]=size(sol.(temp{1}));
rap=f*'lamda';
tempos=zeros(r22,1);
for i=1:r22
    for j=1:r+r1
        rap=subs(rap,temp{j},sol.(temp{j})(i));
    end
if(rap==0)
    tempos(i,1)=rap;
else
    tempos(i,1)=rap./'lamda';
end
rap=f*'lamda';
end
[fmax,~]=max(tempos);
index=find(tempos==fmax);
disp('The max value of function is:')
disp(fmax)
disp('The max value occurs at:')
[siz,~]=size(index);
for i=1:siz
    for j=r1+1:r+r1
        disp(sol.(temp{j})(index(i)))
    end
end
[fmin,~]=min(tempos);
index1=find(tempos==fmin);
disp('The min value of function is:')
disp(fmin)
disp('The min value occurs at:')
[siz,~]=size(index1);
for i=1:siz
    for j=(r1+1):(r+r1)
        disp(sol.(temp{j})(index1(i)))
    end
end

% hessian

hes=hessian(L);
grad=sym(zeros(r,r1));
for i=1:r1
    for j=1:r
        grad(j,i)=diff(c(i,1),var(j,1));
    end
end

emi=zeros(r1,r1);
iden=sym(eye(r));
for i=1:r
        iden(i,i)=iden(i,i).*'lamda';
end
he=hes-iden;
temi=cat(2,he,grad);
grad=grad';
ttemi=cat(2,grad,emi);
emin=cat(1,temi,ttemi);
emind=det(emin);


final=zeros(r22,1);
for i=1:r22
    poh=emind;
     for j=1:r+r1
        poh=subs(poh,temp{j},sol.(temp{j})(i));
     end
     soll=solve(poh);
     [rt,~]=size(soll);
     for kite=1:rt
         if(soll(kite)>=0)
             final(i)=1;
         end
     end
     
end

disp('critical points that give positive roots of determinant of Hessian are:')

if(sum(final==1)==0)
    disp('there are NO critical points that give positive roots of determinant of Hessian')

else
cntr=0;
for i=1:r22
    if(final(i)==1)
        cntr=cntr+1;
        fprintf('Point %d\n',cntr);
        for j=r1+1:r+r1
            disp(sol.(temp{j})(i))
        end
    end
end
end

[pl,~]=size(index);
[pl1,~]=size(index1);
jacob=0;
jacob1=0;
for i=1:pl
    tlp1=index(pl);
    if(final(tlp1)==1)
        jacob=1;
    end
end
end