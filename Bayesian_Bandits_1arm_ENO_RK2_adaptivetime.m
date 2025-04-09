%dtv+(mu+mu*dsv+dqv-mu2)pi+mu2=0
grid=[110;120;130;140;150;160;170];
errors=zeros(size(grid));
for j=1:size(grid,1)
    n=grid(j);
    differences=0;
    error=0;
    numericalminusexact=0;
    r=4;
    %k order of Runge-Kutta
    k=1;
    w=zeros(3*n-2,3*n-2);
    %policy=zeros(n,n,2);%numerical
    policy=zeros(3*n-2,3*n-2,2);
    %exactpolicy=zeros(n,n,2);%HJB equation
    exactpolicy=zeros(3*n-2,3*n-2,2);
    exactBayesianpolicy=zeros(3*n-2,3*n-2,2);%backwards differences
    ds=1/sqrt(n);
    dq=1/n;
    t=domofdep(dq,2,2,1,2,0,1);
    T=size(t,2);
    muhat=((-(2*k*T-2*k):(2*k*T-2*k))/sqrt(n)+1)'./((-(0):(2*k*T-2*k))/n+2);
    %muhat=zeros(r*k*4*n-3,r*k*3*n-2);
    %sigma=(s+a)./(q+b);
    %calculate exact policy
    %exactpolicy(:,:,1)=mu>=1/2*(1-ones(n,1)*(1:n));
    exactpolicy(:,:,1)=((-(n-1):(2*n-2))/n+1)'./((-(n-1):(2*n-2))/n+2)>=1/2;
    exactpolicy(:,:,2)=~exactpolicy(:,:,1);
    %calculate w
    %Z=running total number of cells
    Z=0;
    for i=T-1:-1:1
        dt=t(i+1)-t(i);
        if i==T-1
            %         w=max(mu,1/(2));
            v=dt*max(muhat,1/2);
            %         exactBayesianpolicy(:,:,1)=(mu>=1/(2));
            %         exactBayesianpolicy(:,:,2)=~exactBayesianpolicy(:,:,1);
            %         policy(:,:,1)=(mu>=1/2);
            %         policy(:,:,2)=~policy(:,:,1);
        else
            %         w1=mu(1:i,1:i)+mu(1:i,1:i).*w(2:(i+1),2:(i+1))+(1-mu(1:i,1:i)).*w(1:i,2:(i+1));
            %         w2=1/(2)+w(1:i,1:i);
            %         w=0;
            %         w(1:i,1:i)=max(w1,w2);
            temp=v;
            center=1+(0:2*(k*i-1));
            centers=2*k*(T-1)+1+(-2*(k*i-1):2*(k*i-1));
            [vsp,vsm,vqp]=ENO(v,k*i,T,k,n);
            %vsp=(v(centers+1,centerq)-v(centers,centerq))/ds;
            % vqp=(v(centers,centerq+1)-v(centers,centerq))/dq;
            %g1=v(1:i,1:i)+dt*(mu(1:i,1:i)+mu(1:i,1:i).*(v(2:(i+1),2:(i+1))-v(1:i,2:(i+1)))/ds+(v(1:i,2:(i+1))-v(1:i,1:i))/dq);
            % if i<=50
            %     pause(1)
            % end
            v_ss=(temp(centers+1,center)-2*temp(centers,center)+temp(centers-1,center))/ds^2;
            % v_ss=0;
            g1=v(centers,center)+dt*(muhat(centers,center)+...
                max(muhat(centers,center),0).*vsp+min(muhat(centers,center),0).*vsm+vqp+v_ss);
            %             g1=v(centers,centerq)+dt*(muhat(centers,centerq)+...
            %                 muhat(centers,centerq).*vsp+vqp);
            %g1=dt*mu(1:i,1:i)+mu(1:i,1:i).*v(2:(i+1),2:(i+1))+(1-mu(1:i,1:i)).*v(1:i,2:(i+1));
            g2=v(centers,center)+dt*1/2;
            %v(1:i,1:i)=v(1:i,1:i)+dt*max(g1,1/2);
            v(centers,center)=max(g1,g2);
            % %v(centers,centerq)=max(g1,g2);
            % center=1+(0:2*(k*i-2));
            % centers=2*k*(T-1)+1+(-2*(k*i-2):2*(k*i-2));
            % [tempsp,tempsm,tempqp]=ENO(temp,k*i-1,T,k,n);
            % v_ss=(temp(centers+1,center)-2*temp(centers,center)+temp(centers-1,center))/ds^2;
            % % v_ss=0;
            % g1=temp(centers,center)+dt*(muhat(centers,center)+...
            %     max(muhat(centers,center),0).*tempsp+min(muhat(centers,center),0).*tempsm+tempqp+v_ss);
            % g2=temp(centers,center)+dt*1/2;
            % temp(centers,center)=max(g1,g2);
            % v(centers,center)=1/2*(v(centers,center)+temp(centers,center));
        end
        %     differences=differences+norm(reshape(exactpolicy(1:i,1:i,:)-exactBayesianpolicy(1:i,1:i,:),1,2*i^2),1);
        %         if mod(i,10)==0
        %             surf((2*r*k*n-1):(2*r*k*n-1+i/k),r*k*n:(r*k*n+i/k),max(muhat((2*r*k*n-1):(2*r*k*n-1+i/k),r*k*n:(r*k*n+i/k)),1/2)*(1-(i-k)/(r*k*n))-v((2*r*k*n-1):(2*r*k*n-1+i/k),r*k*n:(r*k*n+i/k)))
        %             pause(1)
        %         end
        %         error=(Z*error+norm(reshape(max(muhat((2*r*k*n-1):(2*r*k*n-1+min(i/k,n)),r*k*n:(r*k*n+min(i/k,n))),1/2)*...
        %             (1-(i-k)/(r*k*n))-v((2*r*k*n-1):(2*r*k*n-1+min(i/k,n)),r*k*n:(r*k*n+min(i/k,n))),...
        %             1,(min(i/k,n)+1)^2),1))/(Z+(min(i/k,n)+1)^2);
        center=1+(0:min(i-1,n));
        centers=2*k*(T-1)+1+(0:min(i-1,n));
        error=(Z*error+norm(reshape(max(muhat(centers,center),1/2)*(1-t(i))-...
            v(centers,center),...
            1,(min(i-1,n)+1)^2),1))/(Z+(min(i-1,n)+1)^2);
        %         error=(Z*error+norm(reshape(max(muhat(1:i,1:i),1/2)*(1-(i-1)/(r*k*n))-...
        %             v(1:i,1:i),...
        %             1,i^2),1))/(Z+i^2);
        Z=Z+(min(i-1,n)+1)^2;
        %     numericalminusexact=numericalminusexact+norm(reshape(policy(n+(-(i-1):2*(i-1)),n+(-(i-1):2*(i-1)),:)-exactpolicy(n+(-(i-1):2*(i-1)),n+(-(i-1):2*(i-1)),:),1,2*(3*i-2)^2),1);
    end
    %12*differences/(n*(n+1)*(2*n+1))
    %error
    %     i=r*k*(n-2);
    %     errors(j)=norm(reshape(max(muhat((2*r*k*n-1):(2*r*k*n-1+min(i/k,n)),r*k*n:(r*k*n+min(i/k,n))),1/2)*...
    %         (1-(i-k)/(r*k*n))-v((2*r*k*n-1):(2*r*k*n-1+min(i/k,n)),r*k*n:(r*k*n+min(i/k,n))),...
    %         1,(min(i/k,n)+1)^2),1)/(min(i/k,n)+1)^2;
    errors(j)=error;
end
disp('yes diffusion, ds=sqrt(n)')
polyfit(log(grid),log(errors),1)
% ans =
% 
%    -1.2548   -5.0948
% grid=[20;40;60;80;100;120;140;160;180;200;220];
% errors =
% 
%    1.0e-03 *
% 
%     0.1482
%     0.0593
%     0.0353
%     0.0245
%     0.0187
%     0.0148
%     0.0123
%     0.0105
%     0.0091
%     0.0081
%     0.0073
%with diffusion RK=1, domofdep(d1,2,2), ENO=2
% ans =
% 
%    -0.5071   -4.3644
% 
% errors
% grid=[20;40;60;80;100];
% errors =
% 
%     0.0028
%     0.0020
%     0.0016
%     0.0014
%     0.0012
% ans =
% 
%    -0.5079   -4.3613
% 
% errors
% grid=[110;120;130;140;150;160;170];
% errors =[0.0012;
%     0.0011;
%     0.0011;
%     0.0010;
%     0.0010;
%     0.0010;
%     0.0009];
function [vsp,vsm,vqp]=ENO(v,i,T,k,n)
ds=1/sqrt(n);
dq=1/n;
center=1+(0:2*(i-1));
left=[1,1+(0:(2*i-3))];
centers=2*k*(T-1)+1+(-2*(i-1):2*(i-1));
rs=v(centers+2,center)-2*v(centers+1,center)+v(centers,center);
ms=v(centers+1,center)-2*v(centers,center)+v(centers-1,center);
ls=v(centers,center)-2*v(centers-1,center)+v(centers-2,center);
choice=(abs(rs)<abs(ms));
%ps/(2*ds^2).*((x-(-(i-2):2*i-1)*ds)*ones(3*i-2,1)).*((x-(-(i-1):2*(i-1))*ds)*ones(3*i-2,1))+(v(-(i-2):2*i-1,-(i-1):2*(i-1))-v(-(i-1):2*(i-1),-(i-1):2*(i-1)))/ds.*((x-(-(i-1):2*(i-1))*ds)*ones(3*i-2,1));
vsright=-rs/(2*ds)+(v(centers+1,center)-v(centers,center))/ds;
vscenter=ms/(2*ds)+(v(centers,center)-v(centers-1,center))/ds;
vsleft=3*ls/(2*ds)+(v(centers-1,center)-v(centers-2,center))/ds;
vsp=choice.*vsright+(~choice).*vscenter;
vsp=vscenter;
% vsp(1,:)=(v(2,center)-v(1,center))/ds;
% vsp=vsright;
choice=(abs(ls)<abs(ms));
vsm=choice.*vsleft+(~choice).*vscenter;
vsm=vscenter;
rq=v(centers,center+2)-2*v(centers,center+1)+v(centers,center);
mq=v(centers,center+1)-2*v(centers,center)+v(centers,left);
choice=(abs(rq)<abs(mq));
%ps/(2*ds^2).*((x-(-(i-2):2*i-1)*ds)*ones(3*i-2,1)).*((x-(-(i-1):2*(i-1))*ds)*ones(3*i-2,1))+(v(-(i-2):2*i-1,-(i-1):2*(i-1))-v(-(i-1):2*(i-1),-(i-1):2*(i-1)))/ds.*((x-(-(i-1):2*(i-1))*ds)*ones(3*i-2,1));
vqright=-rq/(2*dq)+(v(centers,center+1)-v(centers,center))/dq;
vqcenter=mq/(2*dq)+(v(centers,center)-v(centers,left))/dq;
% vqp=vqright;
vqp=choice.*vqright+(~choice).*vqcenter;
vqp(:,1)=(v(centers,2)-v(centers,1))/dq;
% vqp=(v(center,center+1)-v(center,center))/dq;
end
%12*numericalminusexact/(n*(n+1)*(2*n+1))
%calculate indices used for ENO vectorized
% index=zeros(n,n,n,n);
% indexs1p1=zeros(n,n,n,n);
% indexs1m1=zeros(n,n,n,n);
% indexs1m2=zeros(n,n,n,n);
% indexs2p1=zeros(n,n,n,n);
% indexs2m1=zeros(n,n,n,n);
% indexs2m2=zeros(n,n,n,n);
% indexq1m1=zeros(n,n,n,n);
% indexq1p1=zeros(n,n,n,n);
% indexq2m1=zeros(n,n,n,n);
% indexq2p1=zeros(n,n,n,n);
% for j=1:n
%     for k=1:n
%         for l=1:n
%             index(:,j,k,l)=1:n;
%             indexs1p1(:,j,k,l)=2:(n+1);
%             indexs1m1(:,j,k,l)=0:(n-1);
%             indexs1m2(:,j,k,l)=-1:(n-2);
%             indexs2p1(j,:,k,l)=2:(n+1);
%             indexs2m1(j,:,k,l)=0:(n-1);
%             indexs2m2(j,:,k,l)=-1:(n-2);
%             indexq1m1(j,k,:,l)=0:(n-1);
%             indexq1p1(j,k,:,l)=2:(n+1);
%             indexq2m1(j,k,l,:)=0:(n-1);
%             indexq2p1(j,k,l,:)=2:(n+1);
%         end
%     end
% end
% for i=n:-1:1
%     %Upwind from eq (34) in Yuhua's Paper
%     g1= (v(right,:,right,:)-v(:,:,right,:)).*mu1_max/ds+...
%         (v(:,:,right,:)-v(left,:,right,:)).*mu1_min/ds+...
%         (v(:,:,right,:)-v(:,:,:,:))/dq+mu1;
%     g2= (v(:,right,:,right)-v(:,:,:,right)).*mu2_max/ds+...
%         (v(:,:,:,right)-v(:,left,:,right)).*mu2_min/ds+...
%         (v(:,:,right,:)-v(:,:,:,:))/dq+mu2;
%     v=v+dt*max(g1,g2);
%     %Upwind from eq (34) in Yuhua's Paper centered at j
%     %     g1= (v(right,:,:,:)-v(:,:,:,:)).*mu1_max/ds+...
%     %         (v(:,:,:,:)-v(left,:,:,:)).*mu1_min/ds+...
%     %         (v(:,:,right,:)-v(:,:,:,:))/dq+mu1;
%     %     g2= (v(:,right,:,:)-v(:,:,:,:)).*mu2_max/ds+...
%     %         (v(:,:,:,:)-v(:,left,:,:)).*mu2_min/ds+...
%     %         (v(:,:,right,:)-v(:,:,:,:))/dq+mu2;
%     %ENO2
%     %     [v_s1_right,v_s1_left,v_s2_right,v_s2_left,v_q1,v_q2]=eno2vec(v,ds,dq,index,indexs1m1,indexs1m2,indexs1p1,...
%     %         indexs2m1,indexs2m2,indexs2p1,indexq1m1,indexq1p1,indexq2m1,indexq2p1);
%     %     g1= v_s1_right.*mu1_max+...
%     %         v_s1_left.*mu1_min+...
%     %         v_q1+mu1;
%     %     g2= v_s2_right.*mu2_max+...
%     %         v_s2_left.*mu2_min+...
%     %         v_q2+mu2;
%     %      v1=v+dt*max(g1,g2);
%     %     [v_s1_right,v_s1_left,v_s2_right,v_s2_left,v_q1,v_q2]=eno2vec(v1,ds,dq,index,indexs1m1,indexs1m2,indexs1p1,...
%     %         indexs2m1,indexs2m2,indexs2p1,indexq1m1,indexq1p1,indexq2m1,indexq2p1);
%     %     g1= v_s1_right.*mu1_max+...
%     %         v_s1_left.*mu1_min+...
%     %         v_q1+mu1;
%     %     g2= v_s2_right.*mu2_max+...
%     %         v_s2_left.*mu2_min+...
%     %         v_q2+mu2;
%     %     v=1/2*v+1/2*(v1+dt*max(g1,g2));
%     policy(:,:,:,:,i,1)=(g1>=g2);
%     policy(:,:,:,:,i,2)=~(g1>=g2);
% end
% differences=differences+norm(reshape(policy-exactBayesianpolicy,1,2*n^5),1)/(2*n^5)