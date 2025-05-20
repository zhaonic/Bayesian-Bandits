%dtv+mu+mu*dsv+dqv+mu+mu_2*ds_2v+dq_2v=0
grid=[5;10;15;20;30];
% grid=[10];
% grid=[7];
errors=zeros(size(grid));
siginv=1;
a=[0.1;-0.1];
for j=1:size(grid,1)
    n=grid(j);
    differences=0;
    error=0;
    %k order of Runge-Kutta
    k=1;
    ds=1/n;
    dq=1/n;
    Z=0;
    [t,~,~]=domofdep(ds,1,1,1,1,0,1);
    t=0:dq:1;
    fine=120;
    t=0:(1/fine):1;
    T=size(t,2);
    v=zeros(k*n+1,k*n+1,k*n+1);
    sol=zeros(T,T,T);
    s1=(1:T)';
    s2=1:T;
    q1=reshape(1:T,[1,1,T]);
    %r=ratio of numerical grid for ENO-RK compared to "exact" first order
    %grid
    r=fine/(k*n);
    errorarray=[];
    for i=T-1:-1:(T-2*r*k-1)
        s1=(1:T)';
        s2=1:T;
        q1=reshape(1:T,[1,1,T]);
        muhat=repmat(a(1)*s1+a(2)*s2,[1,1,length(q1)])./repmat(1+q1*a(1)^2+(i-q1)*a(2)^2,[length(s1),length(s2),1]);
        % muhat=ones(T,T,T);
        dt=t(i+1)-t(i);
        s1=(1+(0:(i-1)))';
        s2=(1+(0:(i-1)));
        q1=reshape(1+(0:(i-1)),[1,1,i]);
        sol(s1,s2,q1)=sol(s1,s2,q1)+dt*max(a(1)*muhat(s1,s2,q1)+a(1)*muhat(s1,s2,q1).*(sol(s1+1,s2,q1)-sol(s1,s2,q1))/ds+(sol(s1,s2,q1+1)-sol(s1,s2,q1))/dq,...
            a(2)*muhat(s1,s2,q1)+a(2)*muhat(s1,s2,q1).*(sol(s1,s2+1,q1)-sol(s1,s2,q1))/ds);
        % sol(s1,s2,q1)=sol(s1,s2,q1)+dt*muhat(s1,s2,q1);
        count=((T-i)/(k*r));
        if count==round(count) && count>0
            temp=v;
            % s1=(1+(0:r:(i-1+r)))';
            % s2=(1+(0:r:(i-1+r)));
            % q1=reshape(1+(0:r:(i-1+r)),[1,1,(i-1+r)/r+1]);
            % s1r=(1+(0:(i-1+r)/r))';
            % s2r=(1+(0:(i-1+r)/r));
            % q1r=reshape(1+(0:(i-1+r)/r),[1,1,(i-1+r)/r+1]);
            % [vs1p,vs2p,vq1p]=ENO(v,i-1+r,n,r);
            % % vqp=-muhat(index,index).*vsp;
            % temp(s1r,s2r,q1r)=v(s1r,s2r,q1r)+dt*r*max(a(1)*muhat(s1r,s2r,q1r)+a(1)*muhat(s1r,s2r,q1r).*vs1p+vq1p,...
            %     a(2)*muhat(s1r,s2r,q1r)+a(2)*muhat(s1r,s2r,q1r).*vs2p);
            % % temp(s1r,s2r,q1r)=v(s1r,s2r,q1r)+dt*r*k*muhat(s1,s2,q1);
            s1=(1+(0:r:(i-1)))';
            s2=(1+(0:r:(i-1)));
            q1=reshape(1+(0:r:(i-1)),[1,1,(i-1)/r+1]);
            s1r=(1+(0:(i-1)/r))';
            s2r=(1+(0:(i-1)/r));
            q1r=reshape(1+(0:(i-1)/r),[1,1,(i-1)/r+1]);
            [vs1p,vs2p,vq1p]=ENO(temp,i-1,n,r);
            % vqp=-muhat(index,index).*vsp;
            v(s1r,s2r,q1r)=temp(s1r,s2r,q1r)+dt*r*max(a(1)*muhat(s1r,s2r,q1r)+a(1)*muhat(s1r,s2r,q1r).*vs1p+vq1p,...
                a(2)*muhat(s1r,s2r,q1r)+a(2)*muhat(s1r,s2r,q1r).*vs2p);
            % temp(s1r,s2r,q1r)=temp(s1r,s2r,q1r)+dt*r*k*muhat(s1,s2,q1);
            % v(s1r,s2r,q1r)=v(s1r,s2r,q1r)+dt*r*k*muhat(s1,s2,q1);
            % v(s1r,s2r,q1r)=(v(s1r,s2r,q1r)+temp(s1r,s2r,q1r))/2;
            error=(Z*error+norm(reshape(sol(s1,s2,q1)-v(s1r,s2r,q1r),[],1,1),1))/(Z+length(s1)^3);
            var=sol(s1,s2,q1)-v(s1r,s2r,q1r);
            % errorarray=[errorarray,norm(reshape(sol(s1,s2,q1)-v(s1r,s2r,q1r),[],1,1),1)/length(s1)^3];
            Z=Z+length(s1)^3;
        end
        %         error=(Z*error+norm(reshape(max(muhat(1:i,1:i),1/2)*(1-(i-1)/(r*k*n))-...
        %             v(1:i,1:i),...
        %             1,i^2),1))/(Z+i^2);
    end
    %12*differences/(n*(n+1)*(2*n+1))
    %error
    %     i=r*k*(n-2);
    %     errors(j)=norm(reshape(max(muhat((2*r*k*n-1):(2*r*k*n-1+min(i/k,n)),r*k*n:(r*k*n+min(i/k,n))),1/2)*...
    %         (1-(i-k)/(r*k*n))-v((2*r*k*n-1):(2*r*k*n-1+min(i/k,n)),r*k*n:(r*k*n+min(i/k,n))),...
    %         1,(min(i/k,n)+1)^2),1)/(min(i/k,n)+1)^2;
    errors(j)=error;
end
polyfit(log(grid),log(errors),1)
% Linear_Bandit
% 
% ans =
% 
%    -1.0406   -1.4175
% 
% errors
% 
% errors =
% 
%     0.0433
%     0.0230
%     0.0151
%     0.0109
%     0.0066
% domofdep(ds,1,2,1,1,0,1)
% grid
% errors
% Bayesian_Bandits_1storder_sq0_nomax_adapt
%
% ans =
%
%    -1.8776   -6.8830
%
%
% grid =
%
%     10
%     20
%     30
%     40
%     50
%     60
%     70
%     80
%     90
%    100
%
%
% errors =
%
%    1.0e-04 *
%
%     0.1288
%     0.0374
%     0.0182
%     0.0104
%     0.0067
%     0.0047
%     0.0035
%     0.0027
%     0.0021
%     0.0017
%with exact bounds v_s at q:
% grid=[10;20;30;40;50;60;70;80;90;100;110;120];
% ans =
%
%    -1.9479   -4.8945
% errors =
%
%    1.0e-04 *
%
%     0.8165
%     0.2171
%     0.1004
%     0.0578
%     0.0376
%     0.0263
%     0.0195
%     0.0149
%     0.0117
%     0.0095
%     0.0078
%     0.0041
%doing v_s(,j) was blowing up
%after correcting to v_s(,j+e_j)
% grid=[10;20;30;40;50;60;70;80;90;100;110;120];
% ans =
%
%     0.0612  -38.5669
%
% errors
%
% errors =
%
%    1.0e-16 *
%
%     0.2132
%     0.2197
%     0.2108
%     0.2150
%     0.2236
%     0.2246
%     0.2314
%     0.2221
%     0.2354
%     0.2326
%     0.2469
%     0.2557
% adaptive ENO1 RK2
% grid=[10;20;30;40;50;60;70;80;90;100;110;120;130];
% ans =
%
%    -1.0165   -5.9284
%
%
% errors =
%
%    1.0e-03 *
%
%     0.2612
%     0.1252
%     0.0823
%     0.0631
%     0.0499
%     0.0413
%     0.0353
%     0.0307
%     0.0276
%     0.0248
%     0.0224
%     0.0205
%     0.0191

function [vs1p,vs2p,vqp]=ENO(v,i,n,r)
ds=r/n;
dq=r/n;
index=1+(0:(i/r));
order=2;
vs1p=(v(index+1,index,index)-v(index,index,index))/ds;
vs2p=(v(index,index+1,index)-v(index,index,index))/ds;
vqp=(v(index,index,index+1)-v(index,index,index))/dq;
% centers=3*(r*k*n-1)+1+(-3*(i-1):3*(i-1));
% if (i/r)>=2*order
if 1==0
    centers1=1+(order:(i/r-order))';
    centers2=1+(order:(i/r-order));
    centerq=reshape(1+(order:(i/r-order)),1,1,[]);
    allq=reshape(1+(0:(i/r)),1,1,[]);
    alls1=(1+(0:(i/r)))';
    alls2=1+(0:(i/r));
    indexs1(:,:,:,1)=repmat(centers1,1,length(alls2),length(allq));
    indexs1(:,:,:,2)=repmat(centers1+1,1,length(alls2),length(allq));
    indexs2(:,:,:,1)=repmat(centers2,length(alls1),1,length(allq));
    indexs2(:,:,:,2)=repmat(centers2+1,length(alls1),1,length(allq));
    table(:,:,:,1)=v(centers1,alls2,allq);
    table(:,:,:,2)=(v(centers1+1,alls2,allq)-v(centers1,alls2,allq))/ds;
    for ii=2:order
        ip=max(indexs1,[],4)+1;
        im=min(indexs1,[],4)-1;
        left=sub2ind(size(v),reshape(im,[],1),reshape(repmat(alls2,length(centers1),1,length(allq)),[],1),...
            reshape(repmat(allq,length(centers1),length(alls2),1),[],1));
        right=sub2ind(size(v),reshape(ip,[],1),reshape(repmat(alls2,length(centers1),1,length(allq)),[],1),...
            reshape(repmat(allq,length(centers1),length(alls2),1),[],1));
        % tempm=reshape(v(left),size(centers,1),size(allq,2));
        % tempp=reshape(v(right),size(centers,1),size(allq,2));
        tempm=reshape(v(left),length(centers1),length(alls2),[]);
        tempp=reshape(v(right),length(centers1),length(alls2),[]);
        for jj=1:ii
            tempm=(tempm-table(:,:,:,jj))./(im-indexs1(:,:,:,jj))/ds;
            tempp=(tempp-table(:,:,:,jj))./(ip-indexs1(:,:,:,jj))/ds;
        end
        choice=abs(tempp)<abs(tempm);
        table(:,:,:,ii+1)=choice.*tempp+(~choice).*tempm;
        indexs1(:,:,:,ii+1)=choice.*ip+(~choice).*im;
        vs1p(centers1,:,:)=vs1p(centers1,:,:)+table(:,:,:,ii+1).*prod((repmat(indexs1(:,:,:,1),1,1,1,ii-1)-indexs1(:,:,:,2:ii))*ds,4);
    end
    clear table
    table(:,:,:,1)=v(alls1,centers2,allq);
    table(:,:,:,2)=(v(alls1,centers2+1,allq)-v(alls1,centers2,allq))/ds;
    for ii=2:order
        ip=max(indexs2,[],4)+1;
        im=min(indexs2,[],4)-1;
        left=sub2ind(size(v),reshape(repmat(alls1,1,length(centers2),length(allq)),[],1),reshape(im,[],1),...
            reshape(repmat(allq,length(alls1),length(centers2),1),[],1));
        right=sub2ind(size(v),reshape(repmat(alls1,1,length(centers2),length(allq)),[],1),reshape(ip,[],1),...
            reshape(repmat(allq,length(alls1),length(centers2),1),[],1));
        % tempm=reshape(v(left),size(centers,1),size(allq,2));
        % tempp=reshape(v(right),size(centers,1),size(allq,2));
        tempm=reshape(v(left),length(alls1),length(centers2),[]);
        tempp=reshape(v(right),length(alls1),length(centers2),[]);
        for jj=1:ii
            tempm=(tempm-table(:,:,:,jj))./(im-indexs2(:,:,:,jj))/ds;
            tempp=(tempp-table(:,:,:,jj))./(ip-indexs2(:,:,:,jj))/ds;
        end
        choice=abs(tempp)<abs(tempm);
        table(:,:,:,ii+1)=choice.*tempp+(~choice).*tempm;
        indexs2(:,:,:,ii+1)=choice.*ip+(~choice).*im;
        vs2p(:,centers2,:)=vs2p(:,centers2,:)+table(:,:,:,ii+1).*prod((repmat(indexs2(:,:,:,1),1,1,1,ii-1)-indexs2(:,:,:,2:ii))*ds,4);
    end
    indexq(:,:,:,1)=repmat(centerq,length(alls1),length(alls2),1);
    indexq(:,:,:,2)=repmat(centerq+1,length(alls1),length(alls2),1);
    clear table
    table(:,:,:,1)=v(alls1,alls2,centerq);
    table(:,:,:,2)=(v(alls1,alls2,1+centerq)-v(alls1,alls2,centerq))/dq;
    for ii=2:order
        ip=max(indexq,[],4)+1;
        im=min(indexq,[],4)-1;
        left=sub2ind(size(v),reshape(repmat(alls1,1,length(alls2),length(centerq)),[],1),...
            reshape(repmat(alls2,length(alls1),1,length(centerq)),[],1),reshape(im,[],1));
        right=sub2ind(size(v),reshape(repmat(alls1,1,length(alls2),length(centerq)),[],1),...
            reshape(repmat(alls2,length(alls1),1,length(centerq)),[],1),reshape(ip,[],1));
        % tempm=reshape(v(left),size(centers,1),size(allq,2));
        % tempp=reshape(v(right),size(centers,1),size(allq,2));
        tempm=reshape(v(left),length(alls1),length(alls2),[]);
        tempp=reshape(v(right),length(alls1),length(alls2),[]);
        for jj=1:ii
            tempm=(tempm-table(:,:,:,jj))./(im-indexq(:,:,:,jj))/ds;
            tempp=(tempp-table(:,:,:,jj))./(ip-indexq(:,:,:,jj))/ds;
        end
        choice=abs(tempp)<abs(tempm);
        table(:,:,:,ii+1)=choice.*tempp+(~choice).*tempm;
        indexq(:,:,:,ii+1)=choice.*ip+(~choice).*im;
        vqp(:,:,centerq)=vqp(:,:,centerq)+table(:,:,:,ii+1).*prod((repmat(indexq(:,:,:,1),1,1,1,ii-1)-indexq(:,:,:,2:ii))*dq,4);
    end
end
% exactright=(((n+1:n+2)/n+1)'./((0:(n))/n+2))*(1-(i)/(r*n));
% exactdown=(((0:(i-1))/n+1)'./((-2:-1)/n+2))*(1-(i)/(r*n));
% exactup=(((0:n)/n+1)'./((n+1:n+2)/n+2))*(1-(i)/(r*n));
% down=[exactdown(:,2),v(index,1:(i-1))];
% up=[v(:,2:(n+1)),exactup(:,1)];
% up2=[v(:,3:(n+1)),exactup];
% left=[exactleft(2,:);v(1:(i-1),index)];
% left2=[exactleft;v(1:(n-1),:)];
% right=[v(2:(n+1),:);exactright(1,:)];
% right2=[v(3:(n+1),:);exactright];
% vsp=(right-up)/ds;
% vsm=(up-left)/ds;
% vqp=(up-v)/dq;
% rs=v(index2+2,index)-2*v(index2+1,index)+v(index2,index);
% ms=v(index2+1,index)-2*v(index2,index)+v(down,index);
% % ls=up-2*left+left2;
% choice=(abs(rs)<abs(ms));
% %ps/(2*ds^2).*((x-(-(i-2):2*i-1)*ds)*ones(3*i-2,1)).*((x-(-(i-1):2*(i-1))*ds)*ones(3*i-2,1))+(v(-(i-2):2*i-1,-(i-1):2*(i-1))-v(-(i-1):2*(i-1),-(i-1):2*(i-1)))/ds.*((x-(-(i-1):2*(i-1))*ds)*ones(3*i-2,1));
% vsright=-rs/(2*ds)+(v(index2+1,index)-v(index2,index))/ds;
% vscenter=ms/(2*ds)+(v(index2,index)-v(down,index))/ds;
% % vsleft=3*ls/(2*ds)+(left-left2)/ds;
% vsp=zeros(i);
% vsp(index2,index)=choice.*vsright+(~choice).*vscenter;
% vsp(1,index)=(v(2,index)-v(1,index))/ds;
% vsp(i,index)=(v(i+1,index)-v(i,index))/ds;
% % vsp=(right-v)/ds;
% % choice=(abs(ls)<abs(ms));
% uq=v(index,index2+2)-2*v(index,index2+1)+v(index,index2);
% mq=v(index,index2+1)-2*v(index,index2)+v(index,down);
% choice=(abs(uq)<abs(mq));
% % %choice=1;
% % %ps/(2*ds^2).*((x-(-(i-2):2*i-1)*ds)*ones(3*i-2,1)).*((x-(-(i-1):2*(i-1))*ds)*ones(3*i-2,1))+(v(-(i-2):2*i-1,-(i-1):2*(i-1))-v(-(i-1):2*(i-1),-(i-1):2*(i-1)))/ds.*((x-(-(i-1):2*(i-1))*ds)*ones(3*i-2,1));
% vqup=-uq/(2*dq)+(v(index,index2+1)-v(index,index2))/dq;
% vqcenter=mq/(2*dq)+(v(index,index2)-v(index,down))/dq;
% vqp=zeros(i);
% vqp(index,index2)=choice.*vqup+(~choice).*vqcenter;
% vqp(index,1)=(v(index,2)-v(index,1))/dq;
% vqp(index,i)=(v(index,i+1)-v(index,i))/dq;
end