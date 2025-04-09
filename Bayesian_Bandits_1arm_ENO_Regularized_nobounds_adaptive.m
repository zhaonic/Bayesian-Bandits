%v_t+lambda*log(exp(1/lambda*(mu*v_s+v_q+mu))+exp(1/lambda*mu_2))=0
grid=[600];
lambdas=[0.0025];
grid=[200;220]*4;
lambdas=[0.005];
% grid=[120,120,120,120,120,120,120,120,120;
%     ];
% lambdas=[1,1/10,(1/10)^2,0.01*(1/2),0.01*(1/2)^2,0.01*(1/2)^3,0.01*(1/2)^4,0.01*(1/2)^5,0.01*(1/2)^6];
% grid=20;
% lambdas=0.01/2;
errors=zeros(size(grid,1),size(lambdas,2));
order=2;
for jj=1:size(lambdas,2)
    for j=1:size(grid,1)
        n=grid(j,jj);
        error=0;
        %k order of Runge-Kutta
        k=2;
        % v=zeros(2*(r*k*3*n-3)+1,r*k*3*n-2);
        ds=1/n;
        dq=1/n;
        %Doing k=1 doesn't work
        [t,~,~]=domofdep(ds,1,k,1,1,0,1);
        T=size(t,2);
        v=zeros(k*T+1);
        % muhat=((-(3*k*r*n-3):(3*k*r*n-3))/n+1)'./((0:(3*k*r*n-3))/n+2);
        muhat=((-(0):(k*T))/n+1)'./((0:(k*T))/n+2);
        %RHS f corresponding to v=(1-t)max(1/2,muhat^2)
        %     coin=1/2<(muhat.^2);
        %     f=coin.*(-muhat.^2+max(1/2,muhat))+(~coin).*(-1/2+max(1/2,muhat));
        f=0;
        %muhat=zeros(r*k*4*n-3,r*k*3*n-2);
        %Z=running total number of cells
        Z=0;
        for i=T:-1:(2)
            dt=t(i)-t(i-1);
            % if i==r*n*k
            %     big=max(muhat,1/2);
            %     pie=exp(1/lambdas(jj)*(muhat-big))./(exp(1/lambdas(jj)*(1/2-big))+exp(1/lambdas(jj)*(muhat-big)));
            %     v=dt*(pie.*(muhat-lambdas(jj)*log(pie))+(1-pie).*(1/2-lambdas(jj)*log(1-pie))-f);
            %     pause
            % else  
            temp=v;
            center=1+(0:(k*i-1));
            % centers=3*(r*k*n-1)+1+(-3*(i-1):3*(i-1));
            centers=1+(0:(k*i-1));
            [tempsp,tempqp]=ENO(temp,k*i-1,n,lambdas(jj),k);
            % if i<=4
            %     isnan(vsp)
            % end
            %g1=v(1:i,1:i)+dt*(mu(1:i,1:i)+mu(1:i,1:i).*(v(2:(i+1),2:(i+1))-v(1:i,2:(i+1)))/ds+(v(1:i,2:(i+1))-v(1:i,1:i))/dq);
            % g1=v+dt*(pie.*(muhat+...
            %     muhat.*vsp+vqp)+(1-pie)*1/2-lambda*(pie.*log(pie)+(1-pie).*log(1-pie))-f);
            % g1=v+dt*(pie.*(muhat+...
            %     muhat.*(vsm+vsp)/2+(vqp+vqm)/2+muhat.*(vsp-vsm)/2+(vqp-vqm)/2)+(1-pie)*1/2-lambda*(pie.*log(pie)+(1-pie).*log(1-pie))-f);

            % big=max(muhat(centers,center)+max(muhat(centers,center),0).*vsp+min(muhat(centers,center),0).*vsm+vqp,1/2);
            % small=min(muhat(centers,center)+max(muhat(centers,center),0).*vsp+min(muhat(centers,center),0).*vsm+vqp,1/2);
            big=max(muhat(centers,center)+muhat(centers,center).*tempsp+tempqp,1/2);
            small=min(muhat(centers,center)+muhat(centers,center).*tempsp+tempqp,1/2);
            temp(centers,center)=temp(centers,center)+dt*(lambdas(jj)*log(exp((small-big)/lambdas(jj))+1)+big);
            %if you want RK2 uncomment below
            center=1+(0:(k*i-2));
            % centers=3*(r*k*n-1)+1+(-3*(i-2):3*(i-2));
            centers=1+(0:(k*i-2));
            [tempsp,tempqp]=ENO(temp,k*i-2,n,lambdas(jj),k);
            % if i<=4
            %     isnan(g1sp)
            % end
            % g1=g1+dt*(pie.*(muhat+...
            %     muhat.*g1sp+g1qp)+(1-pie)*1/2-lambda*(pie.*log(pie)+(1-pie).*log(1-pie))-f);
            % g1=v+dt*(pie.*(muhat+...
            %     muhat.*(g1sm+g1sp)/2+(g1qp+g1qm)/2+muhat.*(g1sp-g1sm)/2+(g1qp-g1qm)/2)+(1-pie)*1/2-lambda*(pie.*log(pie)+(1-pie).*log(1-pie))-f);
            % big=max(muhat(centers,center)+max(muhat(centers,center),0).*tempsp+min(muhat(centers,center),0).*tempsm+tempqp,1/2);
            % small=min(muhat(centers,center)+max(muhat(centers,center),0).*tempsp+min(muhat(centers,center),0).*tempsm+tempqp,1/2);
            big=max(muhat(centers,center)+muhat(centers,center).*tempsp+tempqp,1/2);
            small=min(muhat(centers,center)+muhat(centers,center).*tempsp+tempqp,1/2);
            temp(centers,center)=temp(centers,center)+dt*(lambdas(jj)*log(exp((small-big)/lambdas(jj))+1)+big);
            v(centers,center)=1/2*(v(centers,center)+temp(centers,center));
            %if you want RK2 uncomment above
            indexs=centers;
            index=center;
            %compared to regularized
            big=max(muhat(indexs,index),1/2);
            pie=exp(1/lambdas(jj)*(muhat(indexs,index)-big))./(exp(1/lambdas(jj)*(muhat(indexs,index)-big))+exp(1/lambdas(jj)*(1/2-big)));
            limit1=isinf(log(pie));
            limit2=isinf(log(1-pie));
            reg1=log(pie);
            reg2=log(1-pie);
            reg1(limit1)=0;
            reg2(limit2)=0;
            sol=(pie.*(muhat(indexs,index)-lambdas(jj)*reg1)+...
                (1-pie).*(1/2-lambdas(jj)*reg2))*(1-t(i-1));
            sol=lambdas(jj)*log(exp(muhat(indexs,indexs)/lambdas(jj))+exp(1/(2*lambdas(jj))))*(1-t(i-1));   
            if i>=2
                % step=(pie.*(muhat(indexs,index)-lambdas(jj)*log(pie))+(1-pie).*(1/2-lambdas(jj)*log(1-pie)))*(1-t(i-1))-v(indexs,index);
                % error=(Z*error+norm(reshape(step,1,[]),1))/(Z+size(index,2)^2);
                error=(Z*error+norm(reshape(sol-v(indexs,index),1,[]),1))/(Z+size(index,2)^2);
                if isinf(error)
                    pause
                end
            end
            %uncomment if you want error surf plots
            % if i==floor(T/2)
            %     figure
            %     surf(sol-v(indexs,index))
            % end
            %compared to unregularized
            % error=(Z*error+norm(reshape(max(muhat(indexs,index),1/2)*(1-t(i-1))-...
            %     v(indexs,index),...
            %     1,[]),1))/(Z+size(index,2)^2);
            Z=Z+size(index,2)^2;
        end
        %12*differences/(n*(n+1)*(2*n+1))
        %error
        %     i=r*k*(n-2);
        %     errors(j)=norm(reshape(max(muhat((2*r*k*n-1):(2*r*k*n-1+min(i/k,n)),r*k*n:(r*k*n+min(i/k,n))),1/2)*...
        %         (1-(i-k)/(r*k*n))-v((2*r*k*n-1):(2*r*k*n-1+min(i/k,n)),r*k*n:(r*k*n+min(i/k,n))),...
        %         1,(min(i/k,n)+1)^2),1)/(min(i/k,n)+1)^2;
        errors(j,jj)=error;
        %graph
        % mesh((pie.*(muhat-lambda*log(pie))+(1-pie).*(1/2-lambda*log(1-pie)))*(1-(i-1)/(r*n)))
        % hold on
        % mesh(v)
        % clf
    end
end
%lambda=1
% grid=[20;40;60;80;100;120;140;160;180;200;220];
% errors=1.0e-05 *[
%     0.7555;
%     0.1302;
%     0.0461;
%     0.0233;
%     0.0166;
%     0.0110;
%     0.008492;
%     0.007530;
%     0.005437;
%     0.005003;
%     0.004046;];
% grid=[20;40;60;80;100;120;140;160;180;200;220]*[1,1,1];
% lambdas=[0.1,0.01];
% errors=1.0e-04 *[
% 
%     0.1851       NaN;
%     0.0740       NaN;
%     0.0350       NaN;
%     0.0267       NaN;
%     0.0184       NaN;
%     0.0154       NaN;
%     0.0126       NaN;
%     0.0120       NaN;
%     0.0099       NaN;
%     0.0090       NaN;
%     0.0090       NaN;];
%only order -1.2889
% grid=[20;40;60;80;100;120;140;160;180;200;220]*[1];
% lambdas=0.01;
% polyfit(log(grid),log(errors),1)
% 
% ans =
% 
%    -1.3851   -6.9964
% errors=
% 1.0e-04 *
% 
%     0.2166
%     0.0471
%     0.0230
%     0.0172
%     0.0138
%     0.0119
%     0.0100
%     0.0079
%     0.0069
%     0.0069
%     0.0067
%Unshifted in q:
% grid=[20;40;60;80;100;120;140;160;180;200;220]*[1,1,1];
% lambdas=[0.02,0.01,0.005];
% errors
% 
% errors =
% 1.0e-04 *
% 
%     0.1917    0.2344    0.2640
%     0.0455    0.0645    0.0827
%     0.0192    0.0291    0.0409
%     0.0106    0.0167    0.0250
%     0.0066    0.0106    0.0167
%     0.0046    0.0074    0.0121
%     0.0033    0.0054    0.0091
%     0.0025    0.0042    0.0071
%     0.0020    0.0033    0.0056
%     0.0016    0.0026    0.0046
%     0.0013    0.0022    0.0038
% polyfit(log(grid(:,1)),log(errors(:,1)),1)
% 
% ans =
% 
%    -2.0780   -4.6433
% polyfit(log(grid(:,2)),log(errors(:,2)),1)
% 
% ans =
% 
%    -1.9609   -4.7399   
% polyfit(log(grid(:,3)),log(errors(:,3)),1)
% 
% ans =
% 
%    -1.7678   -5.1918
% grid=[20;40;60;80;100;120;140;160;180;200;220]*[1];
% lambdas=[0.0025];
% polyfit(log(grid(:,1)),log(errors(:,1)),1)
% 
% ans =
% 
%    -1.5825   -5.7257
% 
% errors
% 
% errors =
% 
%    1.0e-04 *
% 
%     0.2800
%     0.0944
%     0.0500
%     0.0324
%     0.0227
%     0.0171
%     0.0132
%     0.0107
%     0.0087
%     0.0073
%     0.0062
% grid=[20;40;60;80;100;120;140;160;180;200;220]*[2];
% lambdas=[0.0025];
% polyfit(log(grid(:,1)),log(errors(:,1)),1)
% 
% ans =
% 
%    -1.6586   -5.3821
% 
% errors
% 
% errors =
% 
%    1.0e-05 *
% 
%     0.9445
%     0.3238
%     0.1708
%     0.1068
%     0.0730
%     0.0534
%     0.0407
%     0.0321
%     0.0260
%     0.0215
%     0.0180
% grid=[20;40;60;80;100;120;140;160;180;200;220]*[3];
% lambdas=[0.0025];
% polyfit(log(grid(:,1)),log(errors(:,1)),1)
% 
% ans =
% 
%    -1.6981   -5.1806
% errors
% 
% errors =
% 
%    1.0e-05 *
% 
%     0.5004
%     0.1708
%     0.0875
%     0.0534
%     0.0361
%     0.0260
%     0.0195
%     0.0153
%     0.0127
%     0.0106
%     0.0090
% ENO3 RK2
% grid=[20;40;60;80;100;120;140;160;180;200;220]*[1];
% lambdas=[0.0025];
%...
%ENO2 RK2
% grid=[20;40;60;80;100;120;140;160;180;200;220]*[2,4];
% lambdas=[0.01,0.005];
% errors
% 
% errors =
% 
%    1.0e-04 *
% 
%     0.0645    0.0250
%     0.0167    0.0071
%     0.0074    0.0032
%     0.0042    0.0038
%     0.0026    0.0299
%     0.0018    0.0858
%     0.0013    0.2373
%     0.0025    0.1676
%     0.0308    0.2830
%     0.0264    0.1982
%     0.2299    0.2618
% grid=[20;40;60;80;100;120;140;160;180;200;220]*[2];
% lambdas=[0.01];
% errors =1.0e-05 *[0.6446;
%     0.1668;
%     0.0742;
%     0.0415;
%     0.0263;
%     0.0182;
%     0.0134;
%     0.0103;
%     0.0082;
%     0.0067;
%     0.0055];
% grid=[20;40;60;80;100;120;140;160;180;200;220]*[4];
% lambdas=[0.005];
% errors = 1.0e-05 *[0.2503;
%     0.0708;
%     0.0324;
%     0.0185;
%     0.0120;
%     0.0084;
%     0.0064;
%     0.0053;
%     0.0046;
%     0.004517;
%     0.004329]
% grid=[75,150,300,600];
% lambdas=[0.02,0.01,0.005,0.0025];
% errors =1.0e-03 *[0.1246    0.0317    0.0081    0.0021];
% grid=[200;220]*4;
% lambdas=[0.005];
% errors =
% 
%    7.8258e-06
% grid=[200;220]*4;
% lambdas=[0.005];
% errors calculated using the log-exp form of the solution for 
% errors =
% 
%    1.0e-07 *
% 
%     0.4517
%     0.4329
function [vsp,vqp]=ENO(v,i,n,lambda,k)
ds=1/n;
dq=1/n;
index=1+(0:i);
order=2;
vsp=(v(index+1,index)-v(index,index))/ds;
vqp=(v(index,index+1)-v(index,index))/dq;
% centers=3*(r*k*n-1)+1+(-3*(i-1):3*(i-1));
if i>=2*order
    centers=1+(order:(i-order))';
    centerq=1+(order:(i-order));
    allq=1+(0:i);
    alls=(1+(0:i))';
    indexs(:,:,1)=repmat(centers,size(allq));
    indexs(:,:,2)=repmat(centers+1,size(allq));
    table(:,:,1)=v(centers,allq);
    table(:,:,2)=(v(centers+1,allq)-v(centers,allq))/ds;
    % vsp=zeros(size(alls,1),size(allq,2));
    for ii=2:order
        ip=max(indexs,[],3)+1;
        im=min(indexs,[],3)-1;
        left=sub2ind(size(v),reshape(im,[],1), reshape(ones(size(centers))*(allq),[],1));
        right=sub2ind(size(v),reshape(ip,[],1), reshape(ones(size(centers))*(allq),[],1));
        % tempm=reshape(v(left),size(centers,1),size(allq,2));
        % tempp=reshape(v(right),size(centers,1),size(allq,2));
        tempm=reshape(v(left),size(centers,1),[]);
        tempp=reshape(v(right),size(centers,1),[]);
        for jj=1:ii
            tempm=(tempm-table(:,:,jj))./(im-indexs(:,:,jj))/ds;
            tempp=(tempp-table(:,:,jj))./(ip-indexs(:,:,jj))/ds;
        end
        choice=abs(tempp)<abs(tempm);
        table(:,:,ii+1)=choice.*tempp+(~choice).*tempm;
        indexs(:,:,ii+1)=choice.*ip+(~choice).*im;
        vsp(centers,:)=vsp(centers,:)+table(:,:,ii+1).*prod((repmat(indexs(:,:,1),1,1,ii-1)-indexs(:,:,2:ii))*ds,3);
        % if ii==order
        %     discont1=abs(tempp)>10;
        %     discont2=abs(tempm)>10;
        %     temp=vsp(centers,:);
        %     temp(discont1)=0;
        %     temp(discont2)=0;
        %     vsp(centers,:)=temp;
        % end
    end
    indexq(:,:,1)=repmat(centerq,size(alls));
    indexq(:,:,2)=repmat(centerq+1,size(alls));
    clear table
    table(:,:,1)=v(alls,centerq);
    table(:,:,2)=(v(alls,1+centerq)-v(alls,centerq))/dq;
    for ii=2:order
        ip=max(indexq,[],3)+1;
        im=min(indexq,[],3)-1;
        left=sub2ind(size(v),reshape(alls*ones(size(centerq)),[],1),reshape(im,[],1));
        right=sub2ind(size(v),reshape(alls*ones(size(centerq)),[],1),reshape(ip,[],1));
        % tempm=reshape(v(left),size(centers,1),size(allq,2));
        % tempp=reshape(v(right),size(centers,1),size(allq,2));
        tempm=reshape(v(left),[],size(centerq,2));
        tempp=reshape(v(right),[],size(centerq,2));
        for jj=1:ii
            tempm=(tempm-table(:,:,jj))./(im-indexq(:,:,jj))/ds;
            tempp=(tempp-table(:,:,jj))./(ip-indexq(:,:,jj))/ds;
        end
        choice=abs(tempp)<abs(tempm);
        table(:,:,ii+1)=choice.*tempp+(~choice).*tempm;
        indexq(:,:,ii+1)=choice.*ip+(~choice).*im;
        vqp(:,centerq)=vqp(:,centerq)+table(:,:,ii+1).*prod((repmat(indexq(:,:,1),1,1,ii-1)-indexq(:,:,2:ii))*dq,3);
        % if ii==order
        %     discont1=abs(tempp)>10;
        %     discont2=abs(tempm)>10;
        %     temp=vqp(centers,:);
        %     temp(discont1)=0;
        %     temp(discont2)=0;
        %     vqp(centers,:)=temp;
        % end
    end
end
% if i>=2
%     vsp(2,:)=(v(3,1+allq)-v(2,1+allq))/ds;
% end
% vsp=(v(alls+1,allq)-v(alls,allq))/ds;
%%almost ENO3 below
% vsright=(v(alls+1,allq)-v(alls,allq))/ds+(v(alls+2,allq)-2*v(alls+1,allq)+v(alls,allq))/(2*ds^2)*(-ds)+...
%     (v(alls+3,allq)-3*v(alls+2,allq)+3*v(alls+1,allq)-v(alls,allq))/(6*ds^3)*(-ds)*(-2*ds);
% vsleft=(v(alls,allq)-v(max(alls-1,1),allq))/ds+(v(alls+1,allq)-2*v(alls,allq)+v(max(alls-1,1),allq))/(2*ds^2)*ds+...
%     (v(alls+2,allq)-3*v(alls+1,allq)+3*v(alls,allq)-v(max(alls-1,1),allq))/(6*ds^3)*ds*(-ds);
% right=v(alls+3,allq)-3*v(alls+2,allq)+3*v(alls+1,allq)-v(alls,allq);
% left=v(alls+2,allq)-3*v(alls+1,allq)+3*v(alls,allq)-v(max(alls-1,1),allq);
% choice=abs(right)<abs(left);
% vsp=choice.*vsright+(~choice).*vsleft;
% vsp(1,allq)=(v(2,allq)-v(1,allq))/ds;
% vqup=(v(alls,allq+1)-v(alls,allq))/dq+(v(alls,allq+2)-2*v(alls,allq+1)+v(alls,allq))/(2*dq^2)*(-dq)+...
%     (v(alls,allq+3)-3*v(alls,allq+2)+3*v(alls,allq+1)-v(alls,allq))/(6*dq^3)*(-dq)*(-2*dq);
% vqdown=(v(alls,allq)-v(alls,max(allq-1,1)))/dq+(v(alls,allq+1)-2*v(alls,allq)+v(alls,max(allq-1,1)))/(2*dq^2)*ds+...
%     (v(alls,allq+2)-3*v(alls,allq+1)+3*v(alls,allq)-v(alls,max(allq-1,1)))/(6*dq^3)*dq*(-dq);
% up=v(alls,allq+3)-3*v(alls,allq+2)+3*v(alls,allq+1)-v(alls,allq);
% down=v(alls,allq+2)-3*v(alls,allq+1)+3*v(alls,allq)-v(alls,max(allq-1,1));
% choice=abs(up)<abs(down);
% vqp=choice.*vqup+(~choice).*vqdown;
% vqp(alls,1)=(v(alls,2)-v(alls,1))/dq;
%%almost ENO3 above
% down=[1,1+(0:(3*i-4))];
% uq=v(centers,centerq+3)-3*v(centers,centerq+2)+3*v(centers,centerq+1)-v(centers,centerq);
% downq=v(centers,centerq+2)-3*v(centers,centerq+1)+3*v(centers,centerq)-v(centers,down);
% choice=(abs(uq)<abs(downq));
% vqdown=-downq/(6*dq)+(v(centers,centerq+1)-2*v(centers,centerq)+v(centers,down))/(2*dq)+(v(centers,centerq)-v(centers,down))/dq;
% vqup=uq/(3*dq)-(v(centers,centerq+2)-2*v(centers,centerq+1)+v(centers,centerq))/(2*dq)+(v(centers,centerq+1)-v(centers,centerq))/dq;
% vqp=choice.*vqup+(~choice).*vqdown;
% vqp(:,1)=(v(centers,2)-v(centers,1))/dq;
end