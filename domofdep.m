function [t, s, q] = domofdep(dx, eorder, rkstage, alpha, beta, starttime, endtime, C)

   if nargin == 7
      if eorder == 1
         C = 1;
      else
         C = 0.5;
      end
   end

   if length(dx) == 1
      dx(2) = dx(1);
   end
   
   N = 1;
   t(N) = 0;
   while t(N) < abs(endtime-starttime)
      smax = (N+1)*rkstage*eorder*min(dx);
      qmin = 0;
      mumax = (smax+alpha)/(qmin+alpha+beta);

%      dt(N) = min(dx)/4;
      % dt(N) = min(C*min(dx)/(mumax+1),abs(endtime-starttime)-t(N));
      % dt(N) = min(C*min(dx)/(mumax+1),abs(endtime-starttime)-t(N));
      dt(N) = min([C*min(dx)/(mumax+1),min(dx.^2)/2,abs(endtime-starttime)-t(N)]);
%      dt(N) = min(C*min(dx)/(mumax+1),abs(endtime-starttime)-t(N));
%      if eorder == 1 
%         dt(N) = min(C*min(dx)/(mumax+1),abs(endtime-starttime)-t(N));
%      else
%         dt(N) = min(C*min(dx)/sqrt(mumax^2+1),abs(endtime-starttime)-t(N));
%         dt(N) = min(dx)/4;
%      end
      t(N+1) = t(N)+dt(N);

      N = N+1;
   end
   %don't want output to clutter the command window
   % 
   % subplot(221)
   % plot(t)
   % subplot(222)
   % plot(dt)

   % t = flip(t);
   s = [0:(N-1)*rkstage*eorder]*dx(1);
   q = [0:(N-1)*rkstage*eorder]*dx(2);
   % 
   % [min(t) max(t)]
   % subplot(223)
   % plot(t)
   % subplot(224)
   % plot(abs(t(2:end)-t(1:end-1)))
   % 
   % pause(0)
end
