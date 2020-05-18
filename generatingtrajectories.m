%these values must be given!
states=10;
Dmin=5;
Dmax=50;
k=0.2;
steps=500;

Dvec=linspace(Dmin,Dmax,states);
dt=0.1;

%% Transitionmatrices
%The first matrix below is for a state, which can be visited by its two
%nearest neighbours. There is no bias towards Dmin. 
%generate matrix with rates by shifting the row vector kvec
kvec=zeros(1,states);
kvec(1)=-2*k;
kvec(2)=k;
kvec(end)=k;
kalt=zeros(1,states);
kalt(1)=-k;
kalt(2)=k;

K=zeros(states,states);
for i=1:states
    if i==1
       K(i,:)=kalt;
    elseif i==states
       kalt2=-1*circshift(kalt,states-2);
       K(i,:)=kalt2;
    else 
    K(i,:)=circshift(kvec,i-1);
    end
end

% This second matrix is for a system in which the particle tends to states with lower diffusion constants. 
% kvecbig=zeros(1,states);
% kvecbig(1)=-k;
% kvecbig(2)=k;
% kaltbig=zeros(1,states);
% kaltbig(end)=-k;
% Kbig=zeros(states,states);
% for i=1:states
%     if i==1
%         Kbig(i,:)=-1*circshift(kaltbig,2);
%     elseif i==states
%         Kbig(i,:)=kaltbig;
%     else
%     Kbig(i,:)=circshift(kvecbig,i-1);
%     end
% end
% 
% % Symmetric matrix
% ksmall=k/3; %this rate determines the slight upwards diffusion in the state space. 
% kvecsmall=zeros(1,states);
% kvecsmall(1)=-2*ksmall;
% kvecsmall(2)=ksmall;
% kvecsmall(end)=ksmall;
% kaltsmall=zeros(1,states);
% kaltsmall(1)=-ksmall;
% kaltsmall(2)=ksmall;
% 
% Ksmall=zeros(states,states);
% for i=1:states
%     if i==1
%        Ksmall(i,:)=kaltsmall;
%     elseif i==states
%        kaltsmall2=-1*circshift(kaltsmall,states-2);
%        Ksmall(i,:)=kaltsmall2;
%     else 
%     Ksmall(i,:)=circshift(kvecsmall,i-1);
%     end
% end
% 
% K=Ksmall+Kbig;
%%
TRANS=abs(expm(K*dt)); %transition probability matrix

%inverse cumulative probability function
invs=@(z,D) 2*sqrt(D*dt)*erfinv(2*z-1); 

X=[zeros(steps,1), zeros(steps,1)];
dat=zeros(2,steps);
l=zeros(steps,1);
m=rand(steps,1);
u=rand(steps,2); 
for i=1:steps
    %first step
    if i==1
       l(i)=randi(states);
       X(i,1)=invs(u(i,1),Dvec(l(i)));
       X(i,2)=invs(u(i,2),Dvec(l(i)));
    else
        prob=m(i);
        s=l(i-1);
        A=cumsum(TRANS(:,s));
        l(i)=find(A>=prob,1,'first');
        X(i,1)=invs(u(i,1),Dvec(l(i)));
        X(i,2)=invs(u(i,2),Dvec(l(i)));
    end
end

traj=[X(:,1),X(:,2)];
save('randomdata.txt','traj','-ascii')

%uncomment the lines below for a plot of the trajectory
%Z=zeros(steps-1,2);
%for i=1:steps-1
%    if i==1
%    Z(i,:)=traj(i,:);
%    else
%    Z(i,:)=Z(i-1,:)+traj(i,:);
%    end
%end
%plot(Z(:,1),Z(:,2),'-')
%xlabel('x');
%ylabel('y');
