%these values must be given!
states=
Dmin=
Dmax=
k=
steps=
%k1=
%k2=
%k1 and k2 are for the model with bias to lower diffusion states

Dvec=linspace(Dmin,Dmax,states);
dt=

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
% kvecbig(1)=-k2;
% kvecbig(2)=k2;
% kaltbig=zeros(1,states);
% kaltbig(end)=-k2;
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
% kvecsmall(1)=-2*k1;
% kvecsmall(2)=k1;
% kvecsmall(end)=k1;
% kaltsmall=zeros(1,states);
% kaltsmall(1)=-k1;
% kaltsmall(2)=k1;
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


