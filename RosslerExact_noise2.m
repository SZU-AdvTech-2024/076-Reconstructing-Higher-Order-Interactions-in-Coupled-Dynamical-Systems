%%% reseach paper: Reconstruction of complex networks with unknown topology and dynamics
%%% By Federico Battiston, Vito Latora, Massimiliano Nicosia, Giovanni Russo, and Vincenzo Nicosia
%%% Nature Communications -  NCOMMS-13-10612A

%%% As noise strength = 0.05, the reconstruction error is calculated 
%%% for different values of M/H by different methods

clc; clear; close all;

load('ZackaryNet.mat');

%remove the following three lines, if not using parallel computing 
%if isempty(gcp('nocreate'))
%    parpool('local',10);
%end

for i45=1:length(EdgeList)
    A(EdgeList(i45,1),EdgeList(i45,2))=1;
    A(EdgeList(i45,2),EdgeList(i45,1))=1;
end

TriangleList=closedtriangles;

% parameter of the strength of the Gaussian noise
strength = 0.05;

% parameters of the coupling
k=1e-4;
kD=1e-5;

% initial conditions for the oscillators
xoold =(30*rand(N,1)-15)/5;
yoold =(30*rand(N,1)-15)/5;
zoold =(40*rand(N,1)-5)/5;
x0=[xoold; yoold; zoold];

% parameters of integration
dt=0.01;
tmax=100;
T=0:dt:tmax;
options=odeset('abstol',1e-12,'reltol',1e-12);

%range of M/H
imin = 100;
imax = 800;
istep=5;
irange=round(logspace(1,2.9,9));

% initialization of the errors
nvaloriM = length(irange)
err0 = zeros(nvaloriM,1);
err00 = zeros(nvaloriM,1);
err2 = zeros(nvaloriM,1);
err3 = zeros(nvaloriM,1);


% number of unknowns for each node
H=(N-1)+(N-1)*(N-2)/2;

%% 
%change parfor into for, if not using parallel computing
for iMindice = 1:length(irange)
    iM=irange(iMindice);
    disp(iM)
    M = H*iM;
    dt = tmax/M;
    T = 0:dt:tmax;

    % Accurate integration of the equation
    options = odeset('abstol',1e-12,'reltol',1e-12);
    [T,X]=ode45(@(t,x) roessler_hoi(t,x,EdgeList,TriangleList),T,x0,options);
    [raw, col] = size(X);
    for i=1:N
        temp_mean = mean(X(:,i*3-2));
        X(:,i*3-2) = X(:,i*3-2) + temp_mean*strength*(2*rand(raw,1)-1);  %add noise
    end

    nt = length(T);

    I1 = 2:nt-1; % Internal steps
    I2 = 3:nt-2; % Internal internal points

    F = zeros(nt-2,3*N);

    for n=0:nt-1
        pop = X(n+1,:);
        F(n+1,:) = roessler_hoi(n*dt,pop',EdgeList,TriangleList);
        temp_mean = mean(F(n+1,:));
        for i = 1:N
            F(n+1,3*i-2) = F(n+1,3*i-2) + temp_mean*strength*(2*rand-1);
        end
    end
    Ftrue=F;
    
    A00=zeros(H,N);  % reconstructed matrix with OLS
    A2=zeros(H,N);  % reconstructed matrix with NNLS
    A3=zeros(H,N);  % reconstructed matrix with signal lasso
    AA=zeros(H,N);  % "true" matrix
    EdgeList0=[1 1; 2 2];
    TriangleList0=[1 1 1; 2 2 2];
    for n=0:nt-1
        pop = X(n+1,:);
        F(n+1,:) = roessler_hoi(n*dt,pop',EdgeList0,TriangleList0);
        temp_mean = mean(F(n+1,:));
        for i = 1:N
            F(n+1,3*i-2) = F(n+1,3*i-2) + temp_mean*strength*(2*rand-1);
        end
    end

    % Zero order method
    Phi = zeros(nt-2,H);
    
    for i=1:N
        % pairwise interactions j=1:N-1
        Phi(:,1:N-1) = k*(X(I1,[1:i-1 i+1:N])-X(I1,i));
        AA(1:N-1,i)=A(i,[1:i-1 i+1:N]);
        % h.o.i. terms j=N:H
        vtemp=zeros(nt-2,(N-1)*(N-2)/2);
        itemp=1;
        for ii1=[1:i-1 i+1:N]
            for jj1=[1:i-1 i+1:N]
                if jj1>ii1
                    vtemp(:,itemp)=kD*(X(I1,ii1).*X(I1,jj1).^2+X(I1,jj1).*X(I1,ii1).^2-2*X(I1,i).^3);
                    RowIdx = find(ismember(TriangleList, [i ii1 jj1],'rows'));
                    if (~isempty(RowIdx))
                        AA(N-1+itemp,i)=1;
                    end
                    RowIdx = find(ismember(TriangleList, [ii1 i jj1],'rows'));
                    if (~isempty(RowIdx))
                        AA(N-1+itemp,i)=1;
                    end
                    RowIdx = find(ismember(TriangleList, [ii1 jj1 i],'rows'));
                    if (~isempty(RowIdx))
                        AA(N-1+itemp,i)=1;
                    end
                    itemp=itemp+1;
                end
            end
        end
        Phi(:,N:N-1+(N-1)*(N-2)/2)=vtemp;
        
        Yi = Ftrue(I1,i)-F(I1,i);
        Zi0 = lsqminnorm(Phi,Yi,1e-12);  % OLS
        Zi2 = lsqnonneg(Phi,Yi);  % NNLS
        %Zi2 = NNLS_cvx(Yi, Phi, H);
        
        %%%%%signal_lasso
        %w0 = zeros(H,1);
        alpha1=0.001;
        alpha2=0.001;
        max_iters_ = 50;%50000; %50 for a quick check, 50000 for accurate results
        intercept_=0;
        %Zi3 = signal_lasso_3(Phi,Yi,alpha1,alpha2,max_iters_);
        Zi3 = signal_lasso_2(Yi,Phi,H,alpha1,alpha2);
        
        A00(:,i) = Zi0;
        A2(:,i) = Zi2;
        A3(:,i) = Zi3;
    end
    normA = norm(AA,"fro");
    err00(iMindice,1)=norm(AA-A00,"fro")/normA;
    err2(iMindice,1)=norm(AA-A2,"fro")/normA;
    err3(iMindice,1)=norm(AA-A3,"fro")/normA;
    
end

%% draw pictures
%irange = imin:istep:imax;

figure,semilogy(irange,err00,irange,err2,irange,err3);
xlabel('M/H')
ylabel('E')
legend('OLS', 'NNLS','SL')

figure,plot(irange,err00,irange,err2,irange,err3);
xlabel('M/H')
ylabel('E')
legend('OLS', 'NNLS','SL')