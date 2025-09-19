%==========================================================================
%% housekeeping
%==========================================================================
clear variables;close all;userpath('clear');restoredefaultpath;clc;
tic;
rng('default'); % reinitialize the random number generator to its startup configuration
seed =0;
rng(seed,'twister');         % set seed
currdir=pwd;
addpath([currdir,'/data']); % set path to data functions
addpath([currdir,'/helpfunctions']); % set path to helper functions
addpath([currdir,'/helpfunctions/ChrisSimsOptimize']);
addpath([currdir,'/helpfunctions/subroutines']); % set path to helper functions
cd(currdir)


%==========================================================================
%% load the data
%==========================================================================
load('kmData.mat')
% percent change in global oil production, real activity index from Kilian(AER 2009), the log real price of oil, and changes in OECD crude oil inventories

%tstart = find(all((datesq.FEDFUNDS=='1986-03-31'),2)==1); % start estimation with 1986:Q1
%tend   = size(yall,1) - 1;                                % end estimation with 2008:Q3
%YY     = yall(tstart:tend,:);
num    = kmData;

%=========================================================================
%% model setup
%==========================================================================
nlag      = 24;                     % number of lags
nvar      = size(num,2);           % number of endogenous variables
nex       = 12;                     % 12 because of seasonal dummies
m         = nvar*nlag + nex;       % number of exogenous variables
horizon   = 20;                    % maximum horizon for IRFs
horizons  = [0,1,2,3,4,inf];       % horizons upon which sign and zero restrictions can be imposed
NS        = 1 + numel(horizons);   % number of objects in F(THETA) to which we impose sign and zero restrictions: F(THETA)=[A_0;L_{0};L_{1};L_{2};L_{3};L_{4};L_{inf}]
e         = eye(nvar);             % create identity matrix
maxBSQdraws2store  =1e2;
conjugate  = 'none';
iter_show  = 1e4;
fixed_rf   = 0;
prior_only = 0;
label_R = 'kilianmurphy';
prior_type = 'flat';
%==========================================================================
%% Setup info/settings and IDENTIFYING RESTRICTIONS
%==========================================================================
info=SetupInfo(nvar,m, nlag,horizons,@(x)chol(x));


%==========================================================================
%% Mappings
%==========================================================================
fo                 = @(x)f_h(x,info);
fo_inv             = @(x)f_h_inv(x,info);
fo_str2irfs        = @(x)StructuralToIRF(x,info);
fo_str2irfs_inv    = @(x)IRFToStructural(x,info);

%==========================================================================
%% write data in Rubio, Waggoner, and Zha (RES 2010)'s notation
%==========================================================================
% yt(t) A0 = xt(t) Aplus + constant + et(t) for t=1...,T;
% yt(t)    = xt(t) B     + ut(t)            for t=1...,T;
% x(t)     = [constant,yt(t-1), ... ,yt(t-nlag)];
% matrix notation yt = xt*B + ut;
% xt=[yt_{-1} ones(T,1)];
y0bar = mean(num(1:nlag,:),1);% useful for GLP prior
info.y0bar = y0bar;

Tnum = size(num,1);
yt   = num(nlag+1:end,:);
T    = size(yt,1);

%creating SA dummies
x=[eye(11); zeros(1,11)];
X2=[];
for i=1:fix((Tnum-nlag)/12)  %number of years
    X2=[X2;x];
end
[l w] =size(X2);
last=[eye((Tnum-nlag)-l), zeros((Tnum-nlag)-l,11-((Tnum-nlag)-l))];
X2=[X2;last];
X2=[ones(Tnum-nlag,1), X2];


xt = zeros(T,nvar*nlag+nex);
for i=1:nlag
    xt(:,nex+nvar*(i-1)+1:nex+nvar*i) = num((nlag-(i-1)):end-i,:) ;
end



if nex==1
    xt(:,1)=ones(T,1);
end
if nex==12
    xt(:,1:12)=X2;
end


% write data in Zellner (1971, pp 224-227) notation
Y = yt; % T by nvar matrix of observations
X = xt; % T by (nvar*nlag+1) matrix of regressors

Ynd = Y;
Xnd = X;

%% prior for reduced-form parameters


% prior
nnuBar             = 0;
OomegaBarInverse   = zeros(m);
mmuBar             = zeros(m,nvar);  % Psi in the paper
PpsiBar            = zeros(nvar);    % Phi in the paper



%% posterior for reduced-form parameters
nnuTilde            = T +nnuBar;
OomegaTilde         = (X'*X  + OomegaBarInverse)\eye(m);
OomegaTildeInverse  =  X'*X  + OomegaBarInverse;
mmuTilde            = (X'*X  + OomegaBarInverse)\(X'*Y + OomegaBarInverse*mmuBar);
PpsiTilde           = Y'*Y + PpsiBar + mmuBar'*OomegaBarInverse*mmuBar - mmuTilde'*OomegaTildeInverse*mmuTilde;
PpsiTilde           = (PpsiTilde'+PpsiTilde)*0.5;

%% Parameters to fix reduced-form
BHat            = (X'*X)\(X'*Y);
SigmaHat           = (Y-X*BHat)'*(Y-X*BHat)/size(Y,1);
SigmaHat           = (SigmaHat'+SigmaHat)*0.5;


%% useful definitions
% definitios used to store orthogonal-reduced-form draws, volume elements, and unnormalized weights
Bdraws         = cell([maxBSQdraws2store,1]); % reduced-form lag parameters
Sigmadraws     = cell([maxBSQdraws2store,1]); % reduced-form covariance matrices
Qdraws         = cell([maxBSQdraws2store,1]); % orthogonal matrices


% definitions related to IRFs and stability of the coefficients; based on page 12 of Rubio, Waggoner, and Zha (RES 2010)
J      = [e;repmat(zeros(nvar),nlag-1,1)];
A      = cell(nlag,1);
extraF = repmat(zeros(nvar),1,nlag-1);
bigF      = zeros(nlag*nvar,nlag*nvar);
for l=1:nlag-1
    bigF((l-1)*nvar+1:l*nvar,nvar+1:nlag*nvar)=[repmat(zeros(nvar),1,l-1) e repmat(zeros(nvar),1,nlag-(l+1))];
end

% definition to facilitate the draws from B|Sigma
hh              = info.h;
cholOomegaTilde = hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below

if prior_only==1
    cholOomegaBar   = hh(OomegaBar)'; % this matrix is used to draw B|Sigma below
end

%% initialize counters to track the state of the computations
counter = 1;
count   = 0;
tStart = tic;


nshocks=3;

Ss{1,1}= [-1 0 0 0;
    0 -1 0 0;
    0  0 1 0;];


Ss{2,1}=   [1 0 0 0;
    0 1 0 0;
    0  0 1 0;];

Ss{3,1}=   [1 0 0 0;
    0 -1 0 0;
    0  0 1 0;
    0  0 0 1];

while count<maxBSQdraws2store

    if fixed_rf==1


        Sigmadraw     = SigmaHat;
        Bdraw         = BHat ;





    elseif prior_only==1
        %% Draw Sigma
        Sigmadraw     = iwishrnd(PpsiBar,nnuBar);
        cholSigmadraw = hh(Sigmadraw)';
        Bdraw         = kron(cholSigmadraw,cholOomegaBar)*randn(m*nvar,1) + reshape(mmuBar,nvar*m,1);
        Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);
    else

        %% Draw Sigma
        Sigmadraw     = iwishrnd(PpsiTilde,nnuTilde);
        cholSigmadraw = hh(Sigmadraw)';
        Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*nvar,1) + reshape(mmuTilde,nvar*m,1);
        Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);

    end



    %% Draw Q
    Qdraw        = DrawQ(info.nvar);



    %% check if sign restrictions hold

    % build T (m times n), m=structural shocks=3 in this application, and
    % n=variables=4 in this application



    Rr=chol(Sigmadraw)'*Qdraw;


    Tt=zeros(nshocks,nvar);
    for j=1:nshocks

        switch j



            case 1

                for i=1:nvar
                    if all(Ss{j,1}*Rr(:,i)>0)

                        Tt(j,i)=1;
                    elseif all(Ss{j,1}*Rr(:,i)<0)
                        Tt(j,i)=-1;
                    end
                end


            case 2

                for i=1:nvar
                    if all(Ss{j,1}*Rr(:,i)>0) && ((Rr(1,i)/Rr(3,i))<0.0258)
                        Tt(j,i)=1;
                    elseif all(Ss{j,1}*Rr(:,i)<0) && ((Rr(1,i)/Rr(3,i))<0.0258)
                        Tt(j,i)=-1;
                    end
                end

            case 3

                for i=1:nvar
                    if all(Ss{j,1}*Rr(:,i)>0) && ((Rr(1,i)/Rr(3,i))<0.0258)
                        Tt(j,i)=1;
                    elseif all(Ss{j,1}*Rr(:,i)<0) && ((Rr(1,i)/Rr(3,i))<0.0258)
                        Tt(j,i)=-1;
                    end

                end

        end
    end

    checkA1 = any(sum(Tt ~= 0, 1) > 1);

    if checkA1==1
        disp('PROBLEM with Tt')
        return;
    end

    nshocks=3;
    %if rank(Tt)<nshocks
    
    if any(sum(Tt == 0, 2)==size(Tt,2))

    else


        Rstar = zeros(nvar,nvar);
        S = cell(nshocks,1);

        isamples = zeros(nshocks,1);

        % For each row j
        for j = 1:nshocks
      
            S{j} = find(Tt(j,:) == 1 | Tt(j,:) == -1);



            % Get a random index from S_j
            if numel(S{j})>1
                randomIndex = randsample(S{j},1);
            else
                randomIndex = S{j};
            end

            isamples(j,1) = randomIndex;



            if Tt(j,isamples(j))==1
                Rstar(:,j) = Rr(:,isamples(j,1));
            elseif Tt(j,isamples(j))==-1
                Rstar(:,j) = -Rr(:,isamples(j,1));
            end

        end

        % Loop through j = m+1, ..., n
        i_values=zeros(nvar-nshocks,1);
        for j = nshocks+1:nvar
            % Create S_j = {1,...,n} \ {i_1,...,i_j}
            % We need to exclude i_1 through i_(j-1) from the set {1,...,n}
           
            S_j = setdiff(1:nvar, isamples(1:j-1,1));

            % Sample an element uniformly from S_j
            if numel(S_j)>1


                random_index = randsample(S_j,1);

            else
                random_index = S_j;
            end
            sampled_element = random_index;

            % Store the sampled element as i_j
            i_values(j) = sampled_element;

            if rand(1)>0.5
                Rstar(:,j) = Rr(:,i_values(j));
            else
                Rstar(:,j) = -Rr(:,i_values(j));
            end

        end


        count=count+1;

        % store orthogonal reduced-form draws
        Bdraws{count,1}     = Bdraw;
        Sigmadraws{count,1} = Sigmadraw;
        Qdraws{count,1} = (chol(Sigmadraw)')\Rstar;





    end






    if rem(counter,iter_show)==0

        display(['Number of draws = ',num2str(count)])
        display(['Remaining draws = ',num2str(maxBSQdraws2store-(count))])


    end

    if rem(counter,10000)==0

        display(['draws = ',num2str(counter)])


    end

    counter = counter + 1;



end


tElapsed = toc(tStart);


ne = count;


% alternative method for importance weights: imp_w and imp_w_check must be identical
%imp_w_check  = uwcheck/sum(uwcheck);
%max(abs(imp_w-imp_w_check))


%% store draws
L    = zeros(horizon+1,nvar,nvar,count);
cumL    = zeros(horizon+1,nvar,nvar,count);


for s=1:count

 

        Bdraw =     Bdraws{s,1} ;
        Sigmadraw = Sigmadraws{s,1} ;
        Qdraw=Qdraws{s,1};
 
        BSigmaQ = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];



    structpara = f_h_inv(BSigmaQ,info);


    LIRF = IRF_horizons(structpara, nvar, nlag, m, nex, 0:horizon);


 


    for h=0:horizon
        L(h+1,:,:,s) =  LIRF(1+h*nvar:(h+1)*nvar,:);
      
        for i=1:nvar
            cumL(h+1,1:4,i,s)   = sum(L(1:h+1,1:4,i,s),1);
        end

    end



end
telaspsed = toc(tStart); 


stack = cat(3, Sigmadraws{:});        % 10×3×3 array
SigmaMean = mean(stack, 3);  

stack = cat(3, Bdraws{:}); 
BMean = mean(stack, 3);  
cd results

   if ~exist('matfiles', 'dir')
       mkdir('matfiles')
   end


    savefile= ['matfiles/rTalgo_',label_R,'prior_only_',num2str(prior_only),'prior_',prior_type,'fixed_rf_',num2str(fixed_rf),'.mat'];

save(savefile,'L','cumL','horizon','');
cd ..
