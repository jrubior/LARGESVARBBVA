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
if fixed_rf==1
    draw_Sigma=0;
    draw_B=0;
else
    draw_Sigma=1;
    draw_B=1;
end

elliptical_Q = 1; %0 when using Chan, Matthes, Yu, 1 when using elliptical
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
gs_qr    = @(x)qr_unique(x);
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



switch prior_type
    case 'flat'


    %% prior
    nnuBar             = 0;
    OomegaBarInverse   = zeros(m);
    mmuBar             = zeros(m,nvar);  % Psi in the paper
    PpsiBar            = zeros(nvar);    % Phi in the paper

    case 'minnesota'

        % Giannone, Lenza, Primiceri (2015) hyperpriors modes
        mode.llambda = .2;

        % Giannone, Lenza, Primiceri (2015) hyperpriors std deviations
        sd.llambda = .4;

        % Giannone, Lenza, Primiceri (2015) scale and shape of the IG on psi/(d-n-1)
        scalePSI   = 0.02^2;    


         % transform mode and std into Gamma(k,theta) coefficients
        info.priorcoef.llambda= GammaCoef(mode.llambda,sd.llambda,0);
        info.priorcoef.alpha.PSI=scalePSI;
        info.priorcoef.beta.PSI=scalePSI;
        info.pos = [1 4]; % stationary variables


        % bounds for maximization
        MIN.llambda = 0.0001;
        MAX.llambda = 5;



        % initial guess for hyper-parameters
        llambda  = 0.5;
        info.Vc  = 10^7;

        % initial guess for psi (phi in the paper)
        % residual variance of AR(1) for each variable
        SS=zeros(nvar,1);
        for i=1:nvar
            yols=Y(2:end,i);
            xols=[ones(T-1,1),Y(1:end-1,i)];
            bbetaols = (xols'*xols)\(xols'*yols);
            eols = yols-xols*bbetaols;
            SS(i)=eols'*eols/(size(yols,1)-size(xols,2));
        end
        MIN.ppsi  = SS./100;
        MAX.ppsi  = SS.*100;
        ppsi      = SS;
        inllambda = -log((MAX.llambda-llambda)./(llambda-MIN.llambda));
        inppsi    = -log((MAX.ppsi-ppsi)./(ppsi-MIN.ppsi));

        x0 = [inllambda;inppsi];


        f0 = negative_loglike_obj_mp(x0,info,Xnd,Ynd,MIN,MAX);

        crit                              = 1e-16;
        nit                               = 1000;
        H0                                = eye(size(x0,1))*10;

        [fh,xh,~,~,itct,~,~]              = csminwel('negative_loglike_obj_mp',x0,H0,[],crit,nit,info,Xnd,Ynd,MIN,MAX);

        llambda0_xh = MIN.llambda+(MAX.llambda-MIN.llambda)./(1+exp(-xh(1,1)));
        ppsi0_xh    = MIN.ppsi+(MAX.ppsi-MIN.ppsi)./(1+exp(-xh(2:info.nvar+1,1)));%
  
        llambda0 =  llambda0_xh;
        ppsi0    =  ppsi0_xh;


        if itct>1 % we re-set the random seed after using csminwel to guarantee the
            % posterior draws are identical across operating systems.

           

            rng(seed,'twister');
        end


        % set position where variables are in first differences

        % inverse Wishart prior parameters
        nnuBar              = info.nvar + 2;
        PpsiBar             = diag(ppsi0);
        oomegaBar=zeros(info.m,1);
        oomegaBar(1:info.nex,1)=info.Vc;
        for ell=1:info.nlag
            oomegaBar(info.nex+1+(ell-1)*info.nvar:info.nex+ell*info.nvar,1) =  (llambda0^2)*(nnuBar-info.nvar-1)./((ell^2)*ppsi0);
        end
        OomegaBar=diag(oomegaBar);
        OomegaBarInverse=diag(1./oomegaBar);
        mmuBar = zeros(info.m,info.nvar);
        diagmmuBar=ones(info.nvar,1);
        diagmmuBar(info.pos)=0;   % Set to zero the prior mean on the first own lag for variables selected in the vector pos
    
        mmuBar(info.nex+1:info.nvar+info.nex,:)=diag(diagmmuBar);


        Y = Ynd;
        X = Xnd;
        T = size(Y,1);


end





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
info.extraF = extraF;
info.J = J;
info.A = A;
info.bigF = bigF;


% definition to facilitate the draws from B|Sigma
hh              = info.h;
cholOomegaTilde = hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below

if prior_only==1
    cholOomegaBar   = hh(OomegaBar)'; % this matrix is used to draw B|Sigma below
end


Rmean = nan(nvar,nnuTilde);
for i=1:nnuTilde
    Rmean(:,i) = zeros(nvar,1);
end
Rvariance = zeros(nvar*nnuTilde,nvar*nnuTilde);
for i=1:nnuTilde
    Rvariance((i-1)*nvar+1:i*nvar,(i-1)*nvar+1:i*nvar) = inv(PpsiTilde);
end

% initialize Gibbs
Sigmadraw=SigmaHat;
SigmaHatInv = inv(SigmaHat);
% Check if SigmaHatInv is positive definite
% Compute the Cholesky decomposition of SigmaHatInv
L = chol(SigmaHatInv, 'lower');
% L*L' = SigmaHatInv, so R = L
R_base = L;
R = [R_base, zeros(nvar, nnuTilde - nvar)];

init = 0;


info.nshocks=3;

info.Ss{1,1}= [-1 0 0 0;
    0 -1 0 0;
    0  0 1 0;];


info.Ss{2,1}=   [1 0 0 0;
    0 1 0 0;
    0  0 1 0;];

info.Ss{3,1}=   [1 0 0 0;
    0 -1 0 0;
    0  0 1 0;
    0  0 0 1];

Bdraw=BHat;

z_old_R = vec(R);%chol(Rvariance,'lower')*randn(info.nvar*nnuTilde,1);
while init==0


    Xold =randn(info.nvar);
    [Qold,Rold]  = qr(Xold);
    for ii=1:nvar
        if Rold(ii,ii)<0
            Qold(:,ii)=-Qold(:,ii);
        end
    end

    Rr=chol(Sigmadraw)'*Qold;


    Tt=zeros(info.nshocks,nvar);
    for j=1:info.nshocks

        switch j



            case 1

                for i=1:nvar
                    if all(info.Ss{j,1}*Rr(:,i)>0)

                        Tt(j,i)=1;
                    elseif all(info.Ss{j,1}*Rr(:,i)<0)
                        Tt(j,i)=-1;
                    end
                end


            case 2

                for i=1:nvar
                    if all(info.Ss{j,1}*Rr(:,i)>0) && ((Rr(1,i)/Rr(3,i))<0.0258)
                        Tt(j,i)=1;
                    elseif all(info.Ss{j,1}*Rr(:,i)<0) && ((Rr(1,i)/Rr(3,i))<0.0258)
                        Tt(j,i)=-1;
                    end
                end

            case 3

                for i=1:nvar
                    if all(info.Ss{j,1}*Rr(:,i)>0) && ((Rr(1,i)/Rr(3,i))<0.0258)
                        Tt(j,i)=1;
                    elseif all(info.Ss{j,1}*Rr(:,i)<0) && ((Rr(1,i)/Rr(3,i))<0.0258)
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


    %if rank(Tt)<nshocks

    if any(sum(Tt == 0, 2)==size(Tt,2))

    else


        Rstar = zeros(nvar,nvar);
        S = cell(info.nshocks,1);

        isamples = zeros(info.nshocks,1);

        % For each row j
        for j = 1:info.nshocks

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
    
        for j = info.nshocks+1:nvar
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
            isamples(j,1) = sampled_element;

            if rand(1)>0.5
                Rstar(:,j) = Rr(:,isamples(j,1));
            else
                Rstar(:,j) = -Rr(:,isamples(j,1));
            end

        end

        Qdraw=(chol(Sigmadraw)')\Rstar;

        L0old=Rstar;


        %check dynamic sign restrictions
        % dynamic impulse responses
        hSigmadraw = chol(Sigmadraw);
        A0         = hSigmadraw\Qdraw;
        Aplus      = Bdraw*A0;
        for l=1:info.nlag-1
            info.A{l} = Aplus((l-1)*info.nvar+info.nex+1:l*info.nvar+info.nex,1:end);
            info.bigF((l-1)*info.nvar+1:l*info.nvar,1:info.nvar)=info.A{l}/A0;
        end
        A{info.nlag} = Aplus((info.nlag-1)*info.nvar+1:info.nlag*info.nvar,1:end);
        info.bigF((info.nlag-1)*info.nvar+1:info.nlag*info.nvar,:)=[A{info.nlag}/A0 info.extraF];



        activity_index = zeros(12,1);
        oil_price      = zeros(12,1);

        for h=0:1:11

            tmp = (info.J'*((info.bigF')^h)*info.J)*L0old;

            activity_index(h+1) =tmp(2,1);
            oil_price(h+1)      =tmp(3,1);

        end


        init = all(activity_index<0)*all(oil_price>0);




    end
end
% Assuming Qdraw is already defined as an orthogonal matrix of size n×n

% Create a sample upper triangular matrix R
n = size(Qdraw, 1);
Rdraw_check = triu(rand(n, n));  % Random upper triangular matrix

% Compute X
XX_check = Qdraw * Rdraw_check;

% % Verify using QR decomposition
% [Q_check, R_check] = qr(XX_check);
% keyboard
% % Adjust for potential sign differences
% for j = 1:n
%     if sign(Qdraw(1,j)) ~= sign(Q_check(1,j))  % Check first element of each column
%         Q_check(:,j) = -Q_check(:,j);      % Flip signs if different
%         R_check(j,:) = -R_check(j,:);
%     end
% end



%Bdraw=BHat;

%BSigmaQold = [vec(Bdraw);vec(SigmaHat);vec(Qold)];
%SigmaQold  = [vec(SigmaHat);vec(Qold)];



%% initialize counters to track the state of the computations
counter = 1;
count   = 0;
tStart = tic;
record = 1;






switch label_R

    case 'nolabel'

        function_restrictions = @restriction_none;

    case 'kilianmurphy'

        function_restrictions = @restriction_kilianmurphy;

        function_restrictions_Q = @restriction_kilianmurphy_withT;

    otherwise
        disp('choose a valid specification')
        return
end


z_old_X =vec(XX_check);
z_old_B =vec(Bdraw);

posterior_Sigmadraws_accepted = nan(maxBSQdraws2store,1);
tstart = tic;
while record<=maxBSQdraws2store
    record
    if fixed_rf==1
        Sigmadraw     = SigmaHat;
        Bdraw         = BHat ;
    elseif prior_only==1
        %% Draw Sigma
        Sigmadraw     = iwishrnd(PpsiBar,nnuBar);
        cholLSigmadraw = hh(Sigmadraw)';
        Bdraw         = kron(cholLSigmadraw,cholOomegaBar)*randn(m*nvar,1) + reshape(mmuBar,nvar*m,1);
        Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);
    else


        % % if counter_Q==0
        %
        %      %counter_Q = counter_Q +1;
        %      %% Draw Sigma
        %      Sigmadraw     = iwishrnd(PpsiTilde,nnuTilde);
        %      cholSigmadraw = hh(Sigmadraw)';
        %      Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*nvar,1) + reshape(mmuTilde,nvar*m,1);
        %      Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);
        %  %end

    end


    %% setting function restrictions
    function_restrictions_i = @(ag0, ag1, ag2, ag3, ag4, ag5) function_restrictions(ag0, ag1, ag2, ag3, ag4, ag5);


    %% Draw Q

    if elliptical_Q==1
            
        slice.scale_z    = 1;
        slice.mean       = zeros(info.nvar*info.nvar,1);
        slice.chol_cov_z = chol(eye(info.nvar*info.nvar),'lower');

        slice.fcn_lik    = @(z_prop) loglike_Q(z_prop,function_restrictions_i,gs_qr,fo_inv,fo_str2irfs,Bdraw,Sigmadraw,info);
        slice.nobs       = info.nvar*info.nvar;

        % --- actual slice sampling
        % z_old = randn(info.nvar*info.nvar,1);
        lik_old = slice.fcn_lik(z_old_X);


        [z_old_X, ~, n_try] = slice_sampling_v02(slice, z_old_X, lik_old);




        [Qdraw,Rdraw] =  gs_qr(reshape(z_old_X,info.nvar,info.nvar));

    else

        S_prop=0;
        function_restrictions_Q_i = @(ag0, ag1, ag2, ag3, ag4, ag5) function_restrictions_Q(ag0, ag1, ag2, ag3, ag4, ag5);

        while S_prop==0


            [Q_prop,~]=DrawQ(info.nvar);

            S_prop = function_restrictions_Q_i(Bdraw, Sigmadraw, Q_prop,fo_inv,fo_str2irfs,info);

    
        end
       

    end

    %% Draw Sigma
    if draw_Sigma==1
        slice.scale_z    = 1;



        slice.mean       = vec(Rmean);
        slice.chol_cov_z = chol(Rvariance,'lower');
        slice.fcn_lik    = @(z_prop) loglike_Sigma(z_prop,function_restrictions_i,fo_inv,fo_str2irfs,Bdraw,Qdraw,info,nnuTilde,mmuTilde,OomegaTilde);

        slice.nobs       = info.nvar*nnuTilde;

        % --- actual slice sampling

        lik_old = slice.fcn_lik(z_old_R);

        [z_old_R, ~, n_try] = slice_sampling_v02(slice, z_old_R, lik_old);

        R_old = reshape(z_old_R,info.nvar,nnuTilde);
        Sigmadraw=inv(R_old*R_old');

        posterior_Sigmadraws_accepted(record,1)=slice.fcn_lik(z_old_R);


    else
        Sigmadraw=SigmaHat;

    end

    if draw_B==1
        %% Draw B
        %cholSigmadraw = hh(Sigmadraw)';
        %Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*nvar,1) + reshape(mmuTilde,nvar*m,1);
        %Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);

        cholLSigmadraw = hh(Sigmadraw)';
        cholLvarB = kron(cholLSigmadraw,cholOomegaTilde);
        slice.scale_z    = 1;
        slice.mean       = vec(mmuTilde);
        slice.chol_cov_z = cholLvarB;%chol(Rvariance,'lower');

        slice.fcn_lik    = @(z_prop) loglike_B(z_prop,function_restrictions_i,fo_inv,fo_str2irfs,Sigmadraw,Qdraw,info,mmuTilde,OomegaTilde);

        slice.nobs       = info.nvar*info.m;

        % --- actual slice sampling

        lik_old = slice.fcn_lik(z_old_B);

        [z_old_B, ~, n_try] = slice_sampling_v02(slice, z_old_B, lik_old);

        Bdraw = reshape(z_old_B,info.m,info.nvar);

    else
        Bdraw         = BHat ;
    end


    x            = [vec(Bdraw); vec(Sigmadraw); vec(Qdraw)];
    structpara   = f_h_inv(x,info);
    irfpara      = fo_str2irfs(structpara);


    % if counter_Q==1000
    %        counter_Q=0;
    % end



    %% check if sign restrictions hold

    L0 =reshape(irfpara(1:nvar*nvar,:),nvar,nvar);

    BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];




    count=count+1;

    % store orthogonal reduced-form draws
    Bdraws{count,1}     = Bdraw;
    Sigmadraws{count,1} = Sigmadraw;
    Qdraws{count,1} = Qdraw;



    record=record+1;




    if counter==iter_show

        display(['Number of draws = ',num2str(record)])
        display(['Remaining draws = ',num2str(maxBSQdraws2store-(record))])
        counter =0;

    end
    counter = counter + 1;



end

telapsed = toc(tstart);


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


    savefile= ['matfiles/rgibbs_',label_R,'prior_only_',num2str(prior_only),'prior_',prior_type,'fixed_rf_',num2str(fixed_rf),'_ndraws',num2str(maxBSQdraws2store),'.mat'];

save(savefile,'L','cumL','horizon','telapsed');
cd ..
