function [nnuBar,PpsiBar,mmuBar,OomegaBar] = minnesota_prior_with_dummies(model)

%% parameters

llambda = model.llambda;
ppsi = model.ppsi;


Vc   = model.Vc;
m    = model.m;
n    = model.nvar;
p    = model.nlag;
nex  = model.nex;


Ystar       = zeros(m +llambda(3,1)*n,n);
Xstar       = zeros(m +llambda(3,1)*n,m);

nnuBar = size(Ystar,1)-m -n - 1 + (n+1) +2; % the term n+1 comes from combining the likelihood for the dummy observations with an improper prior p(B,Sigma) propto |Sigma|^(-0.5(nvar+1))


% yBar    =   mean(Ypresample)'; % Del Negro and Schorfheide (2010)'s
% sBar    =   std(Ypresample)';  % Del Negro and Schorfheide (2010)'s

sBar=sqrt(ppsi);

% prior on B coefficients
% first lag of B
for i=1:n
    for j=1:n
        if i==j
            Ystar(i,j) = (1/llambda(1,1))*sBar(i,1);
        else
            Ystar(i,j) = 0;
        end
    end
end

for i=1:n
    for ell=1:m
        if i==ell
            Xstar(i,ell) = (1/llambda(1,1))*sBar(i,1);
        else
            Xstar(i,ell) = 0;
        end
    end
end
% prior for other lags of B
for ell=2:p
    for i=1:n

        for j=1:n
            if i==j
                Xstar(n*(ell-1)+i,n*(ell-1)+j) = (1/llambda(1,1))*sBar(i,1)*ell^llambda(2,1);

            else
                Xstar(n*(ell-1)+i,n*(ell-1)+j) = 0;
            end
        end
    end
end

% prior for the constant

for i=1:n

    Xstar(m,m) = 1/sqrt(Vc);

end

% prior for the covariance matrix Sigma
for ii=1:llambda(3,1)
    for i=1:n
        for j=1:n
            if i==j
                Ystar(m + i + n*(ii-1),j) = sBar(i,1);
            else
                Ystar(m + i + n*(ii-1),j) = 0;
            end
        end
    end
end

% % Switch first and last columns
% Xstar_noc = Xstar(:,1:end-nex);
% Xstar_c   = Xstar(:,end);
% Xstar     = [Xstar_c,Xstar_noc];

OomegaBarInverse    = (Xstar'*Xstar);
OomegaBar           = (Xstar'*Xstar)\eye(m);
mmuBar              = OomegaBar*(Xstar'*Ystar);
PpsiBar             = Ystar'*Ystar - mmuBar'*OomegaBarInverse*mmuBar;
PpsiBar = (PpsiBar+PpsiBar')/2;

