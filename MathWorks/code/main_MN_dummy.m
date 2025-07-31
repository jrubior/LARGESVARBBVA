clearvars
close all
clc
% test github

data = fullfile('data.csv');
opts = detectImportOptions(data);
opts = opts.setvartype("date", "datetime");
opts = opts.setvaropts("date", "InputFormat",'yyyy.Q');
tbl = readtimetable(data, opts);
Y = tbl{:,2:end}*100;


%% Model parameters 
nvar      = size(Y,2);       % number of endogenous variables
nlag      = 2;               % number of lags
nex       = 1;               % set equal to 1 if a constant is included; 0 otherwise
m         = nvar*nlag + nex; % number of exogenous variables
T         = size(Y,1) - nlag; % sample
pos       = []; % if there is variables that are stationary.

%% MN Prior Parameters

disp('MN Prior')

Vc=1e4;
lambda1=.1;
lambda3=2;

llambda=zeros(3,1);
llambda(1,1)=1/lambda1;
llambda(2,1)=1;
llambda(3,1)=lambda3/2;

% nnuBar
nnuBar = nvar+2;

diagPphiBar=zeros(nvar,1);
for i=1:nvar
    Mdl = arima(nlag, 0, 0);       % AR(1) model: p=1, d=0, q=0
    EstMdl = estimate(Mdl, Y(:,i),'Display','off');
    diagPphiBar(i)=EstMdl.Variance*(nnuBar-nvar-1);
end
ppsi=diagPphiBar;

model = struct('Vc',Vc,'m',m','nvar',nvar,'nex',nex,'nlag',nlag,'llambda', llambda,'ppsi', ppsi);

[nnuBar,PpsiBar,mmuBar,OomegaBar] = minnesota_prior_with_dummies(model);

save('MNdummies.mat','nnuBar','PpsiBar','mmuBar','OomegaBar') 

