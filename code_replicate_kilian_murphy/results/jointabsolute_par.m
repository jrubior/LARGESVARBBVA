% JOINTABSOLUTE.M
% 
% Additively separable absolute loss function:
% Bayes estimator and joint credible set for any M x nirf matrix IRF, where 
% M is the number of posterior draws and nirf is the number of impulse 
% responses, when Bayes estimator is unique

function [irfabs,credibleset]=jointabsolute_par(IRF)

M=size(IRF,1);
absolute_loss = zeros(M,1);
parfor i=1:M
   for j=1:M
       absolute_loss(i,1) = absolute_loss(i,1)+sum(abs(IRF(j,:)-IRF(i,:)));
   end
end
[~,Imin]=min(absolute_loss);
irfabs = IRF(Imin,:);
[a,I] = sort(absolute_loss);
credibleset = IRF(I(1:ceil(0.68*M)),:);
