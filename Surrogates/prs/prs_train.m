function [prs_model,options] = prs_train(xtr,ytr)
options = srgtsPRSSetOptions(xtr,ytr,1);
prs_model = srgtsPRSFit(options);
return