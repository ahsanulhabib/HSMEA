function [prs_model,options] = prs3_train(xtr,ytr)
options = srgtsPRSSetOptions(xtr,ytr,3);
prs_model = srgtsPRSFit(options);
return