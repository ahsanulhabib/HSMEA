function [prs_model,options] = prs1_train(xtr,ytr)
options = srgtsPRSSetOptions(xtr,ytr,1);
prs_model = srgtsPRSFit(options);
prs_model.surro_type = 'prs1';
return