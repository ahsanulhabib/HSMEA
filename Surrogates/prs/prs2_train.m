function [prs_model,options] = prs2_train(xtr,ytr)
options = srgtsPRSSetOptions(xtr,ytr,2);
prs_model = srgtsPRSFit(options);
prs_model.surro_type = 'prs2';
return