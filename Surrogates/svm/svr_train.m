function [svr_model,options] = svr_train(xtr,ytr)
options = srgtsSVRSetOptions(xtr,ytr);
[svr_model,status] = srgtsSVRFit(options);
svr_model.status = status;
svr_model.surro_type = 'svr';
return