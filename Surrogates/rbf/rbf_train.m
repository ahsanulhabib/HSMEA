function [rbf_model] = rbf_train(xtr,ytr)
options = srgtsRBFSetOptions(xtr,ytr);
rbf_model = srgtsRBFFit(options);
rbf_model.options = options;
rbf_model.surro_type = 'rbf';
return