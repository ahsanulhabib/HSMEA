function [rbnn_model,options] = rbnn_train(xtr,ytr)
options = srgtsRBNNSetOptions(xtr,ytr);
rbnn_model = srgtsRBNNFit(options);
rbnn_model.surro_type = 'rbnn';
return