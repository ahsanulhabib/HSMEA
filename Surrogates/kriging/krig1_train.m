function [krig_model,options] = krig1_train(xtr,ytr)
options = srgtsKRGSetOptions(xtr,ytr,@dace_fit,@dace_regpoly1);
[krig_model,state] = srgtsKRGFit(options);
krig_model.perf = state.KRG_DACEPerf;
return