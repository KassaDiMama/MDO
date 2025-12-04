designVector = DesignVector();
const = Const();
wingDesign = WingDesign(designVector);

W_to = 52390;
W_zf = 46720;

emwet_wrapper(wingDesign,const, "test",W_to,W_zf)