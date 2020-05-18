%% simulate estimated model
clear all


original = csvread('parameters_original_value.csv');
to_be = csvread('parameters_to_be_estimated.csv');
estimated = csvread('estimated_parameters.csv');

param = original;
param(to_be) = estimated;

S = 200;
T = 1500;
parfor s = 1:S
%for s = 1:S   
    [Y_real,noprod,cactual,cdemand1,cdemand2,welfare_cs,welfare_dc1s,welfare_dc2s,welfare_c,welfare_dc1,welfare_dc2,Ygap,rDIV,rDIVK,dincome,dincomek,DIV,DIVK,OCC,wincome,deviation,opt_share,neut_share,pes_share,types_agg,gperiod,Growth,Dek,Dek2,Kdem2,bankruptcy_rate,residualC,YP_e,YP_u,prob_EE,prob_EU,prob_UU,prob_UE,prob_EE_p,prob_EU_p,prob_UU_p,prob_UE_p,prob_DD,prob_DN,prob_NN,prob_ND,prob_DD_k,prob_DN_k,prob_NN_k,prob_ND_k,price,Div_prob,Div_probb,Div_probk,profitsB,Cdivs_p,Kdivs_p,shares,dividendsB,EXPcontrol,Invent,Assets,baryk,Kdes,Kdem,Dexp,Kdivs,Cdivs,valI,unsatisfiedDemandK,EXP,X1,X2,totK, gdp_deflator, Investment,I, consumption, Prod_k,Prod_c, Un, totalDeb, totalDeb_k,public_debt,GB,TA,benefit, wage, desired_consumption,totE,average_interest_rate,impulse_rw_opt,YPE_true,YPE_anchor,YPU_true,YPU_anchor,rwage,Occ_status,Type,et] = Model(s,T, param);
    
    np(s)=noprod;
    Welfare_c(:,s)=welfare_c;
    Welfare_dc1(:,s)=welfare_dc1;
    Welfare_dc2(:,s)=welfare_dc2;
    Welfare_cs(:,s)=welfare_cs;
    Welfare_dc1s(:,s)=welfare_dc1s;
    Welfare_dc2s(:,s)=welfare_dc2s;
    c_ind(:,:,s)=cactual;
    dc1_ind(:,:,s)=cdemand1;
    dc2_ind(:,:,s)=cdemand2;
    Y(:,s)  = Y_real;
    P(:,s)  = price;
    In(:,s)  = Investment;
    brupt(:,s)=bankruptcy_rate;
    Ig(:,s) = I;
    C(:,s)  = consumption;
    DC(:,s) = desired_consumption;
    U(:,s) = Un;
    D(:,s) = totalDeb;
    Dk(:,s) = totalDeb_k;
    W(:,s) = wage;
    YC(:,s) = Prod_c*3;
    YK(:,s) = Prod_k*3;
    bYK(:,s) = baryk;
    elapsed(s) = et;
    Cap1(:,s)=X1;
    Cap2(:,s)=X2;
    K(:,s)=totK;
    G(:,s)=EXP;
    Gc(:,s)=EXPcontrol;
    K_des(:,s)= Kdes;
    K_dem(:,s)= Kdem;
    K_dem2(:,s)=Kdem2;
    Kgap(:,s)=unsatisfiedDemandK;
    D_e(:,s)=Dexp;
    D_ek(:,s)=Dek;
    D_ek2(:,s)=Dek2;
    K_divs(:,:,s)=Kdivs;
    C_divs(:,:,s)=Cdivs;
    K_divsp(:,:,s)=Kdivs_p;
    C_divsp(:,:,s)=Cdivs_p;
    B_divs(:,s)=dividendsB;
    PA(:,s)=Assets;
    Inv(:,s)=Invent;
    share(:,:,s)=shares;
    intervention(s)=gperiod;
    maxY(s)=find(Y_real==max(Y_real(1000:1400)));
    minY(s)=find(Y_real==min(Y_real(1000:1400)));
    Types(:,s)=types_agg;
    Deficit(:,s)=GB;
    Bonds(:,s)=public_debt;
    Taxes(:,s)=TA;
    Benefits(:,s)=benefit;
    dev(:,:,s)=deviation;
    Occ(:,:,s)=Occ_status;
    DD(:,s)=prob_DD;
    Winc(:,s)=wincome;
    Oc(:,s)=OCC;
    Rwage(:,s)=rwage;
    rd(:,s)=rDIV;
    rdk(:,s)=rDIVK;
    Dinc(:,s)=dincome;
    Dinck(:,s)=dincomek;
    div(:,s)=DIV;
    divk(:,s)=DIVK;
    Yg(:,s)=Ygap;
    opt(:,s)=opt_share;
    pes(:,s)=pes_share;
    neut(:,s)=neut_share;
end
