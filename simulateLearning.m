%% simulate estimated model
clear all


original = csvread('parameters_original_value.csv');
to_be = csvread('parameters_to_be_estimated.csv');
estimated = csvread('estimated_parameters.csv');

param = original;
param(to_be) = estimated;
warning('off','MATLAB:nearlySingularMatrix')

S = 1;
T = 5000;
Learn=1;
Act=1;
%parfor s = 1:S
for s = 1:S   
    [Y_real,emp_count,div_count,divk_count,wincome,dincome,dincomek,coeffs,add_w,add_d,add_dk,YP_e,YP_u,YP_d,YP_nd,YP_dk,YP_ndk,gperiods,trans_EU,trans_EE,trans_UE,trans_UU,pub_exp_cr,Growth,prob_EE,prob_EU,prob_UU,prob_UE,prob_DD,prob_DN,prob_NN,prob_ND,prob_DD_k,prob_DN_k,prob_NN_k,prob_ND_k,price,EXPcontrol,Invent,Assets,baryk,valI,actualEXP, gdp_deflator, Investment,I, consumption, Prod_k, Prod_c, Un, totalDeb, totalDeb_k,stock_bonds,GB,TA,G,wages_t, desired_consumption,rwage,et] = learningModel(s,T, param,Learn,Act);
    winc(:,s)=wincome;
    dinc(:,s)=dincome;
    dinck(:,s)=dincomek;
    coef(:,:,s)=coeffs;
    Gov(:,s)=pub_exp_cr;
    DC(:,s)=desired_consumption;
    Y(:,s)=Y_real;
    C(:,s)=consumption;
    Ig(:,s)=I;
    AW(:,s)=add_w;
    AD(:,s)=add_d;
    ADK(:,s)=add_dk;
end
