
function  [Y_real,emp_count,div_count,divk_count,wincome,dincome,dincomek,coeffs,add_w,add_d,add_dk,YP_e,YP_u,YP_d,YP_nd,YP_dk,YP_ndk,gperiods,trans_EU,trans_EE,trans_UE,trans_UU,pub_exp_cr,Growth,prob_EE,prob_EU,prob_UU,prob_UE,prob_DD,prob_DN,prob_NN,prob_ND,prob_DD_k,prob_DN_k,prob_NN_k,prob_ND_k,price,EXPcontrol,Invent,Assets,baryk,valI,actualEXP, gdp_deflator, Investment,I, consumption, Prod_k, Prod_c, Un, totalDeb, totalDeb_k,stock_bonds,GB,TA,G,wages_t, DC,rwage,et] = learningModel(seed, T, par,Learn,Act)
tic 

Bk=1;                               %no. of banks
F=100; %200;%                         %no. of firms/capitalists
W=1000; % 3000;%                        %no. of workers
N = 20;% 50;%                          %no. of capital producing firms
z_c=par(1);                              %no. of aplications in consumption good market
z_k = par(2);                         %no. of aplications in capital good market
z_e = par(3);                            %number of job applications  

rng(seed)

%SET RANDOM NUMBERS
shock_pk = rand(T,N);
shock_p  = rand(T,F);
prob_k = rand(T,F);
permutations_consumption = NaN(T,W+N+F,z_c); %might happen that two are the same...
permutations_un = NaN(T,W,z_e);
permutations_capital = NaN(T,F,z_k);

for tt = 1:T
    for ii = 1:W+N+F
        permutations_consumption(tt,ii,:) = randperm(F,z_c);
    end
    for ii = 1:W
        permutations_un(tt,ii,:) = randperm(F+N,z_e);
    end
    for ii = 1:F
        permutations_capital(tt,ii,:) = randperm(N,z_k);
    end
    
end

seeds_unemployed = randi(2^30,T,1);
seeds_capital = randi(2^30,T,1);
seed_consumption = randi(2^30,T,1);
seeds_surplusk = randi(2^30,T,N);
seeds_surplus = randi(2^30,T,F);

coeffsw=rand(3,1);
Rw=eye(3).*coeffsw*10000000;

coeffsd=rand(3,1);
Rd=eye(3).*coeffsd*10000000;

coeffsdk=rand(3,1);
Rdk=eye(3).*coeffsdk*10000000;

Gov=zeros(1,T);
gprds=randi([700 2500],1,round((2500-700)/10));
for i=1:length(gprds)
    per=gprds(i);
    Gov(per)=max(0,60*normrnd(1,0.3));
end
Gov(695:710)=max(0,60*normrnd(1,0.3,1,16));
nshocks=round((T-2500)/10);
for j=1:nshocks
    per=2500+10*j;
    Gov(per)=max(0,60*normrnd(1,0.3));
end



warning('off','MATLAB:nearlySingularMatrix')

step=1;

YPtest=1;
horizon=400;
periods=400;
beta=0.9999;
tau=2/3;
betasum=1;
for i=1:periods
    betasum=betasum+(beta^i)^(1/tau);
end
burnin=horizon+1;
r_f=0.01;%   0.015;% 0.005          %general refinancing rate
counter=0;
trans_UE=zeros(1,T);
trans_EU=zeros(1,T);
trans_UU=zeros(1,T);
trans_EE=zeros(1,T);
trans_ND=zeros(1,T);
trans_DN=zeros(1,T);
trans_DD=zeros(1,T);
trans_NN=zeros(1,T);
trans_ND_k=zeros(1,T);
trans_DN_k=zeros(1,T);
trans_DD_k=zeros(1,T);
trans_NN_k=zeros(1,T);
prob_EE=zeros(1,T);
prob_EE(1)=0.85;
prob_EU=zeros(1,T);
prob_EU(1)=1-prob_EE(1);
prob_UE=zeros(1,T);
prob_UE(1)=0.5;
prob_UU=zeros(1,T);
prob_UU(1)=1-prob_UE(1);
prob_DD=zeros(1,T);
prob_DD(1)=1;
prob_DN=zeros(1,T);
prob_DN(1)=1-prob_DD(1);
prob_ND=zeros(1,T);
prob_ND(1)=1;
prob_NN=zeros(1,T);
prob_NN(1)=1-prob_ND(1);
prob_DD_k=zeros(1,T);
prob_DD_k(1)=1;
prob_DN_k=zeros(1,T);
prob_DN_k(1)=1-prob_DD_k(1);
prob_ND_k=zeros(1,T);
prob_ND_k(1)=1;
prob_NN_k=zeros(1,T);
prob_NN_k(1)=1-prob_ND_k(1);

rw=zeros(1,T);
rw(1)=1/3;
rwage=zeros(1,T);
rwage(1)=rw(1);
rdiv=zeros(1,T);
rdivk=zeros(1,T);
rdivs=zeros(F,T);
rdivs_k=zeros(N,T);
rd=zeros(1,F);
rdk=zeros(1,N);
yd=zeros(1,F);
ydk=zeros(1,N);
residualC=zeros(1,T);
bonds_real=zeros(1,T);
deficit_real=zeros(1,T);
Ygap=zeros(1,T);

%set by the sensitivity
xi = par(4);                          %memory parameter human wealth
chi = par(5);                            %fraction of wealth devoted to consumption
q_adj = par(6);                            %quantity adjustment parameter
p_adj = par(7);                            %price adjustment parameter    
mu =  par(8);                           %bank's gross mark-up
eta =par(9);                        %capital depreciation
Iprob=par(10);                         %probability of investing
phi =  par(11);                        %bank's leverage parameter
theta=par(12);                         %rate of debt reimbursment
delta = par(13);                        %memory parameter in the capital utilization rate
alpha = par(14);                              %labour productivity
k = par(15);                            %capital productivity
div =par(16);                            %share of dividends
div_B=0.3;
barX=par(17);%0.85;                          %desired capital utilization
inventory_depreciation = par(18);%0.3;           %rate at which capital firms' inventories depreciate
b1 = par(19);   %-15;
b2 = par(20);   %13;                      %Parameters for risk evaluation by banks
b_k1 = par(21); %-5;
b_k2 = par(22); %5 ;

interest_rate = par(23);
subsidy = par(24);
tax_rate = par(27);

wage_update_up = par(28);
wage_update_down = par(29);
u_target = par(30);

%phillips curve
wb=1.5;                                % initial wage rate


%else government will use unemployment subsidies as line below.
bond_interest_rate = interest_rate;   %% ricordare di aggiungere anche questi con liquidità
unemployment_subsidy_init = subsidy;

G=zeros(1,T);                       %government expenditures
TA=zeros(1,T);                      %government income
GB=zeros(1,T);                      %governament budget  GB = TA - G - EXP -->???
EXP = Gov;       % spesa pubblica EROGABILE, update erogata in fondo
actualEXP=EXP;
pub_exp_cr=zeros(1,T+30);
pub_exp_cr(1:T)=EXP;
EXP(1,1)=0;              % inizializzazione, vedi sotto
Gshock=1;
EXPcontrol = zeros(1,T);

pub_exp_c= zeros(1,T);    % fu il valore totale EROGATO alle consumption firms 
tax_rate_d = 0;             %taxes on dividends

exp_c=zeros(F,T);         % valore EROGATO individualmente per updating liquidity
quota_exp_c=zeros(F,T);   % quota singola erogabile impresa per totale
quota_exp_c(:,1)=1/F;     % time 1: occhio se cambi iniziando con t=2
public_dem_c=zeros(F,T);  % domanda pubblica per le consumption

bonds = zeros(1,T);
bonds(1)=0;             
stock_bonds=zeros(1,T);

deficitPil=zeros(1,T);
primary_deficit_pil = zeros(1,T);
primary_GB= zeros(1,T);
dividendsB= zeros(1,T);

average_interest_rate=zeros(1,T);
average_interest_rate_k=zeros(1,T);
unfilledVacancies_k=zeros(1,T);
unfilledVacancies=zeros(1,T);
Occ_status=zeros(1,T,W);
Div_prob=zeros(1,T,F);
Div_probk=zeros(1,T,N);
Div_probb=zeros(1,T);
bankrupt=zeros(1,T,F+N);


%%Bank's Parameters
b = [b1;b2];
b_k = [b_k1;b_k2];
                         
%%Initial conditions
%capital firm
Leff_k = zeros(1,N);        %employees
Y_k =zeros(1,N);
Y_prev_k=3*ones(1,N);       %why is this not initialised to 0?

Y_kd=Y_prev_k;

P_k=3*ones(1,N);
A_k =10*ones(1,N);
liquidity_k = A_k;
De_k=ones(1,N);                 %expected demand of k-Firm
deb_k=zeros(1,N);
price_k=zeros(1,T+1);
price_k(1:2)=mean(P_k);                %capital price index
Q_k=zeros(1,N);                     %sales 
Ftot_k=zeros(1,N);                  %borrowing by k firms

interest_r_k=zeros(1,N);
initialA_k=zeros(1,T);
%firms
value_investments=zeros(1,F);
investment=zeros(1,F);              %pshysical capital acquired in the period
K=10*ones(1,F);
A=10*ones(1,F)+ K*price_k(1);                     %firm equity
liquidity=A-K*price_k(1);           %firm liquid resources

capital_value = K*price_k(1);

P=3*ones(1,F);                        %prices
Y_prev=5*ones(1,F);                   %past production

Yd=Y_prev;

Q=zeros(1,F);                       %sales   
Leff=zeros(1,F);                    %current employees
De=ones(1,F);                       %expected demand --> why is this not set to 5???
deb=zeros(1,F);                     %firm debts

barK=K;                             %long-run K
barYK=Y_prev/k;                     %
x = barYK./barK;                    %this gives initial capacity utilisation of 150? Why would you initialise it like that???
X1=zeros(1,T);                       %Macro-time series for capacity utilisation                        
X2=zeros(1,T);

interest_r=zeros(1,F);
%households
w=zeros(1,W);                       %earned wages
PA=ones(1,W+F+N)+1;                 %household personal assets (saving stock)
Oc=zeros(1,W);                      %employment status: Oc(i)=j --> worker j employed by firm i; if i=0, j is unemployed

totE=zeros(1,T);                    %bank equity time series
E=ones(1,Bk)*3000;                   %bank's equity
loans=0;                            %total loans
totE(1:2)=E;


%macro time series
consumption=zeros(1,T);             %total cunsumption 
price=zeros(1,T+1);
price(1:2)=P(1);                    %consumer price index time series
Un=zeros(1,T);                      %unemployment
wages_t = ones(1,T);                 %wages
dividends=zeros(1,T);               %total dividends

Ftot=zeros(1,F);                    %borrowings
defaults=zeros(1,T);                %number of bankruptcies
profitsB=zeros(T,Bk);               %bank profit time series
unsatisfiedDemand = zeros(1,T);
totK = zeros(1,T);
Investment = NaN(1,T);
inflationRate=zeros(1,T); %%time series of inflation rate
totalDeb=NaN(1,T);
totalDeb_k=NaN(1,T);
defaults_k=zeros(1,T); 

Y_nominal_k = NaN(1,T);
Y_nominal_c = NaN(1,T);
Y_nominal_tot= NaN(1,T);
Y_real =NaN(1,T);
gdp_deflator = NaN(1,T);
I = NaN(1,T);
Prod_c = NaN(1,T);
Prod_k = NaN(1,T);
baryk = NaN(1,T);
Kdes = NaN(1,T);
Kdem = NaN(1,T);
unsatisfiedDemandK = NaN(1,T);
valI=NaN(1,T);
Dexp = NaN(1,T);
Dek=NaN(1,T);
Dek2=NaN(1,T);
YPE=NaN(W,T);
YP_e=zeros(1,T);
YPU=NaN(W,T);
YP_u=zeros(1,T);
YP_d=zeros(1,T);
YP_dk=zeros(1,T);
YP_nd=zeros(1,T);
YP_ndk=zeros(1,T);
YPE(:,1)=rw(1);
YPU(:,1)=rw(1)*subsidy;
Kdivs = NaN(1,T,N);
Cdivs = NaN(1,T,F);
shares = NaN(1,T,F);
Assets = NaN(1,T);
Invent = NaN(1,T);
Kdivs_p=NaN(1,T,N);
Cdivs_p=NaN(1,T,F);
net_money=NaN(1,T);
bankruptcy_rate=NaN(1,T);
P_lower=NaN(1,T);
Kdem2=NaN(1,T);
DC=NaN(1,T);
Growth=NaN(1,T);
credit_mismatch=NaN(1,T);
Occ_status_prev=NaN(T,W);
%%%END OF MODEL INITIALIZATION

% cons_budget = ones(1,F+N+W)*1;
dividends_income_k = zeros(1,N);
dividends_income   = zeros(1,F);
permanent_income = ones(F+W+N,1)./price(1);
%%%HEREAFTER THE MODEL STARTS RUNNING FOR T PERIODS

money(1:2) = sum(PA)+sum(liquidity)+sum(liquidity_k)-sum(deb)-sum(deb_k)+totE(1);  %What is this supposed to be? Why is bank's equity part of "money"? Why is debt subtracted???
learn=Learn;
learn_act=Act;
wincome=zeros(1,T);
dincome=zeros(1,T);
dincomek=zeros(1,T);
emp_count=zeros(1,T);
div_count=zeros(1,T);
divk_count=zeros(1,T);
coeffs=zeros(3,T);
add_w=zeros(1,T);
add_d=zeros(1,T);
add_dk=zeros(1,T);

%%%%%%%%%%%%%%%%%%
 for t=2:T
     
     if learn_act==1 && t>2500 && pub_exp_cr(t+5)>0
         pub_exp_cum=pub_exp_cr(t+5);
         add_w(t)=round(coeffs(1,t-1)*pub_exp_cum)/W;
         add_d(t)=round(coeffs(2,t-1)*pub_exp_cum)/F;
         add_dk(t)=round(coeffs(3,t-1)*pub_exp_cum)/N;
     end
     
     if A<=0                                                             %if all the firms are defaulted, then stop simulation 
        break
     end
    
     EXP(t)=EXP(t)*price(t);
     if Gshock==1
     gexp=transpose(quota_exp_c(:,t-1).*(EXP(t))./P(:));
     else
     gexp=0;
     end
     
    stock=Y_prev-Yd;
    Invent(t) = sum(stock)*3;
        for i=1:F
            if stock(i)<=0 && P(i)>=price(t)
                if Gshock==1 && gexp(i)>0
                    De(i) = Y_prev(i) + (-stock(i))*q_adj + gexp(i);
                else
                    De(i) = Y_prev(i) + (-stock(i))*q_adj;
                end
            elseif stock(i)<=0 && P(i)<price(t)
                P(i)=P(i)*(1+shock_p(t,i)*p_adj);
                if Gshock==1 && gexp(i)>0
                    De(i)=Y_prev(i);
                else
                    De(i)=Y_prev(i);
                end
            elseif stock(i)>0 &&P(i)>price(t)
                P(i)=P(i)*(1-shock_p(t,i)*p_adj);
                if Gshock==1 && gexp(i)>0
                    De(i)= Y_prev(i);
                else
                    De(i)= Y_prev(i);
                end
            elseif stock(i)>0 && P(i)<=price(t)
                if Gshock==1 && gexp(i)>0
                    De(i) = Y_prev(i) - stock(i)*q_adj + gexp(i);
                else
                    De(i) = Y_prev(i) - stock(i)*q_adj;
                end
            end    
            if De(i)<alpha                                                  %in order to hire at least one worker
                De(i)=alpha;    
            end
        end
    Dexp(t)= sum(De)*3;
    
    %CAPITAL PRODUCTION DECISION
    inventory_k=(1-inventory_depreciation)*Y_k;
   
    stock_k = Y_prev_k-Y_kd; 
    for i=1:N
        if stock_k(i)<=0 && P_k(i)>=price_k(t)
             De_k(i) = (Y_prev_k(i) + (-stock_k(i))*q_adj)-inventory_k(i);
        elseif stock_k(i)<=0 && P_k(i)<price_k(t)
            De_k(i)=Y_prev_k(i);
            P_k(i)=P_k(i)*(1+shock_pk(t,i)*p_adj);
            
        elseif stock_k(i)>0 && P_k(i)>price_k(t)
            De_k(i)=Y_prev_k(i); 
            P_k(i)=P_k(i)*(1-shock_pk(t,i)*p_adj);
           
        elseif stock_k(i)>0 && P_k(i)<=price_k(t)
             De_k(i) = (Y_prev_k(i) - (stock_k(i))*q_adj)-inventory_k(i);
        end    
        
        if De_k(i)<alpha                                                  %in order to hire at least one worker
            De_k(i)=alpha;    
        end
    end    
    
   
    
    %% investments

    prob=prob_k(t,:);
    K_dem=zeros(1,F);
    K_des = barYK/barX;
    Kdes(t)= sum(K_des);
    depreciation = K*eta*1/(Iprob);
    K_dem(prob<Iprob)= K_des(prob<Iprob) - K(prob<Iprob)+depreciation(prob<Iprob);
    K_dem(K_dem<0)=0;
    Kdem2(t) = sum(K_dem);
    
    %labour requirement (consumption good)
    Ld=min(ceil(De./alpha), ceil(K.*k/alpha));  
    wages=wb*Ld;   
    p_k=price_k(t);
    
   
    %labour requirement (capital good)
    Ld_k=ceil(De_k./alpha);
    wages_k = wb*Ld_k;

%% CREDIT MARKET OPENS
             
    
    Ftot(:)=0;
    Ftot_k(:)=0;
   
    
    %compute financial gap
    B=wages+K_dem*p_k-liquidity;                                        %financial gap          
    B_k =wages_k-liquidity_k;                                           %financial gap of capital producers
    B(B<0)=0;
    B_k(B_k<0)=0;
    
    lev=(deb+B)./(A+B+deb);                                             %leverage
    lev_k=(deb_k+B_k)./(A_k+B_k+deb_k);
    
    
    loan_applications=find(B>0);                                        %only firms with positive credit demand will apply for an additional loan
    loan_applications_k=find(B_k>0); 
    %evaluate bankruptcy probability and expected survival
    pr(:,1)= exp(b(1)+b(2)*lev)./(1+exp(b(1)+b(2)*lev)); %glmval(b,lev,'logit');  %zeros(F,1);%                           %banks evaluate bankruptcy probability of each firm, given the estimated
    pr_k(:,1)=exp(b_k(1)+b_k(2)*lev_k)./(1+exp(b_k(1)+b_k(2)*lev_k));%glmval(b_k,lev_k,'logit');%zeros(N,1);%                        %parameters (b,b_k)and computed leverage                                                                        
                                                           
    Xi=(1-(1-theta).^(1+1./pr))./(1-(1-theta));                         %this is Xi in the paper
    Xi_k=(1-(1-theta).^(1+1./pr_k))./(1-(1-theta));
    %proposed rate depends on the estimated bankruptcy probability
    proposed_rate=mu*((1+r_f/theta)./Xi - theta)';
    proposed_rate_k=mu*((1+r_f/theta)./Xi_k - theta)';
    
    

    %for each firm the bank computes the maximum loan and gives loans up to the maximum amount 
    for i=loan_applications
        credit=B(i);
        %the bank gives a maximum credit depending
        %on  maximum expected loss
        
        maxL = (phi*totE(t)-pr(i)*deb(i))/pr(i);
        maxL=max(0,maxL); %maxL never negative        
        credit=min(credit,maxL); %%credit given to firm i           

        
        deb_0=deb(i);
        deb(i)=deb(i)+credit;                                           %update firm's debt stock
        loans=loans+credit;                                             %update bank's credit stock
        Ftot(i)=credit;                                                 %record flow of new credit for firm i
        %compute new average interest rate
        if deb(i)>0
            interest_r(i)=(deb_0*interest_r(i)+proposed_rate(i)*credit)/deb(i);
        end
    
    end       
    %weighted average interest rate
     average_interest_rate(t)=sum(interest_r.*deb)/sum(deb);

     
    %mutatis mutandis for capital firms
    for i=loan_applications_k
        credit=B_k(i);
        maxL = (phi*totE(t)-pr_k(i)*deb_k(i))/pr_k(i);
        maxL=max(0,maxL);
        credit=min(credit,maxL);
         
        deb_0=deb_k(i);
        deb_k(i)=deb_k(i)+credit;                                       %update firm's debt stock
        loans=loans+credit;                                             %update bank's credit stock
        Ftot_k(i)=credit;                                               %record flow of new credit for firm i
        if deb_k(i)>0
            interest_r_k(i)=(deb_0*interest_r_k(i)+proposed_rate_k(i)*credit)/deb_k(i);
        end
       
    end 
    average_interest_rate_k(t)=sum(interest_r_k.*deb_k)/sum(deb_k);
    credit_mismatch(t) = sum(B(loan_applications))+sum(B_k(loan_applications_k)) - sum(Ftot)-sum(Ftot_k);

    %%CREDIT MARKET CLOSES    

    %% JOB MARKET OPENS
    Occ_status_prev(t,:)=Oc;
    Occ_status_prev(Occ_status_prev>0)=1;
    
    %determine desired labour and vacancies given available liquidity
    Ld_k=min(Ld_k, (Ftot_k+liquidity_k)/wb);
    Ld_k=floor(Ld_k);
    Ld_k(Ld_k<1)=1;
    vacancies_k=Ld_k-Leff_k;
    surplus_k=find(vacancies_k<0);                                          %firms with too many workers
    
    %%CONSUPTION GOOD

    %%re-define labour demand given available liquidity
    Ld=min(Ld,(Ftot+liquidity)/wb);                                     %demand (stock)     
    Ld=floor(Ld);
    Ld(Ld<1)=1; %%since there is the capital the firms can have positive equity and negative liquidity, in this latter case Ld would be negative, which is impossible
    vacancies=Ld-Leff;                                                  %definitive labour demand (flow)    
    
    %%JOB MARKET OPENS
    surplus=find(vacancies<0);                                          %firms with too many workers
    

    for i=surplus_k
        workforce_k=find(Oc==F+i);
        rng(seeds_surplusk(t,i))%pick all firm i's workers
        f_k=randperm(length(workforce_k));
        f_k=f_k(1:-vacancies_k(i));                                     %take randomly "-vacancies(i)" workers and fire them
        fired_k=workforce_k(f_k);  
        Oc(fired_k)=0;
        w(fired_k)=0;
        Leff_k(i)=Ld_k(i);                                              %update no. of workers
     end


 %firms with excess workforce fire
    for i=surplus   
        workforce=find(Oc==i);  
        rng(seeds_surplus(t,i))%pick all firm i's workers
        f=randperm(length(workforce));
        f=f(1:-vacancies(i));                                           %take randomly "-vacancies(i)" workers and fire them
        fired=workforce(f);  
        Oc(fired)=0;
        w(fired)=0;
        Leff(i)=Ld(i);                                                  %update no. of workers
    end
    
%% UNEMPLOYED WORKERS LOOK FOR A JOB
    
    unemployed=find(Oc==0);
    
    
    rng(seeds_unemployed(t))
    vec = randperm(length(unemployed));
    for un=vec      
        j=unemployed(un);                                               %randomly pick an unemployed worker

        Z_e = permutations_un(t,j,:); 
        flag=1;
        
        while (Oc(j)==0 && flag<=z_e)                              %continue searching until you are unemployed and didn't search at all available firms
            f=Z_e(flag);                                        %pick first best firm; with flag increasing, pick the second best, the third...   
            if f>F %selected firm is a capital firm
                if vacancies_k(f-F)>0                                    %if the selected firm has an open vacancy, take the job
                Oc(j)=f;                                                 %update employed status
                w(j)=wb;                                                 %salary   
                Leff_k(f-F)=Leff_k(f-F)+1;                               %firm's workforce   
                vacancies_k(f-F)=vacancies_k(f-F)-1;                     %firm's vacancies   
                end
            else %selected firm is a consuption firm
                if vacancies(f)>0                                        %if the selected firm has an open vacancy, take the job
                Oc(j)=f;
                w(j)=wb;
                Leff(f)=Leff(f)+1;
                vacancies(f)=vacancies(f)-1;
                end                
            end
            flag=flag+1;                                                %increase counter
        end 
        
    end  
    
    %%JOB MARKET CLOSES
    
    unfilledVacancies_k(t) = sum(vacancies_k(vacancies_k>0)); 
    unfilledVacancies(t) = sum(vacancies(vacancies_k>0));
    
    %% production
    %produce capital
    Dek(t)=sum(De_k);
    Dek2(t)=sum(Leff_k*alpha);
    Y_k=min(De_k,Leff_k*alpha);                              
    Y_prev_k=Y_k;
    Prod_k(t)=sum(Y_prev_k);
  
    Y_k = Y_k+inventory_k;                                   %Y_k is increased by the inventories (capital good is durable)   
    
    %produce consuption
    Yp = min(Leff*alpha, K*k);                               %production frontier given available resources (Leontieff)
    Y=min(De,Yp);                                            %actual production
    
    Ygap(t)=sum(Yp)-sum(Y);
    Y_prev=Y;
    Prod_c(t)=sum(Y_prev);
    X1(t)=sum(Y_prev./sum(Y_prev).*Y_prev./(K*k));
    X2(t)=sum(Y_prev)/(sum(K)*k);
    totK(t)=sum(K);
 
    interests=interest_r.*deb; 
    interests_k=interest_r_k.*deb_k;
    wages=wb*Leff;
    wages_k=wb*Leff_k;

    
    %minimum prices
    
    Pl = 1.01*(wages+interests+Y/k*eta*p_k)./Y_prev;  
    Pl(isinf(Pl)) = 0;
    P(P<Pl) = Pl(P<Pl);

    
    %minimum prices capital
    Pl_k = 1.01*(wages_k+interests_k)./(Y_prev_k);
    Pl_k(isinf(Pl_k)) = 0;
    P_k(P_k<Pl_k) = Pl_k(P_k<Pl_k);
    

  
    %%CAPITAL GOODS MARKET OPENS
    capital_budget=max(0,Ftot+liquidity-Leff*wb);                              %amount of liquidity available to buy capital goods
    %capital_budget=Ftot+liquidity-Leff*wb;                              %amount of liquidity available to buy capital goods
    capdem=min(K_dem,capital_budget./price_k(t));
    Kdem(t) = sum(min(K_dem,capital_budget./price_k(t)));
    capital_demanders=find(K_dem>0);                                    
    investment(:)=0;
    value_investments(:)=0;
    Q_k(:)=0;
    Y_kd=zeros(1,N); 
    
    rng(seeds_capital(t))
    for firm=randperm(length(capital_demanders))
        j=capital_demanders(firm);
        Z = permutations_capital(t,j,:);
        PZ_k=P_k(Z);
        [~,order]=sort(PZ_k);
        flag=1;

        while (capital_budget(j)>0 && K_dem(j)>0 && flag<=z_k) 
            best=Z(order(flag));
            Y_kd(best)=Y_kd(best)+min(capital_budget(j)/P_k(best), K_dem(j));
                if Y_k(best)>0                                                %buy if 'best' firm has still positive stocks 
                    pk=P_k(best);
                    budget=min(capital_budget(j), K_dem(j)*pk);
                    if Y_k(best) > budget/pk           
                        Y_k(best)=Y_k(best)-budget/pk;                   %reduce stocks
                        Q_k(best)=Q_k(best)+budget/pk;                   %update sales
                        K_dem(j)=K_dem(j)-budget/pk;
                        capital_budget(j) = capital_budget(j)-budget;
                        liquidity(j)=liquidity(j)-budget;
                        investment(j)=investment(j)+budget/pk;
                        value_investments(j)= value_investments(j)+budget;
                                                           %j spends all its budget  
                    elseif  Y_k(best) <= budget/pk  
                        K_dem(j)=K_dem(j)-Y_k(best);
                        capital_budget(j)=capital_budget(j) - Y_k(best)*pk;
                        liquidity(j)=liquidity(j)- Y_k(best)*pk;
                        Q_k(best)=Q_k(best)+Y_k(best);    
                        investment(j)=investment(j)+Y_k(best);
                        value_investments(j)= value_investments(j)+Y_k(best)*pk;
                        Y_k(best)=0;    
                    end    
                end
                flag=flag+1;                                                %increase counter
        end  
    end
    
    unsatisfiedDemandK(t) = sum(max(0,capdem-investment));
    valI(t)=sum(value_investments);
    I(t) = sum(investment);

    %%CAPITAL GOOD MARKET CLOSES
    
 
    %%CONSUMPTION GOOD MARKET OPENS
    w(w>0) = wb;
    %taxes on wages
    wn=w*(1-tax_rate);  
    TA(t) = TA(t)+sum(w*tax_rate);
    %workers receive wages

    if YPtest==0 || YPtest==1 && t<500
        PA(1:W)=PA(1:W)+wn;
    end
    
    %%% SET UNEMPLOYMENT SUBSIDY TO MEET THE PUBLIC
    %%% DEFICIT RATIO REQUIREMENT

    unemployed=find(Oc==0);
    
    unemployment_subsidy = unemployment_subsidy_init;

    
    if YPtest==0 || YPtest==1 && t<500
        PA(unemployed)=PA(unemployed)+unemployment_subsidy*wb;
    end
    G(t)=unemployment_subsidy*wb*length(unemployed);
    
    
    workers_income = wn;
    workers_income(unemployed) = unemployment_subsidy*wb;
    
    Occ_status(:,t,:)=Oc;
    Occ_status(Occ_status>0)=1;
    

if YPtest == 1
    for j=1:W
    if Occ_status(1,t,j)==1 && Occ_status_prev(t,j)==0
        trans_UE(t)=trans_UE(t)+1;
    end
    if Occ_status(1,t,j)==0 && Occ_status_prev(t,j)==1
        trans_EU(t)=trans_EU(t)+1;
    end
    if Occ_status(1,t,j)==0 && Occ_status_prev(t,j)==0
        trans_UU(t)=trans_UU(t)+1;
    end
    if Occ_status(1,t,j)==1 && Occ_status_prev(t,j)==1
        trans_EE(t)=trans_EE(t)+1;
    end
    end
    
    if sum(trans_EE(1:t))+sum(trans_EU(1:t))>0
        if t<burnin
        prob_EE(t)=sum(trans_EE(1:t))/(sum(trans_EE(1:t))+sum(trans_EU(1:t)));
        prob_EU(t)=sum(trans_EU(1:t))/(sum(trans_EE(1:t))+sum(trans_EU(1:t)));
        else
        prob_EE(t)=sum(trans_EE((t-(horizon-1)):t))/(sum(trans_EE((t-(horizon-1)):t))+sum(trans_EU((t-(horizon-1)):t)));
        prob_EU(t)=sum(trans_EU((t-(horizon-1)):t))/(sum(trans_EE((t-(horizon-1)):t))+sum(trans_EU((t-(horizon-1)):t)));
        end
    else
    prob_EE(t)=prob_EE(t-1);
    prob_EU(t)=prob_EU(t-1);
    end

    if sum(trans_UE(1:t))+sum(trans_UU(1:t))>0
        if t<burnin
        prob_UU(t)=sum(trans_UU(1:t))/(sum(trans_UE(1:t))+sum(trans_UU(1:t)));
        prob_UE(t)=sum(trans_UE(1:t))/(sum(trans_UE(1:t))+sum(trans_UU(1:t)));
        else
        prob_UU(t)=sum(trans_UU((t-(horizon-1)):t))/(sum(trans_UE((t-(horizon-1)):t))+sum(trans_UU((t-(horizon-1)):t)));
        prob_UE(t)=sum(trans_UE((t-(horizon-1)):t))/(sum(trans_UE((t-(horizon-1)):t))+sum(trans_UU((t-(horizon-1)):t)));   
        end
    else
    prob_UU(t)=prob_UU(t-1);
    prob_UE(t)=prob_UE(t-1);
    end 
    
    
    rwage(t)=(wb*(1-tax_rate))/price(t);
    
    if t < burnin
        rw(t)=sum(rwage(1:t))/t;
        rs=rw(t)*unemployment_subsidy;
    else
        rw(t)=sum(rwage(t-(horizon-1):t))/horizon;
        rs=rw(t)*unemployment_subsidy;
    end
    
    mat=[prob_EE(t) prob_EU(t); prob_UE(t) prob_UU(t)];
    pay=[rw(t); rs];
    permw=pay;
    for i=1:periods
        prod=mat^i;
        permw=permw+prod*pay;
    end
    permwdiff=permw(1)-permw(2);
    if learn_act==1 && t>2500 && pub_exp_cr(t+5)>0
        permw(:)=permw(:)+beta^5*permwdiff*add_w(t);
    end
    if learn_act==1 && t>2500 && pub_exp_cr(t+4)>0
        permw(:)=permw(:)+beta^4*permwdiff*add_w(t-1);
    end
    if learn_act==1 && t>2500 && pub_exp_cr(t+3)>0
        permw(:)=permw(:)+beta^3*permwdiff*add_w(t-2);
    end
    if learn_act==1 && t>2500 && pub_exp_cr(t+2)>0
        permw(:)=permw(:)+beta^2*permwdiff*add_w(t-3);
    end
    if learn_act==1 && t>2500 && pub_exp_cr(t+1)>0
        permw(:)=permw(:)+beta*permwdiff*add_w(t-4);
    end
    if learn_act==1 && t>2500 && pub_exp_cr(t)>0
        permw(:)=permw(:)-permwdiff*add_w(t-5);
    end
    YPE(:,t)=(PA(1:W)/price(t)+permw(1))/betasum;
    YPU(:,t)=(PA(1:W)/price(t)+permw(2))/betasum;
    
    YP_e(t)=permw(1);
    YP_u(t)=permw(2);
    
      
end
    
    if YPtest==0 || YPtest==1 && t<500
        income =  ([workers_income,dividends_income,dividends_income_k]')./price(t);
        permanent_income = permanent_income*xi + (1-xi)*income;
    else
        wi=zeros(1,W);
        for i=1:W
            if Oc(i)>0
                wi(i)=YPE(i,t);
            else
                wi(i)=YPU(i,t);
            end
        end
       income = [wi,yd,ydk]';
       permanent_income = income;
    end
    
    if YPtest==1 && t>=500
        PA(1:W)=PA(1:W)+wn;
        PA(unemployed)=PA(unemployed)+unemployment_subsidy*wb;
    end

    if YPtest==0 || YPtest==1 && t<500
        target = 1*permanent_income' + chi*PA./price(t) ; %0.05
    else
        target = 1*permanent_income';
    end
    
    cons_budget = target.*price(t);
    cons_budget = min(PA,cons_budget);
    PA=PA-cons_budget;        
    consumers=find(cons_budget>0);
    
    DC(t) = (sum(cons_budget)./price(t))*3;

    Q(:)=0;
    Yd=zeros(1,F);
   %search and matching starts

    C = zeros(1,W+F+N);
    
    
    rng(seed_consumption(t))
    vec = randperm(length(consumers));
    for wor=vec     
        j=consumers(wor);                                               %randomly pick a consumer
        Z = permutations_consumption(t,j,:);
        PZ=P(Z);  
        [~,order]=sort(PZ);                                            %sort prices in ascending order
        flag=1;
        
        while (cons_budget(j)>0 && flag<=z_c)                              %continue buying till budget is positive and there are firms available
            best=Z(order(flag));                                        %pick first best firm; with flag increasing, pick the second best, the third...   
            Yd(best)=Yd(best)+cons_budget(j)/P(best);
            if Y(best)>0                                                %buy if 'best' firm has still positive stocks 
                p=P(best);
                if Y(best) > cons_budget(j)/p           
                    Y(best)=Y(best)-cons_budget(j)/p;                   %reduce stocks
                    Q(best)=Q(best)+cons_budget(j)/p; 
                    C(j) = C(j)+cons_budget(j)/p; %update sales
                    consumption(t)=consumption(t)+cons_budget(j)/p;     
                    cons_budget(j)=0;                                   %j spends all its budget  
                elseif  Y(best) <= cons_budget(j)/p  
                    cons_budget(j)=cons_budget(j)- Y(best)*p;   
                    Q(best)=Q(best)+Y(best);
                    C(j) = C(j)+Y(best); 
                    consumption(t)=consumption(t)+Y(best);
                    Y(best)=0;    
                end    
            end
            flag=flag+1;                                                %increase counter
       end 
        
    end    
    
    consumption(t)=consumption(t)*3;
    unsatisfiedDemand(t) = sum(cons_budget);            
    PA=PA+cons_budget;
    
    %%%INTRO SPESA PUBBLICA
    pub_exp_c(t)= 1*EXP(t);
    public_dem_c(:,t)=quota_exp_c(:,t-1).*(pub_exp_c(t)./P');  
    EXPcontrol(t)=sum(public_dem_c(:,t));
    for i=1:F
        %Yd(i)=Yd(i)+public_dem_c(i,t);
        if Y(i)>= public_dem_c(i,t)
           Y(i)=Y(i)-public_dem_c(i,t);   
           Q(i)= Q(i)+public_dem_c(i,t);   
        else 
           Q(i)=Q(i)+Y(i);
           public_dem_c(i,t)=Y(i);              
           Y(i)= 0;      
        end                                     
        %Yd(i)=Yd(i)+public_dem_c(i,t);
    end
    exp_c(:,t)= public_dem_c(:,t).*P';     
    pub_exp_c(t)=sum(exp_c(:,t));
    EXP(t)=pub_exp_c(t);
    actualEXP(t)=pub_exp_c(t)/price(t);
    residualC(t)=sum(Y);
    
%%CONSUMPTION GOOD MARKET CLOSES
    %Capital price index
    if sum(Q_k)>0 
        share_k=(Q_k)/sum(Q_k);
    else
        share_k=ones(1,N)/N;
    end
    
    price_k(t+1)=P_k*share_k';

%%ACCOUNTING

    %capital capacity update
    barYK = delta*barYK + (1-delta)*Y_prev/k;
    baryk(t) = mean(barYK);
    %capital depreciation
    dep = eta*Y_prev/k; %%only the used capital is depreciating
    %capital reduced by depreciation
    K = K - dep;
    %capital increased by bought capital 
    K = K + investment;

    %update capital value in the book
    depreciation_value =  (dep).*capital_value./(K+dep-investment);
    capital_value = capital_value - depreciation_value + value_investments;
    
    %firm revenues
    RIC=P.*Q;   
    %firm update their liquidity, pay wages, interests and installments
    liquidity=liquidity+RIC+Ftot-wages-interests-theta*deb; %%investments already paid
    loans=loans-sum(deb)*theta;
    deb=(1-theta)*deb;
    
    %consumption firm profits
    pi = RIC-wages-interests - depreciation_value;
    %equity law of motion!
    A = A+pi;
    
    %Capital producer accounting
    RIC_k = P_k.*Q_k;
    %update liquidity
    liquidity_k = liquidity_k + RIC_k + Ftot_k-wages_k-interests_k-theta*deb_k;
    %update loans
    loans=loans-sum(deb_k)*theta;
    deb_k=(1-theta)*deb_k;
    %profits
    %inventories missing here.
    pi_k=RIC_k-wages_k-interests_k;
    
    totalDeb(t)=sum(deb);
    totalDeb_k(t)=sum(deb_k);
    
    %dividends
    divstatus_prev=dividends_income;
    divstatus_prev(divstatus_prev>0)=1;
    pospi=find(pi>0);                                                   %pick firms with positive profits
    Div_prob(:,t,:)=pi;
    Div_prob(Div_prob<0)=0;
    Div_prob(Div_prob>0)=1;
    dividends_income(:)=0;
    for i=pospi
        di=div*pi(i);  %%dividends                                              %compute dividends paid by firm i
        divi=di*(1-tax_rate_d); %dividends after taxes
        TA(t)=TA(t)+di*tax_rate_d;
        PA(W+i)=PA(W+i)+divi;                                          %dividends paid to firm i owner
        dividends_income(i) = divi;
        liquidity(i)=liquidity(i)-di;
        A(i)=A(i)-di;
        dividends(t)=dividends(t)+di; %lordi
        pi(i)=pi(i)-divi;
       
    end
    
    divstatus_prev_k=dividends_income_k;
    divstatus_prev_k(divstatus_prev_k>0)=1;
    pospi_k=find(pi_k>0); 
    Div_probk(:,t,:)=pi_k;
    Div_probk(Div_probk<0)=0;
    Div_probk(Div_probk>0)=1;
    dividends_income_k(:)=0;%pick firms with positive profits
    for i=pospi_k
        di=div*pi_k(i);   
        divi=di*(1-tax_rate_d); %dividends after taxes
        TA(t)=TA(t)+di*tax_rate_d;%compute dividends paid by firm i
        PA(W+F+i)=PA(W+F+i)+divi;                                          %dividends paid to firm i owner
        dividends_income_k(i)=divi;
        liquidity_k(i)=liquidity_k(i)-di;
        dividends(t)=dividends(t)+di;
        pi_k(i)=pi_k(i)-divi;
       
    end
    
    
   
    A_k=liquidity_k+Y_k*price_k(t)-deb_k;

    
    %replacement of bankrupted consumption firms 
    piB=0;                                                              %reset bank's profits
       
    
    %% time series (before bankruptcies)
    %inflation rate
    if sum(Q)>0 
        share=(Q)/sum(Q);
    else
        share=ones(1,F)/F;
    end
    P_lower(t)= Pl*share';
    RPI=P*share';                                                       %retail price index
    infla_rate=RPI/price(t);
    %price=[price RPI];
    price(t+1)=RPI;
    inflationRate(t) = infla_rate;
    %unemployment rate
    disocct=(W-length(Oc(Oc>0)))/W;
    Un(t)=disocct;
    
    Y_nominal_k(t) = sum(Y_prev_k)*price_k(t);
    Y_nominal_c(t) = sum(Y_prev)*price(t);
    Y_nominal_tot(t)= Y_nominal_k(t)+Y_nominal_c(t);
    Y_real(t) = sum(Y_prev)*price(1) + sum(Y_prev_k)*price_k(1);
    Growth(t) = (Y_real(t)-Y_real(t-1))/Y_real(t-1);
    gdp_deflator(t) = Y_nominal_tot(t)/ Y_real(t);
     
    Investment(t)=I(t)*price_k(1)+sum(Y_k-inventory_k)*price_k(1)+price(1)*sum(Y); %total investment is investment plus inventory variation
 
    negcash_k=find(A_k<=0);
    negcash=find(A<=0);
    
    NetEq = liquidity-deb;
    Y_prevp=Y_prev(NetEq>0);  
    bankruptcy_rate(t) = (length(negcash_k) + length(negcash))/(F+N);

    %update bankrupted firms!
    for i=negcash                                                       %pick sequentially failed firms
        counter=counter+1;
        bankrupt(1,t,i)=1;
        
        defaults(t)=defaults(t)+1;
        zzz=deb(i);                                                     %take residual debts
        if zzz>0
           
            piB=piB+(liquidity(i)-deb(i));                                               %account for bad debts
            loans=loans-zzz;
        end
        A(i)=PA(W+i)+K(i)*price_k(t+1);                                                   %initialize new firm 
        capital_value(i)=K(i)*price_k(t+1);
        PA(W+i)=0;
        liquidity(i)=A(i)-K(i)*price_k(t+1);
        deb(i)=0;
        P(i)=mean(P);
        targetLev=0.2;
        mxY=((A(i)+targetLev*A(i)/(1-targetLev))*1/wb)*alpha;
        Y_prev(i)=min(trimmean(Y_prevp,10),mxY);
        Yd(i)=Y_prev(i);
        x(i)=Y_prev(i)/k/K(i);
        barK(i)=K(i);
        barYK(i)=Y_prev(i)/k;
        Y(i)=0;
        stock(i)=0;
        interest_r(i)=r_f;
        %%fire workers
        workforce=find(Oc==i);                                          %pick all firm i's workers
        fired=workforce;  
        Oc(fired)=0;
        w(fired)=0;
        Leff(i)=0; 
    end
    
    NetEq_k = A_k;
    if isempty(find(NetEq_k>0, 1))
          warning('all capital firms are bankrupted')
%          display(seed)
        break
        %%keep the same variable of last year        
    else
        Y_prevp_k=Y_prev_k(NetEq_k>0);  
    end
    
    initialA_k(t)=sum(PA(W+F+negcash_k))/length(negcash_k);
    for i=negcash_k                                                       %pick sequentially failed firms
        defaults_k(t)=defaults_k(t)+1;
        counter=counter+1;
        bankrupt(1,t,F+i)=1;
        zzz=deb_k(i);                                                     %take residual debts
        if zzz>0
           
            piB=piB+(liquidity_k(i)-deb_k(i));                                               %account for bad debts
            loans=loans-zzz;
        end
        A_k(i)=PA(W+F+i);                                                   %initialize new firm 
        PA(W+F+i)=0;
        liquidity_k(i)=A_k(i);
        deb_k(i)=0;
        P_k(i)=mean(P_k);
        %maximum initial productin is given by the leverage
        targetLev=0.2;
        mxY=((A_k(i)+targetLev*A_k(i)/(1-targetLev))*1/wb)*alpha;
        
        Y_prev_k(i)=min(trimmean(Y_prevp_k,10),mxY);
        Y_kd(i)=Y_prev_k(i);
        Y_k(i)=0;
        stock_k(i)=0;
        interest_r_k(i)=r_f;
        workforce_k=find(Oc==F+i);                 %pick all firm i's workers
        fired_k=workforce_k;  
        Oc(fired_k)=0;
        w(fired_k)=0;
        Leff_k(i)=0; 
    end
    
    
    %% bank accounting
    piB=piB+sum(interests)+sum(interests_k) + bond_interest_rate*bonds(t-1);                                             %bank profits  

    
    if piB>0 && totE(t)>0
        Div_probb(t)=1;
        dividendsB(t)=div_B*piB;
        piB=(1-div_B)*piB;
        PA(W+1:W+F+N)=PA(W+1:W+F+N)+dividendsB(t)/(F+N);
    else 
         dividendsB(t) = 0;
         Div_probb(t)=0;
    end    
    E=E+piB;                                                            %update bank capital

    %%add bank's dividends to income of capitalists
    Kdivs(:,t,:)=dividends_income_k;
    Cdivs(:,t,:)=dividends_income;
    div_count(t)=sum(dividends_income>0);
    divk_count(t)=sum(dividends_income_k>0);
    dividends_income_k = dividends_income_k + dividendsB(t)/(F+N);
    dividends_income = dividends_income + dividendsB(t)/(F+N);
    Kdivs_p(:,t,:)=dividends_income_k;
    Cdivs_p(:,t,:)=dividends_income;
    divstatus=dividends_income;
    divstatus(divstatus>0)=1;
    divstatus_k=dividends_income_k;
    divstatus_k(divstatus_k>0)=1;
    
    profitsB(t)=piB;

    totE(t+1)=E;
  
  
 GB(t)= TA(t) - G(t)- EXP(t)- bond_interest_rate*bonds(t-1); %JAKOB bonds(t-1) mi impalla tutto ma non potevo
 primary_GB(t) = TA(t) - G(t)- EXP(t);

 stock_bonds(t) =sum(-GB(1:t));
 bonds(t) = max(0,stock_bonds(t));
 bonds_real(t) =stock_bonds(t)/price(t);
 
 quota_exp_c(:,t)= (RIC./sum(RIC));  % quota di mercato calcolata come share entrate
 shares(:,t,:) = quota_exp_c(:,t);
 
 deficitPil(t)= -GB(t)/Y_nominal_tot(t); %% ES POLITICA FISCALE
 deficit_real(t)=-GB(t)/price(t);
 primary_deficit_pil(t) = -primary_GB(t)/Y_nominal_tot(t);

money(t) = sum(PA)+sum(liquidity)+sum(liquidity_k)-sum(deb)-sum(deb_k)+E;
net_money(t) = money(t) - stock_bonds(t);

Assets(t) = sum(PA)/price(t);

if u_target - Un(t)>0
    wb = wb * (1+ wage_update_up * (u_target - Un(t))); 
else
    wb = wb * (1+ wage_update_down * (u_target - Un(t))); 
end

if YPtest==1
    for i=1:F
       if divstatus(i)==1 && divstatus_prev(i)==0
           trans_ND(t)=trans_ND(t)+1;
       end
       if divstatus(i)==1 && divstatus_prev(i)==1
           trans_DD(t)=trans_DD(t)+1;
       end
       if divstatus(i)==0 && divstatus_prev(i)==0
           trans_NN(t)=trans_NN(t)+1;
       end
       if divstatus(i)==0 && divstatus_prev(i)==1
           trans_DN(t)=trans_DN(t)+1;
       end
    end
    
    
    if t<burnin
        if sum(trans_DD(1:t))+sum(trans_DN(1:t))>0
        prob_DD(t)=sum(trans_DD(1:t))/(sum(trans_DD(1:t))+sum(trans_DN(1:t)));
        prob_DN(t)=sum(trans_DN(1:t))/(sum(trans_DD(1:t))+sum(trans_DN(1:t)));
        else
        prob_DD(t)=prob_DD(t-1);
        prob_DN(t)=prob_DN(t-1);
        end
    else
        if sum(trans_DD((t-(horizon-1)):t))+sum(trans_DN((t-(horizon-1)):t))>0
        prob_DD(t)=sum(trans_DD((t-(horizon-1)):t))/(sum(trans_DD((t-(horizon-1)):t))+sum(trans_DN((t-(horizon-1)):t)));
        prob_DN(t)=sum(trans_DN((t-(horizon-1)):t))/(sum(trans_DD((t-(horizon-1)):t))+sum(trans_DN((t-(horizon-1)):t)));
        else
        prob_DD(t)=prob_DD(t-1);
        prob_DN(t)=prob_DN(t-1);
        end
    end

    if t < burnin
        if sum(trans_ND(1:t))+sum(trans_NN(1:t))>0
        prob_ND(t)=sum(trans_ND(1:t))/(sum(trans_ND(1:t))+sum(trans_NN(1:t)));
        prob_NN(t)=sum(trans_NN(1:t))/(sum(trans_ND(1:t))+sum(trans_NN(1:t)));
        else
        prob_ND(t)=prob_ND(t-1);
        prob_NN(t)=prob_NN(t-1);
        end
    else
        if sum(trans_ND((t-(horizon-1)):t))+sum(trans_NN((t-(horizon-1)):t))>0      
        prob_ND(t)=sum(trans_ND((t-(horizon-1)):t))/(sum(trans_ND((t-(horizon-1)):t))+sum(trans_NN((t-(horizon-1)):t)));
        prob_NN(t)=sum(trans_NN((t-(horizon-1)):t))/(sum(trans_ND((t-(horizon-1)):t))+sum(trans_NN((t-(horizon-1)):t)));
        else
        prob_ND(t)=prob_ND(t-1);
        prob_NN(t)=prob_NN(t-1);
        end
    end
    
    rdivs(:,t)=dividends_income./price(t);
    
    if t < burnin
        for i=1:F
            if sum(rdivs(i,1:t)>0)>0
                rd(i)=sum(rdivs(i,1:t))/sum(rdivs(i,1:t)>0);
            else
                rd(i)=0;
            end
        end
    else
        for i=1:F
            rd(i)=sum(rdivs(i,(t-(horizon-1)):t))/sum(rdivs(i,(t-(horizon-1)):t)>0);
        end
            
    end
    
    
        mat=[prob_DD(t) prob_DN(t); prob_ND(t) prob_NN(t)];
        mats=zeros(2,2,periods);
        mats(:,:,1)=mat;
        for i=2:periods
            prod=mat^i;
            mats(:,:,i)=prod;
        end
        for i=1:F
        pay=[rd(i);0];
        permd=pay;
        for j=1:periods
            permd=permd+mats(:,:,j)*pay;
        end
            V_div=permd;
            permddiff=permd(1)-permd(2);
            if learn_act==1 && t>2500 && pub_exp_cr(t+5)>0
            V_div(:)=V_div(:)+beta^5*add_d(t)*permddiff;
            end
            if learn_act==1 && t>2500 && pub_exp_cr(t+4)>0
            V_div(:)=V_div(:)+beta^4*add_d(t-1)*permddiff;
            end
            if learn_act==1 && t>2500 && pub_exp_cr(t+3)>0
            V_div(:)=V_div(:)+beta^3*add_d(t-2)*permddiff;
            end
            if learn_act==1 && t>2500 && pub_exp_cr(t+2)>0
            V_div(:)=V_div(:)+beta^2*add_d(t-3)*permddiff;
            end
            if learn_act==1 && t>2500 && pub_exp_cr(t+1)>0
            V_div(:)=V_div(:)+beta*add_d(t-4)*permddiff;
            end
            if learn_act==1 && t>2500 && pub_exp_cr(t)>0
            V_div(:)=V_div(:)-add_d(t-5)*permddiff;
            end
            if dividends_income(i)>0
                yd(i)=(PA(W+i)/price(t)+V_div(1))/betasum;
            else
                yd(i)=(PA(W+i)/price(t)+V_div(2))/betasum;
            end
            YP_d(t)=YP_d(t)+V_div(1);
            YP_nd(t)=YP_nd(t)+V_div(2);
        end
        YP_d(t)=YP_d(t)/F;
        YP_nd(t)=YP_nd(t)/F;

    for l=1:N
       if divstatus_k(l)==1 && divstatus_prev_k(l)==0
           trans_ND_k(t)=trans_ND_k(t)+1;
       end
       if divstatus_k(l)==1 && divstatus_prev_k(l)==1
           trans_DD_k(t)=trans_DD_k(t)+1;
       end
       if divstatus_k(l)==0 && divstatus_prev_k(l)==0
           trans_NN_k(t)=trans_NN_k(t)+1;
       end
       if divstatus_k(l)==0 && divstatus_prev_k(l)==1
           trans_DN_k(t)=trans_DN_k(t)+1;
       end
    end

    if t<=horizon
        if sum(trans_DD_k(1:t))+sum(trans_DN_k(1:t))>0
        prob_DD_k(t)=sum(trans_DD_k(1:t))/(sum(trans_DD_k(1:t))+sum(trans_DN_k(1:t)));
        prob_DN_k(t)=sum(trans_DN_k(1:t))/(sum(trans_DD_k(1:t))+sum(trans_DN_k(1:t)));
        else
        prob_DD_k(t)=prob_DD_k(t-1);
        prob_DN_k(t)=prob_DN_k(t-1);
        end
    else
        if sum(trans_DD_k((t-(horizon-1)):t))+sum(trans_DN_k((t-(horizon-1)):t))>0
        prob_DD_k(t)=sum(trans_DD_k((t-(horizon-1)):t))/(sum(trans_DD_k((t-(horizon-1)):t))+sum(trans_DN_k((t-(horizon-1)):t)));
        prob_DN_k(t)=sum(trans_DN_k((t-(horizon-1)):t))/(sum(trans_DD_k((t-(horizon-1)):t))+sum(trans_DN_k((t-(horizon-1)):t)));
        else
        prob_DD_k(t)=prob_DD_k(t-1);
        prob_DN_k(t)=prob_DN_k(t-1);
        end
    end

    if t<burnin
        if sum(trans_ND_k(1:t))+sum(trans_NN_k(1:t))>0
        prob_ND_k(t)=sum(trans_ND_k(1:t))/(sum(trans_ND_k(1:t))+sum(trans_NN_k(1:t)));
        prob_NN_k(t)=sum(trans_NN_k(1:t))/(sum(trans_ND_k(1:t))+sum(trans_NN_k(1:t)));
        else
        prob_ND_k(t)=prob_ND_k(t-1);
        prob_NN_k(t)=prob_NN_k(t-1);
        end
    else
        if sum(trans_ND_k((t-(horizon-1)):t))+sum(trans_NN_k((t-(horizon-1)):t))>0
        prob_ND_k(t)=sum(trans_ND_k((t-(horizon-1)):t))/(sum(trans_ND_k((t-(horizon-1)):t))+sum(trans_NN_k((t-(horizon-1)):t)));
        prob_NN_k(t)=sum(trans_NN_k((t-(horizon-1)):t))/(sum(trans_ND_k((t-(horizon-1)):t))+sum(trans_NN_k((t-(horizon-1)):t)));
        else
        prob_ND_k(t)=prob_ND_k(t-1);
        prob_NN_k(t)=prob_NN_k(t-1);
        end
    end
    
    
    rdivs_k(:,t)=dividends_income_k./price(t);
    
      
    if t < burnin
        for l=1:N
            if sum(rdivs_k(l,1:t)>0)>0
                rdk(l)=sum(rdivs_k(l,1:t))/sum(rdivs_k(l,1:t)>0);
            else
                rdk(l)=0;
            end
        end
    else
       for l=1:N
            rdk(l)=sum(rdivs_k(l,(t-(horizon-1)):t))/sum(rdivs_k(l,(t-(horizon-1)):t)>0);
       end 
        
    end
    
        mat=[prob_DD_k(t) prob_DN_k(t); prob_ND_k(t) prob_NN_k(t)];
        mats=zeros(2,2,periods);
        mats(:,:,1)=mat;
        for i=2:periods
            prod=mat^i;
            mats(:,:,i)=prod;
        end
        for l=1:N
        pay=[rdk(l);0];
        permdk=pay;
        for j=1:periods
            permdk=permdk+mats(:,:,j)*pay;
        end
            V_divk=permdk;
            permdkdiff=permdk(1)-permdk(2);
            if learn_act==1 && t>2500 && pub_exp_cr(t+5)>0
            V_divk(:)=V_divk(:)+beta^5*add_dk(t)*permdkdiff;
            end
            if learn_act==1 && t>2500 && pub_exp_cr(t+4)>0
            V_divk(:)=V_divk(:)+beta^4*add_dk(t-1)*permdkdiff;
            end
            if learn_act==1 && t>2500 && pub_exp_cr(t+3)>0
            V_divk(:)=V_divk(:)+beta^3*add_dk(t-2)*permdkdiff;
            end
            if learn_act==1 && t>2500 && pub_exp_cr(t+2)>0
            V_divk(:)=V_divk(:)+beta^2*add_dk(t-3)*permdkdiff;
            end
            if learn_act==1 && t>2500 && pub_exp_cr(t+1)>0
            V_divk(:)=V_divk(:)+beta*add_dk(t-4)*permdkdiff;
            end
            if learn_act==1 && t>2500 && pub_exp_cr(t)>0
            V_divk(:)=V_divk(:)-add_d(t-5)*permdkdiff;
            end
            if dividends_income_k(l)>0
                ydk(l)=(PA(W+F+l)/price(t)+V_divk(1))/betasum;
            else
                ydk(l)=(PA(W+F+l)/price(t)+V_divk(2))/betasum;
            end
            YP_dk(t)=YP_dk(t)+V_divk(1);
            YP_ndk(t)=YP_ndk(t)+V_divk(2);
        end
        YP_dk(t)=YP_dk(t)/N;
        YP_ndk(t)=YP_ndk(t)/N;
    
end

rwage(t)=(wb*(1-tax_rate))/price(t);
rdiv(t)=sum(rdivs(:,t));
rdivk(t)=sum(rdivs_k(:,t));
wincome(t)=sum(workers_income)/price(t);
dincome(t)=sum(dividends_income)/price(t);
dincomek(t)=sum(dividends_income_k)/price(t);
emp_count(t)=sum(Oc>0);

if learn==1
     if t>=700 && pub_exp_cr(t)>0
       
        shocks=[];
        for i=(t-horizon):t
        if pub_exp_cr(i)>0
        shocks=[shocks i-1 i];
        end
        end
        gov=pub_exp_cr(shocks);
        cemp=emp_count(shocks);
        cdiv=div_count(shocks);
        cdivk=divk_count(shocks);
        cons=consumption(shocks);
        invs=rwage(shocks);

        regsw1=transpose(cemp(1,:));
        lyw1=lagmatrix(regsw1,1);
        regsw2=transpose(gov(1,:));
        lyw2=lagmatrix(regsw2,0);
        regsw3=transpose(cons(1,:));
        lyw3=lagmatrix(regsw3,0);
        lyw=[lyw2 lyw3 lyw1];
        regyw=fitlm(lyw,cemp(1,:),'Intercept',false);
        
        regsd1=transpose(cdiv(1,:));
        lyd1=lagmatrix(regsd1,1);
        regsd2=transpose(gov(1,:));
        lyd2=lagmatrix(regsd2,0);
        regsd3=transpose(invs(1,:));
        lyd3=lagmatrix(regsd3,0);
        lyd=[lyd2 lyd3 lyd1];
        regyd=fitlm(lyd,cdiv(1,:),'Intercept',true);
        
        regsdk1=transpose(cdivk(1,:));
        lydk1=lagmatrix(regsdk1,1);
        regsdk2=transpose(gov(1,:));
        lydk2=lagmatrix(regsdk2,0);
        lydk=[lydk2 lydk1];
        regydk=fitlm(lydk,cdivk(1,:),'Intercept',true);

        
        step=step+1;
        coeffs(1,t)=regyw.Coefficients.Estimate(1);
        coeffs(2,t)=regyd.Coefficients.Estimate(2);
        coeffs(3,t)=regydk.Coefficients.Estimate(2);
     else

        coeffs(1,t)=coeffs(1,t-1);
        coeffs(2,t)=coeffs(2,t-1);
        coeffs(3,t)=coeffs(3,t-1);
     end
end


wages_t(t) = wb;

end
gperiods=find(pub_exp_cr>0);


et = toc;
