function [permInc] = projectIncome(probs11,probs12,probs21,probs22,pay1,pay2,periods)

mats=zeros(2,2,periods);
mats(1,1,:)=probs11;
mats(1,2,:)=probs12;
mats(2,1,:)=probs21;
mats(2,2,:)=probs22;
pay=zeros(2,1,periods+1);
pay(1,1,:)=pay1;
pay(2,1,:)=pay2;
perm=pay(:,:,1);
prod=eye(2);
for i=2:periods+1
    prod=prod*mats(:,:,i-1);
    perm=perm+prod*pay(:,:,i);
end
permInc=perm;

end
