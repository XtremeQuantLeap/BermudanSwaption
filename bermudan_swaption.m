function []=bermudan_swaption()

%This program can work for arbitrary no. of factors. You have to specify no. of factors 
%as well as volatility structure for each factor. The volatility structure can be 
%obtained from principal component analysis of correlation matrix and adjusting to 
%calibrated volatilities as done in excellent paper by Rebonato. See my web page for the 
%references(http://www.geocities.com/anan2999). It does not take correlation structure 
%as input. You can also specify CEV constant alpha for skew. Remember changing this constant
%changes effective volatility.

%randn('state',[1541045451;4027226640]) % add a good random number seed here if you wish.
%if you don't matlab will choose its own seed.

delta=.25; %Tenor spacing. usually .25 or .5

P=5000;  % No. of paths, do not try more than 5000 paths unless you are very patient
T_e1=6.0;		%maturity of underlying swap in years(must be an exact multiple of delta)
T_x1=5.75;		%last exercise date of the swaption (must be an exact multiple of delta)
T_s1=3.0;		%lockout date (must be an exact multiple of delta)


T_e=T_e1/delta+1;
T_x=T_x1/delta+1;
T_s=T_s1/delta+1;
N=T_e;

F=2;   % number of factors. If you change this line also change volatility structure appropriately
alpha=1.0;%CEV constant alpha for skew.Remember changing this value changes effective volatility
%It is 1.0 for lognormal model.
k=.1; % strike, fixed coupon
pr_flag=+1; %payer receiver flag; assumes value of +1 for a payer swaption 
%and a value of -1 for a receiver swaption.

           
n_spot=2;
L=repmat(.10,[P,T_e+1]);
vol=zeros([T_e,F]);
for n=1:N
   for f=1:F
      if(f==1)
         vol(n,f)=.15; %volatility of first factor
      end
      if(f==2)
      	 vol(n,f)= (.15-(.009*(n)*.25).^.5); %volatility of second factor
      end
   end
end
%You can add more vaolatility factors in the above line but please also change F accordingly
%drift=repmat(0,[P,F]);
money_market=ones([T_x,P]);
swap=zeros([T_x,P]);
B=ones([P,T_e]);


money_market(2,:)=money_market(1,:).*(1+delta*L(:,1))';
increment=zeros([P,1]);
drift=zeros([P,F]);

for t= 2 : T_x
   
   
  normal_matrix=randn([P,F]);
  drift(:,:)=0;
   for n= t : T_e
      increment(:,1)=0;
     
		%      n
      for f=1:F
         
         drift(:,f)=drift(:,f)+ delta*vol(n-n_spot+1,f).*((L(:,n).^alpha)./(1+delta.*L(:,n))); %  
         
         increment(:,1)=increment(:,1)+vol(n-n_spot+1,f).*(L(:,n).^alpha)./L(:,n)...
            .*(normal_matrix(:,f).*sqrt(delta)-.5.*vol(n-n_spot+1,f).*(L(:,n).^alpha)./L(:,n)...
            .*delta+drift(:,f).*delta);
      end
      
      L(:,n)=L(:,n).*exp(increment(:,1)); 
      L(L(:,n)<.00001,n)=.00001;
   end
   
   B(:,t)=1.0;
   for n=t+1:T_e
      B(:,n)=B(:,n-1)./(1+delta.*L(:,n-1));
   end
   
  
     
   money_market(t+1,:)=money_market(t,:).*(1+delta*L(:,n_spot))';
   
   if((t>= T_s) && (t <=T_x))
      for n=t:(T_e-1) %//the swap leg is determined one date before the end
         
                  swap(t,:)=swap(t,:)+  (B(:,n+1).*(L(:,n)-k).*pr_flag*delta)' ;
      end
   end
   n_spot=n_spot+1;
   
end



value=zeros([P,1]);
stop_rule=repmat(T_x,[P,1]);

value(swap(T_x,:)>0,1) = (swap(T_x,swap(T_x,:)>0))';
coeff=zeros([T_x,6]);

for t=(T_x-1):-1:T_s	
   i=0;
   a=0;
   y=0;
   for p=1:P
      if (swap(t,p)> 0.0)
         i=i+1;
         a(i,1)=1;
         a(i,2)=swap(t,p);
         a(i,3)=swap(t,p)*swap(t,p);
         a(i,4)=money_market(t,p);
         a(i,5)=money_market(t,p)*money_market(t,p);
         a(i,6)=money_market(t,p)*swap(t,p);
                   
         y(i,1)= money_market(t,p)/money_market(stop_rule(p,1),p) * value(p,1);
            
      end
      
   end
   
   
      
   temp=inv(a'*a)\(a'*y);
   coeff(t,:)=temp';
        
   expec_cont_value=zeros([P,1]);
   exer_value=zeros([P,1]);
   
   expec_cont_value(:,1)=(coeff(t,1)+coeff(t,2).*swap(t,:)+coeff(t,3).*swap(t,:)...
      .*swap(t,:)+coeff(t,4).*money_market(t,:)+coeff(t,5).*money_market(t,:)...         
      .*money_market(t,:)+coeff(t,6).*money_market(t,:).*swap(t,:))';

      
exer_value(swap(t,:)>0,1)=(swap(t,swap(t,:)>0))';
    
            
 value((exer_value(:,1)>expec_cont_value(:,1))&(swap(t,:)>0)',1)...
    =exer_value((exer_value(:,1)> expec_cont_value(:,1))&(swap(t,:)>0)',1);
 
 stop_rule((exer_value(:,1)>expec_cont_value(:,1))&(swap(t,:)>0)',1)=t;

                  
end
  
   price=0;
   for p=1:P
      price=price+ (value(p,1)/(money_market(stop_rule(p,1),p)))/P;
            
   end         
   
   
   