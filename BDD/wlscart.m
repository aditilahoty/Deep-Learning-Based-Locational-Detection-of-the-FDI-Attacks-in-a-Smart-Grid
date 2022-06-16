 %% %%%   DESCRIPTION: This program Calculates states of test bus system    
 %% %%%                in Cartesian coordinate systems using conventional   
 %% %%%                power injection flow measurements. This program also 
 %% %%%                performs bad data analysis using maximum normalized  
 %% %%%                residual test.
   
   % Power System State Estimation using Weighted Least Square Method.. 
   
   %nbus represents the number of bus system 
   %nbus = 14 for IEEE 14 bus system 
   
   
   ybus = ybusfunc(nbus); % Get YBus.. 
   zdata = zconv(nbus); % Get  Conventional Measurement data.. 
   [bsh g b] = line_mat_func(nbus); % Get conductance and susceptance matrix  
   type = zdata(:,2);  
   % Type of measurement, 
   % type =1 voltage magnitude p.u 
   % type =2 Voltage phase angle in degree 
   % type =3 Real power injections 
   % type =4 Reactive power injection 
   % type =5 Real power flow           
   % type =6 Reactive power flow       
   z = zdata(:,3); % Measurement values 
   fbus = zdata(:,4); % From bus 
   tbus = zdata(:,5); % To bus 
   Ri = diag(zdata(:,6)); % Measurement Error Covariance matrix 
   e = ones(nbus,1); % Initialize the real part of bus voltages 
   
   f = zeros(nbus,1);% Initialize the imaginary part of bus voltages 
   E = [f;e];  % State Vector comprising of imaginary and real part of voltage 
   G = real(ybus); 
   B = imag(ybus);
   
   ei = find(type == 1); % Index of voltage magnitude measurements.. 
   fi = find(type == 2); % Index of voltage angle measurements.. 
   ppi = find(type == 3); % Index of real power injection measurements.. 
   qi = find(type == 4); % Index of reactive power injection measurements.. 
   pf = find(type == 5); % Index of real power flow measurements.. 
   qf = find(type == 6); % Index of reactive power flow measurements.. 
   Vm=z(ei); 
   Thm=z(fi); 
   z(ei)=Vm.*cosd(Thm); % converting voltage from polar to Cartesian 
   z(fi)=Vm.*sind(Thm);
   
    nei = length(ei); % Number of Voltage measurements(real) 
    nfi = length(fi); % Number of Voltage measurements(imaginary) 
    npi = length(ppi); % Number of Real Power Injection measurements.. 
    nqi = length(qi); % Number of Reactive Power Injection measurements.. 
    npf = length(pf); % Number of Real Power Flow measurements.. 
    nqf = length(qf); % Number of Reactive Power Flow measurements..   
      
    iter = 0; tol = 1; 
    while(tol > 1e-6) 
      
       %Measurement Function, h 
       h1 = e(fbus (ei),1);   % voltage measurement
       h2 = f(fbus (fi),1);   % angle measurement
       h3 = zeros(npi,1);     % real power injection
       h4 = zeros(nqi,1);     % reactive power injection
       h5 = zeros(npf,1);     % real power flow
       h6 = zeros(nqf,1);     % reactive power flow
       
       %Measurement function of power injection 
       for i = 1:npi 
           m = fbus(ppi(i)); 
           for k = 1:nbus 
       %        Real injection 
                h3(i)=h3(i)+(G(m,k)*(e(m)*e(k)+f(m)*f(k))... +B(m,k)*(f(m)*e(k)-e(m)*f(k))); 
       %       Reactive injection        
                h4(i)=h4(i)+(G(m,k)*(f(m)*e(k)-e(m)*f(k))... -B(m,k)*(e(m)*e(k)+f(m)*f(k))); 
            end 
        end 

        %Measurement function of power flow 
        for i = 1:npf 
          m = fbus(pf(i));         
          n = tbus(pf(i));    
   %        Real injection         
           h5(i) =(e(m)^2 + f(m)^2)*g(m,n)...             
                  -(g(m,n)*(e(m)*e(n)+f(m)*f(n))+b(m,n)*(f(m)*e(n)-e(m)*f(n)));      
   %       Reactive injection                
           h6(i) =-g(m,n)*(f(m)*e(n)e(m)*f(n))+b(m,n)*(e(m)*e(n)+f(m)*f(n))...                
                 -(e(m)^2 + f(m)^2)*(b(m,n)+bsh(m,n));     
         end           
          
         h = [h1; h2; h3; h4; h5; h6]; 

         %Residual matrix difference of measurement and the non linear        
         r = z - h;           
  
         % Jacobian..          
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
         %%%%%        Jacobian  Block 1: Derivative of voltage        %%%%%     
         %%%%%            with respect to states                      %%%%%     
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
         H11 = zeros(nei,nbus); % Derivative of e wrt e     
         H12 = zeros(nei,nbus); % Derivative of e wrt f     
         H21 = zeros(nfi,nbus); % Derivative of f wrt e     
         H22 = zeros(nfi,nbus); % Derivative of f wrt f     
         for k = 1:nei          
            H11(k,fbus(k)) = 1;          
            H22(k,fbus(n)) = 1;     
         end      
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
         %%%%%    Jacobian  Block 2: Derivative of Power injection    %%%%%     
         %%%%%            with respect to states                      %%%%%     
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
           H31 = zeros(npi,nbus);  %Derivative of real power injection wrt e     
           H32 = zeros(npi,nbus);  %Derivative of real power injection wrt f     
           H41 = zeros(npi,nbus);  %Derivative of reactive power injection wrt e     
           H42 = zeros(npi,nbus);  %Derivative of reactive power injection wrt f     
           for i = 1:npi         
               m = fbus(ppi(i));         
               for k = 1:(nbus)             
                   if k == m                 
                       for n = 1:nbus                     
                           H31(i,k) = H31(i,k) + (G(m,n)*e(n) - B(m,n)*f(n));                     
                           H32(i,k) = H32(i,k) + (G(m,n)*f(n) + B(m,n)*e(n));                     
                           H41(i,k) = H41(i,k) -G(m,n)*f(n) - B(m,n)*e(n);                     
                           H42(i,k) = H42(i,k) + (G(m,n)*e(n) - B(m,n)*f(n)); 
                        end                 
                        H31(i,k) = H31(i,k) + f(m)*B(m,m) + G(m,m)*e(m);                 
                        H32(i,k) = H32(i,k) - e(m)*B(m,m) + f(m)*G(m,m);                 
                        H41(i,k) = H41(i,k) + f(m)*G(m,m) - e(m)*B(m,m);                 
                        H42(i,k) = H42(i,k) - e(m)*G(m,m) - f(m)*B(m,m);             
                     else                 
                        H31(i,k) = G(m,k)*e(m) + B(m,k)*f(m);                 
                        H32(i,k) =G(m,k)*f(m) - B(m,k)*e(m);                     
                        H41(i,k) = (G(m,k)*f(m) - B(m,k)*e(m));                 
                        H42(i,k) = (-G(m,k)*e(m) - B(m,k)*f(m));            
                      end         
                  end     
              end 
              
          
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
                %%%%%        Jacobian Block 3: Derivative of Power flow       %%%%%     
                %%%%%                 with respect to states                  %%%%%     
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
                H51 = zeros(npf,nbus);     
                H52 = zeros(npf,nbus);    
                H61 = zeros(nqf,nbus);     
                H62 = zeros(nqf,nbus);     
                for i = 1:npf         
                    m = fbus(pf(i));         
                    n = tbus(pf(i));                 
                    H51(i,m) = 2*e(m)*g(m,n) - g(m,n)*e(n) + b(m,n)*f(n);                 
                    H51(i,n) = -g(m,n)*e(m) - b(m,n)*f(m);         
                    H52(i,m) = 2*f(m)*g(m,n) - g(m,n)*f(n) - b(m,n)*e(n);         
                    H52(i,n) = -g(m,n)*f(m) + b(m,n)*e(m);                
                    H61(i,m)=-2*e(m)*(b(m,n)+bsh(m,n))+g(m,n)*f(n)+b(m,n)*e(n);         
                    H61(i,n) = -g(m,n)*f(m) + b(m,n)*e(m);                     
                    H62(i,m)=-2*f(m)*(b(m,n)+bsh(m,n))-g(m,n)*e(n)+b(m,n)*f(n);         
                    H62(i,n) = g(m,n)*e(m) + b(m,n)*f(m);     
                end 
                
                 % Measurement Jacobian, H..     
                 H = [H11 H12;          
                      H21 H22;          
                      H31 H32;          
                      H41 H42;          
                      H51 H52;          
                      H61 H62];           
                  
                  % Gain Matrix, Gm..     
                  Gm = H'*inv(Ri)*H;          
                  
                  %Objective Function..     
                  J = sum(inv(Ri)*r.^2);            
                  
                  %Solving for states iteratively by cholesky factorization and forward
                  %and back substitution. 
                  
                  gm=H'*inv(Ri)*r; 
                  dE=factoriseGbychol(Gm,gm,length(E)); 
                  E = E + dE; 
                  
                  e = E(1:nbus); 
                  f = E(nbus+1:end); 
                  iter = iter + 1; 
                  tol = max(abs(dE)); 
               end 
               
               displayout(E,'a'); % Displaying output in tabular form 
               
               %% 
               %baddata 
               
               O = Ri- H*inv(Gm)*H'; 
               od=diag(O); 
               cm = find(od<=1e-12); 
               r(cm)=0; 
               rN=abs(r)./sqrt(od); 
               displayout(rN,'b');  % Display rN in formatted form







