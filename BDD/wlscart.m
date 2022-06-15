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
           H31 = zeros(npi,nbus);  %Derivative of real power injection wrt e     H32 = zeros(npi,nbus);  %Derivative of real power injection wrt f     H41 = zeros(npi,nbus);  %Derivative of reactive power injection wrt e     H42 = zeros(npi,nbus);  %Derivative of reactive power injection wrt f     for i = 1:npi         m = fbus(ppi(i));         for k = 1:(nbus)             if k == m                 for n = 1:nbus                     H31(i,k) = H31(i,k) + (G(m,n)*e(n) - B(m,n)*f(n));                     H32(i,k) = H32(i,k) + (G(m,n)*f(n) + B(m,n)*e(n));                     H41(i,k) = H41(i,k) -G(m,n)*f(n) - B(m,n)*e(n);                     
H42(i,k) = H42(i,k) + (G(m,n)*e(n) - B(m,n)*f(n)); 
         
         










