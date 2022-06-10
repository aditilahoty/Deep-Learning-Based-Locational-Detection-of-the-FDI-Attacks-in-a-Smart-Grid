function [PdQd,p,q,p_attacked,q_attacked,success] = Attack_generation(Basecase_index)

mpc=Case14PowerFlowData(Basecase_index);
[MVAbase, bus, gen, gencost, branch, f, success, et] = runopf(mpc);
PdQd = [bus([2:6,9:14],3);bus([2:6,9:14],4)];
  
Theta = bus(:,VA)*pi/180;   
V = bus(:,VM);      
for i=1:size(branch,1)
    p(2*i-1,1) = branch(i,PF)/MVAbase;
    p(2*i,1) = branch(i,PT)/MVAbase;
    q(2*i-1,1) = branch(i,QF)/MVAbase;
    q(2*i,1) = branch(i,QT)/MVAbase;
end
I1 = [0.1,0.2,0.3];   
I2 = [2,3,9];
attack_intensity = I1(randi(3)); 
attacked_bus = I2(randi(3));
noise_intensity = 0.01;    
p = p .* (1-noise_intensity+2*noise_intensity*rand(size(p,1),1));
q = q .* (1-noise_intensity+2*noise_intensity*rand(size(q,1),1));

index = [];
for n=1:size(branch,1)
    if (branch(n,1)==attacked_bus || branch(n,2)==attacked_bus)
        index(end+1:end+2) = [2*n-1,2*n];
    end
end

p_delta = zeros(19,1);
p_delta(index,1)=attack_intensity*p(index,1);
p_attacked = p+p_delta;
q_delta = zeros(19,1);q_delta(index,1)=attack_intensity*q(index,1);q_attacked = q+q_delta;

end

