clc
clear


k=0;
Basecase_index = 1;
while(1)

[PdQd,p,q,p_attacked,q_attacked,success] = Attack_generation(Basecase_index);
if success~=1 
    continue;
else
    k = k + 1
    Basecase_index = Basecase_index+1;
    if Basecase_index>100
        Basecase_index = 1;
    end    
    Dataset([2*k-1,2*k],:) = [PdQd',p',q',0;    
                    PdQd',p_attacked',q_attacked',1];             
end
if k>=500000
    break;
end
end
