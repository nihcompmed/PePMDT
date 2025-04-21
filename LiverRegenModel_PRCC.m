function dy = LiverRegenModel_PRCC(t, y, params)
    % Parameters
    M = params(1); ktnf = params(2); kaptnf = params(3); vjak = params(4);
    kmjak = params(5); kapjak = params(6); kprost3 = params(7); vst3 = params(8);
    kmst3 = params(9); kapst3 = params(10); vsoc = params(11); kmsoc = params(12);
    kisoc = params(13); kapsoc = params(14); vie = params(15); kmie = params(16);
    kapie = params(17); kup = params(18); kgf = params(19); kapgf = params(20);
    kapecm = params(21); kdeg = params(22); kecm = params(23); kq = params(24);
    kr = params(25); kp = params(26); kprol = params(27); kreq = params(28);
    thetar = params(29); betar = params(30); kapop = params(31); thetaa = params(32);
    betaa = params(33);
    
    % Differential equations for the model
    dy = zeros(10,1);
    
    % Example equations, based on your provided model equations
    dy(1) = ktnf*M/((y(8)+y(9)+y(10))) - vjak*y(1)/(kmjak+y(1)) - kaptnf*y(1) + (vjak/(kmjak+1) + kaptnf - ktnf*M); % TNF
    dy(2) = vjak*y(1)/(kmjak+y(1)) - kapjak*y(2) + (-vjak/(kmjak+1) + kapjak); % JAK
    dy(3) = vst3*(kprost3^2)*y(2)/(kprost3^2 + kmst3*(1 + y(4)/kisoc)) - vie*y(3)/(y(3)+kmie) - vsoc*y(3)/(y(3)+kmsoc)- kapst3*y(3) + (- vst3*(kprost3^2)/(kprost3^2 + kmst3*(1 + 1/kisoc)) + vie/(1+kmie) + vsoc/(1+kmsoc) + kapst3); % STAT3
    dy(4) = vsoc*y(3)/(y(3)+kmsoc) - kapsoc*y(4) + (-vsoc/(1+kmsoc) + kapsoc); % SOCS3
    dy(5) = -kdeg*y(1)*y(5) - kapecm*y(5) + kapecm; % ECM
    dy(6) = vie*y(3)/(y(3)+kmie) - kapie*y(6) + (- vie/(1 + kmie) + kapie); % IE
    dy(7) = kgf*M/((y(8)+y(9)+y(10))) - kup*y(7)*y(5) - kapgf*y(7) + (kup + kapgf - kgf*M); % GF
    dy(8) = -kq*(y(6)-1)*y(8) + kr*y(5)*y(10) + 0.5*(1 + tanh((thetar-y(7))/betar))*kreq*y(9) - kapop*0.5*(1 + tanh((thetaa - (y(8)+y(9)+y(10)))/betaa))*y(8); % Q
    dy(9) = -kp*(y(7)-1)*y(9) + kq*(y(6)-1)*y(8) - 0.5*(1 + tanh((thetar-y(7))/betar))*kreq*y(9) - kapop*0.5*(1 + tanh((thetaa - (y(8)+y(9)+y(10)))/betaa))*y(9); % P
    dy(10) = -kr*y(5)*y(10) + kp*(y(7)-1)*y(9) + kprol*y(10) - kapop*0.5*(1 + tanh((thetaa - (y(8)+y(9)+y(10)))/betaa))*y(10); % P
end
