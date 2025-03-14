function L_pii_l_round_wire = L_pii_round_wire(mu0, omega,u, l)
    % Calculate the partial self-inductance per unit length with the approximation of 
    % a round wire. Approximation valid for u > 10 and omega = 1.

    % Inputs
    %       mu0: Permeability of free space (H/m)
    %       omega: Ratio of thickness to width 
    %       u: Ration of length to width 
    %       l: length of the conductor (m)
    
    % Outputs
    %       L_pii_l: partial self-inductance per unit lengt with the approximation of 
    %                a round wire (H/m)
    
    tolerance = 1e-1;
    isCloseToOne = abs(omega - 1) < tolerance;

    if u < 10 || ~isCloseToOne
        warning('Check u and omega. HP: u>10 and omega = 1')
    end

    L_pii_l_round_wire = (mu0 / (2*pi)) * (log(4*u) - 1);
    L_pii_l_round_wire = L_pii_l_round_wire / l;
end