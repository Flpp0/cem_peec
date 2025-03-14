function L_pii_l = L_pii(mu0, omega, u, l)
    % Calculate the partial self-inductance per unit length.

    % Inputs
    %       mu0: Permeability of free space (H/m)
    %       omega: Ratio of thickness to width 
    %       u: Ration of length to width 
    %       l: length of the conductor (m)
    
    % Outputs
    %       L_pii_l: partial self-inductance per unit lengt (H/m)

    A1_val = A1(u);
    A2_val = A2(omega);
    A3_val = A3(u, omega);
    A4_val = A4(u, omega);
    A5_val = A5(A4_val, A3_val);
    A6_val = A6(omega, A4_val, A1_val);
    A7_val = A7(u, A4_val, A2_val);

    term1 = (omega^2 / (24 * u)) * (log((1 + A2_val) / omega) - A5_val);
    term2 = (1 / (24 * u * omega)) * (log(omega + A2_val) - A6_val);
    term3 = (omega^2 / (60 * u)) * (A4_val - A3_val);
    term4 = (omega^2 / 24) * (log((u + A3_val) / omega) - A7_val);
    term5 = (omega^2 / (60 * u)) * (omega - A2_val);
    term6 = (1 / (20 * u)) * (A2_val - A4_val);
    term7 = (u / 4) * A5_val - (u^2 / (6 * omega)) * atan(omega / (u * A4_val));
    term8 = u / (4 * omega) * A6_val - (omega / 6) * atan(u / (omega * A4_val));
    term9 = A7_val / 4 - (1 / (6 * omega)) * atan(u * omega / A4_val);
    term10 = (1 / (24 * omega^2)) * (log(u + A1_val) - A7_val) + (u / (20 * omega^2)) * (A1_val - A4_val);
    term11 = (1 / (60 * omega^2 * u)) * (1 - A2_val) + (1 / (60 * omega^2 * u)) * (A4_val - A1_val);
    term12 = (u / 20) * (A3_val - A4_val);
    term13 = (u^3 / (24 * omega^2)) * (log((1 + A1_val) / u) - A5_val);
    term14 = (u^3 / (24 * omega)) * (log((omega + A3_val) / u) - A6_val);
    term15 = (u^3 / (60 * omega^2)) * (A4_val - A1_val + u - A3_val);

    L_pii_l = (2 * mu0 / pi)*(term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10 + term11 + term12 + term13 + term14 + term15);
    L_pii_l = L_pii_l / l;
end