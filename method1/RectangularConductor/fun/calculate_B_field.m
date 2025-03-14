function B = calculate_B_field(x_line, y_line, xi, i_lambda)
    mu0 = 4 * pi * 1e-7;  
    B = zeros(length(x_line), 3);  
    
    for k = 1:length(x_line)
        B_sum = [0, 0, 0]; 
        
        for j = 1:size(xi, 1)
            
            r_vec = [x_line(k) - xi(j, 1), y_line(k) - xi(j, 2), 0];
            distance = norm(r_vec);  
            
            % Avoid singularities
            if distance == 0
                continue;
            end
            
            I_vec = [0, 0, i_lambda(j)];
            
            dB = (mu0 / (2 * pi)) * cross(I_vec, r_vec) / distance^2;
            B_sum = B_sum + dB;  % Sum the contributions
        end
        
        B(k, :) = B_sum;  
    end
end

