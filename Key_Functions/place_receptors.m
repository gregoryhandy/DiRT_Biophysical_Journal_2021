function [rec_center] = place_receptors(cylinder_radius,M,rec_radius, plot_receptors)

max_failures = 100000;

% cylinder_radius = 2*half_R;
rec_center = zeros(M,2);
count = 0;
for i = 1:M
    bad = 1;
    % prevents overlapping receptors
    while bad == 1 && count < max_failures
        bad = 0;
        theta = 2*pi*rand();
        % receptors are completely contains in the PSD
        r = (cylinder_radius-rec_radius)*sqrt(rand());
        rec_center(i,:) = [r*cos(theta) r*sin(theta)];
  
        for j = 1:(i-1)
            if norm(rec_center(i,:) - rec_center(j,:))<2*rec_radius
                bad = 1;
                count = count +1;
                break;
            end 
        end 

    end
end

if count == max_failures
   error_text = ['Unable to prevent receptors from overlapping',newline...
       'Consider making receptors smaller or having less of them'];
   error(error_text); 
end

%% Plot the receptors (if requested)
if plot_receptors == 1
    % domain boarder
    x = [-cylinder_radius:.001:cylinder_radius];
    plot(x,sqrt(cylinder_radius^2-x.^2),'--','linewidth',1.5,'color','k')
    hold on;
    plot(x,-sqrt(cylinder_radius^2-x.^2),'--','linewidth',1.5,'color','k')
    
    % rec boadrers
    x = [-1:.00001:1];
    for i = 1:M
        temp = real(sqrt(rec_radius^2-(x-rec_center(i,1)).^2));
        real_indices = find(temp ~= 0);
    
        plot(x(real_indices),rec_center(i,2)+sqrt(rec_radius^2-(x(real_indices)-rec_center(i,1)).^2),'linewidth',1.5,'color','r')
        plot(x(real_indices),rec_center(i,2)-sqrt(rec_radius^2-(x(real_indices)-rec_center(i,1)).^2),'linewidth',1.5,'color','r')
    end
end

end

