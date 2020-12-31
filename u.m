for i = 1:num_mode
    fprintf('\n')
    for j = 1:length(wi)
            if ~any(strcmp(cs{i,j},'Infeasible')) || ~any(strcmp(cs{i,j},'Failed'))
                Q_bd(i,j) = -0.5 * (real(w(i,j)) / imag(wi(j)));
                fprintf('Q_bd = %f\n',Q_bd(i,j))
            else
                Q_bd(i,j) = -0.5 * (real(w(i,j)) / imag(omega_vals(i)));
                Q_bd(i,j) = Q_ana(i);
                fprintf('Q_bd = Q_ana = %f\n',Q_bd(i,j))
            end
    end
end
Dmax = zeros(num_mode,step);

for i = 1:num_mode
    for j = 2:step
        for k = 1:ND
            if strcmp('Infeasible', cs{i,j}{k}) || strcmp('Failed', cs{i,j}{k})
%             if contains('Infeasible',cs{i,j}{k}) || contains('Failed',cs{i,j}{k})
                Dmax(i,j) = 0;
            elseif strcmp('Solved',cs{i,j}{k})
%             elseif contains('Solved',cs{i,j}{k})
                Dmax(i,j) = 1;
            end
        end
    end
end
% Dmax = zeros(num_mode,step);
% for i = 1:num_mode
%     for j = 1:step
%         for k = 2:ND
%             if isempty(cs{i,j}{k}) == 0
%                 Dmax(i,j) = 0;
%             else
%                 Dmax(i,j) = 1;
%             end
%         end
%     end
% end
Q = zeros(1,num_mode);
Q = Q_ana;
for i = 1:num_mode
    for j = 2:step
            if Dmax(i,j) == 1
                Q(i) = Q_bd(i,j);
                break
            else
                continue
            end
    end
end

for i = 1:num_mode
    if Q(i) < Q_ana(i)
        Q(i) = Q_ana(i);
    end
end

% Dmax = zeros(num_mode,step);
% Q = zeros(1,num_mode);
% for i = 1:num_mode
%     for j = 1:step
%         for k = 3:ND
%             if ~strcmp(cs{i,j}{k},'Infeasible') || ~strcmp(cs{i,j}{k},'Failed')
%                 cs{i,j}{k} = [];
%             end
%         end
%     end
% end
% for i = 1:num_mode
%     for j = 1:step
%         if cellfun('isempty',cs{i,j})
%             Q(i) = max(Q_bd(Q_bd(i,j)~=-Inf));
%         end
%     end
% end



% for i = 1:num_mode
%     fprintf('\n')
%     for j = 1:length(wi)
% %         for k = 1:ND
%             if ~any(strcmp(cs{i,j},'Infeasible')) || ~any(strcmp(cs{i,j},'Failed'))
%                 fprintf('found bd')
%                 Q_bd(i) = -0.5 * (wr(i) / imag(wi(j)));
%                 break
%             else
% %                 Q_bd(i) = - 0.5 * (wr(i) / imag(omega_vals(i)));
%                 fprintf('no bound found qbd = qana')
%                 Q_bd(i) = Q_ana(i);
%                 continue
%             end
% %         end
%     end
% end

% for i = 1:num_mode
%     fprintf('\n')
%     Q_bd(i) = bound(i,cs(i,:),wi,wr,sort(omega_vals),ND);
% end


% for i = 1:num_mode
%     fprintf('\n')
%     for j = 1:length(wi)
%             if ~any(strcmp(cs{i,j},'Infeasible')) || ~any(strcmp(cs{i,j},'Failed'))
%                 Q_bd(i,j) = -0.5 * (real(w(i,j)) / imag(wi(j)));
%                 fprintf('Q_bd = %f\n',Q_bd(i,j))
%             else
%                 Q_bd(i,j) = -0.5 * (real(w(i,j)) / imag(omega_vals(i)));
%                 Q_bd(i,j) = Q_ana(i);
%                 fprintf('Q_bd = Q_ana = %f\n',Q_bd(i,j))
%             end
%     end
% end
  
% Dmax = zeros(num_mode,step);
% Q = zeros(1,num_mode);
% for i = 1:num_mode
%     for j = 1:step
%         for k = 3:ND
%             if cs{i,j}{k} == 'Infeasible' || cs{i,j}{k} == 'Failed'
%                 Dmax(i,j) = 1;
%             else
%                 Dmax(i,j) = 0;
%             end
%         end
%     Q(i) = Q_bd(Dmax(i,j) == 1);
%     end
% end

%plot results
figure
plot(sort(real(omega_vals)),Q_ana)
hold on
plot(sort(real(omega_vals)),Q)
xlabel('Real(Omega)')
ylabel('Qfactor')