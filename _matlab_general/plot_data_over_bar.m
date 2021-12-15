function plot_data_over_bar(labels,obs_vals)
%
obs_vals_orig = obs_vals;

%%
hold on
if size(obs_vals_orig,2)==2
    cols = [0 1 1; 1 0 1];
else
    cols = parula(size(obs_vals_orig,2)+1);
end
bar_dist = linspace(-0.2,0.2,size(obs_vals_orig,2));
bar_width = 1/(size(obs_vals_orig,2)+1);

% For each group
for ii = 1 : size(obs_vals_orig,1)
    for jj = 1:size(obs_vals_orig,2)
        if iscell(obs_vals_orig)
            obs_vals = obs_vals_orig{ii,jj};
        else
            obs_vals = obs_vals_orig(ii,jj);
        end
        
    % Jitter along x axis (to give the cloud appearance)
    x_jitter_obs  = ii + bar_dist(jj)  - (randn(size(obs_vals)) ./ (size(obs_vals_orig,2)*10));
    
    % Draw distributions
%     scatter( x_jitter_perm, perm_vals(:,ii), '.','MarkerEdgeAlpha',0.25,'markerEdgeColor',[0.5 0.5 0.5])
    scatter( x_jitter_obs, obs_vals,20, 'o','MarkerEdgeAlpha',0.5,'MarkerEdgeColor','none','MarkerFaceColor',cols(jj,:),'MarkerFaceAlpha',0.5)
   
    % Draw summary stats
    mean_val = mean( obs_vals);
    std_val = std( obs_vals)/sqrt(length(obs_vals)-1);
    bar( ii +bar_dist(jj) , mean_val, 'FaceColor','none','BarWidth',bar_width)
%     plot( [ii ii], [-std_val std_val]+mean_val,'k')
    errorbar(ii + bar_dist(jj), mean_val, std_val,'k','LineStyle','none','CapSize',5)
%     std_val = std(perm_vals(:,ii));
%     bar( ii + 0.2, mean( perm_vals(:,ii)), 'FaceColor','none','barwidth',0.2)
%     plot( [ii ii] + 0.2, [-std_val std_val]+mean( perm_vals(:,ii)),'k')
%     errorbar(ii+0.2, mean( perm_vals(:,ii)), std_val,'k','LineStyle','none','CapSize',3)
    end
end

set(gca,'xtick',1:numel(labels),'xticklabel',labels)

% xlabel('labels')
% ylabel('% Correct')