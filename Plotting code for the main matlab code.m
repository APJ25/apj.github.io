function plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, fig_num, plot_type, Fr, CopperLoss, PowerLoss, Optimal, Pmax)
    figure(fig_num); clf; hold on; axis square;

    switch plot_type
       case "loss_vs_flux"
    legend_handles = [];  % Store handles for legend
    legend_labels = {};   % Store labels for legend

    for i = 1:numel(core_indices)
        idx = core_indices(i);

        % Plot Series (solid line)
        h_series = plot(Flux.Series{idx}, Pfes{idx}, '-', ...
                        'Color', colors(i,:), 'LineWidth', 2);
        legend_handles(end+1) = h_series;
        legend_labels{end+1} = [core_names{idx} '-S'];

        % Plot Parallel (dashed line)
        h_parallel = plot(Flux.Parallel{idx}, Pfep{idx}, '--', ...
                          'Color', colors(i,:), 'LineWidth', 2);
        legend_handles(end+1) = h_parallel;
        legend_labels{end+1} = [core_names{idx} '-P'];
    end

    xlabel('B (Tesla)');
    ylabel('Core loss (W/m^3)');
    grid on;
    legend(legend_handles, legend_labels, 'Location', 'best');
    xlim([0 0.2]);
    ylim([0 2]);
    
        case "loss_vs_turns"
            legend_entries = strings(1, 2 * numel(core_indices));
for i = 1:numel(core_indices)
    idx = core_indices(i);

    plot(Turns.Series.Primary{idx}, Pfes{idx}, 'Color', colors(i,:), 'LineWidth', 2);
    plot(Turns.Parallel.Primary{idx}, Pfep{idx}, '--', 'Color', colors(i,:), 'LineWidth', 2);

    legend_entries(2*i - 1) = core_names{idx} + " -S";  % Series
    legend_entries(2*i)     = core_names{idx} + " -P";  % Parallel
end

xlabel('Turns');
ylabel('Core Loss (W)');
legend(legend_entries, 'Location', 'best');
xlim([10 50]);
ylim([0 3]);
grid on;

        case "flux_vs_turns"
             legend_handles = [];  
    legend_labels = [];  
            for i = 1:numel(core_indices)
                idx = core_indices(i);
                  h_series = plot(Turns.Series.Primary{idx}, Flux.Series{idx}, '-', ...
                        'Color', colors(i,:), 'LineWidth', 2);
                    legend_handles(end+1) = h_series;
        legend_labels{end+1} = [core_names{idx} '-S'];
                h_parallel = plot(Turns.Parallel.Primary{idx}, Flux.Parallel{idx}, '--', ...
                          'Color', colors(i,:), 'LineWidth', 2);
        legend_handles(end+1) = h_parallel;
        legend_labels{end+1} = [core_names{idx} '-P'];
            end
            xlabel('Turns');
            ylabel('B_{max} (Tesla)');
            % title("Turns vs B_{max} - Series and Parallel");
           legend(legend_handles, legend_labels, 'Location', 'eastoutside');
            xlim([0 200]); ylim([0 0.3]);
            grid on;
        case "compare_core"
            idx = core_indices(1);
            plot(Flux.Series{idx}, Pfes{idx}, 'b', Flux.Parallel{idx}, Pfep{idx}, 'g', 'LineWidth', 3);
            xlabel('B_{max} (Tesla)');
            ylabel('Core loss (W/m^3)');
            title(["Core Losses Comparison for ", core_names{idx}]);
            legend("Series", "Parallel");

        case "layers_vs_turns"
            for i = 1:numel(core_indices)
                idx = core_indices(i);
                plot(Turns.Series.Primary{idx}, Layers.Series.Primary{idx}, 'Color', colors(i,:), 'LineWidth', 2);
                %plot(Turns.Parallel.Primary{idx}, Layers.Parallel.Primary{idx}, '--', 'Color', colors(i,:), 'LineWidth', 2);
            end
            xlabel('Turns');
            ylabel('N Layers');
            
            legend([strcat(core_names(core_indices))], 'Location', 'best');
            grid on;
            axis square;
               case "skeff proximity"
            idx = core_indices(1);

            % Extract all vectors
            x1 = Turns.Series.Primary{idx};
            y1 = Fr.Parallel.Secondary{idx};
            y2 = Fr.Series.Secondary{idx};
            y4 = Fr.Parallel.Primary{idx};
            y5 = Fr.Series.Primary{idx};

            % Determine shortest valid length
            all_lengths = [length(x1), length(y1), length(y2), length(y4), length(y5)];
            minLen = min(all_lengths);

            % Trim all vectors
            x = x1(1:minLen);
            y1 = y1(1:minLen);
            y2 = y2(1:minLen);
            y4 = y4(1:minLen);
            y5 = y5(1:minLen);

            % Compute total losses
            y3 = y1 + y2;
            y6 = y4 + y5;

            % Plot
            plot(x, y1, 'b', 'LineWidth', 3);
            plot(x, y2, 'g', 'LineWidth', 3);
            plot(x, y3, 'r', 'LineWidth', 3);
            plot(x, y4, 'b--', 'LineWidth', 3);
            plot(x, y5, 'g--', 'LineWidth', 3);
            plot(x, y6, 'r--', 'LineWidth', 3);

            xlabel('Number of Turns');
            ylabel('Effect Losses (W)');
            title('Primary and Secondary Skeff-Proeff of Both Configurations');
            legend(...
                'Series Primary', 'Series Secondary', 'Total Series', ...
                'Parallel Primary', 'Parallel Secondary', 'Total Parallel', ...
                'Location', 'best');
            grid on;


        case "copper loss vs turns"
            legend_entries = strings(1, 2 * numel(core_indices));

for i = 1:numel(core_indices)
    idx = core_indices(i);

    % Plot Series (solid line)
    plot(Turns.Series.Primary{idx}, CopperLoss.Series{idx}, 'Color', colors(i,:), 'LineWidth', 2);
    legend_entries(2*i - 1) = core_names{idx} + " -S";  % Series label

    % Plot Parallel (dashed line)
    plot(Turns.Parallel.Primary{idx}, CopperLoss.Parallel{idx}, ':', 'Color', colors(i,:), 'LineWidth', 2);
    legend_entries(2*i) = core_names{idx} + " -P";      % Parallel label
end

xlabel('Turns');
ylabel('Copper Loss (W/m^3)');
title("Turns vs Copper Losses - Series and Parallel");
legend(legend_entries, 'Location', 'best');
xlim([10 150]);
ylim([0 50]);
grid on;

%            %% --- Series Plot ---
% figure;
% hold on;
% legend_entries_series = strings(1, numel(core_indices));
% 
% for i = 1:numel(core_indices)
%     idx = core_indices(i);
%     plot(Turns.Series.Primary{idx}, CopperLoss.Series{idx}, 'Color', colors(i,:), 'LineWidth', 2);
%     legend_entries_series(i) = core_names{idx} ;
% end
% 
% xlabel('Turns');
% ylabel('Copper Loss (W)');
% 
% legend(legend_entries_series, 'Location', 'best');
% xlim([10 150]);
% ylim([0 50]);
% grid on;
% box on;
% axis square;
% 
% %% --- Parallel Plot ---
% figure;
% hold on;
% legend_entries_parallel = strings(1, numel(core_indices));
% 
% for i = 1:numel(core_indices)
%     idx = core_indices(i);
%     plot(Turns.Parallel.Primary{idx}, CopperLoss.Parallel{idx}, 'Color', colors(i,:), 'LineWidth', 2);
%     legend_entries_parallel(i) = core_names{idx} ;
% end
% 
% xlabel('Turns');
% ylabel('Copper Loss (W)');
% 
% legend(legend_entries_parallel, 'Location', 'best');
% xlim([10 150]);
% ylim([0 25]);
% grid on;
% box on;
% axis square;
        case "Total power loss vs flux density series"
    bopt_markers = repmat({'ko'}, 1, numel(core_indices));
    pmax_markers = repmat({'ks'}, 1, numel(core_indices));
    legend_handles = [];
    legend_labels = [];

    for i = 1:numel(core_indices)
        idx = core_indices(i);

        % Plot series total power loss
        h1 = plot(Flux.Series{idx}, PowerLoss.Series{idx}, 'Color', colors(i,:), 'LineWidth', 2);
        legend_handles = [legend_handles, h1];
        legend_labels = [legend_labels, string(sprintf("(%s)_{P_{tot}}", core_names{idx}))];

        % Bopt marker
        Bopt = Optimal.Series.B(idx);
        [~, opt_idx] = min(abs(Flux.Series{idx} - Bopt));
        h2 = plot(Bopt, PowerLoss.Series{idx}(opt_idx), bopt_markers{i}, 'MarkerSize', 6, 'MarkerFaceColor', colors(i,:));
        legend_handles = [legend_handles, h2];
        legend_labels = [legend_labels, string(sprintf("(%s)_{B_{opt}}", core_names{idx}))];

        % Pmax marker
        h3 = plot(Bopt, Pmax(idx), pmax_markers{i}, 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:));
        legend_handles = [legend_handles, h3];
        legend_labels = [legend_labels, string(sprintf("(%s)_{P_{max}}", core_names{idx}))];
    end

    xlim([0, 0.3]);
    ylim([1, 7]);
    xlabel('Flux Density (T)');
    ylabel('Total Losses (W)');
    title('Flux Density vs Total Losses (Series)');
    % legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 8);
    grid on;
    axis square;
         case "Total power loss vs flux density parallel"
    bopt_markers = repmat({'ko'}, 1, numel(core_indices));
    pmax_markers = repmat({'ks'}, 1, numel(core_indices));
    legend_handles = [];
    legend_labels = [];

    for i = 1:numel(core_indices)
        idx = core_indices(i);

        % Plot parallel total power loss
        h1 = plot(Flux.Parallel{idx}, PowerLoss.Parallel{idx}, 'Color', colors(i,:), 'LineWidth', 2);
        legend_handles = [legend_handles, h1];
        legend_labels = [legend_labels, string(sprintf("(%s)_{P_{tot}}", core_names{idx}))];

        % Bopt marker
        Bopt = Optimal.Parallel.B(idx);
        [~, opt_idx] = min(abs(Flux.Parallel{idx} - Bopt));
        h2 = plot(Bopt, PowerLoss.Parallel{idx}(opt_idx), bopt_markers{i}, 'MarkerSize', 6, 'MarkerFaceColor', colors(i,:));
        legend_handles = [legend_handles, h2];
        legend_labels = [legend_labels, string(sprintf("(%s)_{B_{opt}}", core_names{idx}))];

        % Pmax marker
        h3 = plot(Bopt, Pmax(idx), pmax_markers{i}, 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:));
        legend_handles = [legend_handles, h3];
        legend_labels = [legend_labels, string(sprintf("(%s)_{P_{max}}", core_names{idx}))];
    end

    xlim([0, 0.3]);
    ylim([1, 7]);
    xlabel('Flux Density (T)');
    ylabel('Total Losses (W)');
    title('Flux Density vs Total Losses (Parallel)');
    legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 8);
    grid on;
    axis square;
    case "Total power loss vs turns series"
    bopt_markers = repmat({'ko'}, 1, numel(core_indices));
    pmax_markers = repmat({'ks'}, 1, numel(core_indices));
    legend_handles = [];
    legend_labels = [];

    for i = 1:numel(core_indices)
        idx = core_indices(i);

        % Plot series total power loss vs Turns
        h1 = plot(Turns.Series.Primary{idx}, PowerLoss.Series{idx}, 'Color', colors(i,:), 'LineWidth', 2);
        legend_handles = [legend_handles, h1];
        legend_labels = [legend_labels, string(sprintf("(%s)_{P_{tot}}", core_names{idx}))];

        % Nopt marker
        Nopt = Optimal.Series.Np(idx);
        [~, opt_idx] = min(abs(Turns.Series.Primary{idx} - Nopt));
        h2 = plot(Nopt, PowerLoss.Series{idx}(opt_idx), bopt_markers{i}, 'MarkerSize', 6, 'MarkerFaceColor', colors(i,:));
        legend_handles = [legend_handles, h2];
        legend_labels = [legend_labels, string(sprintf("(%s)_{N_{opt}}", core_names{idx}))];

        % Pmax marker
        h3 = plot(Nopt, Pmax(idx), pmax_markers{i}, 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:));
        legend_handles = [legend_handles, h3];
        legend_labels = [legend_labels, string(sprintf("(%s)_{P_{max}}", core_names{idx}))];
    end

    xlim([0, 100]);
    ylim([0, 7]);
    xlabel('Primary Turns (N)');
    ylabel('Total Losses (W)');
    title('Turns vs Total Losses (Series)');
    % legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 8);
    grid on;
    axis square;

        case "Total power loss vs turns parallel"
    bopt_markers = repmat({'ko'}, 1, numel(core_indices));
    pmax_markers = repmat({'ks'}, 1, numel(core_indices));
    legend_handles = [];
    legend_labels = [];

    for i = 1:numel(core_indices)
        idx = core_indices(i);

        % Plot parallel total power loss vs Turns
        h1 = plot(Turns.Parallel.Primary{idx}, PowerLoss.Parallel{idx}, 'Color', colors(i,:), 'LineWidth', 2);
        legend_handles = [legend_handles, h1];
        legend_labels = [legend_labels, string(sprintf("(%s)_{P_{tot}}", core_names{idx}))];

        % Nopt marker
        Nopt = Optimal.Parallel.Np(idx);
        [~, opt_idx] = min(abs(Turns.Parallel.Primary{idx} - Nopt));
        h2 = plot(Nopt, PowerLoss.Parallel{idx}(opt_idx), bopt_markers{i}, 'MarkerSize', 6, 'MarkerFaceColor', colors(i,:));
        legend_handles = [legend_handles, h2];
        legend_labels = [legend_labels, string(sprintf("(%s)_{N_{opt}}", core_names{idx}))];

        % Pmax marker
        h3 = plot(Nopt, Pmax(idx), pmax_markers{i}, 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:));
        legend_handles = [legend_handles, h3];
        legend_labels = [legend_labels, string(sprintf("(%s)_{P_{max}}", core_names{idx}))];
    end

    xlim([0, 100]);
    ylim([0, 7]);
    xlabel('Primary Turns (N)');
    ylabel('Total Losses (W)');
    title('Turns vs Total Losses (Parallel)');
    legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 8);
    grid on;
    axis square;
   case "Total power loss vs turns combined"
    legend_handles = [];
    legend_labels = [];

    for i = 1:numel(core_indices)
        idx = core_indices(i);

        % Check if data exists and is not empty
        if isempty(PowerLoss.Series{idx}) || isempty(PowerLoss.Parallel{idx}) || isempty(Turns.Parallel.Primary{idx})
            continue;
        end

        % Retrieve data
        series_loss = PowerLoss.Series{idx};
        parallel_loss = PowerLoss.Parallel{idx};
        turns = Turns.Parallel.Primary{idx};  % Reference from parallel

        % Ensure same length
        len = min([length(series_loss), length(parallel_loss), length(turns)]);
        total_loss = (series_loss(1:len) + parallel_loss(1:len)) / 2;
        turns = turns(1:len);

        % Plot if valid
        if len > 1
            h = plot(turns, total_loss, 'Color', colors(i,:), 'LineWidth', 2);
            if isgraphics(h)
                legend_handles(end+1) = h;
                legend_labels{end+1} = core_names{idx};
            end
        end
    end

    xlabel('Np (Turns)');
    ylabel('Total Loss (Series + Parallel) (W/m^3)');
    title('Total Transformer Losses (Series + Parallel)');
    xlim([0 200]);
    ylim([0 20]);
    grid on;
    axis square;
    % Only show legend if there's something to show
    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 8);
    end
    case "optimal turns per core"
    optimal_turns = zeros(1, numel(core_indices));
    for i = 1:numel(core_indices)
        idx = core_indices(i);
        % Ensure all vectors have the same length
        len = min([length(PowerLoss.Series{idx}), length(PowerLoss.Parallel{idx}), length(Turns.Parallel.Primary{idx})]);
        if len == 0
            optimal_turns(i) = NaN;
            continue;
        end
        % Compute total loss and find optimal point
        total_loss = (PowerLoss.Series{idx}(1:len) + PowerLoss.Parallel{idx}(1:len)) / 2;
        [~, min_idx] = min(total_loss);
        turns = Turns.Parallel.Primary{idx}(1:len);
        optimal_turns(i) = turns(min_idx);
    end

    core_labels = core_names(core_indices);
    scatter(1:numel(core_indices), optimal_turns, 80, 'filled');
    set(gca, 'XTick', 1:numel(core_indices), 'XTickLabel', core_labels);
    xtickangle(45);
    ylabel('Optimal Turns');
    title('Optimal Number of Turns per Core');
    grid on;
    axis square;
case "turns range errorbar"
    % Define core labels for x-axis (series and parallel for selected cores)
    coreLabels = {'25s', '30s','35s', '40s', '42s', '50s', '55s','60s','ETD39s', ...
                  '25p', '30p','35p', '40p', '42p', '50p', '55p','60p','ETD39p'};

    % Extract min and max of series and parallel turns
    minTurns = [min(Turns.Series.Primary{5}), min(Turns.Series.Primary{6}), min(Turns.Series.Primary{11}), ...
                min(Turns.Series.Primary{7}), min(Turns.Series.Primary{12}), min(Turns.Series.Primary{8}), ...
                min(Turns.Series.Primary{13}), min(Turns.Series.Primary{9}), min(Turns.Series.Primary{10}), ...
                min(Turns.Parallel.Primary{5}), min(Turns.Parallel.Primary{6}), min(Turns.Parallel.Primary{11}), ...
                min(Turns.Parallel.Primary{7}), min(Turns.Parallel.Primary{12}), min(Turns.Parallel.Primary{8}), ...
                min(Turns.Parallel.Primary{13}), min(Turns.Parallel.Primary{9}), min(Turns.Parallel.Primary{10})];

    maxTurns = [max(Turns.Series.Primary{5}), max(Turns.Series.Primary{6}), max(Turns.Series.Primary{11}), ...
                max(Turns.Series.Primary{7}), max(Turns.Series.Primary{12}), max(Turns.Series.Primary{8}), ...
                max(Turns.Series.Primary{13}), max(Turns.Series.Primary{9}), max(Turns.Series.Primary{10}), ...
                max(Turns.Parallel.Primary{5}), max(Turns.Parallel.Primary{6}), max(Turns.Parallel.Primary{11}), ...
                max(Turns.Parallel.Primary{7}), max(Turns.Parallel.Primary{12}), max(Turns.Parallel.Primary{8}), ...
                max(Turns.Parallel.Primary{13}), max(Turns.Parallel.Primary{9}), max(Turns.Parallel.Primary{10})];

    % Midpoints and half-ranges for error bars
    midTurns = (minTurns + maxTurns) / 2;
    errorBars = (maxTurns - minTurns) / 2;

    % Plot
    errorbar(1:numel(coreLabels), midTurns, errorBars, 'k.', 'LineWidth', 1.5, 'MarkerSize', 10);

    % Format axes
    set(gca, 'XTickLabel', coreLabels, 'XTick', 1:numel(coreLabels), 'XTickLabelRotation', 45);
    xlabel('Core Name (s = series, p = parallel)');
    ylabel('Range of Number of Turns');
    
    ylim([0, max(maxTurns) * 1.1]);
    grid on;



    
    end

    
   
end
