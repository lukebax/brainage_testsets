%%

clear

robust_indicator_list = ["no", "yes"];

dataset_list = ["leuven", "oxford"];

for robust_indicator_counter = 1 : length(robust_indicator_list)

    clearvars -except robust_indicator_list dataset_list robust_indicator_counter

    robust_indicator = robust_indicator_list(robust_indicator_counter);

    for dataset_counter = 1 : length(dataset_list)

        clearvars -except robust_indicator_list dataset_list robust_indicator_counter dataset_counter robust_indicator 

        dataset = dataset_list(dataset_counter);

        input_folder = '../data/';

        input_file = strcat('data_', dataset, '.csv');

        output_folder = '../results/';

        T = readtable(strcat(input_folder, input_file));

        %%

        subjectID = string(T.subjectID);

        pma = double(T.pma);

        pmaPredicted = double(T.pmaPredicted);

        num_observations = length(subjectID);

        %%

        level_1 = ones(num_observations, 1);

        level_2 = zeros(num_observations, 1);

        for counter = 1 : num_observations

            ID = subjectID(counter);

            idx = find(strcmp(subjectID(:), ID));

            count = length(idx);

            level_2(idx) = count;

        end

        %%

        max_num_recordings = max(level_2);

        level_3 = zeros(num_observations, 1);


        counter = 1;

        for i = 1 : max_num_recordings

            idx = find(level_2 == i);

            subj_unique = unique(subjectID(idx));

            num_subj_unique = length(subj_unique);

            for j = 1 : num_subj_unique

                k = subj_unique(j);

                l = find(strcmp(subjectID(:), k));

                if i ~= 1
                    counter = counter + 1;
                end

                level_3(l) = counter;

            end

        end

        %%

        level_4 = palm_reindex(level_3,'fixleaves');

        level_4 = level_4(:,2);

        %%

        EB = cat(2, level_1, level_2, level_3, level_4);

        EB(:,1) = EB(:,1) * -1;

        M = (1:size(EB,1))';

        Ptree = palm_tree(EB, M);

        if strcmp(input_file, 'data_leuven.csv')
            ptree_file_name = strcat(output_folder, 'ptree_leuven.dot');
            eb_file_name = strcat(output_folder, 'EB_leuven.csv');
        elseif strcmp(input_file, 'data_oxford.csv')
            ptree_file_name = strcat(output_folder, 'ptree_oxford.dot');
            eb_file_name = strcat(output_folder, 'EB_oxford.csv');
        end

        palm_ptree2dot(Ptree, ptree_file_name);

        writematrix(EB, eb_file_name)

        %% do analysis

        num_resamples = 10000;

        % MAE
        ae = abs(pma - pmaPredicted);
        if strcmp(robust_indicator, "no")
            mae = mean(ae);
        elseif strcmp(robust_indicator, "yes")
            mae = median(ae);
        end

        % R2
        if strcmp(robust_indicator, "no")
            numerator = sum((pma - pmaPredicted).^2);
            denominator = sum((pma - mean(pma)).^2);
            r2 = 1 - (numerator / denominator);
        elseif strcmp(robust_indicator, "yes")
            numerator = median(abs(pma - pmaPredicted));
            denominator = median(abs(pma - median(pma)));
            r2 = 1 - (numerator / denominator).^2;
        end

        % confidence interval (bootstraps)
        if strcmp(robust_indicator, "no")
            mae_ci = bootci(num_resamples, @mean, abs(pma - pmaPredicted));
        elseif strcmp(robust_indicator, "yes")
            mae_ci = bootci(num_resamples, @median, abs(pma - pmaPredicted));
        end

        % p-value (permutations)
        permidx = palm_quickperms(M, EB, num_resamples, true, [], [], []);

        y_perm_all = pmaPredicted(permidx);

        mae_perm_all = zeros(num_resamples, 1);

        for perm_counter = 1 : num_resamples

            y_perm = y_perm_all(:, perm_counter);

            if strcmp(robust_indicator, "no")
                mae_perm = mean(abs(pma - y_perm));
            elseif strcmp(robust_indicator, "yes")
                mae_perm = median(abs(pma - y_perm));
            end

            mae_perm_all(perm_counter) = mae_perm;

        end

        mae_pval = sum((mae >= mae_perm_all)) / num_resamples;


        % store results

        if strcmp(robust_indicator, "no")

            if strcmp(input_file, 'data_leuven.csv')

                results_agePrediction_leuven.numSubjects = length(unique(subjectID));
                results_agePrediction_leuven.numRecordings = length(subjectID);
                results_agePrediction_leuven.r2 = r2;
                results_agePrediction_leuven.mae = mae;
                results_agePrediction_leuven.mae_ci = mae_ci;
                results_agePrediction_leuven.mae_pval = mae_pval;

            elseif strcmp(input_file, 'data_oxford.csv')

                results_agePrediction_oxford.numSubjects = length(unique(subjectID));
                results_agePrediction_oxford.numRecordings = length(subjectID);
                results_agePrediction_oxford.r2 = r2;
                results_agePrediction_oxford.mae = mae;
                results_agePrediction_oxford.mae_ci = mae_ci;
                results_agePrediction_oxford.mae_pval = mae_pval;

            end


        elseif strcmp(robust_indicator, "yes")

            if strcmp(input_file, 'data_leuven.csv')

                results_agePrediction_ROB_leuven.numSubjects = length(unique(subjectID));
                results_agePrediction_ROB_leuven.numRecordings = length(subjectID);
                results_agePrediction_ROB_leuven.r2 = r2;
                results_agePrediction_ROB_leuven.mae = mae;
                results_agePrediction_ROB_leuven.mae_ci = mae_ci;
                results_agePrediction_ROB_leuven.mae_pval = mae_pval;

            elseif strcmp(input_file, 'data_oxford.csv')

                results_agePrediction_ROB_oxford.numSubjects = length(unique(subjectID));
                results_agePrediction_ROB_oxford.numRecordings = length(subjectID);
                results_agePrediction_ROB_oxford.r2 = r2;
                results_agePrediction_ROB_oxford.mae = mae;
                results_agePrediction_ROB_oxford.mae_ci = mae_ci;
                results_agePrediction_ROB_oxford.mae_pval = mae_pval;

            end

        end

        % save results


        if strcmp(robust_indicator, "no")

            if strcmp(input_file, 'data_leuven.csv')
                results_file_name = strcat(output_folder, 'results_agePrediction_leuven.mat');
                save(results_file_name, 'results_agePrediction_leuven')
            elseif strcmp(input_file, 'data_oxford.csv')
                results_file_name = strcat(output_folder, 'results_agePrediction_oxford.mat');
                save(results_file_name, 'results_agePrediction_oxford')

            end

        elseif strcmp(robust_indicator, "yes")

            if strcmp(input_file, 'data_leuven.csv')
                results_file_name = strcat(output_folder, 'results_agePrediction_ROB_leuven.mat');
                save(results_file_name, 'results_agePrediction_ROB_leuven')

            elseif strcmp(input_file, 'data_oxford.csv')
                results_file_name = strcat(output_folder, 'results_agePrediction_ROB_oxford.mat');
                save(results_file_name, 'results_agePrediction_ROB_oxford')

            end

        end

        %% plot results

        axis_limits = [25, 50];

        rgb_normal = [229,39,45] / 255;
        rgb_mild = [69,133,188] / 255;
        rgb_severe = [88,179,85] / 255;
        rgb_oxford = [152,78,163] / 255;
        c = zeros(length(T.group), 3);

        if strcmp(input_file, 'data_leuven.csv')

            for i = 1 : length(T.group)
                test = T.group(i);
                if strcmp(test, 'Normal')
                    c(i,:) = rgb_normal;
                elseif strcmp(test, 'Mild')
                    c(i,:) = rgb_mild;
                elseif strcmp(test, 'Severe')
                    c(i,:) = rgb_severe;
                end
            end

        elseif strcmp(input_file, 'data_oxford.csv')

            for i = 1 : length(pma)
                c(i,:) = rgb_oxford;
            end

        end

        marker_size = 80;
        font_size = 17;
        line_width = 1;

        figure;

        scatter(pma, pmaPredicted, marker_size, c,'filled');
        xlim(axis_limits)
        ylim(axis_limits)
        line(xlim, xlim, 'Color', 'k', 'LineWidth', line_width)
        lsline
        xlabel('PMA (weeks)')
        ylabel('Predicted PMA (weeks)')
        ax = gca;
        ax.FontSize = font_size;
        ax.LineWidth = line_width;

        if strcmp(input_file, 'data_leuven.csv')
            fig_file_name = strcat(output_folder, 'results_agePrediction_scatterPlot_leuven.png');
        elseif strcmp(input_file, 'data_oxford.csv')
            fig_file_name = strcat(output_folder, 'results_agePrediction_scatterPlot_oxford.png');
        end

        saveas(gcf, fig_file_name)

    end

end
