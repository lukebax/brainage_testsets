%% load data

clear

robust_indicator_list = ["no", "yes"];

for robust_indicator_counter = 1 : length(robust_indicator_list)

    clearvars -except robust_indicator_list robust_indicator_counter

    robust_indicator = robust_indicator_list(robust_indicator_counter);

    input_folder = '../data/';

    input_file = 'data_leuven.csv';

    output_folder = '../results/';

    T = readtable(strcat(input_folder, input_file));


    %% remove age-association bias

    pma = T.pma;
    pmaPredicted = T.pmaPredicted;

    % unadjusted values
    predError_unadjusted = pma - pmaPredicted;
    brainAgeGap_unadjusted = abs(predError_unadjusted);

    % adjusted values (remove age-association bias)
    mdl = fitlm(pma, predError_unadjusted, 'RobustOpts','on');
    if strcmp(robust_indicator, "no")
        predError_adjusted = mdl.Residuals.Raw - mean(mdl.Residuals.Raw) + mean(predError_unadjusted);
    elseif strcmp(robust_indicator, "yes")
        predError_adjusted = mdl.Residuals.Raw - median(mdl.Residuals.Raw) + median(predError_unadjusted);
    end

    brainAgeGap_adjusted = abs(predError_adjusted);

    BrainAgeGap_rec = brainAgeGap_adjusted;

    %% get subject-level data from recording-level data

    subjectID = string(T.subjectID);

    bsid_recording = string(T.group);

    subjectID_unique = unique(subjectID);

    num_subjects = length(subjectID_unique);

    BrainAgeGap_subj_all = zeros(num_subjects, 1);

    bsid_subj = strings(num_subjects, 1);

    numRecordings_subj = zeros(num_subjects, 1);

    for subj_counter = 1 : num_subjects

        subj = subjectID_unique(subj_counter);

        subj_idx = find(strcmp(subj, subjectID));

        if strcmp(robust_indicator, "no")
            BrainAgeGap_subj = mean(BrainAgeGap_rec(subj_idx));
        elseif strcmp(robust_indicator, "yes")
            BrainAgeGap_subj = median(BrainAgeGap_rec(subj_idx));
        end

        bsid = unique(bsid_recording(subj_idx));

        numRecordings = length(BrainAgeGap_rec(subj_idx));

        BrainAgeGap_subj_all(subj_counter) = BrainAgeGap_subj;

        bsid_subj(subj_counter) = bsid;

        numRecordings_subj(subj_counter) = numRecordings;

    end

    idx = cat(1, find(strcmp(bsid_subj, 'Normal')), find(strcmp(bsid_subj, 'Mild')), find(strcmp(bsid_subj, 'Severe')));

    BrainAgeGap_subj_all = BrainAgeGap_subj_all(idx);

    Group = bsid_subj(idx);

    Recordings = numRecordings_subj(idx);

    %% do non-parametric stats

    design_matrix_ev1_norm_mean = double(strcmp(Group, 'Normal'));
    design_matrix_ev2_mild_mean = double(strcmp(Group, 'Mild'));
    design_matrix_ev3_sev_mean = double(strcmp(Group, 'Severe'));
    design_matrix_ev4_rec_count = double(Recordings) - mean(double(Recordings));

    design_matrix = cat(2, design_matrix_ev1_norm_mean, design_matrix_ev2_mild_mean, design_matrix_ev3_sev_mean, design_matrix_ev4_rec_count);
    t_contrast_matrix_ttests = [-1,1,0,0; -1,0,1,0; 0,-1,1,0];

    writematrix(BrainAgeGap_subj_all, 'tmp_i.csv');
    writematrix(design_matrix, 'tmp_d.csv');
    writematrix(t_contrast_matrix_ttests, 'tmp_t.csv');

    palm -i tmp_i.csv ...
        -d tmp_d.csv ...
        -t tmp_t.csv ...
        -n 10000 ...
        -ee ...
        -twotail ...
        -corrcon ...
        -cmcx ...
        -saveglm ...
        -o tmp ...
        -quiet


    if strcmp(robust_indicator, "no")

        results_brainAgeGap_leuven.num_subjects = num_subjects;
        results_brainAgeGap_leuven.num_recordings = sum(Recordings);

        results_brainAgeGap_leuven.NORM_num_subjects = sum(strcmp(Group, 'Normal'));
        results_brainAgeGap_leuven.NORM_num_recordings = sum(Recordings(strcmp(Group, 'Normal')));
        results_brainAgeGap_leuven.NORM_mae = mean(BrainAgeGap_subj_all(strcmp(Group, 'Normal')));
        results_brainAgeGap_leuven.MILD_num_subjects = sum(strcmp(Group, 'Mild'));
        results_brainAgeGap_leuven.MILD_num_recordings = sum(Recordings(strcmp(Group, 'Mild')));
        results_brainAgeGap_leuven.MILD_mae = mean(BrainAgeGap_subj_all(strcmp(Group, 'Mild')));
        results_brainAgeGap_leuven.SEV_num_subjects = sum(strcmp(Group, 'Severe'));
        results_brainAgeGap_leuven.SEV_num_recordings = sum(Recordings(strcmp(Group, 'Severe')));
        results_brainAgeGap_leuven.SEV_mae = mean(BrainAgeGap_subj_all(strcmp(Group, 'Severe')));

        results_brainAgeGap_leuven.pairwise_ttest_NORMvsMILD_cope = readmatrix('tmp_dat_cope_c1.csv');
        results_brainAgeGap_leuven.pairwise_ttest_NORMvsMILD_cohen = readmatrix('tmp_dat_cohen_c1.csv');
        results_brainAgeGap_leuven.pairwise_ttest_NORMvsMILD_tstat = readmatrix('tmp_dat_tstat_c1.csv');
        results_brainAgeGap_leuven.pairwise_ttest_NORMvsMILD_tstat_pval = readmatrix('tmp_dat_tstat_cfwep_c1.csv');

        results_brainAgeGap_leuven.pairwise_ttest_NORMvsSEV_cope = readmatrix('tmp_dat_cope_c2.csv');
        results_brainAgeGap_leuven.pairwise_ttest_NORMvsSEV_cohen = readmatrix('tmp_dat_cohen_c2.csv');
        results_brainAgeGap_leuven.pairwise_ttest_NORMvsSEV_tstat = readmatrix('tmp_dat_tstat_c2.csv');
        results_brainAgeGap_leuven.pairwise_ttest_NORMvsSEV_tstat_pval = readmatrix('tmp_dat_tstat_cfwep_c2.csv');

        results_brainAgeGap_leuven.pairwise_ttest_MILDvsSEV_cope = readmatrix('tmp_dat_cope_c3.csv');
        results_brainAgeGap_leuven.pairwise_ttest_MILDvsSEV_cohen = readmatrix('tmp_dat_cohen_c3.csv');
        results_brainAgeGap_leuven.pairwise_ttest_MILDvsSEV_tstat = readmatrix('tmp_dat_tstat_c3.csv');
        results_brainAgeGap_leuven.pairwise_ttest_MILDvsSEV_tstat_pval = readmatrix('tmp_dat_tstat_cfwep_c3.csv');

        results_file_name = strcat(output_folder, 'results_brainAgeGap_leuven.mat');
        save(results_file_name, 'results_brainAgeGap_leuven')

    elseif strcmp(robust_indicator, "yes")

        results_brainAgeGap_ROB_leuven.num_subjects = num_subjects;
        results_brainAgeGap_ROB_leuven.num_recordings = sum(Recordings);

        results_brainAgeGap_ROB_leuven.NORM_num_subjects = sum(strcmp(Group, 'Normal'));
        results_brainAgeGap_ROB_leuven.NORM_num_recordings = sum(Recordings(strcmp(Group, 'Normal')));
        results_brainAgeGap_ROB_leuven.NORM_mae = mean(BrainAgeGap_subj_all(strcmp(Group, 'Normal')));
        results_brainAgeGap_ROB_leuven.MILD_num_subjects = sum(strcmp(Group, 'Mild'));
        results_brainAgeGap_ROB_leuven.MILD_num_recordings = sum(Recordings(strcmp(Group, 'Mild')));
        results_brainAgeGap_ROB_leuven.MILD_mae = mean(BrainAgeGap_subj_all(strcmp(Group, 'Mild')));
        results_brainAgeGap_ROB_leuven.SEV_num_subjects = sum(strcmp(Group, 'Severe'));
        results_brainAgeGap_ROB_leuven.SEV_num_recordings = sum(Recordings(strcmp(Group, 'Severe')));
        results_brainAgeGap_ROB_leuven.SEV_mae = mean(BrainAgeGap_subj_all(strcmp(Group, 'Severe')));

        results_brainAgeGap_ROB_leuven.pairwise_ttest_NORMvsMILD_cope = readmatrix('tmp_dat_cope_c1.csv');
        results_brainAgeGap_ROB_leuven.pairwise_ttest_NORMvsMILD_cohen = readmatrix('tmp_dat_cohen_c1.csv');
        results_brainAgeGap_ROB_leuven.pairwise_ttest_NORMvsMILD_tstat = readmatrix('tmp_dat_tstat_c1.csv');
        results_brainAgeGap_ROB_leuven.pairwise_ttest_NORMvsMILD_tstat_pval = readmatrix('tmp_dat_tstat_cfwep_c1.csv');

        results_brainAgeGap_ROB_leuven.pairwise_ttest_NORMvsSEV_cope = readmatrix('tmp_dat_cope_c2.csv');
        results_brainAgeGap_ROB_leuven.pairwise_ttest_NORMvsSEV_cohen = readmatrix('tmp_dat_cohen_c2.csv');
        results_brainAgeGap_ROB_leuven.pairwise_ttest_NORMvsSEV_tstat = readmatrix('tmp_dat_tstat_c2.csv');
        results_brainAgeGap_ROB_leuven.pairwise_ttest_NORMvsSEV_tstat_pval = readmatrix('tmp_dat_tstat_cfwep_c2.csv');

        results_brainAgeGap_ROB_leuven.pairwise_ttest_MILDvsSEV_cope = readmatrix('tmp_dat_cope_c3.csv');
        results_brainAgeGap_ROB_leuven.pairwise_ttest_MILDvsSEV_cohen = readmatrix('tmp_dat_cohen_c3.csv');
        results_brainAgeGap_ROB_leuven.pairwise_ttest_MILDvsSEV_tstat = readmatrix('tmp_dat_tstat_c3.csv');
        results_brainAgeGap_ROB_leuven.pairwise_ttest_MILDvsSEV_tstat_pval = readmatrix('tmp_dat_tstat_cfwep_c3.csv');

        results_file_name = strcat(output_folder, 'results_brainAgeGap_ROB_leuven.mat');
        save(results_file_name, 'results_brainAgeGap_ROB_leuven')

    end

    delete('tmp*');


    %% write out df for use in R (dabestr)

    dabestr_df = table(Group);

    dabestr_df.BrainAgeGap = BrainAgeGap_subj_all;

    if strcmp(robust_indicator, "no")
        file_name = strcat(input_folder, 'data_dabestrDF_leuven.csv');
    elseif strcmp(robust_indicator, "yes")
        file_name = strcat(input_folder, 'data_dabestrDF_ROB_leuven.csv');
    end

    writetable(dabestr_df, file_name);

end
