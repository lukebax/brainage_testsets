%%

clear

dataset_list = ["leuven", "oxford"];

T_subjectCharacteristics = table;

for dataset_counter = 1 : length(dataset_list)

    dataset = dataset_list(dataset_counter);

    input_folder = '../data/';

    input_file = strcat('data_', dataset, '.csv');

    output_folder = '../results/';

    if ~exist(output_folder, 'dir')
        mkdir(output_folder)
    end

    T = readtable(strcat(input_folder, input_file));

    subjects = string(unique(T.subjectID));

    recordings = string(T.subjectID);

    num_subjects = size(unique(subjects), 1);

    num_recordings = size(T, 1);
    pma_mean = mean(T.pma);
    pma_std = std(T.pma);

    num_rec_per_sub = zeros(num_subjects, 1);
    for sub_counter = 1 : num_subjects
        sub = subjects(sub_counter);
        num_rec_per_sub(sub_counter) = sum(sub == recordings);
    end
    num_rec_per_sub_mean = mean(num_rec_per_sub);
    num_rec_per_sub_std = std(num_rec_per_sub);


    if strcmp(dataset, "leuven")
        bsid2_groups_recordings = string(T.group);
        [~,ia,~] = unique(T.subjectID);
        bsid2_groups_subjects = bsid2_groups_recordings(ia);
        bsid2_normal = sum(strcmp(bsid2_groups_subjects, "Normal"));
        bsid2_mildAbnormal = sum(strcmp(bsid2_groups_subjects, "Mild"));
        bsid2_severeAbnormal = sum(strcmp(bsid2_groups_subjects, "Severe"));
    elseif strcmp(dataset, "oxford")
        bsid2_normal = 0;
        bsid2_mildAbnormal = 0;
        bsid2_severeAbnormal = 0;
    end

    T_subjectCharacteristics.(dataset) = [...
        num_subjects; ...
        bsid2_normal; ...
        bsid2_mildAbnormal; ...
        bsid2_severeAbnormal; ...
        pma_mean; ...
        pma_std; ...
        num_recordings; ...
        num_rec_per_sub_mean; ...
        num_rec_per_sub_std];

end

T_subjectCharacteristics.Properties.RowNames = [...
    "num_subjects"; ...
    "bsid2_normal"; ...
    "bsid2_mildAbnormal"; ...
    "bsid2_severeAbnormal"; ...
    "pma_mean"; ...
    "pma_std"; ...
    "num_recordings"; ...
    "num_rec_per_sub_mean"; ...
    "num_rec_per_sub_std"];

file_name = strcat(output_folder, 'summary_charateristics.csv');

writetable(T_subjectCharacteristics, file_name, "WriteRowNames", true)
