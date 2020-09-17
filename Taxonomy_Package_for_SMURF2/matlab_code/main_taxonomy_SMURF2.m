function main_taxonomy_SMURF2(batch_dir_list, results_filename, taxonomy_path) 

% results_filename = '~/PhD/BACTERIA/Vered/test_for_vered.csv';
% batch_dir_list = {...
%     '/data2/fefuks/Bacterial_Data/Vered/Reco_Data/Example_Batch1'...
%     '/data2/fefuks/Bacterial_Data/Vered/Reco_Data/Example_Batch2'};


% *************************** LOAD DB  ********************
 
% Load the 16S Headers
Header_uni = importdata('Header_uni_smurf2.csv');

% Load the name calls
load(taxonomy_path, 'taxa_name_calls','ranks_to_extract')

% ***************************  PREPARE TAXONOMY ********************
batch_samples_list = {};
for bd = 1:length(batch_dir_list)

    batch_dir = batch_dir_list{bd};    
    disp(['WORKING ON DIR ' num2str(bd) ': ' batch_dir])
    
    directories = dir(batch_dir);
    for dd = 3:length(directories)
        if directories(dd).isdir == 1
            sample_name = directories(dd).name;
            disp(['WORKING ON SAMPLE ' num2str(dd-2) ': ' sample_name])

            readsPath = [batch_dir '/' sample_name];
            save_reconstruction_new_nogroups(readsPath, sample_name, taxa_name_calls, ranks_to_extract, Header_uni)
            batch_samples_list{end+1} = readsPath;
        end
    end
        
end

% *************************** SAVE RECONSTRUCTION FILE ********************
scott_format_newer_func(batch_samples_list,results_filename)
