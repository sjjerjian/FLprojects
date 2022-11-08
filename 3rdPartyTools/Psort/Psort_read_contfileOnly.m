function contfile = Psort_read_contfileOnly(file_fullPath)
% This function is part of PurkinjeSort project
% it reads psort file and returns the psortDataBase

% CF: quick subfunction to grab the filename of the associated continuous
% file (eventually stored in topLevel_data)

file_info = h5info(file_fullPath);
num_slots = length(file_info.Groups.Groups);
slot_name = ['/data/i' num2str(num_slots-1)];
variable_name = 'file_fullPathOriginal';
variable_data = h5read(file_fullPath,[slot_name '/' variable_name]);
contfile = char(variable_data(1:4:end))';
