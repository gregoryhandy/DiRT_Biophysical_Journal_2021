function update_sim_database()

cd('./DiRT_Data')
dir_info = dir;

data_count = 1;
for ii = 1:size(dir_info,1)
    
    DiRT_filename = dir_info(ii).name;
    if contains(DiRT_filename,'DiRT_cyl')
        sim_database{data_count,1} = strcat('./DiRT_Data/',DiRT_filename);
        
        under_loc = strfind(DiRT_filename,'_');
        sim_database{data_count,2} = strcat('./PDE_ODE_Data/PDE_ODE_cyl_',DiRT_filename(under_loc(2)+1:end));
        
        data_count = data_count + 1;
        
    end
end
    
cd ..
save('sim_database','sim_database')
end