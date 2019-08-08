clear
clear
load('grid_search_hr-210')
% load('grid_search_hc-210')
% load('grid_search_lr-210')
% load('grid_search_nl-210')
% load('grid_search_d3-210');y = dx;


all_check_err = [];
test_index = 1:size(A,1);
threshold = 1e-6;
for para_i = 1:length(check)
    
    sys = check{para_i}.sys;
    
    
    allerror = [];
    
    for i = 1:size(sys,2)
        allerror(:,i) = abs(y(test_index) - A(test_index,:)*sys(:,i));
    end
    minerror_pos = [];
    
    for i =1:length(test_index)
        
        min_tmp = find(allerror(i,:)==min(allerror(i,:)));
        minerror_pos(i) = min_tmp(1);
        
    end
    sys_tmp = reshape(sys,size(sys,1)*size(sys,2),1);
    
    nonzeros(para_i) = length(find(abs(sys_tmp)>threshold));
    
    sum_tmp_err(para_i) = norm( min(allerror,[],2));
    
    all_check_err(para_i) = 2*sum_tmp_err(para_i)+2*nonzeros(para_i);
    
end
idx = find(all_check_err==min(all_check_err))
all_check_err(idx)
look = check{idx(1)}