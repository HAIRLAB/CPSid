function [ output_idx alone ] = find_idx( input_idx,threshold)


    k = 1;
    cluster_num = 1;
    cluster_idx = 1;
    cluster{cluster_num } =[];
    seg_start = input_idx(k);

    
    for i = 1 : size(input_idx,2)
        
        
        
        if (input_idx(i) - seg_start) >= i-k   &&  (input_idx(i) - seg_start)<=i-k+threshold  &&  i-k>=0 
            
            cluster{cluster_num } = [ cluster{cluster_num }    input_idx(i)];
            
            
        else
            
            cluster_size(cluster_num ) = size(cluster{cluster_num },2);

            
            cluster_num = cluster_num+1;
            cluster{cluster_num } =[];
            k = i;
            seg_start = input_idx(k);
            cluster{cluster_num } = [ cluster{cluster_num }    input_idx(i)];
            
            
        end
        
    end
    

    
    cluster_size(cluster_num ) = size(cluster{cluster_num },2);
    alone = max(cluster_size);
    if alone<=2
        alone=1;
    end
    
    max_idx = find(cluster_size==max(cluster_size));
    
    output_idx =cluster{max_idx(1)};

end

