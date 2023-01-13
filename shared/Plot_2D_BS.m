function obj = Plot_2D_BS( vertex_symbol_array, k_grid_num, k_point_repeatind, ...
                           vertex_ind_Xcoord, Freq_array,  hax )
                           
% Plot Band Structure of 2D PhC
% vertex_ind_Xcoord is used for locating symbols of BZ_vertex in Xcoord
    
    numofobj = size(Freq_array,1);
    obj = cell(numofobj,1);
    %Set X_coord  
    part_unit = 1.0/k_grid_num;
    X_coord = vertex_ind_Xcoord(1) : vertex_ind_Xcoord(end);
    X_coord = (X_coord-1)*part_unit;
  
    if length(X_coord) ~= length(k_point_repeatind)
        error('length of k_point_repeatind differs from X_coord');
    end

    for i = 1 : length(vertex_symbol_array)
    if vertex_symbol_array{i} == 'G'
        vertex_symbol_array{i} = 'Γ';
    end
    end

    % axes setting
    hold( hax, 'on');
    % Plot Band structure
    for ii = 1 : numofobj
        obj = plot( hax, X_coord, Freq_array(ii,k_point_repeatind), '-', 'Color', ...
                        [0 0.4470 0.7410], 'MarkerSize', 5 );
    end
       
    maxfreq = max(Freq_array(:));
    [minX, maxX] = bounds(X_coord);     
    % Set label
    set(hax,'FontSize'  , 12 );
    set(hax,'XTick'     , X_coord(vertex_ind_Xcoord) );
    set(hax,'XTickLabel', vertex_symbol_array );  
    if ~isscalar(X_coord)
        axis( hax, [ minX, maxX, 0, 1.05*maxfreq ] );
    end
    grid( hax, 'on');
    xlabel(hax,'wave vectors $\mathbf{k}$','interpreter','latex');
    ylabel(hax,'frequency($\omega/2\pi$)','interpreter','latex');
   
end