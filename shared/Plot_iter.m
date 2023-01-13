function Plot_iter ( vertex_symbol_array, k_grid_num, k_point_repeatind, ...
                           vertex_ind_Xcoord, iternum_te, iternum_tm )
                           
% Plot Band Structure of 2D PhC
% vertex_ind_Xcoord is used for locating symbols of BZ_vertex in Xcoord
    

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

    maxfreq = max(iternum_te(:));
    [minX, maxX] = bounds(X_coord); 

    % Plot Band structure
    
    % subplot1
    subplot(2,1,1)
    plot(  X_coord, iternum_te(1,k_point_repeatind), '-', 'Color', ...
                        [0 0.4470 0.7410], 'MarkerSize', 5 );
    % Set label 1
    set(gca, 'FontSize', 12);
    set(gca,'XTick'     , X_coord(vertex_ind_Xcoord) );
    set(gca,'XTickLabel', vertex_symbol_array );  
    if ~isscalar(X_coord)
        axis( gca, [ minX, maxX, 0, 1.05*maxfreq ] );
    end
    grid on
    ylim([20,60])
    xlabel('wave vectors $\mathbf{k}$','interpreter','latex');
    ylabel('iter(pcg) TE','interpreter','latex');

    % subplot2
    subplot(2,1,2)
    plot(  X_coord, iternum_tm(1,k_point_repeatind), '-', 'Color', ...
                        [0 0.4470 0.7410], 'MarkerSize', 5 );

     % Set label 2
    set(gca, 'FontSize', 12);
    set(gca,'XTick'     , X_coord(vertex_ind_Xcoord) );
    set(gca,'XTickLabel', vertex_symbol_array );  
    if ~isscalar(X_coord)
        axis( gca, [ minX, maxX, 0, 1.05*maxfreq ] );
    end
    grid on
    ylim([20,60])
    xlabel('wave vectors $\mathbf{k}$','interpreter','latex');
    ylabel('iter(pcg) TM','interpreter','latex');
       
       
   

end