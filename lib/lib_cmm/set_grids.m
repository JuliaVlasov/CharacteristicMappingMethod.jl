function grid_list = set_grids(params, grid_size_list, grid_names)
    
    grid_list = {};
    for ig = 1:length(grid_size_list)
        
        name = grid_names(ig);
        grid_size = grid_size_list{ig};
        grid.size = grid_size;
        grid.L  = params.L; 
        grid.dx = params.L(1)/grid_size(1);
        grid.dv = params.L(2)/grid_size(2);
        [grid.Xgrid,grid.Vgrid,grid.x,y] = genMeshPoints(params.dom,grid_size);
        [grid.KX, grid.KV]               = genWaveNumbers(params.dom, grid_size);
        grid.kv = grid.KV(:,1)';
        grid.v = y-params.Lv;
        [grid.v_periodic, grid.sigma] = velocity_periodicfication(grid, params.bump_transition_width);
        [Vperiodic,grid.Sigma] = meshgrid(grid.v_periodic,grid.sigma);
        grid.Vperiodic = Vperiodic';
        grid_list.(name) = grid;
    end


end