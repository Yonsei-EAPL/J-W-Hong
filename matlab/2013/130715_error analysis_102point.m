i = 102;
temp = zeros(num_30min_n(i,2),8);
    for j = 1:num_30min_n(i,2)
        if i ==1
            temp(j,1) = data(j,po_u);
            temp(j,2) = data(j,po_v);
            temp(j,3) = data(j,po_w);
            temp(j,4) = data(j,po_Ts);
            temp(j,5) = data(j,po_CO2);
            temp(j,6) = data(j,po_H2O);
            temp(j,7) = data(j,po_cell_tmpr);
            temp(j,8) = data(j,po_cell_prs);
        else
            temp(j,1) = data(num_30min_n(i-1,3)+j,po_u);
            temp(j,2) = data(num_30min_n(i-1,3)+j,po_v);
            temp(j,3) = data(num_30min_n(i-1,3)+j,po_w);
            temp(j,4) = data(num_30min_n(i-1,3)+j,po_Ts);
            temp(j,5) = data(num_30min_n(i-1,3)+j,po_CO2);
            temp(j,6) = data(num_30min_n(i-1,3)+j,po_H2O);
            temp(j,7) = data(num_30min_n(i-1,3)+j,po_cell_tmpr);
            temp(j,8) = data(num_30min_n(i-1,3)+j,po_cell_prs);
        end
    end
    clear j