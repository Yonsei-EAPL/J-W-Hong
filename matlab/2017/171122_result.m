% % % % % le_ctl = nc_varget('ctl.30min.nc','latent_heat');
% % % % % h_ctl = nc_varget('ctl.30min.nc','ftl_gb');
% % % % % rn_ctl = nc_varget('ctl.30min.nc','rad_net');
% % % % 
% % % % 
% % % % %% alnir
% % % % 
% % % % le = zeros(17520,6);
% % % % h = zeros(17520,6);
% % % % rn = zeros(17520,6);
% % % % 
% % % % le(:,1) = nc_varget('alnir_m30.30min.nc','latent_heat');
% % % % h(:,1) = nc_varget('alnir_m30.30min.nc','ftl_gb');
% % % % rn(:,1) = nc_varget('alnir_m30.30min.nc','rad_net');
% % % % 
% % % % le(:,2) = nc_varget('alnir_m20.30min.nc','latent_heat');
% % % % h(:,2) = nc_varget('alnir_m20.30min.nc','ftl_gb');
% % % % rn(:,2) = nc_varget('alnir_m20.30min.nc','rad_net');
% % % % 
% % % % le(:,3) = nc_varget('alnir_m10.30min.nc','latent_heat');
% % % % h(:,3) = nc_varget('alnir_m10.30min.nc','ftl_gb');
% % % % rn(:,3) = nc_varget('alnir_m10.30min.nc','rad_net');
% % % % 
% % % % le(:,4) = nc_varget('alnir_p10.30min.nc','latent_heat');
% % % % h(:,4) = nc_varget('alnir_p10.30min.nc','ftl_gb');
% % % % rn(:,4) = nc_varget('alnir_p10.30min.nc','rad_net');
% % % % 
% % % % le(:,5) = nc_varget('alnir_p20.30min.nc','latent_heat');
% % % % h(:,5) = nc_varget('alnir_p20.30min.nc','ftl_gb');
% % % % rn(:,5) = nc_varget('alnir_p20.30min.nc','rad_net');
% % % % 
% % % % le(:,6) = nc_varget('alnir_p30.30min.nc','latent_heat');
% % % % h(:,6) = nc_varget('alnir_p30.30min.nc','ftl_gb');
% % % % rn(:,6) = nc_varget('alnir_p30.30min.nc','rad_net');
% % % % 
% % % % le_alnir = le;
% % % % h_alnir = h;
% % % % rn_alnir = rn;
% % % % clear le h rn
% % % % 
% % % % 
% % % % 
% % % % %% alpar
% % % % 
% % % % le = zeros(17520,6);
% % % % h = zeros(17520,6);
% % % % rn = zeros(17520,6);
% % % % 
% % % % le(:,1) = nc_varget('alpar_m30.30min.nc','latent_heat');
% % % % h(:,1) = nc_varget('alpar_m30.30min.nc','ftl_gb');
% % % % rn(:,1) = nc_varget('alpar_m30.30min.nc','rad_net');
% % % % 
% % % % le(:,2) = nc_varget('alpar_m20.30min.nc','latent_heat');
% % % % h(:,2) = nc_varget('alpar_m20.30min.nc','ftl_gb');
% % % % rn(:,2) = nc_varget('alpar_m20.30min.nc','rad_net');
% % % % 
% % % % le(:,3) = nc_varget('alpar_m10.30min.nc','latent_heat');
% % % % h(:,3) = nc_varget('alpar_m10.30min.nc','ftl_gb');
% % % % rn(:,3) = nc_varget('alpar_m10.30min.nc','rad_net');
% % % % 
% % % % le(:,4) = nc_varget('alpar_p10.30min.nc','latent_heat');
% % % % h(:,4) = nc_varget('alpar_p10.30min.nc','ftl_gb');
% % % % rn(:,4) = nc_varget('alpar_p10.30min.nc','rad_net');
% % % % 
% % % % le(:,5) = nc_varget('alpar_p20.30min.nc','latent_heat');
% % % % h(:,5) = nc_varget('alpar_p20.30min.nc','ftl_gb');
% % % % rn(:,5) = nc_varget('alpar_p20.30min.nc','rad_net');
% % % % 
% % % % le(:,6) = nc_varget('alpar_p30.30min.nc','latent_heat');
% % % % h(:,6) = nc_varget('alpar_p30.30min.nc','ftl_gb');
% % % % rn(:,6) = nc_varget('alpar_p30.30min.nc','rad_net');
% % % % 
% % % % le_alpar = le;
% % % % h_alpar = h;
% % % % rn_alpar = rn;
% % % % clear le h rn
% % % % 
% % % % 
% % % % %% catch0
% % % % 
% % % % le = zeros(17520,6);
% % % % h = zeros(17520,6);
% % % % rn = zeros(17520,6);
% % % % 
% % % % le(:,1) = nc_varget('catch0_m30.30min.nc','latent_heat');
% % % % h(:,1) = nc_varget('catch0_m30.30min.nc','ftl_gb');
% % % % rn(:,1) = nc_varget('catch0_m30.30min.nc','rad_net');
% % % % 
% % % % le(:,2) = nc_varget('catch0_m20.30min.nc','latent_heat');
% % % % h(:,2) = nc_varget('catch0_m20.30min.nc','ftl_gb');
% % % % rn(:,2) = nc_varget('catch0_m20.30min.nc','rad_net');
% % % % 
% % % % le(:,3) = nc_varget('catch0_m10.30min.nc','latent_heat');
% % % % h(:,3) = nc_varget('catch0_m10.30min.nc','ftl_gb');
% % % % rn(:,3) = nc_varget('catch0_m10.30min.nc','rad_net');
% % % % 
% % % % le(:,4) = nc_varget('catch0_p10.30min.nc','latent_heat');
% % % % h(:,4) = nc_varget('catch0_p10.30min.nc','ftl_gb');
% % % % rn(:,4) = nc_varget('catch0_p10.30min.nc','rad_net');
% % % % 
% % % % le(:,5) = nc_varget('catch0_p20.30min.nc','latent_heat');
% % % % h(:,5) = nc_varget('catch0_p20.30min.nc','ftl_gb');
% % % % rn(:,5) = nc_varget('catch0_p20.30min.nc','rad_net');
% % % % 
% % % % le(:,6) = nc_varget('catch0_p30.30min.nc','latent_heat');
% % % % h(:,6) = nc_varget('catch0_p30.30min.nc','ftl_gb');
% % % % rn(:,6) = nc_varget('catch0_p30.30min.nc','rad_net');
% % % % 
% % % % le_catch0 = le;
% % % % h_catch0 = h;
% % % % rn_catch0 = rn;
% % % % clear le h rn
% % % % 
% % % % %% dcatch
% % % % 
% % % % le = zeros(17520,6);
% % % % h = zeros(17520,6);
% % % % rn = zeros(17520,6);
% % % % 
% % % % le(:,1) = nc_varget('dcatch_m30.30min.nc','latent_heat');
% % % % h(:,1) = nc_varget('dcatch_m30.30min.nc','ftl_gb');
% % % % rn(:,1) = nc_varget('dcatch_m30.30min.nc','rad_net');
% % % % 
% % % % le(:,2) = nc_varget('dcatch_m20.30min.nc','latent_heat');
% % % % h(:,2) = nc_varget('dcatch_m20.30min.nc','ftl_gb');
% % % % rn(:,2) = nc_varget('dcatch_m20.30min.nc','rad_net');
% % % % 
% % % % le(:,3) = nc_varget('dcatch_m10.30min.nc','latent_heat');
% % % % h(:,3) = nc_varget('dcatch_m10.30min.nc','ftl_gb');
% % % % rn(:,3) = nc_varget('dcatch_m10.30min.nc','rad_net');
% % % % 
% % % % le(:,4) = nc_varget('dcatch_p10.30min.nc','latent_heat');
% % % % h(:,4) = nc_varget('dcatch_p10.30min.nc','ftl_gb');
% % % % rn(:,4) = nc_varget('dcatch_p10.30min.nc','rad_net');
% % % % 
% % % % le(:,5) = nc_varget('dcatch_p20.30min.nc','latent_heat');
% % % % h(:,5) = nc_varget('dcatch_p20.30min.nc','ftl_gb');
% % % % rn(:,5) = nc_varget('dcatch_p20.30min.nc','rad_net');
% % % % 
% % % % le(:,6) = nc_varget('dcatch_p30.30min.nc','latent_heat');
% % % % h(:,6) = nc_varget('dcatch_p30.30min.nc','ftl_gb');
% % % % rn(:,6) = nc_varget('dcatch_p30.30min.nc','rad_net');
% % % % 
% % % % le_dcatch = le;
% % % % h_dcatch = h;
% % % % rn_dcatch = rn;
% % % % clear le h rn
% % % % 
% % % % %% dz0v
% % % % 
% % % % le = zeros(17520,6);
% % % % h = zeros(17520,6);
% % % % rn = zeros(17520,6);
% % % % 
% % % % le(:,1) = nc_varget('dz0v_m30.30min.nc','latent_heat');
% % % % h(:,1) = nc_varget('dz0v_m30.30min.nc','ftl_gb');
% % % % rn(:,1) = nc_varget('dz0v_m30.30min.nc','rad_net');
% % % % 
% % % % le(:,2) = nc_varget('dz0v_m20.30min.nc','latent_heat');
% % % % h(:,2) = nc_varget('dz0v_m20.30min.nc','ftl_gb');
% % % % rn(:,2) = nc_varget('dz0v_m20.30min.nc','rad_net');
% % % % 
% % % % le(:,3) = nc_varget('dz0v_m10.30min.nc','latent_heat');
% % % % h(:,3) = nc_varget('dz0v_m10.30min.nc','ftl_gb');
% % % % rn(:,3) = nc_varget('dz0v_m10.30min.nc','rad_net');
% % % % 
% % % % le(:,4) = nc_varget('dz0v_p10.30min.nc','latent_heat');
% % % % h(:,4) = nc_varget('dz0v_p10.30min.nc','ftl_gb');
% % % % rn(:,4) = nc_varget('dz0v_p10.30min.nc','rad_net');
% % % % 
% % % % le(:,5) = nc_varget('dz0v_p20.30min.nc','latent_heat');
% % % % h(:,5) = nc_varget('dz0v_p20.30min.nc','ftl_gb');
% % % % rn(:,5) = nc_varget('dz0v_p20.30min.nc','rad_net');
% % % % 
% % % % le(:,6) = nc_varget('dz0v_p30.30min.nc','latent_heat');
% % % % h(:,6) = nc_varget('dz0v_p30.30min.nc','ftl_gb');
% % % % rn(:,6) = nc_varget('dz0v_p30.30min.nc','rad_net');
% % % % 
% % % % le_dz0v = le;
% % % % h_dz0v = h;
% % % % rn_dz0v = rn;
% % % % clear le h rn
% % % % 
% % % % 
% % % % %% emis
% % % % 
% % % % le = zeros(17520,6);
% % % % h = zeros(17520,6);
% % % % rn = zeros(17520,6);
% % % % 
% % % % le(:,1) = nc_varget('emis_m30.30min.nc','latent_heat');
% % % % h(:,1) = nc_varget('emis_m30.30min.nc','ftl_gb');
% % % % rn(:,1) = nc_varget('emis_m30.30min.nc','rad_net');
% % % % 
% % % % le(:,2) = nc_varget('emis_m20.30min.nc','latent_heat');
% % % % h(:,2) = nc_varget('emis_m20.30min.nc','ftl_gb');
% % % % rn(:,2) = nc_varget('emis_m20.30min.nc','rad_net');
% % % % 
% % % % le(:,3) = nc_varget('emis_m10.30min.nc','latent_heat');
% % % % h(:,3) = nc_varget('emis_m10.30min.nc','ftl_gb');
% % % % rn(:,3) = nc_varget('emis_m10.30min.nc','rad_net');
% % % % 
% % % % le(:,4) = nc_varget('emis_p10.30min.nc','latent_heat');
% % % % h(:,4) = nc_varget('emis_p10.30min.nc','ftl_gb');
% % % % rn(:,4) = nc_varget('emis_p10.30min.nc','rad_net');
% % % % 
% % % % le(:,5) = nc_varget('emis_p20.30min.nc','latent_heat');
% % % % h(:,5) = nc_varget('emis_p20.30min.nc','ftl_gb');
% % % % rn(:,5) = nc_varget('emis_p20.30min.nc','rad_net');
% % % % 
% % % % le(:,6) = nc_varget('emis_p30.30min.nc','latent_heat');
% % % % h(:,6) = nc_varget('emis_p30.30min.nc','ftl_gb');
% % % % rn(:,6) = nc_varget('emis_p30.30min.nc','rad_net');
% % % % 
% % % % le_emis = le;
% % % % h_emis = h;
% % % % rn_emis = rn;
% % % % clear le h rn
% % % % 
% % % % 
% % % % %% glmin
% % % % 
% % % % le = zeros(17520,6);
% % % % h = zeros(17520,6);
% % % % rn = zeros(17520,6);
% % % % 
% % % % le(:,1) = nc_varget('glmin_m30.30min.nc','latent_heat');
% % % % h(:,1) = nc_varget('glmin_m30.30min.nc','ftl_gb');
% % % % rn(:,1) = nc_varget('glmin_m30.30min.nc','rad_net');
% % % % 
% % % % le(:,2) = nc_varget('glmin_m20.30min.nc','latent_heat');
% % % % h(:,2) = nc_varget('glmin_m20.30min.nc','ftl_gb');
% % % % rn(:,2) = nc_varget('glmin_m20.30min.nc','rad_net');
% % % % 
% % % % le(:,3) = nc_varget('glmin_m10.30min.nc','latent_heat');
% % % % h(:,3) = nc_varget('glmin_m10.30min.nc','ftl_gb');
% % % % rn(:,3) = nc_varget('glmin_m10.30min.nc','rad_net');
% % % % 
% % % % le(:,4) = nc_varget('glmin_p10.30min.nc','latent_heat');
% % % % h(:,4) = nc_varget('glmin_p10.30min.nc','ftl_gb');
% % % % rn(:,4) = nc_varget('glmin_p10.30min.nc','rad_net');
% % % % 
% % % % le(:,5) = nc_varget('glmin_p20.30min.nc','latent_heat');
% % % % h(:,5) = nc_varget('glmin_p20.30min.nc','ftl_gb');
% % % % rn(:,5) = nc_varget('glmin_p20.30min.nc','rad_net');
% % % % 
% % % % le(:,6) = nc_varget('glmin_p30.30min.nc','latent_heat');
% % % % h(:,6) = nc_varget('glmin_p30.30min.nc','ftl_gb');
% % % % rn(:,6) = nc_varget('glmin_p30.30min.nc','rad_net');
% % % % 
% % % % le_glmin = le;
% % % % h_glmin = h;
% % % % rn_glmin = rn;
% % % % clear le h rn
% % % % 
% % % % %% lai
% % % % 
% % % % le = zeros(17520,6);
% % % % h = zeros(17520,6);
% % % % rn = zeros(17520,6);
% % % % 
% % % % le(:,1) = nc_varget('lai_m30.30min.nc','latent_heat');
% % % % h(:,1) = nc_varget('lai_m30.30min.nc','ftl_gb');
% % % % rn(:,1) = nc_varget('lai_m30.30min.nc','rad_net');
% % % % 
% % % % le(:,2) = nc_varget('lai_m20.30min.nc','latent_heat');
% % % % h(:,2) = nc_varget('lai_m20.30min.nc','ftl_gb');
% % % % rn(:,2) = nc_varget('lai_m20.30min.nc','rad_net');
% % % % 
% % % % le(:,3) = nc_varget('lai_m10.30min.nc','latent_heat');
% % % % h(:,3) = nc_varget('lai_m10.30min.nc','ftl_gb');
% % % % rn(:,3) = nc_varget('lai_m10.30min.nc','rad_net');
% % % % 
% % % % le(:,4) = nc_varget('lai_p10.30min.nc','latent_heat');
% % % % h(:,4) = nc_varget('lai_p10.30min.nc','ftl_gb');
% % % % rn(:,4) = nc_varget('lai_p10.30min.nc','rad_net');
% % % % 
% % % % le(:,5) = nc_varget('lai_p20.30min.nc','latent_heat');
% % % % h(:,5) = nc_varget('lai_p20.30min.nc','ftl_gb');
% % % % rn(:,5) = nc_varget('lai_p20.30min.nc','rad_net');
% % % % 
% % % % le(:,6) = nc_varget('lai_p30.30min.nc','latent_heat');
% % % % h(:,6) = nc_varget('lai_p30.30min.nc','ftl_gb');
% % % % rn(:,6) = nc_varget('lai_p30.30min.nc','rad_net');
% % % % 
% % % % le_lai = le;
% % % % h_lai = h;
% % % % rn_lai = rn;
% % % % clear le h rn
% % % % 
% % % % 
% % % % %% omnir
% % % % 
% % % % le = zeros(17520,6);
% % % % h = zeros(17520,6);
% % % % rn = zeros(17520,6);
% % % % 
% % % % le(:,1) = nc_varget('omnir_m30.30min.nc','latent_heat');
% % % % h(:,1) = nc_varget('omnir_m30.30min.nc','ftl_gb');
% % % % rn(:,1) = nc_varget('omnir_m30.30min.nc','rad_net');
% % % % 
% % % % le(:,2) = nc_varget('omnir_m20.30min.nc','latent_heat');
% % % % h(:,2) = nc_varget('omnir_m20.30min.nc','ftl_gb');
% % % % rn(:,2) = nc_varget('omnir_m20.30min.nc','rad_net');
% % % % 
% % % % le(:,3) = nc_varget('omnir_m10.30min.nc','latent_heat');
% % % % h(:,3) = nc_varget('omnir_m10.30min.nc','ftl_gb');
% % % % rn(:,3) = nc_varget('omnir_m10.30min.nc','rad_net');
% % % % 
% % % % le(:,4) = nc_varget('omnir_p10.30min.nc','latent_heat');
% % % % h(:,4) = nc_varget('omnir_p10.30min.nc','ftl_gb');
% % % % rn(:,4) = nc_varget('omnir_p10.30min.nc','rad_net');
% % % % 
% % % % le(:,5) = nc_varget('omnir_p20.30min.nc','latent_heat');
% % % % h(:,5) = nc_varget('omnir_p20.30min.nc','ftl_gb');
% % % % rn(:,5) = nc_varget('omnir_p20.30min.nc','rad_net');
% % % % 
% % % % le(:,6) = nc_varget('omnir_p30.30min.nc','latent_heat');
% % % % h(:,6) = nc_varget('omnir_p30.30min.nc','ftl_gb');
% % % % rn(:,6) = nc_varget('omnir_p30.30min.nc','rad_net');
% % % % 
% % % % le_omnir = le;
% % % % h_omnir = h;
% % % % rn_omnir = rn;
% % % % clear le h rn
% % % % 
% % % % %% rootd
% % % % 
% % % % le = zeros(17520,6);
% % % % h = zeros(17520,6);
% % % % rn = zeros(17520,6);
% % % % 
% % % % le(:,1) = nc_varget('rootd_m30.30min.nc','latent_heat');
% % % % h(:,1) = nc_varget('rootd_m30.30min.nc','ftl_gb');
% % % % rn(:,1) = nc_varget('rootd_m30.30min.nc','rad_net');
% % % % 
% % % % le(:,2) = nc_varget('rootd_m20.30min.nc','latent_heat');
% % % % h(:,2) = nc_varget('rootd_m20.30min.nc','ftl_gb');
% % % % rn(:,2) = nc_varget('rootd_m20.30min.nc','rad_net');
% % % % 
% % % % le(:,3) = nc_varget('rootd_m10.30min.nc','latent_heat');
% % % % h(:,3) = nc_varget('rootd_m10.30min.nc','ftl_gb');
% % % % rn(:,3) = nc_varget('rootd_m10.30min.nc','rad_net');
% % % % 
% % % % le(:,4) = nc_varget('rootd_p10.30min.nc','latent_heat');
% % % % h(:,4) = nc_varget('rootd_p10.30min.nc','ftl_gb');
% % % % rn(:,4) = nc_varget('rootd_p10.30min.nc','rad_net');
% % % % 
% % % % le(:,5) = nc_varget('rootd_p20.30min.nc','latent_heat');
% % % % h(:,5) = nc_varget('rootd_p20.30min.nc','ftl_gb');
% % % % rn(:,5) = nc_varget('rootd_p20.30min.nc','rad_net');
% % % % 
% % % % le(:,6) = nc_varget('rootd_p30.30min.nc','latent_heat');
% % % % h(:,6) = nc_varget('rootd_p30.30min.nc','ftl_gb');
% % % % rn(:,6) = nc_varget('rootd_p30.30min.nc','rad_net');
% % % % 
% % % % le_rootd = le;
% % % % h_rootd = h;
% % % % rn_rootd = rn;
% % % % clear le h rn
% % % % 
% % % % 
% % % % %% z0hm
% % % % 
% % % % le = zeros(17520,6);
% % % % h = zeros(17520,6);
% % % % rn = zeros(17520,6);
% % % % 
% % % % le(:,1) = nc_varget('z0hm_m30.30min.nc','latent_heat');
% % % % h(:,1) = nc_varget('z0hm_m30.30min.nc','ftl_gb');
% % % % rn(:,1) = nc_varget('z0hm_m30.30min.nc','rad_net');
% % % % 
% % % % le(:,2) = nc_varget('z0hm_m20.30min.nc','latent_heat');
% % % % h(:,2) = nc_varget('z0hm_m20.30min.nc','ftl_gb');
% % % % rn(:,2) = nc_varget('z0hm_m20.30min.nc','rad_net');
% % % % 
% % % % le(:,3) = nc_varget('z0hm_m10.30min.nc','latent_heat');
% % % % h(:,3) = nc_varget('z0hm_m10.30min.nc','ftl_gb');
% % % % rn(:,3) = nc_varget('z0hm_m10.30min.nc','rad_net');
% % % % 
% % % % le(:,4) = nc_varget('z0hm_p10.30min.nc','latent_heat');
% % % % h(:,4) = nc_varget('z0hm_p10.30min.nc','ftl_gb');
% % % % rn(:,4) = nc_varget('z0hm_p10.30min.nc','rad_net');
% % % % 
% % % % le(:,5) = nc_varget('z0hm_p20.30min.nc','latent_heat');
% % % % h(:,5) = nc_varget('z0hm_p20.30min.nc','ftl_gb');
% % % % rn(:,5) = nc_varget('z0hm_p20.30min.nc','rad_net');
% % % % 
% % % % le(:,6) = nc_varget('z0hm_p30.30min.nc','latent_heat');
% % % % h(:,6) = nc_varget('z0hm_p30.30min.nc','ftl_gb');
% % % % rn(:,6) = nc_varget('z0hm_p30.30min.nc','rad_net');
% % % % 
% % % % le_z0hm = le;
% % % % h_z0hm = h;
% % % % rn_z0hm = rn;
% % % % clear le h rn
% % % % 
% % % % 
% % % alnir = zeros(6,3);
% % % alpar = zeros(6,3);
% % % catch0 = zeros(6,3);
% % % dcatch = zeros(6,3);
% % % dz0v = zeros(6,3);
% % % emis = zeros(6,3);
% % % glmin = zeros(6,3);
% % % lai = zeros(6,3);
% % % omnir = zeros(6,3);
% % % rootd = zeros(6,3);
% % % z0hm = zeros(6,3);
% % % ctl = zeros(1,3);
% % 
% % for i = 1:6
% %     alnir(i,1) = mean(h_alnir(:,i));
% %     alnir(i,2) = mean(le_alnir(:,i));
% %     alnir(i,3) = mean(rn_alnir(:,i));
% % end
% % clear i
% % 
% % for i = 1:6
% %     alpar(i,1) = mean(h_alpar(:,i));
% %     alpar(i,2) = mean(le_alpar(:,i));
% %     alpar(i,3) = mean(rn_alpar(:,i));
% % end
% % clear i
% % 
% % for i = 1:6
% %     catch0(i,1) = mean(h_catch0(:,i));
% %     catch0(i,2) = mean(le_catch0(:,i));
% %     catch0(i,3) = mean(rn_catch0(:,i));
% % end
% % clear i
% % 
% % for i = 1:6
% %     dcatch(i,1) = mean(h_dcatch(:,i));
% %     dcatch(i,2) = mean(le_dcatch(:,i));
% %     dcatch(i,3) = mean(rn_dcatch(:,i));
% % end
% % clear i
% % 
% % for i = 1:6
% %     dz0v(i,1) = mean(h_dz0v(:,i));
% %     dz0v(i,2) = mean(le_dz0v(:,i));
% %     dz0v(i,3) = mean(rn_dz0v(:,i));
% % end
% % clear i
% % 
% % for i = 1:6
% %     emis(i,1) = mean(h_emis(:,i));
% %     emis(i,2) = mean(le_emis(:,i));
% %     emis(i,3) = mean(rn_emis(:,i));
% % end
% % clear i
% % 
% % for i = 1:6
% %     glmin(i,1) = mean(h_glmin(:,i));
% %     glmin(i,2) = mean(le_glmin(:,i));
% %     glmin(i,3) = mean(rn_glmin(:,i));
% % end
% % clear i
% % 
% % for i = 1:6
% %     lai(i,1) = mean(h_lai(:,i));
% %     lai(i,2) = mean(le_lai(:,i));
% %     lai(i,3) = mean(rn_lai(:,i));
% % end
% % clear i
% % 
% % for i = 1:6
% %     omnir(i,1) = mean(h_omnir(:,i));
% %     omnir(i,2) = mean(le_omnir(:,i));
% %     omnir(i,3) = mean(rn_omnir(:,i));
% % end
% % clear i
% % 
% % for i = 1:6
% %     rootd(i,1) = mean(h_rootd(:,i));
% %     rootd(i,2) = mean(le_rootd(:,i));
% %     rootd(i,3) = mean(rn_rootd(:,i));
% % end
% % clear i
% % 
% % for i = 1:6
% %     z0hm(i,1) = mean(h_z0hm(:,i));
% %     z0hm(i,2) = mean(le_z0hm(:,i));
% %     z0hm(i,3) = mean(rn_z0hm(:,i));
% % end
% % clear i
% % 
% % ctl(1,1) = mean(h_ctl);
% % ctl(1,2) = mean(le_ctl);
% % ctl(1,3) = mean(rn_ctl);
% % 
% for i = 1:6
%     for j = 1:3
%         alnir(i,j) = (alnir(i,j)-ctl(1,j))*100/ctl(1,j);
%         alpar(i,j) = (alpar(i,j)-ctl(1,j))*100/ctl(1,j);
%         catch0(i,j) = (catch0(i,j)-ctl(1,j))*100/ctl(1,j);
%         dcatch(i,j) = (dcatch(i,j)-ctl(1,j))*100/ctl(1,j);
%         dz0v(i,j) = (dz0v(i,j)-ctl(1,j))*100/ctl(1,j);
%         emis(i,j) = (emis(i,j)-ctl(1,j))*100/ctl(1,j);
%         glmin(i,j) = (glmin(i,j)-ctl(1,j))*100/ctl(1,j);
%         lai(i,j) = (lai(i,j)-ctl(1,j))*100/ctl(1,j);
%         omnir(i,j) = (omnir(i,j)-ctl(1,j))*100/ctl(1,j);
%         rootd(i,j) = (rootd(i,j)-ctl(1,j))*100/ctl(1,j);
%         z0hm(i,j) = (z0hm(i,j)-ctl(1,j))*100/ctl(1,j);
%     end
% end
% clear i j
% 
% 
% 
