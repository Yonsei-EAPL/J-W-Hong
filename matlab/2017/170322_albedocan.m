%% input description
% 1 : H/W (HWR)
% 2 : albedo road (alb_r)
% 3 : albedo wall (alb_w)
% 4 : zenith angle (in degree)
% 5 : diffused radiation fraction (chi_f)
%% output description (from input 6-8 column)
% 6 : shadowed fraction road (chi_r)
% 7 : shadowed fraction wall (chi_w)
% 8 : canyon albedo (albcan)
%% run information
[size_n, size_var] = size(input);
clear size_var
%% main 
for i = 1:size_n
    HWR = input(i,1);
    alb_r = input(i,2);
    alb_w = input(i,3);
    zenith = input(i,4);
    if zenith ==0
        zenith = 0.1;
    elseif zenith ==90
        zenith = 89.9;
    end
    zenith = zenith/180*pi();
    chi_f = input(i,5);
    
    omega0 = 1/(tan(zenith)*HWR);
    omega0 = min(omega0, 1);
    omega0 = asin(omega0);
    chi_r = (2/pi())*(omega0 - HWR*tan(zenith)*(1 - cos(omega0)));
    chi_w = (1.0 -  chi_r)/(2*HWR);
    
    aa = (1 + HWR^2)^(0.5) - HWR;
    bb = (1 + (1/HWR)^2)^(0.5) - (1/HWR);
    cc = 1 - aa;
    dd = (1 - bb)/2;
    
    % shape factor
    sf(1,1) = 0;
    sf(1,2) = cc;
    sf(1,3) = aa;
    sf(2,1) = dd;
    sf(2,2) = bb;
    sf(2,3) = dd;
    sf(3,1) = aa;
    sf(3,2) = cc;
    sf(3,3) = 0;
    
    % reflections array
    rf(1,1) = 1.0;
    rf(1,2) = -alb_r*cc;
    rf(1,3) = -alb_r*aa;
    rf(2,1) = -alb_w*dd;
    rf(2,2) = (1 - alb_w*bb);
    rf(2,3) = -alb_w*dd;
    rf(3,1) = 0;
    rf(3,2) = 0;
    rf(3,3) = 1;

    % invert matrix
%     n=3;
%     tol = 10^(-6);
%     indxc = zeros(1,3);
%     indxr = zeros(1,3);
%     ipiv = zeros(1,3);
%     dumr = zeros(1,3);
%     dumc = zeros(3,1);
%     for j = 1:3
%         big = 0;
%         for k = 1:3
%             if (ipiv(k)~=1)
%                 for l = 1:3
%                     if (ipiv(l)==0)
%                         if abs(rf(k,l))>=big
%                             big = abs(rf(k,l));
%                             irow = k;
%                             icol = l;
%                         end
%                     end
%                 end
%             end
%         end
%         ipiv(icol) = ipiv(icol) + 1;
%         if (irow~= icol)
%             for l =1:3
%                 dumr(1,l) = rf(irow,l);
%                 rf(irow,l) = rf(icol,l);
%                 rf(icol,l) = dumr(1,l);
%             end
%         end
%         indxr(j) = irow;
%         indxc(j) = icol;
%         pivinv = 1/rf(icol,icol);
%         rf(icol,icol) = 1;
%         for l = 1:3
%             rf(icol,l) = rf(icol,l)*pivinv;
%         end
%         for k =1:3
%             if k~=icol
%                 dum = rf(k,icol);
%                 rf(k,icol)=0;
%                 for l = 1:3
%                     rf(k,l) = rf(k,l) - rf(icol,l)*dum;
%                 end
%             end
%         end
%     end
%     for j=3:-1:1
%         if indxr(j)~=indxc(j)
%             for k = 1:3
%                 dumc(k,1) = rf(k,indxr(j));
%                 rf(k,indxr(j))=rf(k,indxc(j));
%                 rf(k,indxc(j))=dumc(k,1);
%             end
%         end
%     end
    rf = rf^(-1);
%     clear n tol indxc indxr ipiv dumr dumc big j k l pivinv
    
    emis(1) = alb_r*(1 - chi_f)*chi_r;
    emis(2) = alb_w*(1 - chi_f)*chi_w;
    emis(3) = chi_f;

    out = zeros(1,3);
    inn = zeros(1,3);
    for j = 1:3
        for k = 1:3
            out(j) = out(j) + rf(j,k)*emis(k);
        end 
    end
    for j = 1:3
        for k = 1:3
            inn(j) = inn(j) + sf(j,k)*out(k);
        end
    end
    clear j k
    inn(1) = inn(1) + (1-chi_f)*chi_r - out(1);
    inn(2) = inn(2) + (1-chi_f)*chi_w - out(2);
    
    inn(1) = inn(1) + 2*HWR*inn(2);
    if (inn(1)>=0)&&(inn(1)<1)
        albcan = 1- inn(1);
    else
        albcan = -999;
    end
    
    input(i,6) = chi_r;
    input(i,7) = chi_w;
    input(i,8) = albcan;    
   clear HWR alb_r alb_w zenith chi_f omega0 chi_r chi_w albcan aa bb cc dd sf rf emis out inn 
end
clear i size_n


