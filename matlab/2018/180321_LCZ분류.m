% for all
% x=zeros(1,4)
% x(:,1) : site number
% x(:,2) : wind direction
% x(:,3) : empty (zeros)
% x(:,4) : YYYYMMDDHH

for i = 1:length(x)
    if x(i,2)<0
        x(i,3)=888.8;
    elseif x(i,2)>360
        x(i,3)=888.8;
    else
        if x(i,1)==400
            if x(i,2)<270
                x(i,3)=18;
            elseif x(i,2)<315
                x(i,3)=50;
            else
                x(i,3)=52;
            end
        end

        if x(i,1)==401
            if x(i,4)<2016091199
                if x(i,2)<45
                    x(i,3)=69;
                elseif x(i,2)<90
                    x(i,3)=18;
                elseif x(i,2)<135
                    x(i,3)=18;
                elseif x(i,2)<180
                    x(i,3)=18;
                elseif x(i,2)<225
                    x(i,3)=1;
                elseif x(i,2)<270
                    x(i,3)=18;
                elseif x(i,2)<315
                    x(i,3)=18;
                else
                    x(i,3)=69;
                end
            else
                if x(i,2)<45
                    x(i,3)=58;
                elseif x(i,2)<90
                    x(i,3)=58;
                elseif x(i,2)<135
                    x(i,3)=171;
                elseif x(i,2)<180
                    x(i,3)=171;
                elseif x(i,2)<225
                    x(i,3)=174;
                elseif x(i,2)<270
                    x(i,3)=174;
                elseif x(i,2)<315
                    x(i,3)=174;
                else
                    x(i,3)=174;
                end
            end            
        end

        if x(i,1)==402
            if x(i,4)<2004010000
                if x(i,2)<45
                    x(i,3)=57;
                elseif x(i,2)<90
                    x(i,3)=137;
                elseif x(i,2)<135
                    x(i,3)=11;
                elseif x(i,2)<180
                    x(i,3)=1;
                elseif x(i,2)<225
                    x(i,3)=2;
                elseif x(i,2)<270
                    x(i,3)=17;
                elseif x(i,2)<315
                    x(i,3)=1;
                else
                    x(i,3)=1;
                end
            elseif x(i,4)<2008083200
                if x(i,2)<45
                    x(i,3)=57;
                elseif x(i,2)<90
                    x(i,3)=137;
                elseif x(i,2)<135
                    x(i,3)=11;
                elseif x(i,2)<180
                    x(i,3)=1;
                elseif x(i,2)<225
                    x(i,3)=2;
                elseif x(i,2)<270
                    x(i,3)=223;
                elseif x(i,2)<315
                    x(i,3)=223;
                else
                    x(i,3)=223;
                end
            else
                if x(i,2)<45
                    x(i,3)=57;
                elseif x(i,2)<90
                    x(i,3)=137;
                elseif x(i,2)<135
                    x(i,3)=11;
                elseif x(i,2)<180
                    x(i,3)=1;
                elseif x(i,2)<225
                    x(i,3)=2;
                elseif x(i,2)<270
                    x(i,3)=20;
                elseif x(i,2)<315
                    x(i,3)=52;
                else
                    x(i,3)=52;
                end
            end
        end

        if x(i,1)==403
            if (x(i,2)>=180)&&(x(i,2)<225)
                if (x(i,4)>2003010000)&&(x(i,4)<2006010000)
                    x(i,3)=223;
                else
                    x(i,3)=1;
                end
            else
                x(i,3)=1;
            end
        end

        if x(i,1)==404
            if x(i,4)<2005010000
                if x(i,2)<45
                    x(i,3)=154;
                elseif x(i,2)<90
                    x(i,3)=154;
                elseif x(i,2)<135
                    x(i,3)=69;
                elseif x(i,2)<180
                    x(i,3)=21;
                elseif x(i,2)<225
                    x(i,3)=5;
                elseif x(i,2)<270
                    x(i,3)=5;
                elseif x(i,2)<315
                    x(i,3)=154;
                else
                    x(i,3)=154;
                end
            elseif x(i,4)<2007010000
                if x(i,2)<45
                    x(i,3)=154;
                elseif x(i,2)<90
                    x(i,3)=154;
                elseif x(i,2)<135
                    x(i,3)=69;
                elseif x(i,2)<180
                    x(i,3)=21;
                elseif x(i,2)<225
                    x(i,3)=229;
                elseif x(i,2)<270
                    x(i,3)=229;
                elseif x(i,2)<315
                    x(i,3)=154;
                else
                    x(i,3)=154;
                end
            elseif x(i,4)<2017092900
                if x(i,2)<45
                    x(i,3)=154;
                elseif x(i,2)<90
                    x(i,3)=154;
                elseif x(i,2)<135
                    x(i,3)=69;
                elseif x(i,2)<180
                    x(i,3)=21;
                elseif x(i,2)<225
                    x(i,3)=5;
                elseif x(i,2)<270
                    x(i,3)=5;
                elseif x(i,2)<315
                    x(i,3)=154;
                else
                    x(i,3)=154;
                end
            else
                if x(i,2)<45
                    x(i,3)=120;
                elseif x(i,2)<90
                    x(i,3)=120;
                elseif x(i,2)<135
                    x(i,3)=120;
                elseif x(i,2)<180
                    x(i,3)=120;
                elseif x(i,2)<225
                    x(i,3)=120;
                elseif x(i,2)<270
                    x(i,3)=120;
                elseif x(i,2)<315
                    x(i,3)=120;
                else
                    x(i,3)=120;
                end
            end
        end

        if x(i,1)==405
            if x(i,2)<45
                x(i,3)=222;
            elseif x(i,2)<90
                x(i,3)=222;
            elseif x(i,2)<135
                x(i,3)=120;
            elseif x(i,2)<180
                x(i,3)=120;
            elseif x(i,2)<225
                x(i,3)=222;
            elseif x(i,2)<270
                x(i,3)=222;
            elseif x(i,2)<315
                x(i,3)=222;
            else
                x(i,3)=222;
            end
        end

        if x(i,1)==406
            if x(i,2)<45
                x(i,3)=26;
            elseif x(i,2)<90
                x(i,3)=26;
            elseif x(i,2)<135
                x(i,3)=26;
            elseif x(i,2)<180
                x(i,3)=18;
            elseif x(i,2)<225
                x(i,3)=26;
            elseif x(i,2)<270
                x(i,3)=26;
            elseif x(i,2)<315
                x(i,3)=154;
            else
                x(i,3)=154;
            end
        end

        if x(i,1)==407
            if x(i,2)<45
                x(i,3)=125;
            elseif x(i,2)<90
                x(i,3)=205;
            elseif x(i,2)<135
                x(i,3)=205;
            elseif x(i,2)<180
                x(i,3)=69;
            elseif x(i,2)<225
                x(i,3)=52;
            elseif x(i,2)<270
                x(i,3)=69;
            elseif x(i,2)<315
                x(i,3)=171;
            else
                x(i,3)=168;
            end
        end

        if x(i,1)==408
            if x(i,2)<45
                x(i,3)=18;
            elseif x(i,2)<90
                x(i,3)=69;
            elseif x(i,2)<135
                x(i,3)=69;
            elseif x(i,2)<180
                x(i,3)=69;
            elseif x(i,2)<225
                x(i,3)=69;
            elseif x(i,2)<270
                x(i,3)=69;
            elseif x(i,2)<315
                x(i,3)=69;
            else
                x(i,3)=18;
            end
        end

        if x(i,1)==409
            x(i,3)=18;
        end

        if x(i,1)==410
            if x(i,4)<1996070000
                if x(i,2)<45
                    x(i,3)=157;
                elseif x(i,2)<90
                    x(i,3)=157;
                elseif x(i,2)<135
                    x(i,3)=157;
                elseif x(i,2)<180
                    x(i,3)=157;
                elseif x(i,2)<225
                    x(i,3)=149;
                elseif x(i,2)<270
                    x(i,3)=165;
                elseif x(i,2)<315
                    x(i,3)=157;
                else
                    x(i,3)=157;
                end
            elseif x(i,4)>1998070000
                if x(i,2)<45
                    x(i,3)=157;
                elseif x(i,2)<90
                    x(i,3)=157;
                elseif x(i,2)<135
                    x(i,3)=157;
                elseif x(i,2)<180
                    x(i,3)=157;
                elseif x(i,2)<225
                    x(i,3)=219;
                elseif x(i,2)<270
                    x(i,3)=14;
                elseif x(i,2)<315
                    x(i,3)=14;
                else
                    x(i,3)=222;
                end
            else
                if x(i,2)<45
                    x(i,3)=157;
                elseif x(i,2)<90
                    x(i,3)=157;
                elseif x(i,2)<135
                    x(i,3)=157;
                elseif x(i,2)<180
                    x(i,3)=157;
                elseif x(i,2)<225
                    x(i,3)=159;
                elseif x(i,2)<270
                    x(i,3)=223;
                elseif x(i,2)<315
                    x(i,3)=223;
                else
                    x(i,3)=159;
                end
            end
        end

        if x(i,1)==411
            if x(i,4)<2001010000
                x(i,3)=18;
            elseif x(i,4)<2004010000
                if x(i,2)<180
                    x(i,3)=18;
                else
                    x(i,3)=223;
                end
            elseif x(i,4)<2010081599
                if x(i,2)<180
                    x(i,3)=18;
                elseif x(i,2)<315
                    x(i,3)=1;
                else
                    x(i,3)=18;
                end
            else
                if x(i,2)<45
                    x(i,3)=138;
                elseif x(i,2)<90
                    x(i,3)=138;
                elseif x(i,2)<135
                    x(i,3)=154;
                elseif x(i,2)<180
                    x(i,3)=170;
                elseif x(i,2)<225
                    x(i,3)=219;
                elseif x(i,2)<270
                    x(i,3)=171;
                elseif x(i,2)<315
                    x(i,3)=154;
                else
                    x(i,3)=154;
                end
            end
        end

        if x(i,1)==412
            if x(i,2)<45
                x(i,3)=138;
            elseif x(i,2)<90
                x(i,3)=154;
            elseif x(i,2)<135
                x(i,3)=154;
            elseif x(i,2)<180
                x(i,3)=154;
            elseif x(i,2)<225
                x(i,3)=154;
            elseif x(i,2)<270
                x(i,3)=138;
            elseif x(i,2)<315
                x(i,3)=138;
            else
                x(i,3)=154;
            end
        end

        if x(i,1)==413
            if x(i,2)<45
                x(i,3)=143;
            elseif x(i,2)<90
                x(i,3)=143;
            elseif x(i,2)<135
                x(i,3)=143;
            elseif x(i,2)<180
                x(i,3)=143;
            elseif x(i,2)<225
                x(i,3)=223;
            elseif x(i,2)<270
                x(i,3)=223;
            elseif x(i,2)<315
                x(i,3)=175;
            else
                x(i,3)=175;
            end
        end

        if x(i,1)==414
            if x(i,2)<45
                x(i,3)=149;
            elseif x(i,2)<90
                x(i,3)=69;
            elseif x(i,2)<135
                x(i,3)=69;
            elseif x(i,2)<180
                x(i,3)=69;
            elseif x(i,2)<225
                x(i,3)=69;
            elseif x(i,2)<270
                x(i,3)=69;
            elseif x(i,2)<315
                x(i,3)=69;
            else
                x(i,3)=90;
            end
        end

        if x(i,1)==415
            if x(i,4)<1998083200
                if x(i,2)<45
                    x(i,3)=223;
                elseif x(i,2)<90
                    x(i,3)=223;
                elseif x(i,2)<135
                    x(i,3)=223;
                elseif x(i,2)<180
                    x(i,3)=15;
                elseif x(i,2)<225
                    x(i,3)=15;
                elseif x(i,2)<270
                    x(i,3)=223;
                elseif x(i,2)<315
                    x(i,3)=226;
                else
                    x(i,3)=18;
                end
            else
                if x(i,2)<45
                    x(i,3)=165;
                elseif x(i,2)<90
                    x(i,3)=5;
                elseif x(i,2)<135
                    x(i,3)=5;
                elseif x(i,2)<180
                    x(i,3)=5;
                elseif x(i,2)<225
                    x(i,3)=18;
                elseif x(i,2)<270
                    x(i,3)=17;
                elseif x(i,2)<315
                    x(i,3)=1;
                else
                    x(i,3)=5;
                end
            end
        end

        if x(i,1)==416
            x(i,3)=143;
        end

        if x(i,1)==417
            x(i,3)=18;
        end

        if x(i,1)==418
            if x(i,2)<45
                x(i,3)=256;
            elseif x(i,2)<90
                x(i,3)=256;
            elseif x(i,2)<135
                x(i,3)=256;
            elseif x(i,2)<180
                x(i,3)=254;
            elseif x(i,2)<225
                x(i,3)=206;
            elseif x(i,2)<270
                x(i,3)=206;
            elseif x(i,2)<315
                x(i,3)=206;
            else
                x(i,3)=256;
            end
        end

        if x(i,1)==419
            x(i,3)=154;
        end

        if x(i,1)==421
            if x(i,2)<45
                x(i,3)=146;
            elseif x(i,2)<90
                x(i,3)=18;
            elseif x(i,2)<135
                x(i,3)=146;
            elseif x(i,2)<180
                x(i,3)=154;
            elseif x(i,2)<225
                x(i,3)=157;
            elseif x(i,2)<270
                x(i,3)=157;
            elseif x(i,2)<315
                x(i,3)=213;
            else
                x(i,3)=213;
            end
        end

        if x(i,1)==422
            x(i,3)=154;
        end

        if x(i,1)==423
            if x(i,2)<45
                x(i,3)=69;
            elseif x(i,2)<90
                x(i,3)=31;
            elseif x(i,2)<135
                x(i,3)=31;
            elseif x(i,2)<180
                x(i,3)=31;
            elseif x(i,2)<225
                x(i,3)=31;
            elseif x(i,2)<270
                x(i,3)=18;
            elseif x(i,2)<315
                x(i,3)=146;
            else
                x(i,3)=69;
            end
        end

        if x(i,1)==424
            x(i,3)=18;
        end

        if x(i,1)==425
            if x(i,2)<45
                x(i,3)=223;
            elseif x(i,2)<90
                x(i,3)=223;
            elseif x(i,2)<135
                x(i,3)=69;
            elseif x(i,2)<180
                x(i,3)=69;
            elseif x(i,2)<225
                x(i,3)=153;
            elseif x(i,2)<270
                x(i,3)=153;
            elseif x(i,2)<315
                x(i,3)=79;
            else
                x(i,3)=79;
            end
        end

        if x(i,1)==509
            if x(i,4)<2004090500
                if x(i,2)<45
                    x(i,3)=69;
                elseif x(i,2)<90
                    x(i,3)=69;
                elseif x(i,2)<135
                    x(i,3)=69;
                elseif x(i,2)<180
                    x(i,3)=69;
                elseif x(i,2)<225
                    x(i,3)=149;
                elseif x(i,2)<270
                    x(i,3)=149;
                elseif x(i,2)<315
                    x(i,3)=69;
                else
                    x(i,3)=69;
                end
            else
                if x(i,2)<45
                    x(i,3)=69;
                elseif x(i,2)<90
                    x(i,3)=74;
                else
                    x(i,3)=154;
                end
            end   
        end

        if x(i,1)==510
            x(i,3)=18;
        end

        if x(i,1)==889
            if x(i,2)<45
                x(i,3)=174;
            elseif x(i,2)<90
                x(i,3)=174;
            elseif x(i,2)<135
                x(i,3)=174;
            elseif x(i,2)<180
                x(i,3)=205;
            elseif x(i,2)<225
                x(i,3)=205;
            elseif x(i,2)<270
                x(i,3)=205;
            elseif x(i,2)<315
                x(i,3)=171;
            else
                x(i,3)=171;
            end
        end

    end
end
clear i
