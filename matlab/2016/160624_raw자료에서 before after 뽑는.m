% bs = 504001; % 0700-
% be = 540000;
% as = 612001;
% ae = 648000;

bs = 540001; % 0730-
be = bs+36000-1;
as = be+72001;
ae = as+36000-1;

% bs =576001; % 0800-
% be =612000;
% as = 684001;
% ae = 720000;

% bs =612001; % 0830-
% be =648000;
% as = 720001;
% ae = 756000;

% bs = 648001; % 0900-
% be = 684000;
% as = 756001;
% ae = 792000;


data60_0714_b = data60(bs:be,:);
data60_0714_a = data60(as:ae,:);
clear data60
data140_0714_b = data140(bs:be,:);
data140_0714_a = data140(as:ae,:);
clear data140
data300_0714_b = data300(bs:be,:);
data300_0714_a = data300(as:ae,:);
clear data300

clear bs be as ae
