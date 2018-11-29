% data ÇÕÄ¡±â
% 1: site number
% 2: winddirection
% 3: LCZ
% 4: YYYYMMDDHH
% 5: air temperature

load('1995.mat')
data_hour = x;
clear x
load('1996.mat')
data_hour = [data_hour; x];
clear x
load('1997.mat')
data_hour = [data_hour; x];
clear x
load('1998.mat')
data_hour = [data_hour; x];
clear x
load('1999.mat')
data_hour = [data_hour; x];
clear x
load('2000.mat')
data_hour = [data_hour; x];
clear x
load('2001.mat')
data_hour = [data_hour; x];
clear x
load('2002.mat')
data_hour = [data_hour; x];
clear x
load('2003.mat')
data_hour = [data_hour; x];
clear x
load('2004.mat')
data_hour = [data_hour; x];
clear x
load('2005.mat')
data_hour = [data_hour; x];
clear x
load('2006.mat')
data_hour = [data_hour; x];
clear x
load('2007.mat')
data_hour = [data_hour; x];
clear x
load('2008.mat')
data_hour = [data_hour; x];
clear x
load('2009.mat')
data_hour = [data_hour; x];
clear x
load('2010.mat')
data_hour = [data_hour; x];
clear x
load('2011.mat')
data_hour = [data_hour; x];
clear x
load('2012.mat')
data_hour = [data_hour; x];
clear x
load('2013.mat')
data_hour = [data_hour; x];
clear x
load('2014.mat')
data_hour = [data_hour; x];
clear x
load('2015.mat')
data_hour = [data_hour; x];
clear x
load('2016.mat')
data_hour = [data_hour; x];
clear x
load('2017.mat')
data_hour = [data_hour; x];
clear x
save('data_hour.mat','data_hour');


