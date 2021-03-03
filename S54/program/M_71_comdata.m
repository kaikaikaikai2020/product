clear;
%FSTX指数数据，可以做一个爬取方法
addpath(genpath(fullfile(pwd,'jplv7')))
entryZscore=0.1;
symbol = 'FSTX';
%{
data=load('inputDataOHLCDaily_20120517', 'syms', 'tday', 'op', 'hi', 'lo', 'cl');
idx=find(strcmp(symbol, data.syms));
op=data.op(:, idx);
hi=data.hi(:, idx);
lo=data.lo(:, idx);
cl=data.cl(:, idx);
%}
[~,~,data]= xlsread('dataA71.xlsx');
data = data(2:end,:);
tref_num = datenum(data(:,1));
[tref_num,ia] = sort(tref_num);
tref = cellstr(datestr(tref_num,'yyyymmdd'));
tref = cellfun(@str2double,tref);
%id=tref_num<datenum(2012,5,17);
%ia = ia(id);

data = data(ia,:);
cl = cell2mat(data(:,2));
hi = cell2mat(data(:,4));
lo = cell2mat(data(:,5));
op = cell2mat(data(:,3));

data=load('inputDataOHLCDaily_20120517', 'syms', 'tday', 'op', 'hi', 'lo', 'cl');
idx=find(strcmp(symbol, data.syms));
op1=data.op(:, idx);
hi1=data.hi(:, idx);
lo1=data.lo(:, idx);
cl1=data.cl(:, idx);
tday = data.tday(:,idx);

[~,ia,ib] = intersect(tref,tday);
a=[cl(ia),cl1(ib)];