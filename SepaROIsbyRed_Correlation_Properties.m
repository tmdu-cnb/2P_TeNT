

%% ===== kept率の集計（全体 / 高輝度 / 低輝度）=====
inFile  = fullfile('result','roi_intensity_summary_all.csv');
outFile = fullfile('result','kept_percentages_summary.csv');

% ファイル存在チェック
assert(exist(inFile,'file')==2, '%s が見つかりません。先に集計ファイルを作成してください。', inFile);

% 読み込み
T = readtable(inFile);

% 必須列チェック
reqCols = {'status','mean_int'};
assert(all(ismember(reqCols, T.Properties.VariableNames)), ...
    '入力CSVに必要な列がありません（必要: %s）。', strjoin(reqCols, ', '));

% しきい値入力（例: 50）
ans_thr = inputdlg({'Intensity threshold on mean_int (e.g., 50):'}, ...
                   'Kept percentage with intensity gating', 1, {'50'});
thr = str2double(ans_thr{1});
if isnan(thr)
    error('しきい値が数値ではありません。入力を確認してください。');
end

% 型を整える
% status は大小文字無視で "kept" 判定できるように正規化
if iscategorical(T.status), T.status = string(T.status); end
if iscellstr(T.status),     T.status = string(T.status); end
if isstring(T.status)
    status_str = lower(strtrim(T.status));
elseif ischar(T.status)
    status_str = lower(strtrim(string(T.status)));
else
    error('status 列の型を解釈できません。');
end

% mean_int は数値前提（文字なら数値化）
if ~isnumeric(T.mean_int)
    T.mean_int = str2double(string(T.mean_int));
end

% マスク作成
is_kept = (status_str == "kept");
valid_int = ~isnan(T.mean_int);

mask_all  = true(height(T),1);       % 全体
mask_high = valid_int & (T.mean_int >= thr);
mask_low  = valid_int & (T.mean_int <= thr);

% 分母（件数）と kept 件数
N_all   = sum(mask_all);
K_all   = sum(is_kept & mask_all);

N_high  = sum(mask_high);
K_high  = sum(is_kept & mask_high);

N_low   = sum(mask_low);
K_low   = sum(is_kept & mask_low);

% パーセント（分母0のときは NaN）
pct_all  = (K_all  / N_all )*100;
pct_high = (N_high>0) * (K_high / max(N_high,1))*100; if N_high==0, pct_high = NaN; end
pct_low  = (N_low >0) * (K_low  / max(N_low ,1))*100; if N_low ==0, pct_low  = NaN; end

% 結果テーブル
Result = table( ...
    ["overall"; "mean_int>=thr"; "mean_int<=thr"], ...
    [N_all; N_high; N_low], ...
    [K_all; K_high; K_low], ...
    [pct_all; pct_high; pct_low], ...
    'VariableNames', {'category','n_total','n_kept','pct_kept'});

% しきい値を別テーブルで記録（任意）
Meta = table(thr, 'VariableNames', {'threshold'});

% 保存
if ~exist('result','dir'), mkdir result; end
writetable(Result, outFile);
writetable(Meta,   fullfile('result','kept_percentages_threshold.csv'));

% 画面にも概要表示
fprintf('=== kept percentages (threshold=%.3g) ===\n', thr);
disp(Result);
fprintf('保存: %s\n', outFile);

% Correlation始まり


%% === ROIのmean_intを読み込み、F_signal2の各行に対応付け ===
Tmean   = readtable(fullfile('result','roi_intensity_summary.csv')); % cols: roi, area_px, mean_int
roi_vec = Tmean.roi;                      % statの1始まりインデックス
mean_vec= Tmean.mean_int;

% --- 末尾列からROI番号を取り出して double ベクトルに正規化 ---
if istable(F_signal2)
    roi_in_use = F_signal2{:, end};      % table列の"中身"
else
    roi_in_use = F_signal2(:, end);      % 配列をそのまま
end

% roi_in_use を「double の数値ベクトル」に統一（必ず builtin('iscell',...) を使う）
if builtin('iscell', roi_in_use)
    if ~isempty(roi_in_use) && all(cellfun(@isnumeric, roi_in_use))
        roi_in_use = cell2mat(roi_in_use);
    else
        roi_in_use = str2double(string(roi_in_use));
    end
elseif isstring(roi_in_use) || ischar(roi_in_use) || iscategorical(roi_in_use)
    roi_in_use = str2double(string(roi_in_use));
else
    roi_in_use = double(roi_in_use);
end

% rowごとのmean_int（無いROIはNaN）
[tf,pos]      = ismember(roi_in_use, roi_vec);
mean_for_row  = nan(size(roi_in_use));
mean_for_row(tf) = mean_vec(pos(tf));

% 高輝度/低輝度マスク（NaN は除外）
mask_high = ~isnan(mean_for_row) & (mean_for_row >= thr);
mask_low  = ~isnan(mean_for_row) & (mean_for_row <= thr);

%% === ランダムサンプリング（最大300行） ===
n = size(F_signal2, 1);
F_signal2_backup = F_signal2; %#ok<NASGU>
idx_keep = 1:n;
if n > 300
    idx_keep = randperm(n, 300);
end
F_signal2    = F_signal2(idx_keep, :);
mean_for_row = mean_for_row(idx_keep);
mask_high    = mask_high(idx_keep);
mask_low     = mask_low(idx_keep);
n = size(F_signal2, 1); % 更新

% ---- ユーティリティ：相関行列と距離ビン平均を作る関数 ----
function [cormat, Mbin] = compute_corr_and_bins(Fsig, bin_width, bin_max)
    % Fsig: [n x (T+3)] 末尾3列が [x, y, roi] であることを想定（ここでは roi 列は未使用）
    nloc = size(Fsig,1);
    if nloc < 2
        cormat = nan(nloc);
        Mbin   = nan(1, floor(bin_max/bin_width)+1);
        return;
    end
    Tcols = size(Fsig,2) - 3;
    X  = Fsig(:, end-2);
    Y  = Fsig(:, end-1);

    % 事前に基線引きと分散（=二乗和）を計算
    S   = Fsig(:, 1:Tcols);
    S   = S - mean(S,2);            % 行ごとに平均を引く
    A   = sum(S.^2, 2);             % 各行の二乗和（0は定数列）

    % 出力器
    cormat = nan(nloc, nloc, 'like', S);

    % 距離ビン用
    nbin = floor(bin_max/bin_width) + 1;
    Dis_Cor = cell(1, nbin);

    for i = 1:nloc
        s1 = S(i,:);
        a1 = A(i);
        xi = X(i); yi = Y(i);

        for j = 1:nloc
            s2 = S(j,:);
            a2 = A(j);

            if a1==0 || a2==0
                cij = NaN; % 定数列は相関未定義
            else
                cij = sum(s1 .* s2) / sqrt(a1 * a2);
            end
            cormat(i,j) = cij;

            dx  = xi - X(j);
            dy  = yi - Y(j);
            dis = sqrt(dx.^2 + dy.^2);

            bin  = min(floor(dis / bin_width), nbin-1); % 0..nbin-1
            bidx = bin + 1;

            if ~isnan(cij) && (cij ~= 0)
                Dis_Cor{bidx}(end+1,1) = cij; %#ok<AGROW>
            end
        end
    end

    Mbin = nan(1, nbin);
    for k = 1:nbin
        if ~isempty(Dis_Cor{k})
            ck = Dis_Cor{k};
            Mbin(k) = mean(ck); % 既に0は除外済み
        end
    end
end

% ---- 全体 / 高輝度 / 低輝度 でのサンプル数チェック ----
n_overall = n;
n_high    = sum(mask_high);
n_low     = sum(mask_low);

% ---- 出力ディレクトリ ----
if ~exist('result','dir'), mkdir result; end

% 出力パス
p_overall_csv = fullfile('result','Correlation_index_overall.csv');
p_high_csv    = fullfile('result','Correlation_index_high.csv');
p_low_csv     = fullfile('result','Correlation_index_low.csv');

p_overall_svg = fullfile('result','correlation_fig_overall.svg');
p_high_svg    = fullfile('result','correlation_fig_high.svg');
p_low_svg     = fullfile('result','correlation_fig_low.svg');

% ---- 全体 ----
if n_overall >= 2
    [cormat_all, M_all] = compute_corr_and_bins(F_signal2, 5, 800);
    writematrix(M_all, p_overall_csv);
    figure; pcolor(cormat_all); axis ij; axis image; axis off; grid off; shading flat;
    colorbar('location','eastoutside'); caxis([0 1]);
    saveas(gcf, p_overall_svg, 'svg');
else
    % ファイルを作らない／既存があれば削除
    if exist(p_overall_csv,'file')==2, delete(p_overall_csv); end
    if exist(p_overall_svg,'file')==2, delete(p_overall_svg); end
    warning('overall のサンプル数が %d (<2) のため、相関計算・出力をスキップしました（ファイル未作成）。', n_overall);
end

% ---- 高輝度同士 ----
if n_high >= 2
    F_high = F_signal2(mask_high, :);
    [C_high, M_high] = compute_corr_and_bins(F_high, 5, 800);
    writematrix(M_high, p_high_csv);
    figure; pcolor(C_high); axis ij; axis image; axis off; grid off; shading flat;
    colorbar('location','eastoutside'); caxis([0 1]);
    saveas(gcf, p_high_svg, 'svg');
else
    if exist(p_high_csv,'file')==2, delete(p_high_csv); end
    if exist(p_high_svg,'file')==2, delete(p_high_svg); end
    warning('高輝度（>=%.3g）のサンプル数が %d (<2) のため、相関計算・出力をスキップしました（ファイル未作成）。', ...
        thr, n_high);
end

% ---- 低輝度同士 ----
if n_low >= 2
    F_low = F_signal2(mask_low, :);
    [C_low, M_low] = compute_corr_and_bins(F_low, 5, 800);
    writematrix(M_low, p_low_csv);
    figure; pcolor(C_low); axis ij; axis image; axis off; grid off; shading flat;
    colorbar('location','eastoutside'); caxis([0 1]);
    saveas(gcf, p_low_svg, 'svg');
else
    if exist(p_low_csv,'file')==2, delete(p_low_csv); end
    if exist(p_low_svg,'file')==2, delete(p_low_svg); end
    warning('低輝度（<=%.3g）のサンプル数が %d (<2) のため、相関計算・出力をスキップしました（ファイル未作成）。', ...
        thr, n_low);
end


% Calcium properties始まり

%% ==== 基本設定 ====
X          = F_signal2_backup;              % [nROI x nFramesAll]
nROI       = size(X,1);
nFramesAll = size(X,2);
trimTail   = 3;                        % 末尾3列は無視
nFrames    = nFramesAll - trimTail;    % 解析対象フレーム数
Fs         = freq;                     % サンプリング周波数 [Hz]（既存変数を利用）
Total_T    = nFrames / Fs;

maxIndex   = min(5, nROI);

%% ==== しきい値の事前計算 ====
Xmain   = X(:,1:nFrames);                         
mu      = mean(Xmain, 2);
sigma   = std(Xmain, 0, 2);
%Th      = mu + 1.*sigma; 
Th      = mu;
xmax    = max(Xmain, [], 2);
xmin    = min(Xmain, [], 2);
drange  = max(xmax - xmin, eps);                  

MinProm = min(Th, 0.5*drange);


MinH    = Th;

PeakAmp = cell(nROI,1);
location = cell(nROI,1);
wid = cell(nROI,1);
RelativeAmp = cell(nROI,1);

%% ==== 行ごとのピーク検出 ====
for r = 1:nROI
    x = Xmain(r,:);
    mp = max(MinProm(r), 0);                  
    useHeight = xmax(r) > MinH(r);

    if useHeight
        [PeakAmp{r}, location{r}, wid{r}, RelativeAmp{r}] = ...
            findpeaks(x, 'MinPeakProminence', mp, 'MinPeakHeight', MinH(r), 'MinPeakDistance', 4);
    else
        [PeakAmp{r}, location{r}, wid{r}, RelativeAmp{r}] = ...
            findpeaks(x, 'MinPeakProminence', mp, 'MinPeakDistance', 4);
    end
end

all_amp   = cell2mat(PeakAmp(:)');
all_width = cell2mat(wid(:)');

%% ==== 発火頻度 ====
Number_transients = cellfun(@numel, location);              
Fre_all           = Number_transients ./ Total_T;           
NOF               = sum(Fre_all > 0.0033);                  

%% ==== 積分（正値区間の面積） ====
area_All = cell(nROI,1);
all_area = [];

for r = 1:nROI
    sig = X(r,:);
    pos = sig > 0;
    if ~any(pos)
        area_All{r} = [];
        continue
    end
    d      = diff([false, pos, false]);     
    starts = find(d == 1);
    stops  = find(d == -1) - 1;
    A = zeros(1, numel(starts));
    for j = 1:numel(starts)
        A(j) = sum(sig(starts(j):stops(j))); 
    end
    area_All{r} = A;
    all_area    = [all_area, A];
end

%% ==== 60列ごとの最大振幅 ====
step60       = 60;
nInt60       = floor(nFrames / step60);
Amplitude_max_60cols = cell(nROI,1);

for r = 1:nROI
    x  = Xmain(r,:);
    mp = max(MinProm(r), 0);
    mh = MinH(r);
    Am = nan(1, nInt60);
    for k = 1:nInt60
        s = (k-1)*step60 + 1; e = k*step60;
        xi = x(s:e);
        if max(xi) > mh
            pk = findpeaks(xi, 'MinPeakProminence', mp, 'MinPeakHeight', mh, 'MinPeakDistance', 4);
        else
            pk = findpeaks(xi, 'MinPeakProminence', mp, 'MinPeakDistance', 4);
        end
        if ~isempty(pk), Am(k) = max(pk); end
    end
    Amplitude_max_60cols{r} = Am;
end

% ランダム 10 行を抽出して保存
rows10            = randperm(nROI, min(10,nROI));
Amplitude_selected = cell2mat(Amplitude_max_60cols(rows10));
writematrix(Amplitude_selected, 'Amplitude_selected.csv');
movefile('Amplitude_selected.csv', 'result');



% roi_intensity_summary.csv を読み込んで mean_int を取得
Tmean   = readtable(fullfile('result','roi_intensity_summary.csv'));  % 必須: 列名 roi / mean_int
roi_vec = double(Tmean.roi);          
mean_vec= double(Tmean.mean_int);

% この解析で使っている行 → ROI番号（X=F_signal2_2 は末尾列が ROI番号の想定）
roi_in_use = X(:, end);
if builtin('iscell', roi_in_use)
    if ~isempty(roi_in_use) && all(cellfun(@isnumeric, roi_in_use))
        roi_in_use = cell2mat(roi_in_use);
    else
        roi_in_use = str2double(string(roi_in_use));
    end
elseif isstring(roi_in_use) || ischar(roi_in_use) || iscategorical(roi_in_use)
    roi_in_use = str2double(string(roi_in_use));
else
    roi_in_use = double(roi_in_use);
end

% 各行に mean_int を対応付け
[tf,pos] = ismember(roi_in_use, roi_vec);
mean_for_row = nan(size(roi_in_use));
mean_for_row(tf) = mean_vec(pos(tf));

% マスク（高輝度/低輝度）。NaN はどちらにも入れない
mask_high = ~isnan(mean_for_row) & (mean_for_row >= thr);
mask_low  = ~isnan(mean_for_row) & (mean_for_row <= thr);

% インデックスとグループサイズ
idx_high = find(mask_high);
idx_low  = find(mask_low);
N_high   = numel(idx_high);
N_low    = numel(idx_low);
N_overall= nROI;

%% === グループ別データ抽出（元の計算結果を流用） ===
% amp（ピーク高さ）
all_amp_high   = cell2mat(PeakAmp(idx_high)');
all_amp_low    = cell2mat(PeakAmp(idx_low)');
% width
all_width_high = cell2mat(wid(idx_high)');
all_width_low  = cell2mat(wid(idx_low)');
% area
all_area_high  = cell2mat(area_All(idx_high)');
all_area_low   = cell2mat(area_All(idx_low)');
% frequency
Freq_high      = Fre_all(idx_high);
Freq_low       = Fre_all(idx_low);

%% === 出力準備 ===
outdir = 'result';
if ~exist(outdir,'dir'), mkdir(outdir); end

% 便利ムーバ（存在チェックして result へ移動）
safeMove = @(f) (exist(f,'file')==2 && movefile(f, fullfile(outdir, f)));
% 便利デリータ（カレント & result の両方を掃除）
safeDelete = @(f) ( (exist(f,'file')==2 && delete(f)) | ...
                    (exist(fullfile(outdir,f),'file')==2 && delete(fullfile(outdir,f))) );

%% === 全体指標（overall）: N_overall < 2 なら作成しない / 残骸削除 ===
if N_overall >= 2
    writematrix(all_width, 'Width.csv');
    writematrix(Fre_all,   'Freq.csv');
    writematrix(all_amp(:),'all_amp.csv');
    writematrix(all_area(:),'all_area.csv');
    safeMove('Width.csv'); safeMove('Freq.csv');
    safeMove('all_amp.csv'); safeMove('all_area.csv');
else
    % ===== ここを置き換え =====
    files = {'Width.csv','Freq.csv','all_amp.csv','all_area.csv'};
    for ii = 1:numel(files)
        f  = files{ii};
        fp = fullfile(outdir, f);
        if exist(f ,'file')==2, delete(f ); end
        if exist(fp,'file')==2, delete(fp); end
    end
    % ==========================
end


%% === 高輝度グループ: N_high < 2 なら作成しない / 残骸削除 ===
if N_high >= 2
    writematrix(all_amp_high,   'all_amp_high.csv');
    writematrix(all_width_high, 'Width_high.csv');
    writematrix(all_area_high,  'all_area_high.csv');
    writematrix(Freq_high(:),   'Freq_high.csv');
    safeMove('all_amp_high.csv');  safeMove('Width_high.csv');
    safeMove('all_area_high.csv'); safeMove('Freq_high.csv');
else

    % ===== ここを置き換え =====
    files = {'all_amp_high.csv','Width_high.csv','all_area_high.csv','Freq_high.csv'};
    for ii = 1:numel(files)
        f  = files{ii};
        fp = fullfile(outdir, f);
        if exist(f ,'file')==2, delete(f ); end
        if exist(fp,'file')==2, delete(fp); end
    end
    % ==========================
end


%% === 低輝度グループ: N_low < 2 なら作成しない / 残骸削除 ===
if N_low >= 2
    writematrix(all_amp_low,   'all_amp_low.csv');
    writematrix(all_width_low, 'Width_low.csv');
    writematrix(all_area_low,  'all_area_low.csv');
    writematrix(Freq_low(:),   'Freq_low.csv');
    safeMove('all_amp_low.csv');  safeMove('Width_low.csv');
    safeMove('all_area_low.csv'); safeMove('Freq_low.csv');
else
    % ===== ここを置き換え =====
    files = {'all_amp_low.csv','Width_low.csv','all_area_low.csv','Freq_low.csv'};
    for ii = 1:numel(files)
        f  = files{ii};
        fp = fullfile(outdir, f);
        if exist(f ,'file')==2, delete(f ); end
        if exist(fp,'file')==2, delete(fp); end
    end
    % ==========================
end


%% === ログ（gate_summary）は常に作成して残す ===
T_gate = table(thr, N_high, N_low, N_overall, ...
               'VariableNames', {'thr_mean_int','N_high','N_low','N_overall'});
writetable(T_gate, 'gate_summary.csv');
safeMove('gate_summary.csv');

%% === そのほか（任意の図ファイルなどは存在すれば移動） ===
safeMove('F_signal2_3.eps');
safeMove('F_signal2_2.eps');
safeMove('F_signal2_1.eps');
safeMove('presentation.csv');
safeMove('correlation_fig.svg');  % 相関関連を別ブロックで未作成にした場合、ここでは何もしない

% 度数分布開始

% histogram_analysis_grouped.m
% -------------------------------------------------------------------------
% 全体に加えて、roi_intensity_summary.csv の mean_int しきい値で
% 高輝度/低輝度グループの頻度分布も出力
% -------------------------------------------------------------------------
close all; clc;

% 解析結果フォルダ
resultDir = fullfile(pwd, 'result');
if ~isfolder(resultDir)
    error('現在のフォルダ (%s) に result フォルダがありません。先に前段の解析を実行してください。', pwd);
end

%% -------- Default parameters -------------------------------------------
defVals = struct( ...
    'amp',  [1     20     0.5 ], ...   % [binInterval max min]
    'area', [2     40     1   ], ...
    'freq', [0.01 0.2    0.01], ...
    'wid',  [0.5   10     0.25]);
params = [defVals.amp; defVals.area; defVals.freq; defVals.wid];

% %% -------- Single input dialog ------------------------------------------
% prompt = {
%     'Amplitude  –  bin interval', 'Amplitude  –  MAX', 'Amplitude  –  MIN', ...
%     'Area       –  bin interval', 'Area       –  MAX', 'Area       –  MIN', ...
%     'Frequency  –  bin interval', 'Frequency  –  MAX', 'Frequency  –  MIN', ...
%     'Width      –  bin interval', 'Width      –  MAX', 'Width      –  MIN', ...
%     'mean_int threshold (for HIGH/LOW split)'};
% defInput = [ ...
%     arrayfun(@num2str, [defVals.amp defVals.area defVals.freq defVals.wid], 'UniformOutput', false), ...
%     {'50'}];         % mean_int の初期しきい値
% answer = inputdlg(prompt, 'Histogram settings (all features + threshold)', [1 48], defInput);
% 
% if isempty(answer)  % cancelled → keep defaults
%     params = [defVals.amp; defVals.area; defVals.freq; defVals.wid];
%     thr_int = 50;
% else
%     numAns  = cellfun(@str2double, answer);
%     params  = reshape(numAns(1:12), 3, []).';   % each row: [interval max min]
%     thr_int = numAns(13);
% end

%% -------- Overall histograms (存在する場合のみ) -----------------------

filesOverall = struct( ...
    'amp','all_amp.csv', ...
    'area','all_area.csv', ...
    'freq','Freq.csv', ...
    'wid','Width.csv');

hasOverall = all(cellfun(@(f) exist(fullfile(resultDir,f),'file')==2, struct2cell(filesOverall)));

if ~hasOverall
    warning('overall 用の CSV (all_amp / all_area / Freq / Width) が揃っていないため、overall ヒストグラムはスキップします。');
else
    ampData  = readmatrix(fullfile(resultDir,filesOverall.amp));
    areaData = readmatrix(fullfile(resultDir,filesOverall.area));
    freqData = readmatrix(fullfile(resultDir,filesOverall.freq));
    widData  = readmatrix(fullfile(resultDir,filesOverall.wid));
    
    makeHistogram(ampData,  'Amplitude',               fullfile(resultDir,'amp_frequency_distribution.csv'),  params(1,1), params(1,2), params(1,3));
    makeHistogram(areaData, 'Area of Ca2+ transient',  fullfile(resultDir,'area_frequency_distribution.csv'), params(2,1), params(2,2), params(2,3));
    makeHistogram(freqData, 'Frequency (Hz)',          fullfile(resultDir,'freq_frequency_distribution.csv'), params(3,1), params(3,2), params(3,3));
    makeHistogram(widData,  'Width',                   fullfile(resultDir,'width_frequency_distribution.csv'),params(4,1), params(4,2), params(4,3));
end

%% -------- Grouped histograms (HIGH / LOW by mean_int) -------------------
% 前工程で出力したグループ別CSVを使用（無ければスキップ）
filesHigh = struct('amp','all_amp_high.csv','area','all_area_high.csv','freq','Freq_high.csv','wid','Width_high.csv');
filesLow  = struct('amp','all_amp_low.csv', 'area','all_area_low.csv', 'freq','Freq_low.csv', 'wid','Width_low.csv');

hasHigh = all(cellfun(@(f) exist(fullfile(resultDir,f),'file')==2, struct2cell(filesHigh)));
hasLow  = all(cellfun(@(f) exist(fullfile(resultDir,f),'file')==2, struct2cell(filesLow)));

% ログ（参考）
fprintf('mean_int threshold = %.3f\n', thr);
if ~hasHigh, warning('HIGH グループのCSVが揃っていません。前段の解析で N_high < 2 だった可能性があります。'); end
if ~hasLow,  warning('LOW グループのCSVが揃っていません。前段の解析で N_low < 2 だった可能性があります。');  end

if hasHigh
    H_amp  = readmatrix(fullfile(resultDir,filesHigh.amp));
    H_area = readmatrix(fullfile(resultDir,filesHigh.area));
    H_freq = readmatrix(fullfile(resultDir,filesHigh.freq));
    H_wid  = readmatrix(fullfile(resultDir,filesHigh.wid));

    makeHistogram(H_amp,  sprintf('Amplitude  (HIGH: mean\\_int \\ge %.3g)',thr), ...
                  fullfile(resultDir,'amp_frequency_distribution_high.csv'),  params(1,1), params(1,2), params(1,3));
    makeHistogram(H_area, sprintf('Area of Ca2+ transient  (HIGH: mean\\_int \\ge %.3g)',thr), ...
                  fullfile(resultDir,'area_frequency_distribution_high.csv'), params(2,1), params(2,2), params(2,3));
    makeHistogram(H_freq, sprintf('Frequency (Hz)  (HIGH: mean\\_int \\ge %.3g)',thr), ...
                  fullfile(resultDir,'freq_frequency_distribution_high.csv'), params(3,1), params(3,2), params(3,3));
    makeHistogram(H_wid,  sprintf('Width  (HIGH: mean\\_int \\ge %.3g)',thr), ...
                  fullfile(resultDir,'width_frequency_distribution_high.csv'),params(4,1), params(4,2), params(4,3));
end

if hasLow
    L_amp  = readmatrix(fullfile(resultDir,filesLow.amp));
    L_area = readmatrix(fullfile(resultDir,filesLow.area));
    L_freq = readmatrix(fullfile(resultDir,filesLow.freq));
    L_wid  = readmatrix(fullfile(resultDir,filesLow.wid));

    makeHistogram(L_amp,  sprintf('Amplitude  (LOW: mean\\_int \\le %.3g)',thr), ...
                  fullfile(resultDir,'amp_frequency_distribution_low.csv'),  params(1,1), params(1,2), params(1,3));
    makeHistogram(L_area, sprintf('Area of Ca2+ transient  (LOW: mean\\_int \\le %.3g)',thr), ...
                  fullfile(resultDir,'area_frequency_distribution_low.csv'), params(2,1), params(2,2), params(2,3));
    makeHistogram(L_freq, sprintf('Frequency (Hz)  (LOW: mean\\_int \\le %.3g)',thr), ...
                  fullfile(resultDir,'freq_frequency_distribution_low.csv'), params(3,1), params(3,2), params(3,3));
    makeHistogram(L_wid,  sprintf('Width  (LOW: mean\\_int \\le %.3g)',thr), ...
                  fullfile(resultDir,'width_frequency_distribution_low.csv'),params(4,1), params(4,2), params(4,3));
end

%% ====== Local helper function ==========================================
function makeHistogram(data, featureName, csvOut, binInt, maxVal, minVal)
    if isempty(data) || all(~isfinite(data))
        warning('"%s" のデータが空または無効です: %s', featureName, csvOut);
        return;
    end
    data = data(:);
    data = data(isfinite(data));

    % ビン設定と頻度
    binEdges = minVal:binInt:maxVal;
    if numel(binEdges) < 2
        warning('"%s" のビンが不足しています（min/max/binIntervalを見直してください）。', featureName);
        return;
    end
    [counts, edges] = histcounts(data, binEdges);
    relCounts = counts / max(sum(counts),1) * 100;

    % テーブル保存
    tbl = table(edges(1:end-1)', edges(2:end)', counts', relCounts', ...
        'VariableNames', {'Start', 'End', 'Count', 'RelPercent'});
    writetable(tbl, csvOut);

    % 図（相対頻度を散布で）
    figure('Color','w');
    centers = edges(1:end-1) + diff(edges)/2;
    scatter(centers, relCounts, 's', 'MarkerEdgeColor','k', 'MarkerFaceColor','k');
    xlabel(featureName, 'FontSize',12,'FontWeight','bold');
    ylabel('Relative Frequency (%)', 'FontSize',12,'FontWeight','bold');
    xlim([min(binEdges) max(binEdges)]); ylim([0 100]);
    set(gca,'LineWidth',1.2); pbaspect([1 0.6 1]);
    title(sprintf('%s — freq. distribution', featureName), 'Interpreter','none');
end

% move_result_to_root.m
% rootDir で選んだ解析フォルダの中に、カレントの result フォルダを移動する

if ~exist('rootDir','var')
    error('rootDir が workspace にありません。最初のフォルダ選択スクリプトを先に実行してください。');
end

srcResult = fullfile(pwd, 'result');      % 現在のフォルダの result
destResult = fullfile(rootDir, 'result'); % 解析フォルダ内の result

if ~isfolder(srcResult)
    warning('現在のフォルダ (%s) に result フォルダがありません。移動をスキップします。', pwd);
    return;
end

% すでに rootDir 側に result がある場合は、タイムスタンプ付きで退避
if isfolder(destResult)
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    backupDir = fullfile(rootDir, ['result_' timestamp]);
    warning('rootDir 内に既に result があるため、%s にリネームして保存します。', backupDir);
    destResult = backupDir;
end

try
    movefile(srcResult, destResult);
    fprintf('result フォルダを\n  %s\nから\n  %s\nへ移動しました。\n', srcResult, destResult);
catch ME
    warning('result フォルダの移動に失敗しました: %s', ME.message);
end

clear