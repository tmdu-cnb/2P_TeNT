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
Th      = mu + 1.*sigma;                          
xmax    = max(Xmain, [], 2);
xmin    = min(Xmain, [], 2);
drange  = max(xmax - xmin, eps);                  

MinProm = min(Th, 0.5*drange);
MinH    = 0.2*Th;

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

%% ==== mean_int しきい値で ROI を二分（高/低） ====
% しきい値をポップアップで入力（Cancel されたら 50）
ans_thr = inputdlg({'Intensity threshold on mean_int (e.g., 50):'}, ...
                   'Intensity gating for metrics', 1, {'50'});
if isempty(ans_thr)
    thr_int = 50;
else
    thr_int = str2double(ans_thr{1});
    if isnan(thr_int), thr_int = 50; end
end

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
mask_high = ~isnan(mean_for_row) & (mean_for_row >= thr_int);
mask_low  = ~isnan(mean_for_row) & (mean_for_row <= thr_int);

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
T_gate = table(thr_int, N_high, N_low, N_overall, ...
               'VariableNames', {'thr_mean_int','N_high','N_low','N_overall'});
writetable(T_gate, 'gate_summary.csv');
safeMove('gate_summary.csv');

%% === そのほか（任意の図ファイルなどは存在すれば移動） ===
safeMove('F_signal2_3.eps');
safeMove('F_signal2_2.eps');
safeMove('F_signal2_1.eps');
safeMove('presentation.csv');
safeMove('correlation_fig.svg');  % 相関関連を別ブロックで未作成にした場合、ここでは何もしない
