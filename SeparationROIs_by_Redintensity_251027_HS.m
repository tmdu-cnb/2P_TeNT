

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
