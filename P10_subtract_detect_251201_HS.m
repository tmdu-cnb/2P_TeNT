% ===== フォルダ選択 =====
rootDir = uigetdir(pwd, '解析フォルダを選択してください');
if rootDir == 0
    error('フォルダが選択されませんでした。処理を中断します。');
end

disp(['選択フォルダ: ', rootDir]);

% ===== suite2p/plane0/fall.mat を読み込む =====

% まずは plane0 固定で探す
planeDir = fullfile(rootDir, 'suite2p', 'plane0');

% plane0 が無い場合は、suite2p 配下の plane* から1つ探す（保険）
if ~isfolder(planeDir)
    spDir = fullfile(rootDir, 'suite2p');
    planeCandidates = dir(fullfile(spDir, 'plane*'));
    planeCandidates = planeCandidates([planeCandidates.isdir]);

    if isempty(planeCandidates)
        error('suite2p フォルダ内に plane0 も plane* も見つかりませんでした: %s', spDir);
    else
        planeDir = fullfile(spDir, planeCandidates(1).name);
        warning('plane0 が無いため、%s を使用します。', planeDir);
    end
end

fallPath = fullfile(planeDir, 'fall.mat');

if ~isfile(fallPath)
    error('plane フォルダ内に fall.mat が見つかりません: %s', fallPath);
end

disp(['fall.mat を読み込み中... : ', fallPath]);
S = load(fallPath);

% ---- fall.mat の中身をすべて base workspace に展開 ----
fields = fieldnames(S);
for i = 1:numel(fields)
    assignin('base', fields{i}, S.(fields{i}));
end
disp('fall.mat の内容を workspace に展開しました。');


% ===== C=0 / C=1 の tiff を自動検出 =====
tifFiles = dir(fullfile(rootDir, '*.tif'));

if isempty(tifFiles)
    warning('フォルダ内に *.tif が存在しません。');
end

tifFile  = '';
tifFile2 = '';

for k = 1:numel(tifFiles)
    fname = tifFiles(k).name;

    % GCaMP or RFP のチャンネル名が "C=0", "C=1" だと仮定
    if contains(fname, 'C=0') && isempty(tifFile)
        tifFile = fullfile(rootDir, fname);
    end
    if contains(fname, 'C=1') && isempty(tifFile2)
        tifFile2 = fullfile(rootDir, fname);
    end
end

% ===== 結果確認 =====
if isempty(tifFile)
    warning('C=0 を含む tifFile が見つかりませんでした。');
else
    disp(['tifFile (C=0): ', tifFile]);
end

if isempty(tifFile2)
    warning('C=1 を含む tifFile2 が見つかりませんでした。');
else
    disp(['tifFile2 (C=1): ', tifFile2]);
end

disp('--- 全処理完了 ---');

mkdir result % folder name result

% Subtract開始

% % 条件の値を入力するポップアップを表示
% prompt = {'Enter frequency of imaging (freq, e.g., 5):', ...
%           'Enter MINIMUM threshold for npix (ROI size) (e.g., 20):', ...
%           'Enter MAXIMUM threshold for npix (ROI size, 0 to skip) (e.g., 1000):', ...
%           'Enter MINIMUM of y_range / x_range ratio (e.g., 0.5):', ...
%           'Enter MAXIMUM of y_range / x_range ratio (e.g., 2.0):'};
% dlg_title = 'Set Conditions';
% num_lines = 1;
% default_ans = {'1', '20', '1000', '0.5', '2.0'}; % デフォルト値を調整
% answer = inputdlg(prompt, dlg_title, num_lines, default_ans);

% % 入力値を数値に変換
% freq = str2double(answer{1});
% npix_threshold_min = str2double(answer{2});
% npix_threshold_max = str2double(answer{3});
% ratio_min = str2double(answer{4});
% ratio_max = str2double(answer{5});

freq = 1;
npix_threshold_min = 20;
npix_threshold_max = 1000;
ratio_min = 0.5;
ratio_max = 2;

% 入力された値を確認
disp(['Frequency (freq): ', num2str(freq)]);
disp(['npix MINIMUM threshold: ', num2str(npix_threshold_min)]);
disp(['npix MAXIMUM threshold: ', num2str(npix_threshold_max), ' (0=skip)']);
disp(['y_range / x_range ratio range: ', num2str(ratio_min), ' ~ ', num2str(ratio_max)]);


rows_to_delete = [];

% 元のROI番号を追加
for i = 1:numel(stat)
    stat{i}.original_roi_number = i-1; % 元のROI番号をフィールドとして追加
end

% statセル配列の各要素について処理
for i = 1:numel(stat)
    % npixが最小閾値未満の場合、対応するFの行を削除するためにrows_to_deleteに追加
    if stat{i}.npix < npix_threshold_min
        rows_to_delete = [rows_to_delete, i];
    end
    % npixが最大閾値を超える場合も削除（0の場合はスキップ）
    if npix_threshold_max > 0 && stat{i}.npix > npix_threshold_max
        rows_to_delete = [rows_to_delete, i];
    end
end

% 新しい削除条件を追加（min～maxの「間」を削除＝min以下とmax以上を解析対象として残す）
x_range_all = zeros(numel(stat),1);
y_range_all = zeros(numel(stat),1);
ratio_all   = zeros(numel(stat),1);

for i = 1:numel(stat)
    % 各ROIのxpixとypixの範囲を取得
    x_range = max(stat{i}.xpix) - min(stat{i}.xpix);
    y_range = max(stat{i}.ypix) - min(stat{i}.ypix);

    % 0割回避：x_rangeが0なら比率は無限大とみなす（縦長＝解析対象になりやすい）
    if x_range == 0
        ratio = Inf;
    else
        ratio = y_range / x_range;
    end

    x_range_all(i) = x_range;
    y_range_all(i) = y_range;
    ratio_all(i)   = ratio;
end

% 「min以下 or max以上」を“残す”→ その間（min < ratio < max）を削除
mask_between = (ratio_all > ratio_min) & (ratio_all < ratio_max); % これが「間」
% すでに rows_to_delete には npix < 閾値 のROIが入っている前提
rows_to_delete = unique([rows_to_delete, find(mask_between).']);  % 行ベクトルで結合

% 重複削除は unique 済み
% 対応するFの行を削除
F(rows_to_delete, :) = [];

% stat の対応要素を削除（末尾から消すのは元のまま）
for i = numel(rows_to_delete):-1:1
    stat(rows_to_delete(i)) = [];
end


nR_pre = size(F, 1);  % FのROIの数
nF_pre = size(F, 2);  % Fのフレームの数

size2 = size(F(:,1));
size3 = minus(size2(1),1);

for k = 1:size2(1)
  Y{k} = stat{k}.med(1);
end
Y2=Y.';
Y3=double(cell2mat(Y2));
%Y3=vertcat('y',Y2);
%Y4=cell2table(Y3);

for i = 1:size2(1)
  X{i} = stat{i}.med(2); 
end
X2=X.';
X3=double(cell2mat(X2));
%X3=vertcat('x',X2);
%X4=cell2table(X3);

a=(1:size3+1)';
b=zeros(size2(1),1);
%size_x=size(F(1,:));
%size_x2=size_x(2);
%c=1:1:size_x2;
%F1=vertcat(c,F);
F2=horzcat(a,b,F);

% cell_1=iscell(:,1);
%cell_pre3=num2cell(cell_1);
cell_pre1=ones(size2(1),1);
%cell_pre2=num2cell(cell_pre1);
%cell_2=vertcat('is_cell',cell_pre2);
%cell_2_zero=vertcat('is_cell',cell_pre3);
%cell_3_zero=table(cell_2_zero);
%cell_3=table(cell_2);
%F3pre=table(F2);
%F3=horzcat(F3pre,cell_3);
%F3_zero=horzcat(F3pre,cell_3_zero);
F3=horzcat(F2,cell_pre1);
%F3cell=num2cell(F3);
%F3cell_zero=num2cell(F3_zero);
F4=horzcat(X3,Y3);

size10 = size(F4, 1);
%F5=table(F4);
F6=horzcat(F3,F4);

% 列を削除するためのインデックスを作成
rows_to_remove2 = false(size(F,1), 1);  % 削除するROIを識別する論理配列

% for j = 1:nR_pre % 各ROIについて
% first_100_mean = mean(F6(j, 3:52), 'omitnan');
%     last_100_mean = mean(F6(j, end-52:end-3), 'omitnan');
% 
%     if first_100_mean >= 5 * last_100_mean
%         rows_to_remove2(j) = true;  % このROIを削除
%     end
% end

% 条件に合致する列を削除
F_filtered = F6(~rows_to_remove2, :);
F=F_filtered(:, 3:nF_pre+2);

nR_pre = size(F, 1);  % 新しいFのROIの数
nF_pre = size(F, 2);  % 新しいFのフレームの数

for i=1:nR_pre
min_F(i,1)=min(F(i,1:nF_pre));
end

min_F2=repmat(min_F,1,nF_pre);
F0_pre=F-min_F2;  %すべてのフレームで最も小さいシグナル値を引く。

%緩やかな変動をまっすぐにする
Frame_Size=50*freq;

for i = 1:nR_pre % 各ROIについて
    signal_pre = F0_pre(i, 1:nF_pre);

    for k = 1:fix(nF_pre / Frame_Size)  % 各開始フレームについて
        T_pre = mean(signal_pre((k-1) * Frame_Size + 1 : k * Frame_Size)) + std(signal_pre((k-1) * Frame_Size + 1 : k * Frame_Size));  % 現在のフレームからFsizeフレーム分の信号の閾値
        T2_pre = repmat(T_pre, 1, Frame_Size);
        signal_1 = signal_pre((k-1) * Frame_Size + 1 : k * Frame_Size);
        signal_2 = signal_1;
        signal_1(signal_1 > T2_pre) = NaN;
        T3 = mean(signal_1, 'omitnan');
        T4 = repmat(T3, 1, Frame_Size);
        signal_3 = signal_2 - T4;
        F_pre(i, (k-1) * Frame_Size + 1 : k * Frame_Size) = signal_3;
    end

    % 最後の残りのフレームを処理
    if mod(nF_pre, Frame_Size) > 0
        start_idx = fix(nF_pre / Frame_Size) * Frame_Size + 1;
        end_idx = nF_pre;
        T_pre = mean(signal_pre(start_idx:end_idx)) + std(signal_pre(start_idx:end_idx));
        T2_pre = repmat(T_pre, 1, end_idx - start_idx + 1);
        signal_1 = signal_pre(start_idx:end_idx);
        signal_2 = signal_1;
        signal_1(signal_1 > T2_pre) = NaN;
        T3 = mean(signal_1, 'omitnan');
        T4 = repmat(T3, 1, end_idx - start_idx + 1);
        signal_3 = signal_2 - T4;
        F_pre(i, start_idx:end_idx) = signal_3;
    end
end

F_pre2 = F_pre + min_F2(:, 1:size(F_pre, 2));


F6=horzcat(F_filtered(:,1:2), F_pre2,F_filtered(:,nF_pre+3:nF_pre+5));

% ★ここから追加: F6行に対応する original ROI 番号ベクタを用意
% いまの stat は F/npix/ratio フィルタ後のものなので、各行→stat{i} が1対1対応
orig_roi_vec = zeros(size(F6,1),1);
for ii = 1:size(F6,1)
    % Suite2Pの初期index-1（あなたが先に stat{i}.original_roi_number = i-1 と付与済み）
    orig_roi_vec(ii) = stat{ii}.original_roi_number;
end
% 削除ログ用の蓄積配列（step, original_roi_number）
removed_log_step   = strings(0,1);
removed_log_roi    = [];

% for k=1:10
% 
%     subplot(10,1,k);
%     axis([0 nF_pre -0.5 10])
%     plot(F_pre(k,1:nF_pre))
%     grid off
%     axis off
% 
% end
% filename_head = 'wave';
% filename = strcat( filename_head, num2str(i) ); 
% saveas(gcf,filename,'jpg')
% close(gcf);


%Detect開始



% ポップアップで閾値の設定値を入力させる


prompt = { ...
    'Threshold for the fluorescent intensity value (Check fluorescent intensity in Suite2P) (e.g.40):', ...
    'Threshold for the number of sd(T) (e.g.,1.5):', ...
    'Threshold for the number of sd(F) (e.g.,5):', ...
    'Threshold for removing low ROIs (例: 0.2):', ...
    'Number of the last frame for traces (例: 1000):', ...
    'Extra margin for the max-vs-mean test (e.g.,15):', ...
    'ROI番号をカンマ区切りで入力 (例: 5,12,23, 空欄=ランダム10個):' ...  % ← ROI選択を追加
};
dlg_title = 'Set Conditions';
num_lines = 1;
default_values = {'10', '1', '1', '0', '180', '2', ''}; % ← ROI番号の空欄を追加
answer = inputdlg(prompt, dlg_title, num_lines, default_values);

% 入力値を数値として取得
threshold_mean = str2double(answer{1});
std_multiplier_T = str2double(answer{2});
std_multiplier_F = str2double(answer{3});
roi_remove_ratio = str2double(answer{4});
Num_frames = str2double(answer{5});
max_margin = str2double(answer{6});  % ← 追加: 余白（+15 など）
roi_input = answer{7};  % ← ROI番号の入力（文字列）


% 必要に応じて結果を表示
disp('入力された閾値:');
disp(['平均値の閾値: ', num2str(threshold_mean)]);
disp(['標準偏差1倍の閾値: ', num2str(std_multiplier_T)]);
disp(['標準偏差5倍の閾値: ', num2str(std_multiplier_F)]);
disp(['ROI削除割合: ', num2str(roi_remove_ratio)]);
if isempty(roi_input)
    disp('波形表示: ランダム10個のROIを選択');
else
    disp(['波形表示ROI番号: ', roi_input]);
end
disp('処理完了');



% 元の配列のサイズを取得
nR = size(F6, 1);  % FneuのROIの数
nF = size(F6, 2);  % Fneuのフレームの数

% 各行の平均値と最大値を計算（列3～nF-3を評価対象に）
mean_values = mean(F6(:, 3:nF-3), 2);
max_values  = max(F6(:, 3:nF-3), [], 2);
min_values  = min(F6(:, 3:nF-3), [], 2);

% ①平均でのしきい値判定（従来）
rows_mean_fail = mean_values <= threshold_mean;

% ②最大値が（行の平均*max_marging）を超えない行は削除
rows_max_fail  = max_values <= (min_values*max_margin);

% ③いずれかに該当する行を削除
rows_to_remove = rows_mean_fail | rows_max_fail;

% ★ここから追加: この段階で落ちる original ROI を記録
removed_idx_1   = find(rows_to_remove);
removed_roi_1   = orig_roi_vec(removed_idx_1);
removed_log_step = [removed_log_step; repmat("mean_or_max_filter", numel(removed_roi_1), 1)];
removed_log_roi  = [removed_log_roi;  removed_roi_1];
% 残す側に対応ベクタを更新
orig_roi_vec = orig_roi_vec(~rows_to_remove);

F6(rows_to_remove, :) = [];

% **追加: F6の状態を確認**
F6_check = F6;


% 初期設定とサイズ計算
n_f = size(F6, 2) - 3;
n_R = size(F6, 1);
file_2 = F6(:, 3:end);
n_f2 = size(file_2, 2);
file_3 = file_2(:, 1:n_f2 - 3);
% TとT2の計算
T = mean(file_3, 2) + std_multiplier_T * std(file_3, 0, 2);
T2 = repmat(T, 1, n_f2 - 3);
% NaNで置き換え
file_3(file_3 > T2) = NaN;
% M_file_3, S_file_3, FBackを計算
M_file_3 = mean(file_3, 2, 'omitnan');
S_file_3 = std(file_3, 0, 2, 'omitnan');
FBack = M_file_3 + S_file_3;
% 閾値を設定して行をフィルタリング
F_Thre = M_file_3 + std_multiplier_F * S_file_3;
rows_to_keep = any(file_2(:, 1:n_f2-3) > F_Thre, 2);  % rows_to_keep を計算

% ★ここから追加: keep=0 で落ちる original ROI を記録
removed_idx_2    = find(~rows_to_keep);
removed_roi_2    = orig_roi_vec(removed_idx_2);
removed_log_step = [removed_log_step; repmat("background_screen", numel(removed_roi_2), 1)];
removed_log_roi  = [removed_log_roi;  removed_roi_2];
% 残す側に対応ベクタを更新
orig_roi_vec = orig_roi_vec(rows_to_keep);
% ★ここまで追加

% F6とFBackをフィルタリング
F6 = F6(rows_to_keep, :);  % F6をフィルタリング
FBack = FBack(rows_to_keep, :);  % FBackをフィルタリング
% F_signalの計算（bsxfunを使用してバックグラウンド値を引く）
F_signal = bsxfun(@minus, F6(:, 4:end-3), FBack);  % FBackをbsxfunで引く
F_signal = [F_signal, F6(:, 1)];  % ROI番号を最終列に保持

% **追加: F_signalの状態を確認**
F_signal_check = F_signal;

% F_signal内の負の値を0に置き換え
F_signal(F_signal < 0) = 0;
% F_signal2の計算（F_signalの末尾列を追加）
F_signal2_data = horzcat(F_signal(:, 1:end-1), F6(:, end-2:end));  % 信号データに末尾列を追加
F_signal2 = [F_signal2_data, F_signal(:, end)];  % ROI番号を最終列に保持

% **追加: F_signal2の状態を確認**
F_signal2_check = F_signal2;

% % サイズ情報を取得
% n_R3 = size(F_signal2, 1);
% n_f_minus_2 = n_f2 - 2;
% % 各行について処理
% for i = 1:n_R3
%    signal = F_signal2(i, :);  % 現在の行のデータを取得
%    % 条件1: 0より大きい値が1列のみで前後の列が0の場合
%    for j = 2:n_f_minus_2
%        if signal(j) > 0 && signal(j-1) <= 0 && signal(j+1) <= 0
%            signal(j) = 0;
%        end
%    end
%    % 条件2: 連続する2列の値がどちらも0より大きい値で前後の列が0の場合
%    for j = 2:n_f_minus_2-1
%        if signal(j) > 0 && signal(j+1) > 0 && signal(j-1) <= 0 && signal(j+2) <= 0
%            signal(j) = 0;
%            signal(j+1) = 0;
%        end
%    end
%    % 処理した結果をF_signal2に反映
%    F_signal2(i, :) = signal;
% end
% F_signal2_check_3 = F_signal2;
% 
% % widthの平均値が低いROIを最低から10%まで削除
% nR1_3 = size(F_signal2, 1);
% w_size = round(nR1_3 * roi_remove_ratio);
% A = zeros(nR1_3, 1); % Aを初期化
% for i = 1:nR1_3  % 各ROIについて
%    signal3 = F_signal2(i, 1:end-3);
%    [~, ~, width] = findpeaks(signal3);
%    A(i, 1) = mean(width);
% end
% [~, sortedIndices] = sort(A, 'ascend');
% rowsToDelete = sortedIndices(1:w_size);  % 該当列の最も低い1番目からw_size番目の値のインデックスを取得
% % 行を削除
% F_signal2(rowsToDelete, :) = [];
% 
% % **追加: F_signal2の状態を確認（削除後）**
% F_signal2_check2 = F_signal2;

% 任意の数だけすべて列が0の行を作成する数を指定
num_zero_columns = 0; % ここで任意の数を指定

% ★ROI番号の処理：入力があればそのROIを選択、なければランダム
selected_roi_numbers = [];
if ~isempty(roi_input) && ~isempty(strtrim(roi_input))
    % カンマ区切りの文字列を数値配列に変換
    selected_roi_numbers = str2num(roi_input); %#ok<ST2NM>
end

% Figure作成
min_value = min(Num_frames, nF_pre); % 120とnF_preの小さい方を選択
ranges = {[1:min_value], [1:min_value], [1:min_value]}; % 任意のframe数の範囲を指定
%ranges = {[1:min_value], [121:min_value+120], [181:min_value+180]}; % 任意のframe数の範囲を指定
titles = {'Columns 1', 'Columns 2', 'Columns 3'};
file_names = {'F_signal2_1.eps', 'F_signal2_2.eps', 'F_signal2_3.eps'};
total_lines = 10; % 合計の行数を計算
colors = lines(total_lines); % 必要な行数の異なる色を取得
for k = 1:3
   range = ranges{k};
   figure;
   hold on;
  
   % ROI番号が指定されている場合
   if ~isempty(selected_roi_numbers)
       % F_signal2内でこれらのROI番号を探す（末尾列がROI番号）
       if builtin('iscell', F_signal2)
           all_roi_numbers = cell2mat(F_signal2(:, end));
       else
           all_roi_numbers = F_signal2(:, end);
       end
       randomRows = [];
       for roi = selected_roi_numbers
           idx = find(all_roi_numbers == roi, 1);
           if ~isempty(idx)
               randomRows = [randomRows, idx]; %#ok<AGROW>
           end
       end
       if isempty(randomRows)
           warning('指定されたROI番号が見つかりません。ランダムに選択します。');
           if size(F_signal2, 1) < 10
               randomRows = 1:size(F_signal2, 1);
           else
               randomRows = randperm(size(F_signal2, 1), 10-num_zero_columns);
           end
       end
   % ROI番号が指定されていない場合はランダム
   elseif size(F_signal2, 1) < 10
       randomRows = 1:size(F_signal2, 1);
   else
       % ランダムに10行選択
       randomRows = randperm(size(F_signal2, 1), 10-num_zero_columns);
   end
   % 最大値を取得
   max_value = max(F_signal2(randomRows, range), [], 2);
   offset = 10; % 各グラフをオフセットする量
   % ランダムに選ばれた行と0の列の行を混合するためのインデックス
   allRows = [randomRows, zeros(1, num_zero_columns)]; % 0の行を追加
   allRows = allRows(randperm(length(allRows))); % シャッフル
   for i = 1:length(allRows)
       if allRows(i) == 0
           % 0の行の場合
           plot(zeros(1, length(range)) + (i-1) * offset, 'Color', colors(i, :));
          
           % 縦軸のスケール追加（スケールバーの縦線のみ）
           scale_value = 1; % すべて0の行のスケールは1とする
           plot([0 0], [0 scale_value] + (i-1) * offset, 'Color', colors(i, :), 'LineWidth', 1); % スケールバーの縦線
       else
           % ランダムに選ばれた行の場合
           % 各行の信号を最大値で正規化し、値をオフセットにスケール
           plot(F_signal2(allRows(i), range) + (i-1) * offset, 'Color', colors(i, :));
          
           % 縦軸のスケール追加（スケールバーの縦線のみ）
           scale_value = max_value(find(randomRows == allRows(i))) * 0.2; % 最大値の20%
           plot([0 0], [0 scale_value] + (i-1) * offset, 'Color', colors(i, :), 'LineWidth', 1); % スケールバーの縦線
           
           % ★ROI番号を図の中に表示
           % F_signal2の最後の列からROI番号を取得
           if builtin('iscell', F_signal2)
               roi_num = F_signal2{allRows(i), end};
           else
               roi_num = F_signal2(allRows(i), end);
           end
           % 波形の左側にROI番号を表示（白背景で見やすく）
         % 左側に大きく余白をとって番号を配置
x_pos = -round(length(range) * 0.15);  % 左に 5% の余白（必要なら 0.1 にするともっと左へ）

text(x_pos, (i-1) * offset + offset/2, sprintf('%d', roi_num), ...
    'Color', colors(i, :), 'FontSize', 12, 'FontWeight', 'bold', ...
    'BackgroundColor', 'none', 'EdgeColor', 'none');
       end
   end
  
   hold off;
   % X軸を表示、Y軸は非表示に設定
   set(gca, 'XColor', 'k', 'YColor', 'none');
   xlabel('Frame Number', 'FontSize', 10, 'FontWeight', 'bold');
   % 余白を追加（左右に5%ずつ）
   margin = length(range) * 0.05;
   xlim([1-margin length(range)+margin]);
   title(sprintf('Traces from F\\_signal2 (Frames: %d-%d)', range(1), range(end)), 'FontSize', 12);
  
   % EPSファイルとして保存
   saveas(gcf, file_names{k}, 'epsc');
end
n_R4 = size(F_signal2, 1);
% for n = 1:n_R4
%     F_2(n, 1:n_f2-3) = F_signal2(n, 1:n_f2-3) / max(F_signal2(n, 1:n_f2-3));
% end







% F_signal2をcell配列に変換
F_signal3 = num2cell(F_signal2);
% 必要な列を抽出
is_cell = vertcat('is_cell', F_signal3(:, n_f2-3));  % is_cell列
is_cellY = vertcat('y', F_signal3(:, n_f2-1));      % y座標列
is_cellX = vertcat('x', F_signal3(:, n_f2-2));        % x座標列
% 必要なデータを結合
d = horzcat(is_cell, is_cellY, is_cellX);  % is_cell, x, yを結合
d2 = table(d);  % テーブル化
% F_signal3からROI番号列を除去し、残りのデータを整理
F_signal3_cleaned = F_signal3(:, 2:n_f2-3);  % 信号データ部分のみ抽出
% ROI番号を最後の列に移動（1回のみ）

ROI_numbers = F_signal3(:, end);  % ROI番号列を抽出
F_signal4 = horzcat(F_signal3_cleaned, ROI_numbers);  % 信号データの最後にROI番号を追加

% ヘッダーを作成（列数を F_signal4 に合わせる）
header = [strcat("Frame_", string(1:(size(F_signal4, 2) - 1))), "ROI_number"];  % ヘッダーを自動調整
F_signal4_with_header = vertcat(header, F_signal4);  % ヘッダーを追加

% % % ヘッダーを作成
% header = [strcat("Frame_", string(1:(n_f2-3))), "ROI_number"];  % ヘッダー（フレーム列 + ROI番号）
% F_signal4_with_header = vertcat(header, F_signal4);  % ヘッダーを追加
% プレゼンテーション用のテーブルを作成
F_table = array2table(F_signal4_with_header);  % 信号データとROI番号を含むデータをテーブル化
F_signal5 = [F_table, d2];  % is_cell, x, y列を結合
% CSVファイルに書き出し
writetable(F_signal5, 'presentation.csv');

% ★ここから追加: F6→F_signal5 の間で削除された original ROI 番号ログを保存
% ログをテーブル化
T_removed = table(removed_log_step, removed_log_roi, ...
                  'VariableNames', {'step', 'original_ROI_number'});
% 1) 段階別の明細
% writetable(T_removed, fullfile('result','removed_original_ROIs_F6_to_F_signal5_detail.csv'));
% 2) 重複なしのユニークリスト（番号のみ）
unique_removed = unique(removed_log_roi);
T_unique = table(unique_removed, 'VariableNames', {'original_ROI_number'});
writetable(T_unique, fullfile('result','removed_original_ROIs_F6_to_F_signal5_unique.csv'));
% ★ここまで追加

% ROI_numbersとoriginal_ROI_numberを連結したデータを作成して格納

% ROI_numbersに対応するoriginal_roi_numberを取得
original_ROI_number_data = cell(length(ROI_numbers), 2);  % 2列のセル配列を初期化

for i = 1:length(ROI_numbers)
    % 現在のROI番号を取得
    roi_num = ROI_numbers{i};

    % 1列目: ROI番号
    original_ROI_number_data{i, 1} = roi_num;

    % 2列目: statから対応するoriginal_roi_numberを取得
    if isnumeric(roi_num) && ~isempty(stat{roi_num}) && isfield(stat{roi_num}, 'original_roi_number')
        original_ROI_number_data{i, 2} = stat{roi_num}.original_roi_number;
    else
        original_ROI_number_data{i, 2} = NaN;  % 該当するデータがない場合はNaNを格納
    end
end

% 結果をoriginal_ROI_numberに代入
original_ROI_number = original_ROI_number_data;

% ヘッダーを追加
header = {'ROI_number', 'Suite2P_ROI_number'};  % ヘッダーを定義

% ヘッダーをoriginal_ROI_numberの上に追加
original_ROI_number = vertcat(header, original_ROI_number);







