

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
dlg_title = 'Condition settings';
num_lines = 1;
default_values = {'10', '1', '2', '0', '100', '1.2', ''}; % ← ROI番号の空欄を追加
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

% サイズ情報を取得
n_R3 = size(F_signal2, 1);
n_f_minus_2 = n_f2 - 2;
% 各行について処理
for i = 1:n_R3
   signal = F_signal2(i, :);  % 現在の行のデータを取得
   % 条件1: 0より大きい値が1列のみで前後の列が0の場合
   for j = 2:n_f_minus_2
       if signal(j) > 0 && signal(j-1) <= 0 && signal(j+1) <= 0
           signal(j) = 0;
       end
   end
   % 条件2: 連続する2列の値がどちらも0より大きい値で前後の列が0の場合
   for j = 2:n_f_minus_2-1
       if signal(j) > 0 && signal(j+1) > 0 && signal(j-1) <= 0 && signal(j+2) <= 0
           signal(j) = 0;
           signal(j+1) = 0;
       end
   end
   % 処理した結果をF_signal2に反映
   F_signal2(i, :) = signal;
end
F_signal2_check_3 = F_signal2;

% widthの平均値が低いROIを最低から10%まで削除
nR1_3 = size(F_signal2, 1);
w_size = round(nR1_3 * roi_remove_ratio);
A = zeros(nR1_3, 1); % Aを初期化
for i = 1:nR1_3  % 各ROIについて
   signal3 = F_signal2(i, 1:end-3);
   [~, ~, width] = findpeaks(signal3);
   A(i, 1) = mean(width);
end
[~, sortedIndices] = sort(A, 'ascend');
rowsToDelete = sortedIndices(1:w_size);  % 該当列の最も低い1番目からw_size番目の値のインデックスを取得
% 行を削除
F_signal2(rowsToDelete, :) = [];

% **追加: F_signal2の状態を確認（削除後）**
F_signal2_check2 = F_signal2;

% 任意の数だけすべて列が0の行を作成する数を指定
num_zero_columns = 0; % ここで任意の数を指定

% ★ROI番号の処理：入力があればそのROIを選択、なければランダム
selected_roi_numbers = [];
if ~isempty(roi_input) && ~isempty(strtrim(roi_input))
    % カンマ区切りの文字列を数値配列に変換
    selected_roi_numbers = str2num(roi_input); %#ok<ST2NM>
end

% Figure作成
min_value = min(120, nF_pre); % 120とnF_preの小さい方を選択
ranges = {[1:Num_frames], [1:Num_frames], [1:Num_frames]}; % 任意のframe数の範囲を指定
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
           % 波形の左側にROI番号を表示（番号のみ、枠なし）
           text(-3, (i-1) * offset + offset/2, sprintf('%d', roi_num), ...
               'Color', colors(i, :), 'FontSize', 4, 'FontWeight', 'bold', ...
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








