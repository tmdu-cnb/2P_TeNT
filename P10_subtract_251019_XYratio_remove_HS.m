mkdir result % folder name result

% 条件の値を入力するポップアップを表示
prompt = {'Enter frequency of imaging (freq, e.g., 5):', ...
          'Enter threshold for npix (ROI size) (e.g., 2.5 * 2.5 * 3.14):', ...
          'Enter MINIMUM of y_range / x_range ratio (e.g., 2.5 for dendrite):', ...
          'Enter MAXIMUM of y_range / x_range ratio (e.g., 200 for dendrite):'};
dlg_title = 'Set Conditions';
num_lines = 1;
default_ans = {'1', '20', '0.5', '2.0'}; % デフォルト値を調整
answer = inputdlg(prompt, dlg_title, num_lines, default_ans);

% 入力値を数値に変換
freq = str2double(answer{1});
npix_threshold = str2double(answer{2});
ratio_min = str2double(answer{3});
ratio_max = str2double(answer{4});

% 入力された値を確認
disp(['Frequency (freq): ', num2str(freq)]);
disp(['npix threshold: ', num2str(npix_threshold)]);
disp(['y_range / x_range ratio range: ', num2str(ratio_min), ' ~ ', num2str(ratio_max)]);


rows_to_delete = [];

% 元のROI番号を追加
for i = 1:numel(stat)
    stat{i}.original_roi_number = i-1; % 元のROI番号をフィールドとして追加
end

% statセル配列の各要素について処理
for i = 1:numel(stat)
    % npixが50未満の場合、対応するFの行を削除するためにrows_to_deleteに追加
    if stat{i}.npix < npix_threshold
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
Frame_Size=20*freq;

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




