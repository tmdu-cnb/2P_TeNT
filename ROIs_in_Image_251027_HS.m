% ===== 白黒ベース画像 + カラーROI + ROI番号ラベル =====
tifFile = 'AVG_C2-済TeNT1_Anterior_15min_2.tif';
% for GCaMP
tifFile2 = 'AVG_TeNT1_Anterior_15min_1.nd2 - C=0.tif';
% for RFP

img0 = imread(tifFile);

% グレースケール化
if ndims(img0) == 3
    img = rgb2gray(img0);
else
    img = img0;
end
[H, W] = size(img);

% CSVからROI番号を読み込み（末尾から4列目想定）
csvFile = 'presentation.csv';
data = readmatrix(csvFile);
roi_numbers = data(3:end, end-3);
roi_numbers = roi_numbers(~isnan(roi_numbers));
roi_numbers = unique(round(roi_numbers(:)'));   % 念のため丸め & 重複排除

% ★CSVが original(0始まり) の場合は true → 1始まりへ変換
useOriginalIndex = false;   % ←必要に応じて true に
if useOriginalIndex
    roi_numbers = roi_numbers + 1;
end

% ★ROI番号指定で強調表示＋波形データ表示機能を追加（図を描く前に聞く）
prompt_highlight = {'強調表示したいROI番号をカンマ区切りで入力 (例: 5,12,23, 空欄=なし):'};
answer_highlight = inputdlg(prompt_highlight, 'Highlight Specific ROIs (Optional)', 1, {''});

highlight_rois = [];
if ~isempty(answer_highlight) && ~isempty(strtrim(answer_highlight{1}))
    highlight_rois = str2num(answer_highlight{1}); %#ok<ST2NM>
end

% 表示（白黒固定）
figure('Color','w');
imshow(img, [], 'InitialMagnification', 'fit');
colormap(gray(256));
axis image; hold on;
set(gca,'YDir','reverse');   % y軸を画像の行方向（上→下で増加）に合わせる

% ROI色
num_rois = numel(roi_numbers);
colors = lines(num_rois);

% ROIごとに描画 & ラベリング
for i = 1:num_rois
    r = roi_numbers(i);
    if r < 1 || r > numel(stat)
        warning('ROI %d は stat の範囲外です（スキップ）。', r);
        continue;
    end

    x = double(stat{r}.xpix(:));
    y = double(stat{r}.ypix(:));

    % 範囲クリップ（はみ出し防止）
    x = min(max(x, 0.5), W+0.5);
    y = min(max(y, 0.5), H+0.5);

    if numel(x) < 3
        warning('ROI %d の頂点数が少なすぎます（スキップ）。', r);
        continue;
    end

    % 外周ポリゴン（なめらかめ）
    k  = boundary(x, y, 0.8);
    px = x(k); py = y(k);

    % ★ROIの形を強調：輪郭を太く、塗りは薄く
    fill(px, py, colors(i,:), ...
        'EdgeColor', colors(i,:), 'LineWidth', 3, 'FaceAlpha', 0.15);

    % 重心をマスクから推定して番号描画（小さめ＆透明背景で控えめに）
    mask = poly2mask(px, py, H, W);
    statsC = regionprops(mask, 'Centroid');
    if ~isempty(statsC)
        c = statsC(1).Centroid;   % [x y]
    else
        c = [mean(px) mean(py)];  % フォールバック
    end

    % ROI番号を小さく控えめに表示
    text(c(1), c(2), sprintf('%d', r), ...
        'Color','w','FontSize',7,'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'BackgroundColor', 'none', 'EdgeColor', 'none');
end

% 全体が見えるように軸固定
xlim([0.5, W+0.5]); ylim([0.5, H+0.5]);

% ★強調表示する ROI がある場合の処理
if ~isempty(highlight_rois)
        % 強調表示するROI
        for i = 1:numel(highlight_rois)
            r = highlight_rois(i);
            if r < 1 || r > numel(stat)
                warning('ROI %d は stat の範囲外です（スキップ）。', r);
                continue;
            end
            
            x = double(stat{r}.xpix(:));
            y = double(stat{r}.ypix(:));
            x = min(max(x, 0.5), W+0.5);
            y = min(max(y, 0.5), H+0.5);
            
            if numel(x) < 3
                continue;
            end
            
            % 外周ポリゴン
            k = boundary(x, y, 0.8);
            px = x(k); py = y(k);
            
            % ★太い線＋塗りつぶしで強調表示（テキストなし）
            fill(px, py, 'r', 'EdgeColor', 'r', 'LineWidth', 5, 'FaceAlpha', 0.3);
        end
        
        % ★波形データを表示
        if exist(csvFile, 'file') == 2
            try
                % presentation.csvからデータを読み込み
                csv_data = readmatrix(csvFile);
                csv_roi_col = csv_data(3:end, end-3);  % ROI番号列
                
                % 各強調ROIの波形を表示
                figure('Name', 'Highlighted ROI Traces', 'Position', [100, 100, 1200, 800]);
                num_plots = numel(highlight_rois);
                colors_trace = lines(num_plots);
                
                for i = 1:num_plots
                    roi_to_find = highlight_rois(i);
                    
                    % CSVからROIのデータを検索
                    roi_idx = find(csv_roi_col == roi_to_find);
                    
                    if isempty(roi_idx)
                        warning('ROI %d の波形データが見つかりません。', roi_to_find);
                        continue;
                    end
                    
                    % 波形データを取得（2列目から末尾-3まで）
                    trace_data = csv_data(2+roi_idx, 2:end-3);
                    
                    % サブプロット
                    subplot(num_plots, 1, i);
                    plot(trace_data, 'Color', colors_trace(i,:), 'LineWidth', 2);
                    title(sprintf('ROI #%d - Calcium Trace', roi_to_find), 'FontSize', 14, 'FontWeight', 'bold');
                    xlabel('Frame Number', 'FontSize', 12);
                    ylabel('Fluorescence Intensity', 'FontSize', 12);
                    grid on;
                    xlim([1, length(trace_data)]);
                end
                
                disp('強調表示したROIの波形データを表示しました。');
            catch ME
                warning('波形データの読み込みに失敗しました: %s', ME.message);
            end
        else
            warning('presentation.csv が見つからないため、波形データは表示できません。');
        end
end

% 保存
outFile = 'output_image_with_ROI_colored_G.tif';
exportgraphics(gca, outFile, 'Resolution', 300);
dest = 'result'; if ~exist(dest,'dir'), mkdir(dest); end
movefile(outFile, fullfile(dest, outFile));

% 画像読み込み
img = imread(tifFile2);
if ndims(img) == 3, img = rgb2gray(img); end
img = double(img);
[H, W] = size(img);

% stat チェック
assert(exist('stat','var')==1, 'Workspace に stat がありません。');

% iscell の衝突を回避してアクセサ作成
if builtin('iscell', stat)
    stat_at = @(i) stat{i};
    nStat   = numel(stat);
else
    stat_at = @(i) stat(i);
    nStat   = numel(stat);
end

% ROI 番号の取得（CSVが読めなければ全ROI）
roi_numbers = [];
if exist(csvFile,'file') == 2
    try
        data = readmatrix(csvFile);
        tmp = data(3:end, end-3);
        tmp = tmp(~isnan(tmp));
        roi_numbers = unique(tmp(:)');
    catch
        warning('ROI番号CSVの読み込みに失敗。全ROIを処理します。');
    end
end
if isempty(roi_numbers)
    roi_numbers = 1:nStat;
end

% 可視化準備
figure('Color','w');
imagesc(img); axis image off; colormap gray; set(gca,'YDir','reverse'); hold on;
colors = lines(numel(roi_numbers));

% 出力用
res = struct('roi', num2cell(roi_numbers), ...
             'area_px', [], 'mean_int', []);

% ROIごとにマスク化＆計測
for k = 1:numel(roi_numbers)
    r = roi_numbers(k);
    S = stat_at(r);

    x = double(S.xpix(:));
    y = double(S.ypix(:));

    % 画像範囲にクリップ（はみ出し対策）
        % 画像範囲にクリップ（はみ出し対策）
    x = min(max(x, 0.5), W+0.5);
    y = min(max(y, 0.5), H+0.5);

    % ==== 計測用マスク（従来どおり）====
    mask = poly2mask(x, y, H, W);
    vals = img(mask);

    if isempty(vals)
        warning('ROI %d のマスクが空です。頂点順や座標系を確認してください。', r);
        res(k).area_px = 0;
        res(k).mean_int = NaN;
    else
        res(k).area_px  = nnz(mask);
        res(k).mean_int = mean(vals);
    end

    % ==== 表示用：輪郭だけ（外周をなめらか化して線描画）====
    % 外周ポリゴン（なめらかめ）
    try
        kk = boundary(x, y, 0.8);   % 0<shrink<1 小さいほど外周が滑らか
        px = x(kk); py = y(kk);
    catch
        % boundary が使えない場合のフォールバック
        px = x;    py = y;
    end

    % 輪郭線のみ（塗り無し）
    plot([px; px(1)], [py; py(1)], '-', ...
         'Color', colors(k,:), 'LineWidth', 1.2);

    % ==== ラベル用の重心（表示用外周→マスク化、空なら計測用マスク）====
    mask_disp = poly2mask(px, py, H, W);
    if ~any(mask_disp,'all'), mask_disp = mask; end

    statsC = regionprops(mask_disp, 'Centroid');
    if ~isempty(statsC)
        c = statsC(1).Centroid;   % [x y]
    else
        c = [mean(px) mean(py)];  % 最後のフォールバック
    end

    % ==== ROI番号ラベル（二重描画で視認性UP）====
    text(c(1), c(2), sprintf('%d', r), ...
        'Color', 'k', 'FontSize', 11, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    text(c(1), c(2), sprintf('%d', r), ...
        'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
% === ここまで置換 ===

% 保存
outImg = 'output_image_with_ROI_colored_R.tif';
exportgraphics(gca, outImg, 'Resolution', 300);
dest = 'result'; if ~exist(dest,'dir'), mkdir(dest); end
movefile(outImg, fullfile(dest, outImg));

T = struct2table(res);
writetable(T, fullfile(dest, 'roi_intensity_summary.csv'));
fprintf('完了: %s に結果CSVを保存しました。\n', fullfile(dest, 'roi_intensity_summary.csv'));




%% 追加: 除外された original ROI 番号の輝度も測定する
removed_csv = fullfile('result','removed_original_ROIs_F6_to_F_signal5_unique.csv');
if exist(removed_csv,'file') ~= 2
    warning('除外ROIリスト %s が見つかりません。除外ROIの輝度測定はスキップします。', removed_csv);
else
    % 0始まりの original ROI 番号を読み込み
    T_removed = readtable(removed_csv);
    if ~ismember('original_ROI_number', T_removed.Properties.VariableNames)
        warning('除外ROI CSVに original_ROI_number 列がありません。スキップします。');
    else
       removed_orig = T_removed.original_ROI_number;

% --- 型に依存せず数値ベクトルへ変換（ビルトイン iscell を強制使用）---
if builtin('iscell', removed_orig)
    % セル配列。文字→数値、数値セル→そのまま数値化
    if all(cellfun(@(x) ischar(x) || isstring(x), removed_orig))
        removed_orig = str2double(string(removed_orig));
    else
        removed_orig = cellfun(@double, removed_orig);
    end
elseif isstring(removed_orig) || ischar(removed_orig) || iscategorical(removed_orig)
    removed_orig = str2double(string(removed_orig));
else
    removed_orig = double(removed_orig);
end

% --- NaNを除去して一意化（0-basedの original ROI 番号）---
removed_orig = removed_orig(~isnan(removed_orig));
removed_orig = unique(removed_orig(:)');


        % original_roi_number (0-based) -> stat index (1-based) の対応表
        origvec = nan(1, nStat);
        for ii = 1:nStat
            Sii = stat_at(ii);
            if isfield(Sii, 'original_roi_number')
                origvec(ii) = double(Sii.original_roi_number); % 0-based
            end
        end

        % 対応する stat インデックスを探す（見つからないものは捨てる）
        removed_idx = [];
        for v = removed_orig
            hit = find(origvec==v, 1, 'first');
            if ~isempty(hit)
                removed_idx(end+1) = hit; %#ok<AGROW>
            else
                warning('original ROI %d は現在の stat で見つかりませんでした。', v);
            end
        end
        removed_idx = unique(removed_idx);

        % 測定（mean/area のみ）
        res_removed = struct('roi', num2cell(removed_idx), ...
                             'original_roi_number', [], ...
                             'area_px', [], 'mean_int', []);
        for k = 1:numel(removed_idx)
            r = removed_idx(k);
            S = stat_at(r);

            x = double(S.xpix(:));
            y = double(S.ypix(:));
            x = min(max(x, 0.5), W+0.5);
            y = min(max(y, 0.5), H+0.5);

           % x,y の基本チェック
if ~isfield(S,'xpix') || ~isfield(S,'ypix') || numel(S.xpix) < 3
    warning('除外ROI %d (orig=%s) に有効な xpix/ypix がありません。', r, string(S.original_roi_number));
    res_removed(k).area_px  = 0;
    res_removed(k).mean_int = NaN;
    res_removed(k).original_roi_number = double(getfield(S, 'original_roi_number')); %#ok<GFLD>
    continue;
end

x = double(S.xpix(:));
y = double(S.ypix(:));
x = min(max(x, 0.5), W+0.5);
y = min(max(y, 0.5), H+0.5);

% 1) まずは多角形としてマスク化
mask = poly2mask(x, y, H, W);

% 2) 空なら「内部点群→散布→充填」でフォールバック
if ~any(mask, 'all')
    xi = min(max(round(x), 1), W);
    yi = min(max(round(y), 1), H);
    mask = false(H, W);
    mask(sub2ind([H, W], yi, xi)) = true;
    mask = imfill(mask, 'holes');  % 必要なら bwconvhull(mask,'union') も可
end

vals = img(mask);
res_removed(k).original_roi_number = double(S.original_roi_number); % 0-based
if isempty(vals)
    warning('除外ROI %d (orig=%d) のマスクが空です。', r, res_removed(k).original_roi_number);
    res_removed(k).area_px  = 0;
    res_removed(k).mean_int = NaN;
else
    res_removed(k).area_px  = nnz(mask);
    res_removed(k).mean_int = mean(vals);
end

        end

        % 保存（除外ROIのみ）
        T_removed_sum = struct2table(res_removed);
        writetable(T_removed_sum, fullfile(dest, 'roi_intensity_summary_removed.csv'));

        % 保持ROIテーブルに status を付与して結合（mean/area のみ）
        T_kept = T;
        T_kept.status = repmat("kept", height(T_kept), 1);

        % 除外は列を合わせる（roi/area_px/mean_int）+status
        T_removed_sum2 = T_removed_sum(:, {'roi','area_px','mean_int'});
        T_removed_sum2.status = repmat("removed", height(T_removed_sum2), 1);

        T_all = [T_kept; T_removed_sum2];
        writetable(T_all, fullfile(dest, 'roi_intensity_summary_all.csv'));

        fprintf('除外ROIの輝度測定完了：\n  - %s\n  - %s\n', ...
            fullfile(dest,'roi_intensity_summary_removed.csv'), ...
            fullfile(dest,'roi_intensity_summary_all.csv'));
    end
end