function ROIs_in_Image_Interactive()
% ROIs_in_Image_Interactive - ROI番号を指定してリアルタイムで画像と波形を表示するGUI
%
% 必要なもの（caller / base workspace に存在していること）:
%   - stat   : Suite2P の ROI 構造体 / セル配列
%   - img_G  : GCaMP画像 (2D or 3D, 512x512 double など)
%   - img    : RFP画像   (2D or 3D)
%   - presentation.csv : F_signal5 を書き出した CSV

% ========================
% 0. workspace から必須変数を取得
% ========================

% --- stat ---
if evalin('caller', 'exist(''stat'', ''var'')') ~= 1
    error('Workspace に stat がありません。Suite2P の fall.mat から stat をロードしてください。');
end
stat_data = evalin('caller', 'stat');

% iscell の衝突を避けつつアクセサ作成
if builtin('iscell', stat_data)
    stat_at = @(i) stat_data{i};
    nStat   = numel(stat_data);
else
    stat_at = @(i) stat_data(i);
    nStat   = numel(stat_data);
end

% --- GCaMP画像 (img_G) ---
if evalin('caller', 'exist(''img_G'', ''var'')') ~= 1
    error('Workspace に img_G がありません。img_G を作成してから実行してください。');
end
img_G = evalin('caller', 'img_G');

% --- RFP画像 (img) ---
if evalin('caller', 'exist(''img'', ''var'')') ~= 1
    error('Workspace に img がありません。tif を読み込んだ img を作成してから実行してください。');
end
img_RFP_raw = evalin('caller', 'img');

% --- CSV ファイル名 ---
csvFile = 'presentation.csv';

% ========================
% 1. 画像の整形
% ========================

% GCaMP (img_G) → グレースケール
if ndims(img_G) == 3
    img_gcamp = rgb2gray(img_G);
else
    img_gcamp = img_G;
end
img_gcamp = double(img_gcamp);
[H, W] = size(img_gcamp);

% RFP (img) → グレースケール
if ndims(img_RFP_raw) == 3
    img_rfp = rgb2gray(img_RFP_raw);
else
    img_rfp = img_RFP_raw;
end
img_rfp = double(img_rfp);
[H_rfp, W_rfp] = size(img_rfp);

% ========================
% 2. CSV から ROI 番号を取得
% ========================

all_roi_numbers = [];
csv_data = [];
csv_roi_col = [];

if exist(csvFile,'file') == 2
    try
        csv_data = readmatrix(csvFile);
        % 3行目以降, 末尾から4列目に ROI番号がある前提
        tmp = csv_data(3:end, end-3);
        tmp = tmp(~isnan(tmp));
        all_roi_numbers = unique(tmp(:)');
        csv_roi_col = tmp;
    catch ME
        warning('CSV %s の読み込みに失敗: %s', csvFile, ME.message);
    end
else
    warning('CSVファイル %s が見つかりません。波形表示には CSV が必要です。', csvFile);
end

if isempty(all_roi_numbers)
    all_roi_numbers = 1:nStat;
    fprintf('CSVが使えないため，全ROI（%d個）を表示対象とします。\n', nStat);
end

% ========================
% 3. ROI 形状を前計算（高速化ポイント）
% ========================

roiInfo = [];   % struct 配列
cnt = 0;

for i = 1:numel(all_roi_numbers)
    r = all_roi_numbers(i);
    if r < 1 || r > nStat
        continue;
    end

    S = stat_at(r);
    if ~isfield(S,'xpix') || ~isfield(S,'ypix')
        continue;
    end

    x = double(S.xpix(:));
    y = double(S.ypix(:));

    if numel(x) < 3 || numel(y) < 3
        continue;
    end

    % --- GCaMP 用 ---
    xg = min(max(x, 0.5), W+0.5);
    yg = min(max(y, 0.5), H+0.5);

    try
        kg = boundary(xg, yg, 0.8);
        pxg = xg(kg);
        pyg = yg(kg);
    catch
        pxg = xg;
        pyg = yg;
    end

    if numel(pxg) < 3
        continue;
    end

    mask_g = poly2mask(pxg, pyg, H, W);
    statsC_g = regionprops(mask_g, 'Centroid');
    if ~isempty(statsC_g)
        cg = statsC_g(1).Centroid;
    else
        cg = [mean(pxg), mean(pyg)];
    end

    % --- RFP 用 ---
    xr = min(max(x, 0.5), W_rfp+0.5);
    yr = min(max(y, 0.5), H_rfp+0.5);

    try
        kr = boundary(xr, yr, 0.8);
        pxr = xr(kr);
        pyr = yr(kr);
    catch
        pxr = xr;
        pyr = yr;
    end

    if numel(pxr) >= 3
        mask_r = poly2mask(pxr, pyr, H_rfp, W_rfp);
        statsC_r = regionprops(mask_r, 'Centroid');
        if ~isempty(statsC_r)
            cr = statsC_r(1).Centroid;
        else
            cr = [mean(pxr), mean(pyr)];
        end
    else
        pxr = [];
        pyr = [];
        cr  = [NaN, NaN];
    end

    % struct に格納
    cnt = cnt + 1;
    roiInfo(cnt).roi        = r;      %#ok<*AGROW>
    roiInfo(cnt).px_gcamp   = pxg;
    roiInfo(cnt).py_gcamp   = pyg;
    roiInfo(cnt).cx_gcamp   = cg(1);
    roiInfo(cnt).cy_gcamp   = cg(2);

    roiInfo(cnt).px_rfp     = pxr;
    roiInfo(cnt).py_rfp     = pyr;
    roiInfo(cnt).cx_rfp     = cr(1);
    roiInfo(cnt).cy_rfp     = cr(2);
end

if isempty(roiInfo)
    error('有効な ROI が見つかりませんでした。stat / ROI番号を確認してください。');
end

num_rois = numel(roiInfo);

% ========================
% 4. GUI 作成
% ========================

fig = figure('Name', 'ROI Interactive Viewer', ...
             'Position', [100, 100, 1800, 800], ...
             'Color', 'w', ...
             'CloseRequestFcn', @closeFigure);

% 画像表示用 axes
ax_image_gcamp = axes('Parent', fig, ...
                      'Position', [0.02, 0.1, 0.38, 0.85], ...
                      'Units', 'normalized');
ax_image_rfp   = axes('Parent', fig, ...
                      'Position', [0.42, 0.1, 0.38, 0.85], ...
                      'Units', 'normalized');
axis(ax_image_gcamp, 'off');
axis(ax_image_rfp,   'off');

% コントロールパネル
ctrlPanel = uipanel('Parent', fig, ...
                    'Units', 'normalized', ...
                    'Position', [0.82, 0.1, 0.17, 0.8], ...
                    'BackgroundColor', 'w', ...
                    'Title', 'Controls', ...
                    'FontSize', 10);

% ROI 入力ラベル
uicontrol('Parent', ctrlPanel, ...
          'Style', 'text', ...
          'String', 'ROI番号（カンマ区切り）:', ...
          'Units', 'normalized', ...
          'Position', [0.05, 0.85, 0.9, 0.08], ...
          'HorizontalAlignment', 'left', ...
          'FontSize', 11, ...
          'BackgroundColor', 'w');

% ROI 入力欄
edit_roi = uicontrol('Parent', ctrlPanel, ...
                     'Style', 'edit', ...
                     'String', '', ...
                     'Units', 'normalized', ...
                     'Position', [0.05, 0.75, 0.9, 0.08], ...
                     'FontSize', 11, ...
                     'BackgroundColor', 'w', ...
                     'HorizontalAlignment', 'left');

% 「表示」ボタン
uicontrol('Parent', ctrlPanel, ...
          'Style', 'pushbutton', ...
          'String', '表示', ...
          'Units', 'normalized', ...
          'Position', [0.05, 0.64, 0.42, 0.09], ...
          'FontSize', 12, ...
          'FontWeight', 'bold', ...
          'Callback', @(~,~) updateDisplay());

% 「クリア」ボタン
uicontrol('Parent', ctrlPanel, ...
          'Style', 'pushbutton', ...
          'String', 'クリア', ...
          'Units', 'normalized', ...
          'Position', [0.53, 0.64, 0.42, 0.09], ...
          'FontSize', 12, ...
          'Callback', @(~,~) clearDisplay());

% 「波形を表示」ボタン
uicontrol('Parent', ctrlPanel, ...
          'Style', 'pushbutton', ...
          'String', '波形を表示', ...
          'Units', 'normalized', ...
          'Position', [0.05, 0.52, 0.9, 0.09], ...
          'FontSize', 12, ...
          'Callback', @(~,~) showTraces());

% ステータス表示
if ~isempty(csv_data)
    status_msg = 'ROI番号を入力して「表示」をクリック';
else
    status_msg = 'CSVがないため、波形表示は利用できません。';
end

txt_status = uicontrol('Parent', ctrlPanel, ...
                       'Style', 'text', ...
                       'String', status_msg, ...
                       'Units', 'normalized', ...
                       'Position', [0.05, 0.30, 0.9, 0.20], ...
                       'FontSize', 10, ...
                       'ForegroundColor', [0.4, 0.4, 0.4], ...
                       'BackgroundColor', 'w', ...
                       'HorizontalAlignment', 'left');

% UserData に全部突っ込む
ud = struct();
ud.img_gcamp   = img_gcamp;
ud.img_rfp     = img_rfp;
ud.H           = H;
ud.W           = W;
ud.H_rfp       = H_rfp;
ud.W_rfp       = W_rfp;
ud.stat_at     = stat_at;
ud.nStat       = nStat;
ud.roiInfo     = roiInfo;
ud.csvFile     = csvFile;
ud.csv_data    = csv_data;
ud.csv_roi_col = csv_roi_col;
ud.highlight_rois = [];
ud.edit_roi    = edit_roi;
ud.txt_status  = txt_status;
ud.ax_image_gcamp = ax_image_gcamp;
ud.ax_image_rfp   = ax_image_rfp;

set(fig, 'UserData', ud);

% 初期描画
redrawImage();

% ========================
% 5. ネスト関数
% ========================

    function updateDisplay()
        ud = get(fig, 'UserData');
        roi_input = strtrim(get(ud.edit_roi, 'String'));

        if isempty(roi_input)
            ud.highlight_rois = [];
            set(ud.txt_status, 'String', '強調表示なし');
        else
            try
                hr = str2num(roi_input); %#ok<ST2NM>
                if isempty(hr)
                    set(ud.txt_status, 'String', '無効な入力です');
                    return;
                end
                ud.highlight_rois = hr;
                set(ud.txt_status, 'String', sprintf('ROI %s を強調表示', num2str(hr)));
            catch
                set(ud.txt_status, 'String', '入力エラー');
                return;
            end
        end

        set(fig, 'UserData', ud);
        redrawImage();
    end

    function clearDisplay()
        ud = get(fig, 'UserData');
        set(ud.edit_roi, 'String', '');
        ud.highlight_rois = [];
        set(ud.txt_status, 'String', '強調表示をクリアしました');
        set(fig, 'UserData', ud);
        redrawImage();
    end

    function redrawImage()
        ud = get(fig, 'UserData');
        axG = ud.ax_image_gcamp;
        axR = ud.ax_image_rfp;

        img_gcamp = ud.img_gcamp;
        img_rfp   = ud.img_rfp;
        roiInfo   = ud.roiInfo;
        H         = ud.H;
        W         = ud.W;
        H_rfp     = ud.H_rfp;
        W_rfp     = ud.W_rfp;
        highlight_rois = ud.highlight_rois;

        % ==== GCaMP ====
        cla(axG);
        imagesc(axG, img_gcamp);
        colormap(axG, gray(256));
        axis(axG, 'image');
        set(axG, 'YDir', 'reverse');
        title(axG, 'GCaMP', 'FontSize', 12, 'FontWeight', 'bold');
        hold(axG, 'on');

        % ==== RFP ====
        cla(axR);
        imagesc(axR, img_rfp);
        colormap(axR, gray(256));
        axis(axR, 'image');
        set(axR, 'YDir', 'reverse');
        title(axR, 'RFP', 'FontSize', 12, 'FontWeight', 'bold');
        hold(axR, 'on');

        % ==== ROI 描画 ====
        num_rois = numel(roiInfo);
        colors = lines(num_rois);

        for k = 1:num_rois
            info = roiInfo(k);
            r    = info.roi;

            is_highlight = ~isempty(highlight_rois) && any(highlight_rois == r);

            % --- GCaMP ---
            pxg = info.px_gcamp;
            pyg = info.py_gcamp;
            if numel(pxg) >= 3
                if is_highlight
                    patch(axG, pxg, pyg, 'r', ...
                          'EdgeColor', 'r', 'LineWidth', 2, 'FaceAlpha', 0);
                else
                    patch(axG, pxg, pyg, colors(k,:), ...
                          'EdgeColor', colors(k,:), 'LineWidth', 1.2, 'FaceAlpha', 0);
                end

                cxg = info.cx_gcamp; cyg = info.cy_gcamp;
                if ~isnan(cxg)
                    text(axG, cxg, cyg, sprintf('%d', r), ...
                        'Color', 'w', 'FontSize', 7, 'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'BackgroundColor', 'none', 'EdgeColor', 'none');
                end
            end

            % --- RFP ---
            pxr = info.px_rfp;
            pyr = info.py_rfp;
            if numel(pxr) >= 3
                if is_highlight
                    patch(axR, pxr, pyr, 'r', ...
                          'EdgeColor', 'r', 'LineWidth', 2, 'FaceAlpha', 0);
                else
                    patch(axR, pxr, pyr, colors(k,:), ...
                          'EdgeColor', colors(k,:), 'LineWidth', 1.2, 'FaceAlpha', 0);
                end

                cxr = info.cx_rfp; cyr = info.cy_rfp;
                if ~isnan(cxr)
                    text(axR, cxr, cyr, sprintf('%d', r), ...
                        'Color', 'w', 'FontSize', 7, 'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'BackgroundColor', 'none', 'EdgeColor', 'none');
                end
            end
        end

        xlim(axG, [0.5, W+0.5]); ylim(axG, [0.5, H+0.5]);
        xlim(axR, [0.5, W_rfp+0.5]); ylim(axR, [0.5, H_rfp+0.5]);

        drawnow;
    end

    function showTraces()
        ud = get(fig, 'UserData');
        txt_status = ud.txt_status;
        highlight_rois = ud.highlight_rois;
        csv_data = ud.csv_data;
        csv_roi_col = ud.csv_roi_col;

        if isempty(highlight_rois)
            set(txt_status, 'String', 'ROI番号を入力してください');
            return;
        end
        if isempty(csv_data)
            set(txt_status, 'String', 'CSVが読み込まれていないため、波形表示はできません。');
            return;
        end

        try
            fig_trace = figure('Name', 'ROI Traces', ...
                               'Position', [200, 200, 1000, 600], ...
                               'Color', 'w');

            hr = highlight_rois(:)';
            num_traces = numel(hr);
            colors_trace = lines(num_traces);
            traces_shown = 0;

            for idx = 1:num_traces
                roi = hr(idx);
                row_idx = find(csv_roi_col == roi, 1);
                if isempty(row_idx)
                    warning('ROI %d の波形データが見つかりません', roi);
                    continue;
                end

                traces_shown = traces_shown + 1;

                % 3行目が row_idx=1 に対応 → 実際の行は row_idx+2
                % 最後の4列はメタ情報 → それ以前を波形とみなす
                trace_row = row_idx + 2;
                trace_data = csv_data(trace_row, 1:end-4);

                subplot(num_traces, 1, traces_shown, 'Parent', fig_trace);
                plot(trace_data, 'Color', colors_trace(idx,:), 'LineWidth', 1.5);
                title(sprintf('ROI %d', roi), 'FontSize', 11, 'FontWeight', 'bold');
                ylabel('Signal');
                grid on;
                if traces_shown == num_traces
                    xlabel('Frame Number');
                end
            end

            if traces_shown > 0
                set(txt_status, 'String', sprintf('波形データを表示しました（%d個）', traces_shown));
            else
                set(txt_status, 'String', '選択したROIの波形データが見つかりませんでした');
            end

        catch ME
            warning('波形データの表示に失敗: %s', ME.message);
            set(txt_status, 'String', sprintf('波形表示中にエラー: %s', ME.message));
        end
    end

    function closeFigure(~, ~)
        delete(fig);
    end

end
