function ROIs_in_Image_Interactive()
% ROIs_in_Image_Interactive - ROI番号を指定してリアルタイムで画像と波形を表示するGUI
% 
% 使用方法:
%   ROIs_in_Image_Interactive()
%
% GUI操作:
%   - ROI番号をカンマ区切りで入力（例: 5,12,23）
%   - 「表示」ボタンをクリックして画像を更新
%   - 「クリア」ボタンで強調表示を解除
%   - 「波形を表示」ボタンで波形データを表示

% ===== 初期化 =====
% ファイルパスの設定
tifFile = 'AldoCP71a.tif';  % for GCaMP
tifFile2 = 'AVG_TeNT1_Anterior_15min_1.nd2 - C=0.tif';  % for RFP
csvFile = 'presentation.csv';

% stat チェック
if evalin('caller', 'exist(''stat'', ''var'')') ~= 1
    error('Workspace に stat がありません。stat を読み込んでください。');
end

% stat をワークスペースから取得（関数名の衝突を回避）
stat_data = evalin('caller', 'stat');

% iscell の衝突を回避してアクセサ作成
if builtin('iscell', stat_data)
    stat_at = @(i) stat_data{i};
    nStat   = numel(stat_data);
else
    stat_at = @(i) stat_data(i);
    nStat   = numel(stat_data);
end

% GCaMP画像の読み込み
if exist(tifFile, 'file') ~= 2
    warning('画像ファイル %s が見つかりません。ファイルを選択してください。', tifFile);
    [file, path] = uigetfile({'*.tif;*.tiff', 'TIFF files (*.tif, *.tiff)'; ...
                              '*.png', 'PNG files (*.png)'; ...
                              '*.jpg;*.jpeg', 'JPEG files (*.jpg, *.jpeg)'; ...
                              '*.*', 'All files'}, ...
                             'GCaMP画像ファイルを選択');
    if isequal(file, 0)
        error('画像ファイルが選択されませんでした。処理を終了します。');
    end
    tifFile = fullfile(path, file);
end

try
    img0 = imread(tifFile);
catch ME
    warning('画像ファイルの読み込みに失敗しました: %s', ME.message);
    [file, path] = uigetfile({'*.tif;*.tiff', 'TIFF files (*.tif, *.tiff)'; ...
                              '*.png', 'PNG files (*.png)'; ...
                              '*.jpg;*.jpeg', 'JPEG files (*.jpg, *.jpeg)'; ...
                              '*.*', 'All files'}, ...
                             'GCaMP画像ファイルを再選択');
    if isequal(file, 0)
        error('画像ファイルが選択されませんでした。処理を終了します。');
    end
    tifFile = fullfile(path, file);
    img0 = imread(tifFile);
end

% GCaMP画像のグレースケール化
if ndims(img0) == 3
    img_gcamp = rgb2gray(img0);
else
    img_gcamp = img0;
end
[H, W] = size(img_gcamp);

% RFP画像の読み込み
img_rfp = [];
if exist(tifFile2, 'file') == 2
    % ファイルが存在する場合は読み込みを試みる
    try
        img0_rfp = imread(tifFile2);
        if ndims(img0_rfp) == 3
            img_rfp = rgb2gray(img0_rfp);
        else
            img_rfp = img0_rfp;
        end
    catch ME
        warning('RFP画像ファイル %s の読み込みに失敗: %s', tifFile2, ME.message);
        warning('ファイルを選択してください。');
        % 読み込み失敗時もファイル選択ダイアログを表示
        [file2, path2] = uigetfile({'*.tif;*.tiff', 'TIFF files (*.tif, *.tiff)'; ...
                                    '*.png', 'PNG files (*.png)'; ...
                                    '*.jpg;*.jpeg', 'JPEG files (*.jpg, *.jpeg)'; ...
                                    '*.*', 'All files'}, ...
                                   'RFP画像ファイルを選択（キャンセルでスキップ）');
        if ~isequal(file2, 0)
            tifFile2 = fullfile(path2, file2);
            try
                img0_rfp = imread(tifFile2);
                if ndims(img0_rfp) == 3
                    img_rfp = rgb2gray(img0_rfp);
                else
                    img_rfp = img0_rfp;
                end
            catch ME2
                warning('RFP画像の読み込みに失敗: %s', ME2.message);
                img_rfp = [];
            end
        end
    end
else
    % ファイルが存在しない場合は警告を表示してファイル選択ダイアログを表示
    warning('RFP画像ファイル %s が見つかりません。ファイルを選択してください。', tifFile2);
    [file2, path2] = uigetfile({'*.tif;*.tiff', 'TIFF files (*.tif, *.tiff)'; ...
                                '*.png', 'PNG files (*.png)'; ...
                                '*.jpg;*.jpeg', 'JPEG files (*.jpg, *.jpeg)'; ...
                                '*.*', 'All files'}, ...
                               'RFP画像ファイルを選択（キャンセルでスキップ）');
    if ~isequal(file2, 0)
        tifFile2 = fullfile(path2, file2);
        try
            img0_rfp = imread(tifFile2);
            if ndims(img0_rfp) == 3
                img_rfp = rgb2gray(img0_rfp);
            else
                img_rfp = img0_rfp;
            end
        catch ME
            warning('RFP画像の読み込みに失敗: %s', ME.message);
            img_rfp = [];
        end
    end
end

% CSVからROI番号を読み込み
all_roi_numbers = [];
if exist(csvFile,'file') == 2
    try
        data = readmatrix(csvFile);
        tmp = data(3:end, end-3);
        tmp = tmp(~isnan(tmp));
        all_roi_numbers = unique(tmp(:)');
    catch ME
        warning('ROI番号CSVの読み込みに失敗: %s', ME.message);
    end
else
    % CSVファイルが見つからない場合、選択ダイアログを表示
    warning('CSVファイル %s が見つかりません。ファイルを選択してください。', csvFile);
    [file, path] = uigetfile({'*.csv', 'CSV files (*.csv)'; '*.*', 'All files'}, ...
                             'presentation.csv を選択（キャンセルで全ROIを表示）');
    if ~isequal(file, 0)
        csvFile = fullfile(path, file);
        try
            data = readmatrix(csvFile);
            tmp = data(3:end, end-3);
            tmp = tmp(~isnan(tmp));
            all_roi_numbers = unique(tmp(:)');
        catch ME
            warning('CSVファイルの読み込みに失敗: %s', ME.message);
        end
    end
end
if isempty(all_roi_numbers)
    all_roi_numbers = 1:nStat;
    fprintf('全ROI（%d個）を表示します。\n', nStat);
end

% ===== GUIの作成 =====
% メインFigureウィンドウ
fig = figure('Name', 'ROI Interactive Viewer', ...
             'Position', [100, 100, 1800, 800], ...
             'Color', 'w', ...
             'CloseRequestFcn', @closeFigure);

% 画像表示用のaxes（左側: GCaMP、右側: RFP）
% すべて正規化単位で配置してウィンドウサイズ変更に対応

% コントロールパネルの位置パラメータを先に定義
if ~isempty(img_rfp)
    % RFPがある場合：コントロールパネルは左から70%の位置から開始
    control_x_start = 0.70;
    control_width = 0.28;  % コントロールパネルの幅（右端まで使う、最大0.98まで）
    % 右端がウィンドウを超えないように調整
    if control_x_start + control_width > 0.98
        control_width = 0.98 - control_x_start;
    end
else
    % RFPがない場合：コントロールパネルは左から58%の位置から開始
    control_x_start = 0.58;
    control_width = 0.38;  % コントロールパネルの幅（右端まで使う、最大0.98まで）
    % 右端がウィンドウを超えないように調整
    if control_x_start + control_width > 0.98
        control_width = 0.98 - control_x_start;
    end
end

if ~isempty(img_rfp)
    % RFPがある場合：2つのaxesを並べて配置
    % 左から2%、幅32%、間隔2%で配置（コントロールパネルとの間に余白を確保）
    ax_image_gcamp = axes('Parent', fig, ...
                          'Position', [0.02, 0.1, 0.32, 0.85], ...
                          'Units', 'normalized', ...
                          'Clipping', 'on', ...
                          'XLimMode', 'manual', ...
                          'YLimMode', 'manual');
    ax_image_rfp = axes('Parent', fig, ...
                        'Position', [0.36, 0.1, 0.32, 0.85], ...
                        'Units', 'normalized', ...
                        'Clipping', 'on', ...
                        'XLimMode', 'manual', ...
                        'YLimMode', 'manual');
    axis(ax_image_gcamp, 'off');
    axis(ax_image_rfp, 'off');
    ax_image = ax_image_gcamp;  % 互換性のため
else
    % RFPがない場合：GCaMPのみ
    ax_image_gcamp = axes('Parent', fig, ...
                          'Position', [0.02, 0.1, 0.55, 0.85], ...
                          'Units', 'normalized', ...
                          'Clipping', 'on', ...
                          'XLimMode', 'manual', ...
                          'YLimMode', 'manual');
    ax_image_rfp = [];
    axis(ax_image_gcamp, 'off');
    ax_image = ax_image_gcamp;
end

% コントロールパネル（右側）- uipanelとして作成して確実に表示
[control_panel, edit_roi, txt_status] = createControlPanel(fig, control_x_start, control_width, csvFile);

% パネルを前面に表示
try
    uistack(control_panel, 'top');
catch ME
    warning('uistackで前面表示に失敗: %s', ME.message);
end
drawnow;

% RFP画像のサイズを取得（存在する場合）
if ~isempty(img_rfp)
    [H_rfp, W_rfp] = size(img_rfp);
else
    H_rfp = H;
    W_rfp = W;
end

% データをFigureのUserDataに保存
set(fig, 'UserData', struct('img_gcamp', img_gcamp, ...
                            'img_rfp', img_rfp, ...
                            'H', H, ...
                            'W', W, ...
                            'H_rfp', H_rfp, ...
                            'W_rfp', W_rfp, ...
                            'stat_at', stat_at, ...
                            'nStat', nStat, ...
                            'all_roi_numbers', all_roi_numbers, ...
                            'csvFile', csvFile, ...
                            'highlight_rois', [], ...
                            'edit_roi', edit_roi, ...
                            'txt_status', txt_status, ...
                            'ax_image_gcamp', ax_image_gcamp, ...
                            'ax_image_rfp', ax_image_rfp));

% 初期画像を表示
redrawImage();

% ===== コールバック関数 =====

    function updateDisplay()
        % データを取得
        ud = get(fig, 'UserData');
        edit_roi = ud.edit_roi;
        txt_status = ud.txt_status;
        
        % 入力されたROI番号を取得
        roi_input = get(edit_roi, 'String');
        roi_input = strtrim(roi_input);
        
        % ROI番号をパース
        if isempty(roi_input)
            highlight_rois = [];
            set(txt_status, 'String', '強調表示なし');
        else
            try
                highlight_rois = str2num(roi_input); %#ok<ST2NM>
                if isempty(highlight_rois)
                    set(txt_status, 'String', '無効な入力です');
                    return;
                end
                set(txt_status, 'String', sprintf('ROI %s を強調表示', num2str(highlight_rois)));
            catch
                set(txt_status, 'String', '入力エラー');
                return;
            end
        end
        
        % UserDataを更新
        ud.highlight_rois = highlight_rois;
        set(fig, 'UserData', ud);
        
        % 画像を再描画
        redrawImage();
    end

    function clearDisplay()
        % データを取得
        ud = get(fig, 'UserData');
        edit_roi = ud.edit_roi;
        txt_status = ud.txt_status;
        
        % 入力フィールドをクリア
        set(edit_roi, 'String', '');
        
        % UserDataを更新
        ud.highlight_rois = [];
        set(fig, 'UserData', ud);
        
        % ステータスを更新
        set(txt_status, 'String', '強調表示をクリアしました');
        
        % 画像を再描画
        redrawImage();
    end

    function redrawImage()
        % データを取得
        ud = get(fig, 'UserData');
        ax_image_gcamp = ud.ax_image_gcamp;
        ax_image_rfp = ud.ax_image_rfp;
        img_gcamp = ud.img_gcamp;
        img_rfp = ud.img_rfp;
        stat_at = ud.stat_at;
        all_roi_numbers = ud.all_roi_numbers;
        highlight_rois = ud.highlight_rois;
        
        % GCaMP画像の描画
        cla(ax_image_gcamp);
        imshow(img_gcamp, [], 'Parent', ax_image_gcamp);
        colormap(ax_image_gcamp, gray(256));
        axis(ax_image_gcamp, 'image');
        hold(ax_image_gcamp, 'on');
        set(ax_image_gcamp, 'YDir', 'reverse');
        [H, W] = size(img_gcamp);
        
        % タイトルを追加
        title(ax_image_gcamp, 'GCaMP', 'FontSize', 12, 'FontWeight', 'bold');
        
        % RFP画像の描画（存在する場合）
        if ~isempty(img_rfp) && ~isempty(ax_image_rfp)
            cla(ax_image_rfp);
            imshow(img_rfp, [], 'Parent', ax_image_rfp);
            colormap(ax_image_rfp, gray(256));
            axis(ax_image_rfp, 'image');
            hold(ax_image_rfp, 'on');
            set(ax_image_rfp, 'YDir', 'reverse');
            [H_rfp, W_rfp] = size(img_rfp);
            
            % タイトルを追加
            title(ax_image_rfp, 'RFP', 'FontSize', 12, 'FontWeight', 'bold');
        end
        
        % 全てのROIを描画
        num_rois = numel(all_roi_numbers);
        colors = lines(num_rois);
        
        for i = 1:num_rois
            r = all_roi_numbers(i);
            if r < 1 || r > ud.nStat
                continue;
            end
            
            % 強調表示対象かチェック
            is_highlighted = ~isempty(highlight_rois) && any(highlight_rois == r);
            
            try
                S = stat_at(r);
                x = double(S.xpix(:));
                y = double(S.ypix(:));
                
                % GCaMP画像にROIを描画
                x_gcamp = min(max(x, 0.5), W+0.5);
                y_gcamp = min(max(y, 0.5), H+0.5);
                
                if numel(x_gcamp) >= 3
                    k = boundary(x_gcamp, y_gcamp, 0.8);
                    px_gcamp = x_gcamp(k); py_gcamp = y_gcamp(k);
                    
                    if is_highlighted
                        fill(ax_image_gcamp, px_gcamp, py_gcamp, 'r', ...
                             'EdgeColor', 'r', 'LineWidth', 2, 'FaceAlpha', 0.3);
                    else
                        fill(ax_image_gcamp, px_gcamp, py_gcamp, colors(i,:), ...
                             'EdgeColor', colors(i,:), 'LineWidth', 2, 'FaceAlpha', 0.15);
                        
                        mask = poly2mask(px_gcamp, py_gcamp, H, W);
                        statsC = regionprops(mask, 'Centroid');
                        if ~isempty(statsC)
                            c = statsC(1).Centroid;
                        else
                            c = [mean(px_gcamp) mean(py_gcamp)];
                        end
                        
                        text(ax_image_gcamp, c(1), c(2), sprintf('%d', r), ...
                             'Color', 'w', 'FontSize', 7, 'FontWeight', 'bold', ...
                             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                             'BackgroundColor', 'none', 'EdgeColor', 'none');
                    end
                end
                
                % RFP画像にROIを描画（存在する場合）
                if ~isempty(img_rfp) && ~isempty(ax_image_rfp)
                    x_rfp = min(max(x, 0.5), W_rfp+0.5);
                    y_rfp = min(max(y, 0.5), H_rfp+0.5);
                    
                    if numel(x_rfp) >= 3
                        k_rfp = boundary(x_rfp, y_rfp, 0.8);
                        px_rfp = x_rfp(k_rfp); py_rfp = y_rfp(k_rfp);
                        
                        if is_highlighted
                            fill(ax_image_rfp, px_rfp, py_rfp, 'r', ...
                                 'EdgeColor', 'r', 'LineWidth', 2, 'FaceAlpha', 0.3);
                        else
                            fill(ax_image_rfp, px_rfp, py_rfp, colors(i,:), ...
                                 'EdgeColor', colors(i,:), 'LineWidth', 2, 'FaceAlpha', 0.15);
                            
                            mask_rfp = poly2mask(px_rfp, py_rfp, H_rfp, W_rfp);
                            statsC_rfp = regionprops(mask_rfp, 'Centroid');
                            if ~isempty(statsC_rfp)
                                c_rfp = statsC_rfp(1).Centroid;
                            else
                                c_rfp = [mean(px_rfp) mean(py_rfp)];
                            end
                            
                            text(ax_image_rfp, c_rfp(1), c_rfp(2), sprintf('%d', r), ...
                                 'Color', 'w', 'FontSize', 7, 'FontWeight', 'bold', ...
                                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                                 'BackgroundColor', 'none', 'EdgeColor', 'none');
                        end
                    end
                end
            catch ME
                warning('ROI %d の描画に失敗: %s', r, ME.message);
                continue;
            end
        end
        
        % 軸を固定
        xlim(ax_image_gcamp, [0.5, W+0.5]);
        ylim(ax_image_gcamp, [0.5, H+0.5]);
        
        if ~isempty(img_rfp) && ~isempty(ax_image_rfp)
            xlim(ax_image_rfp, [0.5, W_rfp+0.5]);
            ylim(ax_image_rfp, [0.5, H_rfp+0.5]);
        end
        
        drawnow;
    end

    function showTraces()
        % データを取得
        ud = get(fig, 'UserData');
        txt_status = ud.txt_status;
        highlight_rois = ud.highlight_rois;
        csvFile = ud.csvFile;
        
        if isempty(highlight_rois)
            set(txt_status, 'String', 'ROI番号を入力してください');
            return;
        end
        
        % 波形データを表示
        csvFile_current = csvFile;
        if exist(csvFile_current, 'file') ~= 2
            % CSVファイルが見つからない場合、選択ダイアログを表示
            [file, path] = uigetfile({'*.csv', 'CSV files (*.csv)'; '*.*', 'All files'}, ...
                                     'presentation.csv を選択');
            if isequal(file, 0)
                set(txt_status, 'String', 'CSVファイルが選択されませんでした');
                return;
            end
            csvFile_current = fullfile(path, file);
            % UserDataのcsvFileも更新
            ud.csvFile = csvFile_current;
            set(fig, 'UserData', ud);
        end
        
        try
            csv_data = readmatrix(csvFile_current);
            csv_roi_col = csv_data(3:end, end-3);
            
            % 波形データ用のFigure
            fig_trace = figure('Name', 'ROI Traces', ...
                               'Position', [200, 200, 1000, 600], ...
                               'Color', 'w');
            
            num_traces = numel(highlight_rois);
            colors_trace = lines(num_traces);
            traces_shown = 0;
            
            for idx = 1:num_traces
                roi = highlight_rois(idx);
                row_idx = find(csv_roi_col == roi, 1);
                
                if isempty(row_idx)
                    warning('ROI %d の波形データが見つかりません', roi);
                    continue;
                end
                
                traces_shown = traces_shown + 1;
                
                % 波形データを取得（最後の列以外）
                trace_data = csv_data(row_idx + 2, 1:end-4);
                
                % サブプロットで表示
                subplot(num_traces, 1, traces_shown);
                plot(trace_data, 'Color', colors_trace(idx, :), 'LineWidth', 1.5);
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
            warning('波形データの読み込みに失敗: %s', ME.message);
            set(txt_status, 'String', sprintf('波形データの表示に失敗: %s', ME.message));
        end
    end

    function closeFigure(~, ~)
        % Figureを閉じる
        delete(fig);
    end

    % 初期表示
    redrawImage();
    
    % ===== ネスト関数 =====
    
    function [control_panel, edit_roi, txt_status] = createControlPanel(parent_fig, x_start, width, csv_file)
        % createControlPanel - コントロールパネルを作成
        
        % コントロールパネルを作成
        control_panel = uipanel('Parent', parent_fig, ...
                                'Title', '', ...
                                'Position', [x_start, 0.1, width, 0.85], ...
                                'Units', 'normalized', ...
                                'BackgroundColor', [0.98, 0.98, 0.98], ...
                                'BorderType', 'none', ...
                                'Visible', 'on');
        
        % 固定ピクセル位置で要素を配置（元の動作していたコード）
        % ROI番号入力フィールドのラベル
        uicontrol('Parent', control_panel, ...
                  'Style', 'text', ...
                  'String', 'ROI番号（カンマ区切り）:', ...
                  'Position', [10, 380, 280, 25], ...
                  'Units', 'pixels', ...
                  'HorizontalAlignment', 'left', ...
                  'FontSize', 11, ...
                  'BackgroundColor', [0.95, 0.95, 0.95], ...
                  'Visible', 'on');
        
        % ROI番号入力フィールド
        edit_roi = uicontrol('Parent', control_panel, ...
                             'Style', 'edit', ...
                             'String', '', ...
                             'Position', [10, 350, 280, 30], ...
                             'Units', 'pixels', ...
                             'FontSize', 11, ...
                             'BackgroundColor', 'w', ...
                             'HorizontalAlignment', 'left', ...
                             'Visible', 'on');
        
        % 表示ボタン
        btn_width = 135;  % ボタンの幅（ピクセル単位）
        uicontrol('Parent', control_panel, ...
                  'Style', 'pushbutton', ...
                  'String', '表示', ...
                  'Position', [10, 300, btn_width, 40], ...
                  'Units', 'pixels', ...
                  'FontSize', 12, ...
                  'FontWeight', 'bold', ...
                  'Visible', 'on', ...
                  'Callback', @(~,~) updateDisplay());
        
        % クリアボタン
        uicontrol('Parent', control_panel, ...
                  'Style', 'pushbutton', ...
                  'String', 'クリア', ...
                  'Position', [155, 300, btn_width, 40], ...
                  'Units', 'pixels', ...
                  'FontSize', 12, ...
                  'Visible', 'on', ...
                  'Callback', @(~,~) clearDisplay());
        
        % 波形表示ボタン
        uicontrol('Parent', control_panel, ...
                  'Style', 'pushbutton', ...
                  'String', '波形を表示', ...
                  'Position', [10, 240, 280, 40], ...
                  'Units', 'pixels', ...
                  'FontSize', 12, ...
                  'Visible', 'on', ...
                  'Callback', @(~,~) showTraces());
        
        % ステータステキスト
        if exist(csv_file, 'file') == 2
            status_msg = 'ROI番号を入力して「表示」をクリック';
        else
            status_msg = 'CSVファイルが見つかりません。波形表示にはCSVが必要です。';
        end
        
        txt_status = uicontrol('Parent', control_panel, ...
                               'Style', 'text', ...
                               'String', status_msg, ...
                               'Position', [10, 150, 280, 60], ...
                               'Units', 'pixels', ...
                               'FontSize', 10, ...
                               'ForegroundColor', [0.2, 0.2, 0.2], ...
                               'BackgroundColor', [0.95, 0.95, 0.95], ...
                               'HorizontalAlignment', 'left', ...
                               'Visible', 'on');
    end  % createControlPanel
    
end  % ROIs_in_Image_Interactive

