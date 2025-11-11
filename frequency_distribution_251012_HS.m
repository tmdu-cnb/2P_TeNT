% histogram_analysis_grouped.m
% -------------------------------------------------------------------------
% 全体に加えて、roi_intensity_summary.csv の mean_int しきい値で
% 高輝度/低輝度グループの頻度分布も出力
% -------------------------------------------------------------------------
close all; clc;
if ~exist('result','dir'), mkdir result; end

%% -------- Default parameters -------------------------------------------
defVals = struct( ...
    'amp',  [1     20     0.5 ], ...   % [binInterval max min]
    'area', [2     40     1   ], ...
    'freq', [0.01 0.2    0.01], ...
    'wid',  [0.5   10     0.25]);

%% -------- Single input dialog ------------------------------------------
prompt = {
    'Amplitude  –  bin interval', 'Amplitude  –  MAX', 'Amplitude  –  MIN', ...
    'Area       –  bin interval', 'Area       –  MAX', 'Area       –  MIN', ...
    'Frequency  –  bin interval', 'Frequency  –  MAX', 'Frequency  –  MIN', ...
    'Width      –  bin interval', 'Width      –  MAX', 'Width      –  MIN', ...
    'mean\_int threshold (for HIGH/LOW split)'};
defInput = [ ...
    arrayfun(@num2str, [defVals.amp defVals.area defVals.freq defVals.wid], 'UniformOutput', false), ...
    {'50'}];         % mean_int の初期しきい値
answer = inputdlg(prompt, 'Histogram settings (all features + threshold)', [1 48], defInput);

if isempty(answer)  % cancelled → keep defaults
    params = [defVals.amp; defVals.area; defVals.freq; defVals.wid];
    thr_int = 50;
else
    numAns  = cellfun(@str2double, answer);
    params  = reshape(numAns(1:12), 3, []).';   % each row: [interval max min]
    thr_int = numAns(13);
end

%% -------- Overall histograms (existing) ---------------------------------
ampData  = readmatrix(fullfile('result','all_amp.csv'));
areaData = readmatrix(fullfile('result','all_area.csv'));
freqData = readmatrix(fullfile('result','Freq.csv'));
widData  = readmatrix(fullfile('result','Width.csv'));

makeHistogram(ampData,  'Amplitude',               fullfile('result','amp_frequency_distribution.csv'),  params(1,1), params(1,2), params(1,3));
makeHistogram(areaData, 'Area of Ca2+ transient',  fullfile('result','area_frequency_distribution.csv'), params(2,1), params(2,2), params(2,3));
makeHistogram(freqData, 'Frequency (Hz)',          fullfile('result','freq_frequency_distribution.csv'), params(3,1), params(3,2), params(3,3));
makeHistogram(widData,  'Width',                   fullfile('result','width_frequency_distribution.csv'),params(4,1), params(4,2), params(4,3));

%% -------- Grouped histograms (HIGH / LOW by mean_int) -------------------
% 前工程で出力したグループ別CSVを使用（無ければスキップ）
filesHigh = struct('amp','all_amp_high.csv','area','all_area_high.csv','freq','Freq_high.csv','wid','Width_high.csv');
filesLow  = struct('amp','all_amp_low.csv', 'area','all_area_low.csv', 'freq','Freq_low.csv', 'wid','Width_low.csv');

hasHigh = all(cellfun(@(f) exist(fullfile('result',f),'file')==2, struct2cell(filesHigh)));
hasLow  = all(cellfun(@(f) exist(fullfile('result',f),'file')==2, struct2cell(filesLow)));

% ログ（参考）
fprintf('mean_int threshold = %.3f\n', thr_int);
if ~hasHigh, warning('HIGH グループのCSVが揃っていません。前段の解析を実行してください。'); end
if ~hasLow,  warning('LOW グループのCSVが揃っていません。前段の解析を実行してください。');  end

if hasHigh
    H_amp  = readmatrix(fullfile('result',filesHigh.amp));
    H_area = readmatrix(fullfile('result',filesHigh.area));
    H_freq = readmatrix(fullfile('result',filesHigh.freq));
    H_wid  = readmatrix(fullfile('result',filesHigh.wid));

    makeHistogram(H_amp,  sprintf('Amplitude  (HIGH: mean\\_int \\ge %.3g)',thr_int), ...
                  fullfile('result','amp_frequency_distribution_high.csv'),  params(1,1), params(1,2), params(1,3));
    makeHistogram(H_area, sprintf('Area of Ca2+ transient  (HIGH: mean\\_int \\ge %.3g)',thr_int), ...
                  fullfile('result','area_frequency_distribution_high.csv'), params(2,1), params(2,2), params(2,3));
    makeHistogram(H_freq, sprintf('Frequency (Hz)  (HIGH: mean\\_int \\ge %.3g)',thr_int), ...
                  fullfile('result','freq_frequency_distribution_high.csv'), params(3,1), params(3,2), params(3,3));
    makeHistogram(H_wid,  sprintf('Width  (HIGH: mean\\_int \\ge %.3g)',thr_int), ...
                  fullfile('result','width_frequency_distribution_high.csv'),params(4,1), params(4,2), params(4,3));
end

if hasLow
    L_amp  = readmatrix(fullfile('result',filesLow.amp));
    L_area = readmatrix(fullfile('result',filesLow.area));
    L_freq = readmatrix(fullfile('result',filesLow.freq));
    L_wid  = readmatrix(fullfile('result',filesLow.wid));

    makeHistogram(L_amp,  sprintf('Amplitude  (LOW: mean\\_int \\le %.3g)',thr_int), ...
                  fullfile('result','amp_frequency_distribution_low.csv'),  params(1,1), params(1,2), params(1,3));
    makeHistogram(L_area, sprintf('Area of Ca2+ transient  (LOW: mean\\_int \\le %.3g)',thr_int), ...
                  fullfile('result','area_frequency_distribution_low.csv'), params(2,1), params(2,2), params(2,3));
    makeHistogram(L_freq, sprintf('Frequency (Hz)  (LOW: mean\\_int \\le %.3g)',thr_int), ...
                  fullfile('result','freq_frequency_distribution_low.csv'), params(3,1), params(3,2), params(3,3));
    makeHistogram(L_wid,  sprintf('Width  (LOW: mean\\_int \\le %.3g)',thr_int), ...
                  fullfile('result','width_frequency_distribution_low.csv'),params(4,1), params(4,2), params(4,3));
end

%% ====== Local helper functions (keep at end of script) =================
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
