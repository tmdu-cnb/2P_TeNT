%% === しきい値入力（mean_int） ===
ans_thr = inputdlg({'Intensity threshold on mean_int (e.g., 50):'}, ...
                   'Intensity gating for correlation', 1, {'50'});
if isempty(ans_thr)
    thr_int = 50; % Cancel されたらデフォルト
else
    thr_int = str2double(ans_thr{1});
    if isnan(thr_int), thr_int = 50; end
end

%% === ROIのmean_intを読み込み、F_signal2の各行に対応付け ===
Tmean   = readtable(fullfile('result','roi_intensity_summary.csv')); % cols: roi, area_px, mean_int
roi_vec = Tmean.roi;                      % statの1始まりインデックス
mean_vec= Tmean.mean_int;

% --- 末尾列からROI番号を取り出して double ベクトルに正規化 ---
if istable(F_signal2)
    roi_in_use = F_signal2{:, end};      % table列の“中身”
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
mask_high = ~isnan(mean_for_row) & (mean_for_row >= thr_int);
mask_low  = ~isnan(mean_for_row) & (mean_for_row <= thr_int);

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
        thr_int, n_high);
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
        thr_int, n_low);
end
