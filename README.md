# 2P_TeNT - Two-Photon Calcium Imaging Analysis Pipeline

二光子カルシウムイメージング解析用のMATLABスクリプト集

## ファイル一覧

### メイン処理スクリプト
- **P10_subtract_251019_XYratio_remove_HS.m** - 前処理とXY比補正・ノイズ除去
- **P10_detectROI__251019_HS.m** - ROI検出とフィルタリング（メインスクリプト）
- **Calcium_properties_251019_HS.m** - カルシウムシグナルの特性解析（振幅、頻度、面積など）
- **ROIs_in_Image_251027_HS.m** - ROIの可視化と輝度測定

### 解析スクリプト
- **Correlation_2501019_HS.m** - ROI間の相関解析
- **frequency_distribution_251012_HS.m** - 頻度分布解析
- **SeparationROIs_by_Redintensity_251027_HS.m** - 赤色蛍光強度によるROI分類

## 主な機能

### データ処理
- Suite2Pで検出されたROIのフィルタリング
- バックグラウンド除去とシグナル正規化
- カルシウムトランジェントの自動検出
- ノイズROIの除去（平均値、最大値、ピーク幅ベース）

### 解析
- ピーク振幅、頻度、持続時間の定量化
- ROI間の相関解析
- 輝度ベースのROI分類（高輝度/低輝度）
- 時系列波形の可視化

### 可視化
- ROIの空間配置の可視化（カラーマップ付き）
- 個別ROIの波形トレース
- 統計データのエクスポート（CSV形式）

## 使用方法

1. Suite2Pで解析したデータを用意
2. `P10_detectROI__251019_HS.m` を実行してROIをフィルタリング
3. `Calcium_properties_251019_HS.m` でカルシウム特性を解析
4. `ROIs_in_Image_251027_HS.m` でROIを可視化

## 出力ファイル

- `presentation.csv` - フィルタリング後のROIデータ
- `roi_intensity_summary.csv` - ROI輝度情報
- `removed_original_ROIs_F6_to_F_signal5_unique.csv` - 除外されたROIリスト
- 各種統計CSV（振幅、頻度、幅など）
- EPS形式の波形図

## 必要な環境

- MATLAB (Signal Processing Toolbox, Image Processing Toolbox)
- Suite2P出力データ（stat.mat, F.npy, Fneu.npyなど）

## ライセンス

MIT License

