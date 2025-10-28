# 取出样本坐标
coords <- RPCA$dtt1_T0$biplot[, c("PC1", "PC2","PC3")]

# 计算每个点到总体中心的距离
center <- colMeans(coords)
dist_to_center <- sqrt(rowSums((coords - center)^2))

# 找出距离最大的样本（偏移点）
outliers <- RPCA$dtt1_T0$biplot[order(dist_to_center, decreasing = TRUE), ]
strange_sample <- head(outliers[, c("SampleID", "RUN", "Project", "Batch_ID","TimePoint")], 15)

# 假设 strange 数据框已经存在（包含 SampleID 列）
strange_ids <- strange_sample$SampleID
strange_bath <- 
# 从 Report 中筛出这些样本的测序深度信息
depth_info <- Report %>%
  filter(`sample-id` %in% strange_ids) %>%
  select(`sample-id`, input, filtered, denoised, non.chimeric = `non-chimeric`)

meta_info <- meta %>%
  mutate(Global_ID = paste0("np",Global_ID)) %>%
  filter(Global_ID %in% strange_ids)

R60 <- meta %>%
  mutate(Global_ID = paste0("np",Global_ID)) %>%
  filter(RUN == "060")

!R60$Global_ID %in% meta_info$Global_ID

differ_id <- R60$Global_ID[!R60$Global_ID %in% meta_info$Global_ID]
print(differ_id)
# 查看前几行
head(depth_info, 50)


# ---- 取对象 ----
fit <- res_asv_time_proj$fit   # 改成你当前的对象

# ---- 自动辨别矩阵朝向 ----
feat_names <- unique(fit$results$feature)
F <- fit$fitted
R <- fit$residuals
if (nrow(F) == length(feat_names)) { F <- t(F); R <- t(R) }   # 转置：行=样本,列=特征
stopifnot(ncol(F) == length(feat_names), ncol(R) == length(feat_names))

# ---- ① 每个特征的诊断指标 ----
message("▶ 计算特征级残差诊断指标...")

# 为防止内存溢出，这里做矢量化计算
mean_fit <- colMeans(F, na.rm = TRUE)
sd_resid <- apply(R, 2, sd, na.rm = TRUE)
rmse     <- sqrt(colMeans(R^2, na.rm = TRUE))

# 使用 tryCatch 防止部分特征导致崩溃
safe_cor <- function(x, y) {
  res <- try(suppressWarnings(cor(abs(x), y, use = "pairwise")), silent = TRUE)
  if (inherits(res, "try-error") || is.na(res)) return(NA_real_) else return(res)
}

cor_abs <- sapply(seq_len(ncol(R)), function(j) safe_cor(R[, j], F[, j]))

# 异方差检验: 随机抽样部分特征减少负担
message("▶ 进行异方差性检验 (Levene-like)...")
max_feats <- min(300, ncol(R))  # 控制计算量
sample_feats <- sample(seq_len(ncol(R)), max_feats)

p_hetero <- rep(NA_real_, ncol(R))
for (j in sample_feats) {
  dat <- data.frame(resid = abs(R[, j]), fit = F[, j])
  # 小样本安全检查
  if (nrow(dat) > 10 && sd(dat$fit, na.rm = TRUE) > 0) {
    fit_lm <- try(lm(resid ~ fit, data = dat), silent = TRUE)
    if (!inherits(fit_lm, "try-error")) {
      s <- summary(fit_lm)
      p_hetero[j] <- s$coefficients[2, 4]
    }
  }
}

q_hetero <- p.adjust(p_hetero, "BH")
diag_df <- data.frame(
  feature = feat_names, mean_fit, sd_resid, rmse, cor_abs,
  p_hetero, q_hetero,
  hetero_flag = ifelse(!is.na(q_hetero) & q_hetero < 0.05 & abs(cor_abs) > 0.2, TRUE, FALSE)
) %>% arrange(desc(hetero_flag), desc(rmse))

# ---- ② 图1：每个特征的残差波动 vs 拟合均值 ----
message("▶ 绘制图1：每特征残差SD vs 平均拟合值")
ggplot(diag_df, aes(mean_fit, sd_resid)) +
  geom_point(alpha = .6, size = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "red") +
  labs(
    title = "Per-feature residual SD vs mean fitted",
    x = "Mean fitted value",
    y = "Residual SD"
  ) +
  theme_minimal()

# ---- ③ 图2：整体残差趋势图（采样 + GAM）----
message("▶ 绘制图2：全体 pooled 残差 vs 拟合值 (抽样)")
pooled <- data.frame(fitted = as.vector(F), residual = as.vector(R))
n_sample <- min(20000, nrow(pooled))
set.seed(42)
pooled_sub <- pooled[sample(seq_len(nrow(pooled)), n_sample), ]

ggplot(pooled_sub, aes(fitted, residual)) +
  geom_bin2d(bins = 60) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "red") +
  labs(
    title = "Pooled residuals vs fitted (sampled, all features)",
    subtitle = paste("Sampled", n_sample, "points"),
    x = "Fitted value",
    y = "Residual"
  ) +
  theme_minimal()

# ---- ④ 查看最可疑特征 ----
message("▶ 输出最可疑的特征（异方差 + 高RMSE）")
head(diag_df %>% filter(hetero_flag) %>%
       select(feature, rmse, cor_abs, q_hetero), 10)

suspect_features <- diag_df %>%
  mutate(
    flag_rmse = rmse > quantile(rmse, 0.95, na.rm = TRUE),
    flag_cor  = abs(cor_abs) > 0.2,
    flag_hetero = !is.na(q_hetero) & q_hetero < 0.1
  ) %>%
  mutate(
    suspicious = flag_rmse + flag_cor + flag_hetero,
    suspect_flag = suspicious >= 2   # 同时满足≥2项
  ) %>%
  arrange(desc(suspicious), desc(rmse))
table(suspect_features$suspect_flag)

both_bad <- merged_asv_M$all_with_both %>%
  filter(Remarkable) %>%                                # 只取 both_remarkable 的特征
  left_join(suspect_features %>% select(feature, suspect_flag), by="feature") %>%
  mutate(model_quality = ifelse(suspect_flag, "Suspect", "Reliable"))

table(both_bad$model_quality)

ggplot(both_bad, aes(x = log2FC72, y = log2FC_real, color = model_quality)) +
  +     geom_point(size=2) +
  +     geom_hline(yintercept = 0, linetype="dashed") +
  +     geom_vline(xintercept = 0, linetype="dashed") +
  +     scale_color_manual(values=c("Reliable"="steelblue", "Suspect"="red")) +
  +     labs(title="Both-Remarkable features: Model reliability overlay",
             +          x="Model-estimated log2FC72", y="Observed log2FC_real") +
  +     theme_minimal()

