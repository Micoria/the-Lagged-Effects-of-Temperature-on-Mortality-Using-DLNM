### --- 0. 设置工作环境和加载包 ---
setwd("/Users/micoria/Documents/research_projects/DLNM")
library(dplyr)
library(lubridate)
library(zoo)
library(dlnm)
library(splines)
library(ggplot2)

### --- 1. 读取气温数据并预处理 ---
df <- read.csv("data/amd_temp.csv")
df$Date <- as.Date(df$Date, format = "%d-%m-%Y")
df <- df %>% arrange(Date)

# 计算 Temp_Mean
df$Temp.Max <- as.numeric(df$Temp.Max)
df$Temp.Min <- as.numeric(df$Temp.Min)
df$Temp_Mean <- (df$Temp.Max + df$Temp.Min) / 2

# 筛选 2013–2018 财政年数据
df <- df %>%
  filter(Date >= as.Date("2013-04-01") & Date <= as.Date("2018-03-31")) %>%
  select(Date, Temp_Mean)

# 加入财政年度
df <- df %>%
  mutate(
    Year = year(Date),
    Month = month(Date),
    Fiscal_Year = ifelse(Month >= 4,
                         paste0(Year, "-", Year + 1),
                         paste0(Year - 1, "-", Year))
  )

### --- 2. 导入死亡数据并模拟每日死亡数 ---
death <- data.frame(
  Fiscal_Year = c("2013-2014", "2014-2015", "2015-2016", "2016-2017", "2017-2018"),
  Total_Deaths = c(44769, 47638, 47346, 50951, 50427)  # 单位已是“人”
)

df <- df %>%
  left_join(death, by = "Fiscal_Year") %>%
  group_by(Fiscal_Year) %>%
  mutate(
    Days_in_Year = n(),
    Baseline_Deaths_Per_Day = Total_Deaths / Days_in_Year,
    MMT = 27.5,
    Risk_Modifier = 1 + 0.03 * abs(Temp_Mean - MMT),
    Expected_Deaths = Baseline_Deaths_Per_Day * Risk_Modifier
  ) %>%
  ungroup()

set.seed(42)
df$Deaths <- rpois(n = nrow(df), lambda = df$Expected_Deaths)

### --- 3. 构建 DLNM 模型 ---
cb_temp <- crossbasis(
  df$Temp_Mean,
  lag = 21,
  argvar = list(fun = "ns", df = 4),
  arglag = list(fun = "ns", df = 3)
)

model <- glm(
  Deaths ~ cb_temp + ns(as.numeric(Date), df = 8 * 5),
  family = quasipoisson(),
  data = df
)

### --- 4. 生成 cumulative RR（修复！）---
pred_cum <- crosspred(
  cb_temp, model,
  cen = 27.5,
  at = quantile(df$Temp_Mean, probs = seq(0.05, 0.95, by = 0.01)),
  cumul = TRUE  # ✅ 关键参数
)

### --- 5. 可视化并保存图像 ---

# 5a. 保存 contour 图
pdf("DLNM_contour_plot.pdf", width = 8, height = 5)
plot(pred_cum, "contour", 
     xlab = "Temperature (°C)", ylab = "Lag (days)",
     key.title = title("RR"), main = "Temperature–Lag–Mortality Risk Surface")
dev.off()

# 5b. 保存 cumulative RR 图
pdf("DLNM_cumulative_plot.pdf", width = 7, height = 5)
plot(pred_cum, "overall", xlab = "Temperature (°C)", ylab = "Cumulative RR",
     main = "Cumulative Exposure–Response")
dev.off()


pred_cum$predvar

# 找最接近 32°C 的温度点
closest_temp <- pred_cum$predvar[which.min(abs(pred_cum$predvar - 32))]

# 保存滞后响应图
png("DLNM_lag_response_32C.png", width = 800, height = 600)
plot(pred_cum, "slices", var = closest_temp,
     main = paste0("Lag–response at ", round(closest_temp, 2), "°C"),
     xlab = "Lag (days)", ylab = "RR")
dev.off()


### --- 6. 导出 cumulative RR 表格 ---
rr_table <- data.frame(
  Temp = pred_cum$predvar,
  RR = pred_cum$allRRfit,
  LowCI = pred_cum$allRRlow,
  HighCI = pred_cum$allRRhigh
)

write.csv(rr_table, "DLNM_Cumulative_RR_Table.csv", row.names = FALSE)

### --- 7. 归因死亡数（可选）---
### --- 7. 归因风险分析 ---

# 重新运行预测，只对唯一温度值生成 RR
pred_all <- crosspred(cb_temp, model, cen = 27.5,
                      at = unique(df$Temp_Mean), cumul = TRUE)

# 构建映射表：每个 Temp_Mean -> cumulative RR
rr_map <- data.frame(
  Temp_Mean = pred_all$predvar,
  RR = pred_all$allRRfit
)

# 把 RR 映射回原始 df
df <- df %>% left_join(rr_map, by = "Temp_Mean")

# 计算 AF 和归因死亡
df <- df %>%
  mutate(
    AF = 1 - (1 / RR),
    Attr_Deaths = Deaths * AF
  )

# 汇总归因结果
total_deaths <- sum(df$Deaths, na.rm = TRUE)
total_attr <- sum(df$Attr_Deaths, na.rm = TRUE)
AF_total <- total_attr / total_deaths

cat("\n======= 归因风险分析结果 =======\n")
cat("总死亡人数：", round(total_deaths), "\n")
cat("归因于非最适温度的死亡人数：", round(total_attr), "\n")
cat("归因分数（Attributable Fraction）：", round(AF_total * 100, 2), "%\n")

# 导出结果
write.csv(df[, c("Date", "Temp_Mean", "Deaths", "AF", "Attr_Deaths")],
          "DLNM_Attributable_Deaths.csv", row.names = FALSE)
