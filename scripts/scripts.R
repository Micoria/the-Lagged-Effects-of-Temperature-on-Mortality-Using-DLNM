
setwd("/Users/micoria/Documents/research_projects/DLNM")
library(dplyr)
library(lubridate)
library(zoo)
library(dlnm)
library(splines)
library(ggplot2)


df <- read.csv("data/amd_temp.csv")
df$Date <- as.Date(df$Date, format = "%d-%m-%Y")
df <- df %>% arrange(Date)


df$Temp.Max <- as.numeric(df$Temp.Max)
df$Temp.Min <- as.numeric(df$Temp.Min)
df$Temp_Mean <- (df$Temp.Max + df$Temp.Min) / 2


df <- df %>%
  filter(Date >= as.Date("2013-04-01") & Date <= as.Date("2018-03-31")) %>%
  select(Date, Temp_Mean)


df <- df %>%
  mutate(
    Year = year(Date),
    Month = month(Date),
    Fiscal_Year = ifelse(Month >= 4,
                         paste0(Year, "-", Year + 1),
                         paste0(Year - 1, "-", Year))
  )


death <- data.frame(
  Fiscal_Year = c("2013-2014", "2014-2015", "2015-2016", "2016-2017", "2017-2018"),
  Total_Deaths = c(44769, 47638, 47346, 50951, 50427)  
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


pred_cum <- crosspred(
  cb_temp, model,
  cen = 27.5,
  at = quantile(df$Temp_Mean, probs = seq(0.05, 0.95, by = 0.01)),
  cumul = TRUE 
)


pdf("DLNM_contour_plot.pdf", width = 8, height = 5)
plot(pred_cum, "contour", 
     xlab = "Temperature (°C)", ylab = "Lag (days)",
     key.title = title("RR"), main = "Temperature–Lag–Mortality Risk Surface")
dev.off()


pdf("DLNM_cumulative_plot.pdf", width = 7, height = 5)
plot(pred_cum, "overall", xlab = "Temperature (°C)", ylab = "Cumulative RR",
     main = "Cumulative Exposure–Response")
dev.off()


pred_cum$predvar

closest_temp <- pred_cum$predvar[which.min(abs(pred_cum$predvar - 32))]

png("DLNM_lag_response_32C.png", width = 800, height = 600)
plot(pred_cum, "slices", var = closest_temp,
     main = paste0("Lag–response at ", round(closest_temp, 2), "°C"),
     xlab = "Lag (days)", ylab = "RR")
dev.off()


rr_table <- data.frame(
  Temp = pred_cum$predvar,
  RR = pred_cum$allRRfit,
  LowCI = pred_cum$allRRlow,
  HighCI = pred_cum$allRRhigh
)

write.csv(rr_table, "DLNM_Cumulative_RR_Table.csv", row.names = FALSE)



pred_all <- crosspred(cb_temp, model, cen = 27.5,
                      at = unique(df$Temp_Mean), cumul = TRUE)

rr_map <- data.frame(
  Temp_Mean = pred_all$predvar,
  RR = pred_all$allRRfit
)

df <- df %>% left_join(rr_map, by = "Temp_Mean")

df <- df %>%
  mutate(
    AF = 1 - (1 / RR),
    Attr_Deaths = Deaths * AF
  )

total_deaths <- sum(df$Deaths, na.rm = TRUE)
total_attr <- sum(df$Attr_Deaths, na.rm = TRUE)
AF_total <- total_attr / total_deaths

cat("\n======= Attribution risk analysis results =======\n")
cat("Alive：", round(total_deaths), "\n")
cat("Estimated deaths due to non-optimal temperatures：", round(total_attr), "\n")
cat("Attributable fraction（Attributable Fraction）：", round(AF_total * 100, 2), "%\n")

write.csv(df[, c("Date", "Temp_Mean", "Deaths", "AF", "Attr_Deaths")],
          "DLNM_Attributable_Deaths.csv", row.names = FALSE)
