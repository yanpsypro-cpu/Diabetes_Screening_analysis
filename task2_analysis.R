# Task 2: Comparison with findings from the previous screening programme

# Questions addressed:
#   1. Is the diabetes proportion different from 14%?
#   2. Is the mean BMI different from 28 kg/m^2?
#   3. Is the mean age different from 55 years?
#   4. Is smoking associated with diabetes diagnosis?
#   5. Can age, sex or BMI help explain the findings?
#
# a one-sample proportion z test, one-sample and Welch two sample t tests,
# Pearson chi-squared tests, odds and risk ratios for 2 by 2 tables (Epi::twoby2),
# and stratified 2 by 2 tables.

options(stringsAsFactors = FALSE)

input_candidates <- c(
  "assessment1/SHDS_Assessment1_data.csv",
  "SHDS_Assessment1_data.csv"
)
input_file <- input_candidates[file.exists(input_candidates)][1]

if (is.na(input_file)) {
  input_file <- file.choose()
}

# Save outputs in a folder next to the data file
output_dir <- file.path(dirname(input_file), "task2_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(input_file, check.names = FALSE)

# Remove the row identifier supplied as the first column.
if (ncol(dat) == 10 && names(dat)[1] %in% c("", "X", "...1", "X.1")) {
  dat <- dat[-1]
}

required_variables <- c("Diabetes_012", "BMI", "Age", "Smoker")
missing_variables <- setdiff(required_variables, names(dat))
if (length(missing_variables) > 0) {
  stop("Missing expected variables: ", paste(missing_variables, collapse = ", "))
}

if (!all(dat$Diabetes_012 %in% c(0, 1)) || !all(dat$Smoker %in% c(0, 1))) {
  stop("Diabetes_012 and Smoker must both be coded 0 or 1.")
}

alpha <- 0.05

# 2.1 Diabetes proportion compared with 14%

n_total <- nrow(dat)
n_diabetes <- sum(dat$Diabetes_012 == 1)
p_observed <- n_diabetes / n_total
p_null <- 0.14

# Check the normal approximation under the null hypothesis.
expected_diabetes_null <- n_total * p_null
expected_no_diabetes_null <- n_total * (1 - p_null)

se_proportion_null <- sqrt(p_null * (1 - p_null) / n_total)
z_proportion <- (p_observed - p_null) / se_proportion_null
p_value_proportion <- 2 * pnorm(-abs(z_proportion))

# A 95% confidence interval for the observed proportion.
se_proportion_observed <- sqrt(p_observed * (1 - p_observed) / n_total)
ci_proportion <- p_observed + c(-1, 1) * qnorm(1 - alpha / 2) * se_proportion_observed

# R's proportion test is included as a check of the manual calculation.
proportion_test_r <- prop.test(
  x = n_diabetes,
  n = n_total,
  p = p_null,
  alternative = "two.sided",
  correct = FALSE
)

proportion_results <- data.frame(
  n = n_total,
  diabetes_count = n_diabetes,
  observed_proportion = p_observed,
  null_proportion = p_null,
  expected_diabetes_under_null = expected_diabetes_null,
  expected_no_diabetes_under_null = expected_no_diabetes_null,
  standard_error_under_null = se_proportion_null,
  z_statistic = z_proportion,
  p_value = p_value_proportion,
  ci_95_lower = ci_proportion[1],
  ci_95_upper = ci_proportion[2]
)

write.csv(
  proportion_results,
  file.path(output_dir, "task2_diabetes_proportion_test.csv"),
  row.names = FALSE
)

# 2.2 Mean BMI compared with 28 kg/m^2

bmi_null <- 28
bmi_test <- t.test(
  dat$BMI,
  mu = bmi_null,
  alternative = "two.sided",
  conf.level = 0.95
)

bmi_results <- data.frame(
  variable = "BMI",
  unit = "kg/m^2",
  n = sum(!is.na(dat$BMI)),
  sample_mean = mean(dat$BMI, na.rm = TRUE),
  sample_sd = sd(dat$BMI, na.rm = TRUE),
  null_mean = bmi_null,
  mean_difference = mean(dat$BMI, na.rm = TRUE) - bmi_null,
  standard_error = sd(dat$BMI, na.rm = TRUE) / sqrt(sum(!is.na(dat$BMI))),
  t_statistic = unname(bmi_test$statistic),
  degrees_of_freedom = unname(bmi_test$parameter),
  p_value = bmi_test$p.value,
  ci_95_lower = bmi_test$conf.int[1],
  ci_95_upper = bmi_test$conf.int[2]
)

# 2.3 Mean age compared with 55 years

age_null <- 55
age_test <- t.test(
  dat$Age,
  mu = age_null,
  alternative = "two.sided",
  conf.level = 0.95
)

age_results <- data.frame(
  variable = "Age",
  unit = "years",
  n = sum(!is.na(dat$Age)),
  sample_mean = mean(dat$Age, na.rm = TRUE),
  sample_sd = sd(dat$Age, na.rm = TRUE),
  null_mean = age_null,
  mean_difference = mean(dat$Age, na.rm = TRUE) - age_null,
  standard_error = sd(dat$Age, na.rm = TRUE) / sqrt(sum(!is.na(dat$Age))),
  t_statistic = unname(age_test$statistic),
  degrees_of_freedom = unname(age_test$parameter),
  p_value = age_test$p.value,
  ci_95_lower = age_test$conf.int[1],
  ci_95_upper = age_test$conf.int[2]
)

mean_test_results <- rbind(bmi_results, age_results)
write.csv(
  mean_test_results,
  file.path(output_dir, "task2_bmi_age_tests.csv"),
  row.names = FALSE
)

# 2.4 Relationship between smoking and diabetes

smoking_table <- table(
  Smoker = factor(dat$Smoker, levels = c(0, 1), labels = c("Non-smoker", "Smoker")),
  Diabetes = factor(dat$Diabetes_012, levels = c(0, 1), labels = c("No diabetes", "Diabetes"))
)

smoking_row_percent <- prop.table(smoking_table, margin = 1) * 100

# The course uses the Pearson chi-squared test without Yates' correction.
smoking_chisq <- chisq.test(smoking_table, correct = FALSE)

# Direction and size of the association, expressed using the observed odds.
# Convert the counts to double precision before multiplication. This prevents
# integer overflow in a large table.
smoking_odds_ratio <-
  (as.numeric(smoking_table["Smoker", "Diabetes"]) *
     as.numeric(smoking_table["Non-smoker", "No diabetes"])) /
  (as.numeric(smoking_table["Smoker", "No diabetes"]) *
     as.numeric(smoking_table["Non-smoker", "Diabetes"]))

smoking_results <- data.frame(
  smoking_status = rownames(smoking_table),
  no_diabetes_n = as.integer(smoking_table[, "No diabetes"]),
  diabetes_n = as.integer(smoking_table[, "Diabetes"]),
  total_n = rowSums(smoking_table),
  diabetes_percent = as.numeric(smoking_row_percent[, "Diabetes"])
)

write.csv(
  smoking_results,
  file.path(output_dir, "task2_smoking_diabetes_summary.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame.matrix(smoking_chisq$expected),
  file.path(output_dir, "task2_smoking_diabetes_expected_counts.csv")
)

# 2.5 Effect size, confounding and possible explanations

library(Epi)

# 2 by 2 table ordered for twoby2 (smoker and diabetes first);
# stored as double to avoid integer overflow with large counts
make_2x2 <- function(d) {
  tab <- table(
    Smoker = factor(d$Smoker, levels = c(1, 0), labels = c("Smoker", "Non-smoker")),
    Diabetes = factor(d$Diabetes_012, levels = c(1, 0), labels = c("Diabetes", "No diabetes"))
  )
  storage.mode(tab) <- "double"
  tab
}

# (a) Risk ratio and odds ratio with 95% CIs
smoking_measures <- twoby2(make_2x2(dat), print = FALSE)$measures[1:2, ]
dimnames(smoking_measures) <- list(c("Risk ratio", "Odds ratio"),
                                   c("estimate", "ci_95_lower", "ci_95_upper"))

# (b) Age and sex of smokers compared with non-smokers
age_by_smoking_test <- t.test(Age ~ Smoker, data = dat)
sex_by_smoking_table <- table(Sex = dat$Sex, Smoker = dat$Smoker)
sex_by_smoking_test <- chisq.test(sex_by_smoking_table, correct = FALSE)

# Age groups for stratification and the trend test
dat$age_group <- cut(dat$Age, breaks = c(18, 45, 60, 75, 93),
                     right = FALSE, include.lowest = TRUE,
                     labels = c("18-44", "45-59", "60-74", "75-93"))

# (c) Smoking and diabetes within age groups and within sex
stratum_or <- function(d) {
  m <- twoby2(make_2x2(d), print = FALSE)$measures
  c(n = nrow(d), odds_ratio = m[2, 1], ci_95_lower = m[2, 2], ci_95_upper = m[2, 3])
}
or_by_age <- t(sapply(split(dat, dat$age_group), stratum_or))
or_by_sex <- t(sapply(split(dat, factor(dat$Sex, labels = c("Female", "Male"))), stratum_or))

# (d) Diabetes proportion by age group
age_diabetes_table <- table(dat$age_group, dat$Diabetes_012)
diabetes_by_age <- data.frame(
  age_group = rownames(age_diabetes_table),
  n = as.numeric(rowSums(age_diabetes_table)),
  diabetes_n = as.numeric(age_diabetes_table[, "1"]),
  diabetes_percent = 100 * as.numeric(prop.table(age_diabetes_table, 1)[, "1"])
)
# (e) BMI by diabetes status
bmi_by_diabetes_test <- t.test(BMI ~ Diabetes_012, data = dat)

write.csv(smoking_measures, file.path(output_dir, "task2_smoking_rr_or.csv"))
write.csv(or_by_age, file.path(output_dir, "task2_or_by_age.csv"))
write.csv(or_by_sex, file.path(output_dir, "task2_or_by_sex.csv"))
write.csv(diabetes_by_age, file.path(output_dir, "task2_diabetes_by_age.csv"), row.names = FALSE)

# Figures
# Figure 1: sample estimates and the values reported by the previous programme.
png(
  file.path(output_dir, "figure_5_sample_estimates_vs_previous_programme.png"),
  width = 1800,
  height = 700,
  res = 180
)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

plot(
  p_observed * 100, 1,
  xlim = range(c(ci_proportion * 100, p_null * 100)),
  ylim = c(0.7, 1.3),
  yaxt = "n",
  ylab = "",
  xlab = "Diabetes (%)",
  main = "Diabetes proportion",
  pch = 19,
  col = "#2C6E9B",
  cex = 1.3
)
segments(ci_proportion[1] * 100, 1, ci_proportion[2] * 100, 1, col = "#2C6E9B", lwd = 3)
abline(v = p_null * 100, col = "#B44A4A", lty = 2, lwd = 2)

plot(
  mean(dat$BMI), 1,
  xlim = range(c(bmi_test$conf.int, bmi_null)),
  ylim = c(0.7, 1.3),
  yaxt = "n",
  ylab = "",
  xlab = expression("BMI (kg/m"^2*")"),
  main = "Mean BMI",
  pch = 19,
  col = "#2C6E9B",
  cex = 1.3
)
segments(bmi_test$conf.int[1], 1, bmi_test$conf.int[2], 1, col = "#2C6E9B", lwd = 3)
abline(v = bmi_null, col = "#B44A4A", lty = 2, lwd = 2)

plot(
  mean(dat$Age), 1,
  xlim = range(c(age_test$conf.int, age_null)),
  ylim = c(0.7, 1.3),
  yaxt = "n",
  ylab = "",
  xlab = "Age (years)",
  main = "Mean age",
  pch = 19,
  col = "#2C6E9B",
  cex = 1.3
)
segments(age_test$conf.int[1], 1, age_test$conf.int[2], 1, col = "#2C6E9B", lwd = 3)
abline(v = age_null, col = "#B44A4A", lty = 2, lwd = 2)

dev.off()

# Figure 2: observed diabetes percentage by smoking status.
png(
  file.path(output_dir, "figure_6_diabetes_by_smoking_status.png"),
  width = 1250,
  height = 900,
  res = 180
)
par(mar = c(5, 5, 3, 1))
diabetes_by_smoking <- smoking_results$diabetes_percent
bar_positions <- barplot(
  diabetes_by_smoking,
  names.arg = smoking_results$smoking_status,
  col = c("grey70", "#2A9D8F"),
  border = "white",
  ylim = c(0, max(diabetes_by_smoking) * 1.22),
  ylab = "Diabetes diagnosis (%)",
  xlab = "Smoking status",
  main = "Diabetes diagnosis by smoking status"
)
text(
  x = bar_positions,
  y = diabetes_by_smoking,
  labels = paste0(round(diabetes_by_smoking, 1), "%"),
  pos = 3
)
dev.off()

# Figure 7: odds ratio for smoking, overall and within strata
or_plot <- rbind(smoking_measures["Odds ratio", ], or_by_age[, -1], or_by_sex[, -1])
rownames(or_plot) <- c("All visits", "Age 18 to 44", "Age 45 to 59", "Age 60 to 74",
                       "Age 75 or over", "Female", "Male")
or_y <- c(9, 7:4, 2:1)
or_col <- c("#264653", rep("#2A9D8F", 4), rep("#E76F51", 2))  # overall, age, sex

png(file.path(output_dir, "figure_7_smoking_or_by_stratum.png"), width = 1800, height = 800, res = 180)
par(mar = c(4.5, 8, 1, 1))
plot(or_plot[, 1], or_y, type = "n", log = "x", xlim = c(0.95, 1.7), ylim = c(0.5, 9.5),
     yaxt = "n", ylab = "", xlab = "Odds ratio for diabetes, smokers vs non-smokers (log scale)")
abline(v = 1, lty = 2, col = "grey55")
segments(or_plot[, 2], or_y, or_plot[, 3], or_y, col = or_col, lwd = 3)
points(or_plot[, 1], or_y, pch = 21, bg = or_col, col = "white", cex = 1.6)
axis(2, at = or_y, labels = rownames(or_plot), las = 1, tick = FALSE)
dev.off()

# Figure 8: diabetes proportion by age group with 95% CIs
p_age <- diabetes_by_age$diabetes_n / diabetes_by_age$n
half_width_age <- qnorm(1 - alpha / 2) * sqrt(p_age * (1 - p_age) / diabetes_by_age$n)
share_of_visits <- 100 * diabetes_by_age$n / sum(diabetes_by_age$n)
bar_col <- ifelse(p_age > p_null, "#2A9D8F", "grey78")  # highlight groups above 14%

png(file.path(output_dir, "figure_8_diabetes_by_age_group.png"), width = 1800, height = 800, res = 180)
par(mar = c(4.5, 4.5, 1, 1))
age_bars <- barplot(100 * p_age, names.arg = c("18 to 44", "45 to 59", "60 to 74", "75 or over"),
                    col = bar_col, border = NA, ylim = c(0, 30),
                    xlab = "Age group (years)", ylab = "Visits with diabetes (%)")
arrows(age_bars, 100 * (p_age - half_width_age), age_bars, 100 * (p_age + half_width_age),
       angle = 90, code = 3, length = 0.04, col = "grey25")
abline(h = 100 * p_null, lty = 2, lwd = 2, col = "#E76F51")
text(age_bars, 100 * p_age, paste0(round(share_of_visits, 1), "% of visits"),
     pos = 3, cex = 0.85, col = "grey25")
legend("topleft", legend = "Previous programme (14%)", lty = 2, lwd = 2,
       col = "#E76F51", bty = "n")
dev.off()

# Text log for review before writing the report

format_p_value <- function(p) {
  if (p < 0.001) {
    format(p, scientific = TRUE, digits = 3)
  } else {
    format(round(p, 3), nsmall = 3)
  }
}

sink(file.path(output_dir, "task2_results_log.txt"))

cat("TASK 2: COMPARISON WITH THE PREVIOUS SCREENING PROGRAMME\n\n")

cat("1. DIABETES PROPORTION COMPARED WITH 14%\n")
cat("H0: the population proportion with diabetes is 0.14\n")
cat("HA: the population proportion with diabetes is not 0.14\n")
cat("Observed count:", n_diabetes, "out of", n_total, "\n")
cat("Observed percentage:", round(100 * p_observed, 3), "%\n")
cat("Expected counts under H0:", expected_diabetes_null, "and", expected_no_diabetes_null, "\n")
cat("z statistic:", round(z_proportion, 4), "\n")
cat("p value:", format_p_value(p_value_proportion), "\n")
cat("95% CI for the observed percentage:",
    round(100 * ci_proportion[1], 3), "to", round(100 * ci_proportion[2], 3), "%\n\n")

cat("R prop.test output used as a calculation check\n")
print(proportion_test_r)

cat("\n2. MEAN BMI COMPARED WITH 28 KG/M^2\n")
cat("H0: the population mean BMI is 28 kg/m^2\n")
cat("HA: the population mean BMI is not 28 kg/m^2\n")
print(bmi_test)
cat("Observed mean difference:", round(bmi_results$mean_difference, 4), "kg/m^2\n\n")

cat("3. MEAN AGE COMPARED WITH 55 YEARS\n")
cat("H0: the population mean age is 55 years\n")
cat("HA: the population mean age is not 55 years\n")
print(age_test)
cat("Observed mean difference:", round(age_results$mean_difference, 4), "years\n\n")

cat("4. SMOKING AND DIABETES\n")
cat("H0: smoking status and diabetes diagnosis are independent\n")
cat("HA: smoking status and diabetes diagnosis are associated\n\n")
cat("Observed counts\n")
print(smoking_table)
cat("\nRow percentages\n")
print(round(smoking_row_percent, 2))
cat("\nExpected counts\n")
print(round(smoking_chisq$expected, 2))
cat("\nPearson chi-squared test without continuity correction\n")
print(smoking_chisq)
cat("\nOdds ratio for diabetes, smokers compared with non-smokers:",
    round(smoking_odds_ratio, 4), "\n")

cat("\n5. EFFECT SIZE, CONFOUNDING AND POSSIBLE EXPLANATIONS\n")
cat("\n(a) Smokers compared with non-smokers\n")
print(round(smoking_measures, 4))
cat("\n(b) Age by smoking status\n")
print(age_by_smoking_test)
cat("Sex (0 = female, 1 = male) by smoking status (0 = non-smoker, 1 = smoker)\n")
print(sex_by_smoking_table)
cat("\nColumn percentages\n")
print(round(100 * prop.table(sex_by_smoking_table, 2), 2))
print(sex_by_smoking_test)
cat("\n(c) Odds ratio for smoking within age groups\n")
print(round(or_by_age, 3))
cat("\nOdds ratio for smoking within sex\n")
print(round(or_by_sex, 3))
cat("\n(d) Diabetes proportion by age group\n")
print(diabetes_by_age, row.names = FALSE, digits = 4)
cat("\n(e) BMI by diabetes status\n")
print(bmi_by_diabetes_test)


sink()

message("Task 2 outputs saved to: ", normalizePath(output_dir))
