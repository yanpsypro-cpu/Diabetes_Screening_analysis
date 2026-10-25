# Task 1: Exploratory and descriptive analysis only


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
output_dir <- file.path(dirname(input_file), "task1_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(input_file, check.names = FALSE)

# remove it.
if (ncol(dat) == 10 && names(dat)[1] %in% c("", "X", "...1", "X.1")) {
  dat <- dat[-1]
}

required_variables <- c(
  "Diabetes_012", "HighChol", "CholCheck", "BMI", "Smoker",
  "Stroke", "Sex", "Age", "MentHlth"
)
missing_variables <- setdiff(required_variables, names(dat))
if (length(missing_variables) > 0) {
  stop("Missing expected variables: ", paste(missing_variables, collapse = ", "))
}

continuous_variables <- c("BMI", "Age", "MentHlth")
binary_variables <- c("Diabetes_012", "HighChol", "CholCheck", "Smoker", "Stroke", "Sex")

binary_labels <- c(
  Diabetes_012 = "Diabetes diagnosis",
  HighChol = "High cholesterol",
  CholCheck = "Cholesterol checked",
  Smoker = "Smoker",
  Stroke = "History of stroke",
  Sex = "Male"
)

continuous_summary <- function(x) {
  c(
    n = sum(!is.na(x)),
    missing = sum(is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q1 = unname(quantile(x, 0.25, na.rm = TRUE)),
    q3 = unname(quantile(x, 0.75, na.rm = TRUE)),
    iqr = IQR(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

continuous_table <- data.frame(
  variable = continuous_variables,
  rbind(
    continuous_summary(dat$BMI),
    continuous_summary(dat$Age),
    continuous_summary(dat$MentHlth)
  ),
  row.names = NULL,
  check.names = FALSE
)

categorical_table <- do.call(rbind, lapply(binary_variables, function(variable) {
  counts <- table(factor(dat[[variable]], levels = c(0, 1)), useNA = "no")
  data.frame(
    variable = variable,
    label = unname(binary_labels[variable]),
    category = c("No / 0", "Yes / 1"),
    n = as.integer(counts),
    percent = round(100 * as.integer(counts) / sum(counts), 2),
    row.names = NULL
  )
}))

variable_check <- data.frame(
  variable = names(dat),
  type = ifelse(names(dat) %in% continuous_variables, "Continuous or count", "Binary categorical"),
  n_missing = vapply(dat, function(x) sum(is.na(x)), numeric(1)),
  n_unique = vapply(dat, function(x) length(unique(x)), numeric(1)),
  row.names = NULL
)

# Share of BMI values recorded as integers
bmi_integer_share <- mean(dat$BMI %% 1 == 0)
# Summaries by diabetes status
by_diabetes <- split(dat, dat$Diabetes_012)
continuous_by_diabetes <- do.call(rbind, lapply(names(by_diabetes), function(g) {
  data.frame(diabetes = g, variable = continuous_variables,
             t(sapply(by_diabetes[[g]][continuous_variables], continuous_summary)),
             row.names = NULL)
}))
binary_by_diabetes <- sapply(by_diabetes, function(d)
  100 * colMeans(d[setdiff(binary_variables, "Diabetes_012")] == 1))

write.csv(continuous_by_diabetes, file.path(output_dir, "task1_continuous_by_diabetes.csv"), row.names = FALSE)
write.csv(binary_by_diabetes, file.path(output_dir, "task1_binary_by_diabetes.csv"))

write.csv(continuous_table, file.path(output_dir, "task1_continuous_summary.csv"), row.names = FALSE)
write.csv(categorical_table, file.path(output_dir, "task1_categorical_summary.csv"), row.names = FALSE)
write.csv(variable_check, file.path(output_dir, "task1_data_check.csv"), row.names = FALSE)

png(file.path(output_dir, "figure_1_bmi_age_distributions.png"), width = 1800, height = 850, res = 180)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
hist(dat$BMI, breaks = 40, col = "grey80", border = "white",
     main = "BMI distribution", xlab = expression("BMI (kg/m"^2*")"), ylab = "Number of visits")
abline(v = mean(dat$BMI), col = "firebrick", lwd = 2, lty = 2)
abline(v = median(dat$BMI), col = "navy", lwd = 2, lty = 3)
legend("topright", legend = c("Mean", "Median"), col = c("firebrick", "navy"),
       lty = c(2, 3), lwd = 2, bty = "n")

hist(dat$Age, breaks = 40, col = "grey80", border = "white",
     main = "Age distribution", xlab = "Age (years)", ylab = "Number of visits")
abline(v = mean(dat$Age), col = "firebrick", lwd = 2, lty = 2)
abline(v = median(dat$Age), col = "navy", lwd = 2, lty = 3)
legend("topleft", legend = c("Mean", "Median"), col = c("firebrick", "navy"),
       lty = c(2, 3), lwd = 2, bty = "n")
dev.off()

png(file.path(output_dir, "figure_2_qq_plots_bmi_age.png"), width = 1800, height = 850, res = 180)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
qqnorm(dat$BMI, main = "Normal Q-Q plot: BMI", pch = 16, cex = 0.35)
qqline(dat$BMI, col = "firebrick", lwd = 2)
qqnorm(dat$Age, main = "Normal Q-Q plot: Age", pch = 16, cex = 0.35)
qqline(dat$Age, col = "firebrick", lwd = 2)
dev.off()

positive_percent <- sapply(binary_variables, function(variable) mean(dat[[variable]] == 1, na.rm = TRUE) * 100)
png(file.path(output_dir, "figure_3_binary_variable_prevalence.png"), width = 1700, height = 1000, res = 180)
par(mar = c(9, 5, 3, 1))
bar_positions_binary <- barplot(positive_percent,
        names.arg = unname(binary_labels[binary_variables]),
        las = 2,
        col = ifelse(binary_variables == "Diabetes_012", "#2A9D8F", "grey70"),
        border = "grey30",
        ylim = c(0, max(positive_percent) * 1.18),
        ylab = "Percentage of visits (%)",
        main = "Binary characteristics of screening visits")
text(x = bar_positions_binary, y = positive_percent,
     labels = paste0(round(positive_percent, 1), "%"), pos = 3, cex = 0.8)
dev.off()

png(file.path(output_dir, "figure_4_mental_health_days.png"), width = 1200, height = 850, res = 180)
par(mar = c(4.5, 4.5, 3, 1))
mental_counts <- table(factor(dat$MentHlth, levels = 0:30))
barplot(mental_counts,
        names.arg = 0:30,
        col = "grey70",
        border = "white",
        xlab = "Days of poor mental health in the last 30 days",
        ylab = "Number of visits",
        main = "Self-reported poor mental health days")
dev.off()

# Figure 4b: BMI and age by diabetes status
png(file.path(output_dir, "figure_4b_bmi_age_by_diabetes.png"), width = 1800, height = 800, res = 180)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 1, 1))
diabetes_group <- factor(dat$Diabetes_012, levels = c(0, 1), labels = c("No diabetes", "Diabetes"))
box_col <- c("grey78", "#2A9D8F")
# outline = FALSE hides outlying points, which are too numerous to display
boxplot(dat$BMI ~ diabetes_group, outline = FALSE, col = box_col, border = "grey25",
        xlab = "", ylab = expression("BMI (kg/m"^2*")"))
boxplot(dat$Age ~ diabetes_group, outline = FALSE, col = box_col, border = "grey25",
        xlab = "", ylab = "Age (years)")
dev.off()

sink(file.path(output_dir, "task1_results_log.txt"))
cat("TASK 1: DESCRIPTIVE ANALYSIS\n\n")
cat("Number of visits:", nrow(dat), "\n\n")
cat("DATA CHECK\n")
print(variable_check, row.names = FALSE)
cat("\nCONTINUOUS VARIABLE SUMMARY\n")
print(continuous_table, row.names = FALSE, digits = 3)
cat("\nShare of BMI values recorded as integers:", round(bmi_integer_share, 4), "\n")
cat("\nCATEGORICAL VARIABLE SUMMARY\n")
print(categorical_table, row.names = FALSE)
sink()

message("Task 1 outputs saved to: ", normalizePath(output_dir))
