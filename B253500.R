#install.packages("readxl")
#install.packages("tidyverse")
#install.packages("reshape2")
#install.packages("caret")
#install.packages("randomForest")
#install.packages("glmnet")
#install.packages("car")
# You may need to first install packages before you can run libraries if so please remove # from above and run codes

# Please kindly note for that the names in datasets were changed for avoidance of errors and comfortable coding. 
# Further values in Smoker and Sex columns were changed from numerical values (1,2) to corresponding string values male, female for Sex and smoker, non-smoker for Smoker columns.
# In regression models values in above columns were changed back to represent initial values
# Reminder that in Smoker column 1 represents smoker, 2 represents non-smoker, while in Sex column 1 represents male and 2 represents female

library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(caret)
library(randomForest)
library(glmnet)
library(car)
library(mgcv)


# Reading the necessary libraries for data

biomarkers <- read_excel("/Users/ak/Desktop/Study/biomarkers.xlsx")
covariates <- read_excel("/Users/ak/Desktop/Study/covariates.xlsx")
# Reading excel files


## Data exploration

# Initial observation
biomarkers
covariates

# Checking the data formats for both datasets
str(biomarkers)
str(covariates)

# Getting summary statistics for both dataset
summary(biomarkers)
summary(covariates)


## Data cleaning and featuring

# Checking if there is any missing values in data
anyNA(biomarkers)
anyNA(covariates)

# Finding the row and column of missing values in the covariates dataset
missing_indices <- which(is.na(covariates), arr.ind = TRUE)
print(missing_indices)

# Calculating the mean of the Vas-12months column, excluding NA values
mean_vas_12months <- mean(covariates$`Vas-12months`, na.rm = TRUE)

# Replacing NA values in the Vas-12months column with the calculated mean
covariates$`Vas-12months`[is.na(covariates$`Vas-12months`)] <- mean_vas_12months

# Getting the summary statistics for covariates to make sure that values afer replacement are not significantly different 
summary(covariates[,6])

# Checking for duplicate rows in covariates
duplicate_rows_covariates <- duplicated(covariates)
covariates[duplicate_rows_covariates, ]

# Checking for duplicate rows in biomarkers
duplicate_rows_biomarkers <- duplicated(biomarkers)
biomarkers[duplicate_rows_biomarkers, ]

# Changing the column names to make sure minus, plus or equal sign don't present there for avoidance of errors in coding and for my own comfort
colnames(covariates) <- c("PatientID", "Age", "Sex", "Smoker", "VAS_at_inclusion", "VAS_12months")
colnames(biomarkers) <- c("Biomarker", "IL_8", "VEGF_A", "OPG", "TGF_beta_1", "IL_6", "CXCL9", "CXCL1", "IL_18", "CSF_1", "PatientID")

# Converting Sex column values 1 and 2, for male and female respectively
covariates$Sex <- ifelse(covariates$Sex == 1, "male", "female")

# Converting Smoker column values 1 and 2, for smoker and non_smoker
covariates$Smoker <- ifelse(covariates$Smoker == 1, "yes", "no")

# Printing updated dataframe
print(covariates)

# Creating new column in biomarkers as PatientID and taking value before "-" in Biomarker column as ID for furtue merging and working with data
biomarkers$PatientID <- as.numeric(sub("-.*", "", biomarkers$Biomarker))


## Data visualisation

# Boxplot comparing VAS at inclusion between male and female
ggplot(covariates, aes(x = Sex, y = VAS_at_inclusion, fill = Sex)) +
  geom_boxplot() +
  labs(x = "Sex", y = "VAS at inclusion", title = "Comparison of VAS at inclusion")

# Boxplot comparing VAS at 12 months between male and female
ggplot(covariates, aes(x = Sex, y = VAS_12months, fill = Sex)) +
  geom_boxplot() +
  labs(x = "Sex", y = "VAS at 12 months", title = "Comparison of VAS at 12 months")

# Boxplot comparing VAS at inclusion  between smokers and non-smokers
ggplot(covariates, aes(x = Smoker, y = VAS_at_inclusion, fill = Smoker)) +
  geom_boxplot() +
  labs(x = "Smoker", y = "VAS at Inclusion", title = "Comparison of VAS at inclusion")

# Boxplot comparing VAS at 12 months between smokers and non-smokers
ggplot(covariates, aes(x = Smoker, y = VAS_12months, fill = Smoker)) +
  geom_boxplot() +
  labs(x = "Smoker", y = "VAS at 12 months", title = "Comparison of VAS at 12 months")

# Scatterplot for comparision age and VAS levels
plot(covariates$Age, covariates$VAS_at_inclusion, 
     xlab = "Age", ylab = "VAS at Inclusion", 
     main = "Scatterplot of age vs VAS levels",
     pch = 16, col = "blue")
points(covariates$Age, covariates$VAS_12months, 
       pch = 17, col = "red")
legend("topright", legend = c("VAS at Inclusion", "VAS at 12 months"), 
       col = c("blue", "red"), pch = c(16, 17), cex = 0.8)

# Pearson correlations
correlation_vas_inclusion <- cor(covariates$Age, covariates$VAS_at_inclusion)
correlation_vas_12months <- cor(covariates$Age, covariates$VAS_12months)

# Printing correlations
print(paste("Correlation between age and VAS at inclusion:", correlation_vas_inclusion))
print(paste("Correlation between age and VAS at 12 months:", correlation_vas_12months))



## Data merging
# biomarkets and covariates dataset were merged by PatientID column
combined_data <- merge(biomarkers, covariates, by = "PatientID")
merged_df <- combined_data[order(combined_data$PatientID), ]
str(merged_df)

# Further 12 months dataset was featured as this study focused on biomarkers value at 12 months and this dataset further splitted into smokers and non-smokers to compares values belonging to each of them
twelve_months_data <- filter(merged_df, grepl("12months", Biomarker))
smokers <- twelve_months_data[twelve_months_data$Smoker == 'yes', ]
non_smokers <- twelve_months_data[twelve_months_data$Smoker == 'no', ]
str(smokers)
str(non_smokers)

# Observation of the split data
smokers
non_smokers

# Creating biomarkers list
biomarkers <- colnames(twelve_months_data)[grep("IL_8|VEGF_A|OPG|TGF_beta_1|IL_6|CXCL9|CXCL1|IL_18|CSF_1", colnames(twelve_months_data))]

# Creating a function to plot histograms for each subset with a geom density
plot_histograms_with_density <- function(data, subset_name) {
  melted_data <- melt(data[, c(biomarkers)], variable.name = "Biomarker", value.name = "Value")
  
  p <- ggplot(melted_data, aes(x=Value)) + 
    geom_histogram(aes(y=..density..), bins=30, alpha=0.7, fill="blue") + 
    geom_density(colour="red", fill="red", alpha=0.2) +
    facet_wrap(~Biomarker, scales="free") + 
    labs(title=paste("Biomarkers distribution for", subset_name, ""), x="Value", y="Density") +
    theme_minimal() +
    theme(legend.position="none") 
  
  print(p)
}

# Ploting histograms for smokers with density
plot_histograms_with_density(smokers, "smokers")

# Plotng histograms for non-smokers with density
plot_histograms_with_density(non_smokers, "non-smokers")


## Checking which hypothesis test to apply

# Normality test
# Initialising an empty list to store the results
normality_results <- list()

for (biomarker in biomarkers) {
  sw_smoker <- shapiro.test(smokers[[biomarker]])
  sw_non_smoker <- shapiro.test(non_smokers[[biomarker]])
  normality_results[[biomarker]] <- list("Smoker" = sw_smoker$p.value,
                                         "Non-Smoker" = sw_non_smoker$p.value)
}
# Creating a loop through each biomarker and perform Shapiro-Wilk test
normality_results_df <- as.data.frame(do.call(rbind, normality_results))
rownames(normality_results_df) <- biomarkers

# Converting the results list to a data frame for better viewing
normality_results_df <- t(normality_results_df)
colnames(normality_results_df) <- biomarkers

# Transposing for similar structure 

# Printing the results
normality_results_df

# Checking the rows and columns size for each smokers and non-smokers dataset to identify hypothesis test to be applied
dim(smokers)
dim(non_smokers)

## Variances check

check_variances <- function(smokers, non_smokers, biomarkers) {
  for (biomarker in biomarkers) {
    smokers_values <- smokers[[biomarker]]
    non_smokers_values <- non_smokers[[biomarker]]
    
    # Combine values into a single vector and create a group vector
    combined_values <- c(smokers_values, non_smokers_values)
    group_vector <- factor(c(rep('smokers', length(smokers_values)), rep('non_smokers', length(non_smokers_values))))
    
    # Performing Levene's test for equality of variances
    levene_result <- leveneTest(combined_values, group_vector)
    
    # Printing results of Levene's test
    cat('Levene\'s Test for equality of variances in', biomarker, ':\n')
    print(levene_result)
    cat('\n---------------------------------\n')
  }
}
check_variances(smokers, non_smokers, biomarkers)


## HYPOTHESIS TEST

compare_biomarkers_all <- function(smokers, non_smokers, biomarkers) {
  for (biomarker in biomarkers) {
    smokers_values <- smokers[[biomarker]]
    non_smokers_values <- non_smokers[[biomarker]]
    
    # Performing Student's t-test (assuming equal variances)
    t_test_result <- t.test(smokers_values, non_smokers_values, var.equal = TRUE)
    
    # Performing Mann-Whitney U test
    mw_test_result <- wilcox.test(smokers_values, non_smokers_values)
    
    # Printing results
    cat('Results for', biomarker, ':\n')
    cat('Student\'s T-test (assuming equal variances):\n')
    print(t_test_result)
    
    cat('\nMann-Whitney U test:\n')
    print(mw_test_result)
    cat('\n---------------------------------\n')
  }
}

compare_biomarkers_all(smokers, non_smokers, biomarkers)


## Hypothesis tests with Bonferroni Correction
compare_biomarkers_all_bonferroni <- function(smokers, non_smokers, biomarkers) {
  alpha <- 0.05 # Original significance level
  n <- length(biomarkers) # Number of tests
  alpha_adjusted <- alpha / n # Adjusted significance level for Bonferroni correction
  
  for (biomarker in biomarkers) {
    smokers_values <- smokers[[biomarker]]
    non_smokers_values <- non_smokers[[biomarker]]
    
    # T-test
    t_test_result <- t.test(smokers_values, non_smokers_values)
    
    # Mann-Whitney U test
    mw_test_result <- wilcox.test(smokers_values, non_smokers_values)
    
    cat('Results for', biomarker, ':\n')
    cat('T-test p-value:', t_test_result$p.value, '- Bonferroni adjusted significance level:', alpha_adjusted, '\n')
    if (t_test_result$p.value < alpha_adjusted) {
      cat('T-test result is significant at the Bonferroni-adjusted level\n')
    } else {
      cat('T-test result is not significant at the Bonferroni-adjusted level\n')
    }
    
    cat('Mann-Whitney U test p-value:', mw_test_result$p.value, '- Bonferroni adjusted significance level:', alpha_adjusted, '\n')
    if (mw_test_result$p.value < alpha_adjusted) {
      cat('Mann-Whitney U test result is significant at the Bonferroni-adjusted level\n')
    } else {
      cat('Mann-Whitney U test result is not significant at the Bonferroni-adjusted level\n')
    }
    
    cat('\n---------------------------------\n')
  }
}

compare_biomarkers_all_bonferroni(smokers, non_smokers, biomarkers)

# Type I error calculation

# Significance level for individual tests
alpha <- 0.05

# Number of independent tests
n_tests <- 9

# Calculating the probability of NOT making a Type I error in all tests
probability_no_type_I_error_all_tests <- (1 - alpha)^n_tests

# Calculating the probability of making at least one Type I error
probability_at_least_one_type_I_error <- 1 - probability_no_type_I_error_all_tests

print(probability_at_least_one_type_I_error)


## REGRESSION MODELS
merged_df$Sex <- ifelse(merged_df$Sex == "male", 1, 2)
merged_df$Smoker <- ifelse(merged_df$Smoker == "yes", 1, 2)
# Converting values in Sex and Smoker columns to factor back ti initial values with appropriate labels for regressions models

view(merged_df)
# Its for observation purposes

## MULTIPLE LINEAR REGRESSION
# Featuring new dataframe for regression purpose
regression_data <- filter(merged_df, grepl("0weeks", Biomarker))

set.seed(123) # For reproducibility
trainIndex <- createDataPartition(regression_data$VAS_12months, p = .8, list = FALSE, times = 1)
trainData <- regression_data[trainIndex,]
testData <- regression_data[-trainIndex,]
# Split the data into training and testing sets

model <- lm(VAS_12months ~ IL_8 + VEGF_A + OPG + TGF_beta_1 + IL_6 + CXCL9 + CXCL1 + IL_18 + CSF_1 + Age + Sex + Smoker + VAS_at_inclusion, data = trainData)
# Fit the linear regression model with variables

# View the summary of the model to interpret the fitted parameters
summary(model)

# Predicting the VAS_12months for the test data
predictions <- predict(model, newdata = testData)

# Comparing predictions with the actual VAS_12months values
comparison <- data.frame(Actual = testData$VAS_12months, Predicted = predictions)

# Calculating Mean Squared Error (MSE) for the test set
mse <- mean((comparison$Predicted - comparison$Actual)^2)
print(paste("Mean Squared Error (MSE) on Test Set:", mse))

# Calculating R-squared for the test set
TSS <- sum((testData$VAS_12months - mean(testData$VAS_12months))^2)
RSS <- sum((testData$VAS_12months - predictions)^2)
R_squared <- 1 - (RSS / TSS)
print(paste("R-squared on Test Set:", R_squared))

# Getting the coefficients table
model_summary <- summary(model)  # 'model' is a linear regression model object
coef_table <- model_summary$coefficients  # Extracting the coefficients table
print(coef_table[, "Estimate"], digits = 10)

# Reshaping the data, excluding PatientID and Biomarker columns because they are not required for the model
long_data <- pivot_longer(regression_data, 
                          cols = -c(VAS_12months, PatientID, Biomarker, Sex, Smoker), 
                          names_to = "Variable", 
                          values_to = "Value")

# Ploting with correlation curve (linear regression line)
ggplot(long_data, aes(x = Value, y = VAS_12months)) +
  geom_point(aes(color = Variable), alpha = 0.5) + # Scatterplot points
  geom_smooth(method = "lm", color = "blue", se = FALSE) + # Linear regression line
  facet_wrap(~Variable, scales = "free", ncol = 3) + # Adjusted for better layout
  theme_minimal() +
  labs(title = "Scatterplot of biomarkers vs VAS_12months",
       x = "Variable Value",
       y = "VAS_12months") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Lasso and Ridge regressions
# Preparing the model matrix for training data (excluding the intercept)
x_train <- model.matrix(VAS_12months ~ IL_8 + VEGF_A + OPG + TGF_beta_1 + IL_6 + CXCL9 + CXCL1 + IL_18 + CSF_1 + Age + Sex + Smoker + VAS_at_inclusion - 1, data = trainData)
y_train <- trainData$VAS_12months

# Preparing the model matrix for test data (similarly, excluding the intercept)
x_test <- model.matrix(VAS_12months ~ IL_8 + VEGF_A + OPG + TGF_beta_1 + IL_6 + CXCL9 + CXCL1 + IL_18 + CSF_1 + Age + Sex + Smoker + VAS_at_inclusion - 1, data = testData)
y_test <- testData$VAS_12months

set.seed(123)

# Lasso Regression
cv_lasso <- cv.glmnet(x_train, y_train, alpha = 1, type.measure = "mse")
lasso_model <- glmnet(x_train, y_train, alpha = 1, lambda = cv_lasso$lambda.min)
predictions_lasso <- predict(lasso_model, s = cv_lasso$lambda.min, newx = x_test)
mse_lasso <- mean((predictions_lasso - y_test)^2)
print(paste("Lasso MSE on Test Set:", mse_lasso))

# Ridge Regression
cv_ridge <- cv.glmnet(x_train, y_train, alpha = 0, type.measure = "mse")
ridge_model <- glmnet(x_train, y_train, alpha = 0, lambda = cv_ridge$lambda.min)
predictions_ridge <- predict(ridge_model, s = cv_ridge$lambda.min, newx = x_test)
mse_ridge <- mean((predictions_ridge - y_test)^2)
print(paste("Ridge MSE on Test Set:", mse_ridge))

# Calculating R-squared for Lasso Regression
ss_res_lasso <- sum((predictions_lasso - y_test)^2)
ss_tot_lasso <- sum((y_test - mean(y_test))^2)
r_squared_lasso <- 1 - (ss_res_lasso / ss_tot_lasso)
print(paste("Lasso R-squared on Test Set:", r_squared_lasso))

# Calculating R-squared for Ridge Regression
ss_res_ridge <- sum((predictions_ridge - y_test)^2)
ss_tot_ridge <- sum((y_test - mean(y_test))^2)
r_squared_ridge <- 1 - (ss_res_ridge / ss_tot_ridge)
print(paste("Ridge R-squared on Test Set:", r_squared_ridge))


## Random forest regression
set.seed(123) # For reproducibility
rf_model <- randomForest(VAS_12months ~ IL_8 + VEGF_A + OPG + TGF_beta_1 + IL_6 + CXCL9 + CXCL1 + IL_18 + CSF_1 + Age + Sex + Smoker + VAS_at_inclusion, data = trainData)

# Printing the model summary to see the overall model performance on the training set
print(rf_model)

# Predictnig the VAS_12months for the test data
predictions_rf <- predict(rf_model, newdata = testData)

# Calculating Mean Squared Error (MSE) for the Random Forest predictions
mse_rf <- mean((predictions_rf - testData$VAS_12months)^2)
print(paste("Random Forest MSE on Test Set:", mse_rf))

# Calculating R-squared for the Random Forest predictions
ss_res_rf <- sum((predictions_rf - y_test)^2)
ss_tot_rf <- sum((y_test - mean(y_test))^2)
r_squared_rf <- 1 - (ss_res_rf / ss_tot_rf)
print(paste("Random Forest R-squared on Test Set:", r_squared_rf))

## Polynominal

# Fit the polynomial regression model
poly_model <- lm(VAS_12months ~ poly(IL_8, degree = 2) + poly(VEGF_A, degree = 2) + poly(OPG, degree = 2) + poly(TGF_beta_1, degree = 2) + poly(IL_6, degree = 2) + poly(CXCL9, degree = 2) + poly(CXCL1, degree = 2) + poly(IL_18, degree = 2) + poly(CSF_1, degree = 2) + poly(Age, degree = 2) + Sex + Smoker + VAS_at_inclusion, data = trainData)

# View the summary of the polynomial regression model
summary(poly_model)

# Predicting the VAS_12months for the test data using the polynomial regression model
predictions_poly <- predict(poly_model, newdata = testData)

# Calculating Mean Squared Error (MSE) for the polynomial regression predictions
mse_poly <- mean((predictions_poly - testData$VAS_12months)^2)
print(paste("Polynomial Regression MSE on Test Set:", mse_poly))

# Calculating R-squared for the polynomial regression predictions
ss_res_poly <- sum((predictions_poly - testData$VAS_12months)^2)
ss_tot_poly <- sum((testData$VAS_12months - mean(testData$VAS_12months))^2)
r_squared_poly <- 1 - (ss_res_poly / ss_tot_poly)
print(paste("Polynomial Regression R-squared on Test Set:", r_squared_poly))