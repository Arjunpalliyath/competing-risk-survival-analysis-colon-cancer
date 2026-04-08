# Load required packages
library(nnet)
library(survival)
library(riskRegression)
library(prodlim)
library(fastDummies)
library(dplyr)



library(readxl)
coloncancerinterval <- read_excel("coloncancer.xlsx")


# Step 1: Load your data
# Replace with your path if not already in the environment
# Step 2: Add ID if not present
coloncancerinterval$id <- 1:nrow(coloncancerinterval)

# Step 3: Censor at 143 months if needed
coloncancerinterval$status <- coloncancerinterval$status * (coloncancerinterval$months <= 143)
coloncancerinterval$months <- pmin(coloncancerinterval$months, 143)

# Step 4: Create long format data (train)
train_creator_colon <- function(data) {
  N <- nrow(data)
  intervals <- 15
  breaks <- seq(0, 143, length.out = intervals + 1)
  
  data$interval <- cut(data$months, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  data$survival <- cut(data$months, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  n.times <- data$interval
  data_long <- data[rep(seq_len(N), times = n.times), ]
  
  for (i in unique(data_long$id)) {
    n_length <- length(data_long$interval[data_long$id == i])
    data_long$interval[data_long$id == i] <- 1:n_length
  }
  
  data_long$status_long <- 0
  for (i in 1:nrow(data_long)) {
    if (data_long$status[i] != 0 &&
        data_long$survival[i] == data_long$interval[i])
      data_long$status_long[i] <- data_long$status[i]
  }
  
  events <- fastDummies::dummy_cols(as.factor(data_long$status_long))
  colnames(events) <- gsub(".data", "event", colnames(events))
  data_long <- data.frame(data_long, events[, 2:4])  # event_0, event_1, event_2
  data_long$interval_scaled <- scale(data_long$interval)
  
  return(data_long)
}

# Step 5: Split into train-test (e.g., 70% train)
set.seed(42)
train_ids <- sample(coloncancerinterval$id, size = 0.7 * nrow(coloncancerinterval))
train_data <- coloncancerinterval[coloncancerinterval$id %in% train_ids, ]
test_data <- coloncancerinterval[!coloncancerinterval$id %in% train_ids, ]

# Step 6: Create long format
train_long <- train_creator_colon(train_data)
test_long <- train_creator_colon(test_data)  # can reuse same function

# Step 7: Create x and y matrices
# ➤ Select your predictors + interval_scaled
predictors <- c("age","sex","race",	"Marital","location",	"grade","size","histologic","tstage", "interval_scaled")  # adjust as needed age	sex	race	Marital	location	grade	size	histologic	tstage

train_x <- train_long[, predictors]
train_y <- train_long[, c("event_0", "event_1", "event_2")]
test_x <- test_long[, predictors]
test_y <- test_long[, c("event_0", "event_1", "event_2")]

# Drop factor levels if needed
train_x <- model.matrix(~., train_x)[, -1]
test_x <- model.matrix(~., test_x)[, -1]

# Step 8: Fit the model using `nnet`
set.seed(42)
model_planncr <- nnet(x = train_x, y = train_y,
                      size = 50, decay = 0.001,
                      maxit = 2000, trace = TRUE, softmax = TRUE)

# Step 9: Predict probabilities
pred_probs <- data.frame(predict(model_planncr, test_x, type = "raw"))
colnames(pred_probs) <- c("prob_0", "prob_1", "prob_2")
pred_probs$id <- test_long$id
pred_probs$interval <- test_long$interval
pred_probs$survival <- test_long$survival

# Step 10: Cumulative incidence functions
groups <- split(pred_probs, pred_probs$id)






# Number of intervals you used (plus 1 because we cbind(0, ...))
max_intervals <- 16

risks1 <- lapply(groups, function(x) {
  r <- cumsum(cumprod(1 - (x$prob_1 + x$prob_2)) * x$prob_1)
  length(r) <- max_intervals - 1  # pad with NA if shorter
  r[is.na(r)] <- 0
  return(r)
})

risks2 <- lapply(groups, function(x) {
  r <- cumsum(cumprod(1 - (x$prob_1 + x$prob_2)) * x$prob_2)
  length(r) <- max_intervals - 1
  r[is.na(r)] <- 0
  return(r)
})

# Now this will work without warnings
risks1_mat <- cbind(0, do.call("rbind", risks1))
risks2_mat <- cbind(0, do.call("rbind", risks2))



# Step 1: Convert target months to interval indices
get_interval_index <- function(months, max_month = 143, intervals = 15) {
  breaks <- seq(0, max_month, length.out = intervals + 1)
  sapply(months, function(m) which.min(abs(breaks - m)))
}

target_months <- c(12, 36, 60,120)
target_indices <- get_interval_index(target_months, max_month = 143, intervals = 15)

# Step 2: Subset risks1_mat and risks2_mat to only include those 3 columns
risk1_sub <- risks1_mat[, target_indices]
risk2_sub <- risks2_mat[, target_indices]

# Step 3: Ensure test data is in proper format
df_test <- data.frame(time = test_data$months,
                      status = test_data$status)

# Step 4: Evaluate model performance for cause 1 and 2
score1 <- riskRegression::Score(
  list(planncr = risk1_sub),
  formula = Hist(time, status) ~ 1,
  data = df_test,
  times = target_months,
  cause = 1,
  metrics = c("Brier", "AUC"),
  conf.int = T
)

score2 <- riskRegression::Score(
  list(planncr = risk2_sub),
  formula = Hist(time, status) ~ 1,
  data = df_test,
  times = target_months,
  cause = 2,
  metrics = c("Brier", "AUC"),
  conf.int = T
)

# Step 5: View results
score1$Brier$score
score1$AUC$score

score2$Brier$score
score2$AUC$score









#######################################################fine and gray model





library(riskRegression)
library(data.table)
library(prodlim)

# Make sure these variables are all present in your dataset
predictors <- c("age","sex","race","Marital","location","grade","size","histologic","tstage")

# Fit Fine and Gray model for cause 1 (event of interest)
fg.fit1 <- FGR(Hist(months, status) ~ age+sex+race+Marital+location+grade+size+histologic+tstage,
               data = train_data, cause = 1)

# Fit Fine and Gray model for cause 2 (competing event)
fg.fit2 <- FGR(Hist(months, status) ~ age+sex+race+Marital+location+grade+size+histologic+tstage ,
               data = train_data, cause = 2)


# Convert test data to data.table
test_data_dt <- as.data.table(test_data)

# Target time points (in months)
target_months <- c(12, 36, 60,120)

# Predict risk (cumulative incidence) for cause 1 and 2
risk_fg1 <- predictRisk(fg.fit1, newdata = test_data_dt, times = target_months)
risk_fg2 <- predictRisk(fg.fit2, newdata = test_data_dt, times = target_months)

# Create test dataset in required format
df_test <- data.frame(time = test_data$months,
                      status = test_data$status)

# Evaluate for Cause 1
score_fg1 <- Score(list(FG = risk_fg1),
                   formula = Hist(time, status) ~ 1,
                   data = df_test,
                   times = target_months,
                   cause = 1,
                   metrics = c("Brier", "AUC"),
                   conf.int = T)

# Evaluate for Cause 2
score_fg2 <- Score(list(FG = risk_fg2),
                   formula = Hist(time, status) ~ 1,
                   data = df_test,
                   times = target_months,
                   cause = 2,
                   metrics = c("Brier", "AUC"),
                   conf.int = T)





# Brier Score for cause 1 and 2
score_fg1$Brier$score
score_fg2$Brier$score

# AUC for cause 1 and 2
score_fg1$AUC$score
score_fg2$AUC$score

#######################################









#############################################at12 months


# ---------------------
# Step 1: Set evaluation time
# ---------------------
eval_time <- 12

# ---------------------
# Step 2: Define true class based on time and status
# ---------------------
get_true_event <- function(time, status, eval_time) {
  if (time <= eval_time && status == 1) return(1)
  if (time <= eval_time && status == 2) return(2)
  return(0)
}

df_test$true_class <- mapply(get_true_event, df_test$time, df_test$status, MoreArgs = list(eval_time = eval_time))

# ---------------------
# Step 3: Get predictions at 12 months from PLANNCR
# ---------------------
# Use the correct index for 12 months
interval_index_12 <- get_interval_index(12, max_month = 143, intervals = 15)

# Extract risks at 12 months
pred1 <- risks1_mat[, interval_index_12]
pred2 <- risks2_mat[, interval_index_12]

# Predicted class based on PLANNCR risks at 12 months
pred_class_planncr <- ifelse(pred1 >= 0.5 & pred1 > pred2, 1,
                             ifelse(pred2 >= 0.5 & pred2 > pred1, 2, 0))

# ---------------------
# Step 4: Get predictions at 12 months from Fine-Gray
# ---------------------
# Use the correct column from risk_fg1 and risk_fg2 (for time = 12 months)
# Assuming target_months = c(12, 36, 60), so column 1 = 12 months

pred_class_fg <- ifelse(risk_fg1[, 1] >= 0.5 & risk_fg1[, 1] > risk_fg2[, 1], 1,
                        ifelse(risk_fg2[, 1] >= 0.5 & risk_fg2[, 1] > risk_fg1[, 1], 2, 0))

# ---------------------
# Step 5: Generate confusion matrices
# ---------------------
library(caret)

# For PLANNCR
cm_planncr <- confusionMatrix(factor(pred_class_planncr),
                              factor(df_test$true_class),
                              mode = "everything")
print(cm_planncr)

# For Fine-Gray
cm_fg <- confusionMatrix(factor(pred_class_fg),
                         factor(df_test$true_class),
                         mode = "everything")
print(cm_fg)


table(train_long$status)
table(test_long$status)

summary(train_data$months)
summary(test_data$months)
View(train_data)
