# competing-risk-survival-analysis-colon-cancer
Comparative study of Fine-Gray model and PLANNCR for predicting colon cancer mortality using competing risk survival analysis
# 📊 Competing Risk Survival Analysis for Colon Cancer Prediction

## 🧠 Overview
This project presents a comparative analysis of statistical and machine learning approaches for predicting colon cancer mortality in the presence of competing risks.

Two models are implemented and compared:
- Fine-Gray subdistribution hazard model  
- Neural network-based discrete-time survival model (PLANNCR approach)  

The analysis includes time-dependent evaluation metrics and bootstrapping for internal validation.

---

## 🎯 Objective
To compare the predictive performance of:
- Fine-Gray model  
- Neural Network model  

for estimating cause-specific cumulative incidence in colon cancer patients.

---

## 📂 Dataset
- Source: Secondary dataset (SEER-based published study)  
- Follow-up: Up to 143 months  

### Outcome variable:
- 0 → Censored  
- 1 → Colon cancer death  
- 2 → Death due to other causes  

### Predictors:
- Age  
- Sex  
- Race  
- Marital status  
- Tumor location  
- Grade  
- Tumor size  
- Histologic type  
- T-stage  

---

## ⚙️ Methodology

### 1. Data Preprocessing
- Censoring applied at 143 months  
- Time converted into discrete intervals  
- Data transformed into long format for neural network  
- One-hot encoding for multi-class outcomes  

---

### 2. Neural Network Model (PLANNCR)
- Implemented using `nnet` package  
- Multi-class classification framework  
- Softmax output layer  
- Hidden layer size: 50  

Predicted probabilities:
- No event  
- Event 1 (cancer death)  
- Event 2 (competing risk)  

Cumulative incidence functions (CIF) derived from predicted probabilities.

---

### 3. Fine-Gray Model
- Implemented using `riskRegression::FGR`  
- Separate models for:
  - Cause 1 (cancer death)  
  - Cause 2 (competing event)  

Direct estimation of cumulative incidence functions (CIF).

---

## 📈 Model Evaluation

### Time Points:
- 12 months  
- 36 months  
- 60 months  
- 120 months  

### Metrics:
- Brier Score (prediction error)  
- Time-dependent AUC (discrimination)  
- Confusion Matrix (at 12 months)  

---

## 🔁 Bootstrapping (Internal Validation)

Bootstrapping was used to improve reliability of model evaluation:

- Multiple resampled datasets generated  
- Models refitted on each bootstrap sample  
- Performance metrics averaged across samples  

### Purpose:
- Reduce overfitting  
- Assess model stability  
- Provide robust performance estimates  

---

## 📊 Results Summary

- Both models demonstrated ability to estimate cumulative incidence under competing risks  
- Neural network captured non-linear relationships effectively  
- Fine-Gray model provided stable and interpretable estimates  
- Bootstrapping confirmed consistency of performance metrics  

---

## 🛠️ Tools & Packages

- R  
- survival  
- riskRegression  
- nnet  
- prodlim  
- fastDummies  
- dplyr  

---


