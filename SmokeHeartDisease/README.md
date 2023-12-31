# Final Project: An Investigation into the Relationship between Smoking and Heart Disease

For this study, I chose the Stroke Prediction Dataset from Kaggle. The raw data are publicly accessible via <https://www.kaggle.com/datasets/fedesoriano/stroke-prediction-dataset>. <br>
This dataset consists of 5,110 patient records (276 presenting with heart disease and 3,834 not presenting with heart disease) and contains 12 features. In my study, the outcome of interest is the presence of heart disease, and the exposure is the smoking status.

Before analysis, I conducted a sanity check on the original data to select the analytical sample. I found that the dataset had some outliers and missing data. For example, there was one patient of gender ‘other’, and some patients had BMI greater than 50, which were likely to be mistakenly recorded. Moreover, some children was recorded as ‘smokes’ or ‘formerly smoked’. For the validity of my analytical data, I excludes the patient of gender ‘other’, and the patients with BMI larger than 50 (this is an empirically chosen number). To dichotomize the exposure, I excluded the patients with smoking status ‘formerly smoked’ or ‘unknown’. Then, ‘smokes’ and ‘never smoked’ were encoded as 1 and 0, respectively. Because the exposure was longitudinal, and heart diseases observed in children are usually congenital (not caused by smoking), I also excluded the patients of age smaller than 18, so that my sample consisted of only adults. After sample selection, the final analytical data included 2,398 patients, with 135 presenting of heart disease. <br>

For mediation analysis, I generated a new column ```cholsterol``` which records the blood cholesterol level (mmol/L) of patients. Please note that these data are FICTITIOUS and are only for educational purpose. If you are conducting a scientific research, please do not use this datatset.

The processed dataset is uploaded to ```smoke_heart.CSV```.
