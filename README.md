# DNN-Boost
Somatic mutation identification of tumor-only whole-exome sequencing data using Deep Neural Network and XGBoost.

Note: 
Before running this tutorial, please install the following pacakages: pandas, matplotlib, Tensorflow, Sklearn, imblearn, and xgboost packages.
The DNN-Boost was developed using Python 3.6.9 and the packages with the following versions:
-	Tensorflow: 2.1.0
-	Numpy: 1.18.1
-	Pandas: 1.0.3
-	Xgboost: 1.4.1

This repository is a machine learning model approach to classify variants into two classes, which are somatic mutation variants and germline variants. Deep Neural Network (DNN) was implemented as a classifier and the XGBoost method was used for feature selection. The aim of using these methods is to create a somatic mutation classifier system that is trained with tumor-normal paired data which then the trained model will be employed to classify somatic mutation variants from tumor-only data. 

The tumor-normal paired data of pancreatic cancer whole-exome sequencing was obtained from NCBI Sequence Read Archive (SRA) data portal (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA386668). The tumor-only data was acquired from Therragen Etex whole-exome sequencing of pancreatic cancer.




