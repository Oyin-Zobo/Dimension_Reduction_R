# Dimension_Reduction_R

- Glass data 
need space
The study of classification of types of glass is motivated by criminological investigations. At the scene
of a crime, the glass left can be used as evidence... if it is correctly identified.
The data set we consider consists of 213 unique glass samples labeled as one of six class categories1:
type description
1 building windows float processed
2 building windows non-float processed
3 vehicle windows float processed
5 containers
6 tableware
7 headlamps
There are nine predictors, including the refractive index and percentages of the following eight elements
found in the glass: Na (Sodium), Mg (Magnesium), Al (Aluminum), Si (Silicon), K (Potassium), Ca
(Calcium), Ba (Barium), and Fe (Iron).
The data is available here: http://archive.ics.uci.edu/ml/datasets/Glass+Identification and is also available in the mlbench package as the dataset Glass.

- Use the data in the file “FB-metrics.csv” for the following problem.
The data is related to Facebook posts published during 2014 on the Facebook page of a renowned
cosmetics brand. This dataset contains 495 rows and 19 features: the first 7 are known prior to post
publication (e.g., current total page likes, post month, etc.) and last 11 features used for evaluating post
impact (e.g., total number of impressions, etc.).
Additionally, please see the paper “Predicting social media performance metrics and evaluation of the
impact on brand building: A data mining approach” by Moro et al. (2016) for more information on the
original study.
In this task, we want to explore the 11 evaluation features using PCA and t-SNE for dimension reduction.
A description of the 11 features are found in Table 2 of the Moro et al. (2016) (Note: I removed “total
interactions” as it is simply the sum of other features.)
(a) Use PCA to analyze the 11 evaluation features. Provide visualizations, interpretations,
and comments as appropriate.
(b)Use t-SNE from the Rtsne package in R to explore 2 or 3-dimensional representations
of the data. Can you find a visualization you find interesting?
For both parts above, consider using different colors to highlight factors associated with the 7 input
features, e.g., paid vs. not-paid, category, type, month, etc. I am interested in what you can find. Do
not be afraid to play with the data, slice and dice as you see fit.
