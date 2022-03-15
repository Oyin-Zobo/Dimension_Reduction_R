#install.packages("mlbench")
library(mlbench)
#install.packages("Rtsne")
library("Rtsne")
library(tidyverse)
#install.packages("caret")
library(caret)
#install.packages("MASS")
library(MASS)
#install.packages("magrittr")
#library(magrittr)
#install.packages("devtools")
#library(devtools)
#install.packages("ggbiplot")
#library(ggbiplot)
#install_github("vqv/ggbiplot")
library("Seurat")
library("FactoMineR")
#install.packages("factoextra")
library(factoextra)
library(MASS)
library(umap)


data("Glass")
str(Glass)
#Glass = matrix(Glass)

Glass [duplicated(Glass),]  
#to find the duplicate row 
Glass_Unique = Glass [!duplicated(Glass),]
#remove the duplicate row 
Glass_matrix = as.matrix(Glass_Unique[,-10])
#store the new dataframe as a matrix 

corMat = cor(Glass_matrix)
#the correlation matrix for the predictor variables
#in the Glass data set. 
#this shows the collinearity in the data and to see highly colinear data. 
#corMat

#aii) 
eigen(corMat)
#the eigen valuesand values for the glass data set. SHown below 
#aiii)
pca = prcomp(Glass_matrix, scale= T)
#the pca would project the dimensions of the dataset into smaller dimension.
#this preserves the majority of the inofrmation 
#from the data and also eases visualization of data relationship

#aiv)
#The eigenvalues of the correlation matrix 
#are the squares of the standard deviation i.e variances
#of the principal components themselves are same as eigenvectors
#of the correlation matrix. Even though the signs 
#maybe opposite as is the case here. 

#av)
#pca
pca_columns = pca$rotation

PC1 = pca_columns[,1]
PC2 = pca_columns[,2]
PC1%*%PC2
# two vectors are orthogonal when the inner product is 0. The inner products of 
#PC1 and PC2 was approximately zero whcih means they are orthogonal 

#1b 
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50))
#The diagram shows the number of dimensions that the data set was reduced to
#and the variance explained by each of the dimensions. 
#The first dimension explains the most varinace (27.9%). The first 5 dimensions 
#explains 89.3% variance of the original data. 

biplot(pca)
#this does not show a good represenatatio of the data due to the large number of
#individual parameters. Color would also make the plot look better.

fviz_pca_var(pca, col.var = "black", repel = TRUE)
#variable correlation plot.The distance of the varibale from the center of the circle
#represent the quality of the variables on the factor map. 
#Variables that are far away are well represented i.e with longer arrow length. 
#i.e the the longer the arrow lines, the more influence 
#the variable has on the principal component. 

fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
#The cos2 values are used to estimate the quality of the representation
#The closer a variable is to the circle of correlations, the better its 
#representation on the factor map 
#(and the more important it is to interpret these components)
#Variables that are closed to the center of the plot are 
#less important for the first components.
#variables with high correlations are colored red, blue for the low effects 
#and orange for average

#pca
#the shows all the pCs sorted accoring to the importance 
#and weight of each of the variables. 


glass.pca <- PCA(Glass_matrix, scale.unit = TRUE)

glass.desc <- dimdesc(glass.pca)
# Description of dimension 1
glass.desc$Dim.1
#the first PC represents 27.89% of the variance of the data. 
#the table shows the % of each of the variable 
#explained by PC1. tHE RI variable has the largest contribution to PC1
#from the biplot plot, the variable with high 
#contribution increase or decrease across the horizontal section. 


glass.desc$Dim.2
#the second PC represents 22.87% of the variance of the data. 
#the table shows the % of each of the variable explained by PC1. 
#tHE Mg variable has the largest contribution to PC1
#from the bilpot plot,  the variable with high 
#contribution increase or decrease across the vertical section. 



#Glass_Unique$Type
type = as.factor(Glass_Unique$Type)

fviz_pca_ind(glass.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = Glass_Unique$Type, # color by groups
             palette = c("red", "blue", "green", "yellow", "purple", "orange"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

#this plot shows the types of glass and the groups of each of the glass type. 
#The red represents the type 1, blue color represents Type 2, green color 
#represents Type 3, Yellow color represents Type 5, Purple color represents 6,
#Orange color represents Type 7 glass. 

#The glass type and cluster overlap on the PC1 and PC2 which shows that the 
#glass type share similar properties. 


rbind(
  SD = sqrt(pca$sd^2),
  per_var = (pca$sd^2/sum(pca$sd^2))*100,
  per_var_cum = (cumsum(pca$sd^2)/sum(pca$sd^2))*100)
#the dimension can be reduced to 5 dimensions
#this tables shows the standard deviation represented by 
#each of the PC and the standard deviation 
#decreases as the PC number increases.
#This is also represented in the percentage variance explained by each of the 
#principal comp. 9 PCs explain 100 variation is the data 
#but 5 PCs represents 89.3% which is very good representation. 

#biii
#yes, The dimension can be reduced to 5 PCs and 
#the still preserve 89.3% of the data. 

#1c
preproc.param = Glass %>% preProcess(method = c("center", "scale"))
#preprocess does mean centering and scaling and 
#it is not affected by factor level. It ignores them 

#transform the data using the estimated parameters
transformed  = preproc.param %>% predict(Glass)

#transformed
#fit the model 

lda.model  = lda(Type ~., data = transformed)
#lda.model
#The purpose of the linear discriminant analysis is to find 
#combination of the variables that give 
#best possible separation between groups (glass Type) in our data set.
#group probability means the initial data proportion 



Coefficients_of_linear_discriminants = lda.model$scaling
#the coefficient that shows the proportion of each of the
#variables represented by the LDA, to get the . 
Coefficients_of_linear_discriminants  

one_lda = Coefficients_of_linear_discriminants[,1]
one_lda
#the first LDA1, it contains 0.8145 of the data. 
#The linear discriminant function from the result in above is
#0.94*RI+1.944Na+1.068*Mg+1.667*Al+1.8989*Si+1.024*K+
#1.432*Ca+1.15*Ba-0.0498*Fe

 
#(the variable returned by the lda() function) is the percentage 
#separation achieved by each discriminant function. 
#LDA1 independently achieve a separation of 0.8145. 

two_lda = Coefficients_of_linear_discriminants[,2]
two_lda
#LDA2 independently achieve a separation of 0.1169. 

lda.model_values = predict(lda.model)
#linear regression model for the hitogram 

par(mar=c(1, 1, 1, 1))
ldahist(lda.model_values$x[,1], g = type)
#From the graph above, we have histogram from LD1, 
#the type 1,2,3 have some overlap between them 
#and we can see that the separation between 1, 2, 3 and the other three Species 
#is quite small  with some overlap. On the contrary, there is a certain amount  
#of overlapping between type 5,6,7. We already said that the tight
#percentage of separation archived by LD1 is 81.45%, that is there is some 
# clear separation from the histogram above. 

#Now, we can try to do the same for LD2.

ldahist(lda.model_values$x[,2], g = type)

#the distance between the separation between the type is very small and we can 
#see the overlapping between all the types of glasses with no clear separation. 
#this is expected due to the smaller separation of 11.69% achieved by LDA. 
#The LDA does not provide a good separtion between the glass types. 

#tsne_out = Rtsne(Glass[,-10], perplexity=5, theta = 0.0, 
                 #max_iter = 2000, num_threads = 6)

#df = data.frame(x=tsne_out$Y[,1], y= tsne_out$Y[,2], type = Glass[,10])
#ggplot(data=df, aes(x=x, y=y, group=type, color=type))+geom_point()



#2
FB_metric = read.csv(file = 'FB-metrics.csv', 
                   header = TRUE, sep = ",", 
                   quote = "\"" , dec = ".", 
                   fill =TRUE, comment.char = "")
#FB_metric
#summary(FB_metric)

FB_metric_new = FB_metric[,8:18]
#extract the 11 columns needed for the analysis 

#str(FB_metric_new)
#the summary of the each columns 
FB_pca = prcomp(FB_metric_new, scale= T)
#model the pca on the 11 variable dataframe to reduce the dimension. 
summary(FB_pca)
#gives the summary. the first three PCs give a cumulative of 84%

rbind(
  SD = sqrt(FB_pca$sd^2),
  per_var = (FB_pca$sd^2/sum(FB_pca$sd^2))*100,
  per_var_cum = (cumsum(FB_pca$sd^2)/sum(FB_pca$sd^2))*100)

heatmap(cor(FB_metric_new))
#shows the closely correlated variables in the data. Total impression and  
#reach on people who liked the page and engagged the post are highly correlated. 
corMat_FB = cor(FB_metric_new)
#corMat_FB

fviz_eig(FB_pca, addlabels = TRUE, ylim = c(0, 70))
#the first three dimensions represents 84% of the data. 
#4 dimensions represents 89.7%. The data can be well represented with 
#4 PCs 

fviz_pca_var(FB_pca, col.var = "black", repel = TRUE)
#the bilpot shows hows the variables are represented on the PCA biplot. 
#The longer length arrow have a higher weight/effect on the PC 1 and PC2. 
#the PC1 has more effect on the the representation. This is shown on the plot 
#as more of the variables increase towards the horizontal direction. 

fviz_cos2(FB_pca, choice = "var", axes = 1:2)
#visualize the cos2 to represent the quality of 
#the variables on factor map. square cosine, squared coordinates. 
# The number of engaged users has the highest effect and impressions 
#on people has the lowest effect on the PC1 and PC2. 

library("corrplot")
var = get_pca_var(FB_pca)
corrplot(var$cos2, is.corr=FALSE)


fviz_pca_var(FB_pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
#The cos2 values are used to estimate the quality of the representation
#The closer a variable is to the circle of correlations, the better its 
#representation on the factor map 
#(and the more important it is to interpret these components)
#Variables that are closed to the center of the plot are 
#less important for the first components.
#variables with high correlations are colored red, blue for the low effects 
#and orange for average. The varibales with the highest imapact is the 
#post consumers and the engaged users. 

#scatter plot matrix for the PCS

# Contributions of individual variables to PC1
fviz_contrib(FB_pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(FB_pca, choice = "var", axes = 2, top = 10)
# Contributions of variables to PC3
fviz_contrib(FB_pca, choice = "var", axes = 3, top = 10)
# Contributions of variables to PC4
fviz_contrib(FB_pca, choice = "var", axes = 4, top = 10)


fviz_pca_biplot(FB_pca, 
                col.ind = FB_metric$Type, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Type") 
#the data is classified and colored with regards to the Type of media 
#link - blue, Photo - yellow, Status - grey, Video -red. The clusters can be 
#seen to overlap on the two dimensions. 
#This shows they either have close properties 
#or the PCA does not do a good job with separation. This biplot is for the PC1 
#and PC2 

fviz_pca_biplot(FB_pca, 
                col.ind = as.factor(FB_metric$Paid), palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Paid")

#the data is classifed accoring to whether it was paid or unpaid. This shows 
#the paid ad more largr reach and an increase in Likes, comments, 
#total impressions, 
#and other variables i.e arrows on the biplot. This biplot is for the PC1 and
#PC2 
fviz_pca_biplot(FB_pca, 
                col.ind = FB_metric$Type, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Type", axes = c(2,3)) 
#this plot is for PC 2 and PC 3. the separation in this section is less clear as 
#this two dimensions have lower representation of the data. 

fviz_pca_biplot(FB_pca, 
                col.ind = FB_metric$Type, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Type", axes = c(3,4)) 
#this plot is for PC 3 and PC 4. the separation in this section is less clear as 
#this two dimensions have lower representation of the data. 



FB_PCA = PCA(FB_metric_new)
FB.desc <- dimdesc(FB_PCA)
# Description of dimension 1
FB.desc$Dim.1
#dimension description is used to identify the 
#most significantly associated variables with a given principal component. 
#the most associated in this case of Dim1  is Engaged user with weight of 0.876
FB.desc$Dim.2
#the most associated in this case of Dim1  is total 
#impressions with a weight of 0.435 with a positive correlation. 
FB.desc$Dim.3

#pca for individuals 
ind  = get_pca_ind(FB_pca)

# Coordinates of individuals
#head(ind$coord)
# Quality of individuals
#head(ind$cos2)
# Contributions of individuals
#head(ind$contrib)
#the individuals cannot be plotted as they are very many and would not have a good 
#interpretation. 

fviz_pca_ind(FB_pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = FB_metric$Type, # color by groups
             palette = c("red", "blue", "green", "yellow"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

#these shows the various type of post and the colors. The types have a large 
#representation because of the high representation of PC 1 and PC2 which are 53.7%
#and 15.6% respectively. A clear distinction cannot be seen as they are overlapping
#but PCA is not a very good clustering method but it is for unsupervised learning.
# a good dimension reduction is observed, 
fviz_pca_ind(FB_pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = FB_metric$Type, # color by groups
             palette = c("red", "blue", "green", "yellow"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             axes = c(2,3)
)
#project on dimensions 2 and 3. The groups have smaller circles here which show
#smaller imapact of the post type on these PCs this is because they have smaller 
#representations of 15.6% and 14.7%. They have less percentage of the total data. 

fviz_pca_ind(FB_pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = FB_metric$Type, # color by groups
             palette = c("red", "blue", "green", "yellow"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             axes = c(3,4)
)
#project on dimensions 3 and 4. They have lesser percentage of the total columns 
#as they only represent 14.7% and 5.7% repectively. 

fviz_pca_biplot(FB_pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

#THE BILPOT PLOT is unable to read due to the large number of varibales 
#an individual that is on the same side of a 
#given variable has a high value for this variable;
#an individual that is on the opposite side of a 
#given variable has a low value for this variable.


#b
tsne_out = Rtsne(FB_metric_new[,1:11], perplexity=50, 
                 theta = 0.0, max_iter = 5000, num_threads = 6)

#tsne_out
cols = rainbow(4)
plot(tsne_out$Y, t='n')
text(tsne_out$Y, labels=(FB_metric[,3]), col=cols[FB_metric[,3]])
#the shows the category of the data colored according to the number of the 
#category. The tsne do not show good clustering for the catergory type. 


str(FB_metric)
#head(tsne_out)
df = data.frame(x=tsne_out$Y[,1], y= tsne_out$Y[,2], type = FB_metric[,2])
df_1 = data.frame(x=tsne_out$Y[,1], y= tsne_out$Y[,2], 
                  type = as.factor(FB_metric[,3]))
#the category of the data 
df_2 = data.frame(x=tsne_out$Y[,1], y= tsne_out$Y[,2], 
                  type = as.factor(FB_metric[,7]))
#is the ad paid or unpaid? 
df_3 = data.frame(x=tsne_out$Y[,1], y= tsne_out$Y[,2], 
                  type = as.factor(FB_metric[,4]))
#the month of the post 
ggplot(data=df, aes(x=x, y=y, group=type, color=type))+geom_point()
#the type of the post was used to cluster the data. 
#The T-SNE didnt do a good job
#with the clusters 
ggplot(data=df_1, aes(x=x, y=y, group=type, color=type))+geom_point()
##the category of the post was used to cluster the data. 
#The T-SNE didnt do a good job
#with the clusters 
ggplot(data=df_2, aes(x=x, y=y, group=type, color=type))+geom_point()
#the type of the paid or unpaid option was used to cluster the data. 
#The T-SNE didnt do a good job with the clusters 
ggplot(data=df_3, aes(x=x, y=y, group=type, color=type))+geom_point()
#the type of the month was used to cluster the data. 
#The T-SNE didnt do a good job with the clusters 

#install.packages("umap")

#reducing clustering for the purpose of clustering is 
#where the difference between tsne and umap starts to show. 
FB_umap = umap(FB_metric_new, n_components = 2, 
               n_neighbors = 50, min_dist = 0.01)

head(FB_umap$layout)

fb.umap = umap(FB_metric[,8:18])
fb_umap_df = data.frame(x=fb.umap$layout[,1], y = fb.umap$layout[,2], 
                        Type = FB_metric$Type, 
                        Category = as.factor(FB_metric$Category), 
                        Post.Month = as.factor(FB_metric$Post.Month),
                        Paid = as.factor(FB_metric$Paid))

ggplot(data = fb_umap_df, aes(x=x, y = y, color = Type))+geom_point()
#the UMAP tries to cluster using the type of post. It however does not show good 
#clusters for the data. The photo post is shown to dominate the clusters as the 
#biggest group. 
ggplot(data = fb_umap_df, aes(x=x, y = y, color = Category))+geom_point()
#UMAP tries to cluster the data with the post category. There is no clear cluster 
#separation and no district clusters are seen. 
ggplot(data = fb_umap_df, aes(x=x, y = y, color = Post.Month))+geom_point()
#UMAP tries to cluster the data with the post month There is no clear cluster 
#separation and no district clusters are seen. 
ggplot(data = fb_umap_df, aes(x=x, y = y, color = Paid))+geom_point()
#UMAP tries to cluster the data with whether they are paid or not.  There is no clear cluster 
#separation and no district clusters are seen. 

#the UMAP parameter optimization as done to achieve better results. There was 
#however not very obvious results. 
#the UMAP shows a better representation of clusters than Tsne. Even though the 
#clusters were not very clear, the structure looked better with the UMAP. 
#Parameter optimization was done for TSNE and UMAP, but the results still dont 
#show very good trends and results. 