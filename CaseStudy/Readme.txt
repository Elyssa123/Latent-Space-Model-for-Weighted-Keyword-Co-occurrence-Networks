keyword_adj_20core.rda: the 20-core of the keyword co-occurrence network

-The data is stored in the form of a tensor, with the size of 1279*1279*10. Here, 10 represents the number of time slices, and 1279 stands for the number of distinct keywords throughout the entire time period. 

-Each element on the tensor represents the co-occurrence count of two keywords in a specific time period. If a keyword doesn't exist in a particular time period, the values in the corresponding row and column for that time slice are empty.

CaseStudy-code.Rmd：The R code for the case study. 

CaseStudy-code.html: The returned results for CaseStudy-code.Rmd. 

beta-alpha.png, line_plot_latent_space.png, Z.png: The returned figures for CaseStudy-code.Rmd.