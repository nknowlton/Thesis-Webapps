---
output: 
  html_document: 
    theme: cerulean
    toc: yes
---
Welcome to the Multi Omics Viewer Tutorial
=======================
**This is still a Work in Progress**

The Multi Omics Viewer Tutorial provides an interactive way to examine the data relationships within a TCGA dataset (currently Metastatic Melanoma and Primary Breast Invasive Carcinoma). 
Relationships are based upon a linear regression of RNAseq gene expression vs Methylation, Copy Number, and Genetic Variants.
The following tutorial will walk through an example of each tab showing how the options change the presentation and offering advice on how to interpret what is being displayed.

## Tabbed Browsing

To navigate around the use of `tabs' is essential.  As you are reading this page now, congratulations, you have already use the tab at least once! There are 6 tabs to choose from; the currently selected tab will be slightly darker and have a thin green line beside the selection, as the Tutorial Tab does in the image below.

![](images/tabselect.png)


## Data Set

The Data Set tab is the default tab when the application starts and loads sets the data import options for the rest of the tabs. Starting in the upper left hand corner, the "Select from" drop down allows the selection of two TCGA datasets, namely Metastatic Skin Cutaneous Melanoma and Breast Invasive Carcinoma. Moving to the lower right hand corner, the "Select Copy Number Flavour" from the drop down allows the selection of Gene Segment Variants and Gistic. Gene Segment Variants (the default) is a continuous estimate of the copy number. This is contrasted to the Gistic estimate which bins copy number into -2,-1,0,1,2. These numbers correlate with loss of 2 copies of a gene (-2) to gaining 2 copies (+2). The upper right hand corner has the "Include Variant Information?" drop down which allows selection of No, Yes- Present/Absent, Yes- Factors. Since the tool is based upon a linear model any missing values in any variable would remove a gene completely. When the option of "No" is selected the overall number of genes possible to interrogate is at its maximum. When the option of "Yes- Present/Absent" is chosen, the Genetic Variants are modelled as binary with no regard to type or how many are observed in a given tumour. When the option "Yes- Factors" there are 11 variant types modelled each modelled as binary. In both "Yes" options the minimum number of mutations for a gene (Present/Absent) or mutation type (Factors) is 3 to allow estimation of the beta value from the linear model. The bottom right corner has the option "Select the association model" and allows the selection of Relative Importance and Correlation. Given that the Methylation, Copy Numbers and Variants are assumed to be associated, some more strongly than others, a relative importance measure was calculated to indicate the amount of variability explained by each of the variables using the RELIMPO package. Correlation is the Pearson correlation of every variable vs the RNAseq expression in a univariate association. Finally, the blue outlined box labelled "Database Information" shows a message from the R backend as to whether the database was loaded correctly and how many samples and genes are available for a given combination. Please wait for this value to display before advancing to the next tab.

![](images/dataset-all.jpg)

## Select Gene

After selecting a dataset, whether to include variant information, Copy Number Flavour and Association model you will be presented with a table similar to the one below:

<div style="width:300px; height=200px">
![](images/genetab.jpg)
</div>

The top of the tab has an orange help tab that can be expanded to remind you what is computed in this tab. The Select Gene of Interest table provides a list of all genes available for interrogation along with the respective HUGO gene symbol. As the data integrates Methylation data from the 450K Illumina Methylation Arrays there are many cg sites can map to a single gene. The current table presents the one with the most RNAseq variance explained, identified under the Meth Site Column. Other sites related to the selected gene can be examined on the next tab. The Total column is present only if the Relative Importance Association was selected and represents the overall amount of variance in RNAseq is represented by all terms in the current model (selected on the previous tab). The Copy Number, Methylation, and Variant (if selected) columns are sortable from low to high and high to low by click on their title. Once sorted a single gene symbol can be selected by clicking on the symbol anywhere in the row. Once clicked the name of the gene selected will appear in under the blue Data Status box in the upper right hand column. A gene symbol can be searched directly using the Search box. Please wait for the gene you have selected to appear in the blue box before continuing to the next tab.

## Single Gene Plots

The Single Gene Plots is the graphical heart of the application.Here you are presented with a table of all matching Methylation sites for the selected gene, a graph of Methylation level vs Copy Number, a graph of Methylation vs Expression, and a graph of Expression vs Copy Number.

![](images/sgplot-part1.jpg)

![](images/sgplot-part2.jpg)

Starting with the Methylation Site table, the symbol column is the Gene Name associated with a given Methylation site. In the example above, the cg site is associated with both TP53 and the p53 antisense transcript WRAP53. The location column is a concatenated list of two predicted locations for the cg site separated with an '_'. In the example, the cg06317056 site is associated with a 5'UTR and a CG Island. The other columns function as in the Select Gene tab allowing the sorting here of characteristics within a HUGO symbol.

Moving to the Meth vs Copy Number coloured by Expression plot we can see that there is a strong negative relationship between Methylation site cg06317056 and Copy Number. The coloured points are the expression values with low expression represented as green and high expression represented as red. In this example, the low expression values are associated with high Methylation values and vice versa. The bars at the side of the plot are the marginal histograms for the relevant terms. 

The two plots at the bottom of the image are more two way relationships of the data. The bottom left is Methylation vs RNAseq expression. A spline curve with 95% confidence intervals is plotted. In this example the relationship is almost non existent. This agrees with with an explained variance for the Methylation site selected of 4.07%. Moving on to the Expression vs Copy Number plot we can see a strong positive relationship between RNAseq Expression and Copy Number. 

## Variant Plots

The relationship between genetic variants and expression can vary based upon the type of mutation. In an effort to allow these relationships to be more easily visualized the following two plots were created.

![](images/varplots.jpg)

If 'Yes-Factors' was selected on the Include Variant Information page then the plots will be available. The tan dots are every tumour sample and how they relate to the RNAseq Expression and Copy Number.  When a variant is selected using the check box on top of the plot the indicated variants are highlighted in purple. In the example Missense Mutation is selected on the left plot.  
Not visible in the example are the 11 different types of mutations that are recorded in the datasets. These include: Missense, Silent, Nonsense, Splice Site, RNA, Frame Shift Insertion, Frame Shift Deletion, In Frame Insertion,In Frame Deletion, Nonstop, and Translation Start Site. The option for each is only presented if a given mutation type is present. The one line table at the bottom shows the amount of variation explained by each mutation type. The sum of these mutations is the 'Variant' column.


## Computational Model Limitations

This webapp relies on a linear model of RNAseq vs Methylation + Copy Number + Variant. To ensure that the model is always estimable there are some restrictions placed on the terms. The variant term must have at least 3 unique mutations or the variant term is not estimated. This rule applies for both the Present/Absent and the Factors selections. If there are less than 3 of a given mutation type they are still available for plotting on the Variant Plot tab, but their variance explained is NA.

## Contact us:

Just to say Hi or provide feedback please contact us.

### Via LinkedIN

[Nick Knowlton](https://nz.linkedin.com/in/nsknowlton)

[Cris Print](https://www.linkedin.com/in/cristin-print-47b72b17)

### Via Email
[Nick Knowlton](mailto:n.knowlton@auckland.ac.nz ?subject=Multi Omics Viewer: Shiny)

[Cris Print](mailto:c.print@auckland.ac.nz ?subject=Multi Omics Viewer: Shiny)


