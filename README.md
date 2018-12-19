# Sematics of Statistical models

This analysis is described in a paper submitted to *Geographical Analysis* in October 2018. The code is in the file `Semantics_code_data_git.R` and this loads data from this repository.

Please contact Lex Comber [a.comber@leeds.ac.uk](a.comber@leeds.ac.uk) if you have any questions.

## Paper title: The forgotten semantics of regression modelling in Geography
Alexis Comber<sup>1</sup>, Paul Harris<sup>2</sup>, Yihe Lü <sup>3</sup>, Lianhai Wu <sup>2</sup> and Pete Atkinson<sup>4</sup> 

<sup>1</sup>School of Geography, University of Leeds, Leeds, UK LS2 9JT\
<sup>2</sup>Rothamsted Research, North Wyke, Okehampton, Devon, UK EX20 2SB\
<sup>3</sup>Chinese Academy of Sciences, Beijing, 100085, China\
<sup>4</sup>Faculty of Science and Technology, Lancaster University, UK

## Abstract
This paper is concerned with the semantics associated with the statistical analysis of spatial data. It takes the simplest case of the prediction of *y* as a function of *x*, in which predicted y is always an approximation of *y* and is always, and can only ever be, a function of *x*, and illustrates a number of core issues using ‘synthetic’ remote sensing and ‘real’ soils case studies. Specifically, the outputs of regression models and therefore the meaning of predicted *y*, are shown to vary due to 1) choices about data: specification of *x* (which covariates to include), the support of *x* (measurement scales, granularity), the measurement of *x* and the error of *x*, and 2) choices about the model including its functional form and the method of model identification. Some of these issues are more widely recognised than others. The case studies illustrate the effects of data and model choices and their impacts on model outputs. The study provides definition to the multiple ways in which prediction from regression may be affected and shows how regression prediction and inference are affected by data and model choices. The paper invites researchers to pause and consider the implications of predicted y being nothing more than a scaled version of a single covariate, inheriting the same spatial correlation, and argues that it is naïve to ignore this. 

## Acknowledgements
This research was supported by the China-UK bilateral collaborative research on critical zone science (the Natural Environment Research Council Newton Fund NE/N007433/1, the National Natural Science Foundation of China NO. 41571130083) and the National Key Research and Development Program of China (No. 2016YFC0501601). All of the data preparation, analyses and mappings were undertaken in R 3.5.1, the open source software. 
