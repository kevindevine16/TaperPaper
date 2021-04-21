# TaperPaper

Code used to produce figures for the formation of air gaps in continuous steel casting. (See https://doi.org/10.1093/imamat/hxab003)

Code uses finite difference methods to solve PDE for a thermal model. 

The location of the air-gap initial formation is determined using fzero on a function which also uses finite difference methods.

The model is then compared to data from Kelly <i> et al. </i> https://doi.org/10.1007/BF02645486 which is imported in the code.
