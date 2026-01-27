# CLIC

![Project Status](http://www.repostatus.org/badges/latest/active.svg)

### Overview

The **CLIC** package provides the tool for feature selection in the context of integrating unpaired scRNA-seq and scATAC-seq data. 
Unlike traditional feature selection methods that measure solely the gene's expression variability, CLIC prioritizes genes that show 
high empirical correlation between expression and nearby accessibility as "high confidence links". 

**CLIC scores**, used as a feature's confidence score is computed using diverse single-cell multiome data from ENCODE.
The script used for computing CLIC scores is available [here](https://github.com/oldvalley49/ENCODE_score)

For more details, please refer to our paper describing the method.

(Preprint link pending)

### Download

Please use GitHub repo to download the most updated package.

```R
devtools::install_github("oldvalley49/CLIC")
```

### User Manual

Please view the vigentte for instructions on how to use the package:
[Introduction to CLIC package](https://htmlpreview.github.io/?https://raw.githubusercontent.com/oldvalley49/CLIC/main/doc/Introduction.html)

### Citation
(Pending)

### License
This software is licensed under the [MIT License](./LICENSE.txt).

### Contact
Should you encounter any bugs or have any suggestions, please feel free to contact Tomoya Furutani at tfuruta1@jh.edu or open an issue on GitHub