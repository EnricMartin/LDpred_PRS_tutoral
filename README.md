## A Tutorial for LDpred 

### 1.	Intro: Why you would need this tutorial

**LDpred is one commonly used method for computing Polygenic Risk Scores.** It accounts for the LD information from an external reference panel, thereby improving the prediction accuracy of the PRS compared to the traditional Clumping + Threshold method[^1]. LDpred has been cited and applied in large amounts of studies.

**However, when you search for the potential tutorials at Bing, the results are lacking.** Also, the first search result is guidance for performing LDpred-2 using R[^2], but according to my practical experiences, performing LDpred-2 using R presents a huge challenge to your computer's memory and disk space. During the process of calculating the LD matrix, a huge .sbk file would be produced and kill the program.

**What's more important is that the author of LDpred believes much in our ability to learn its usage, thereby the description of LDpred on the GitHub website is also very simple.**

The above content constitutes the major reason for producing this tutorial.

### 2.	What you should prepare

#### 2.1	Data Preparation

**A.	The post-QCed GWAS summary statistic**

- How do we perform quality control on GWAS summary data? I recommend this GitHub page [^3].
- What is the required format for the GWAS summary data? The author has answered this question on the Q & A page on GitHub[^4]. 

**B.	The reference panel for calculating LD [PLINK binary type data]**

- The reference panel should be from the same genetic ancestry as the individuals in the GWAS.

- **<u>I have uploaded one 1000G Eastern Asia reference panel after basic quality control</u>**

  ```shell
  ./plink --bfile 1000geas --maf 0.01 --mind 0.01 --geno 0.01 --make-bed --out 1000geas.QC
  ```

**C.	Individual genotype data for computing the scores [PLINK binary type data]**

- Standard quality control should be conducted [^5].

**D.	(Optional) HapMap3/HapMap3+ variants**

- It is recommended to use a set of HapMap3/HapMap3+ variants for LDpred, because they provide a good coverage of the genome and are generally well-imputed and available in most studies.

- Download HapMap3+ variants using R:

  ```R
  info <- readRDS(runonce::download_file(
    "https://figshare.com/ndownloader/files/37802721",
    dir = "tmp-data", fname = "map_hm3_plus.rds"))
  ```

- **<u>I have uploaded the HapMap3+ variants (named HapMap.rds) to this repository.</u>**

#### 2.2	Environment requirements

To run the Ldpred successfully, you should prepare the following Python environment:

* Python 3.8
* h5py 
* scipy 1.7.0 
* plinkio
* hickle
* ldpred 1.0.10

**<u>Note!  Not the newer the better.</u>**

```shell
conda create -n ldpred python=3.8.10
source activate ldpred
python -m venv myenv
source myenv/bin/activate
pip --default-timeout=1000 install h5py
pip --default-timeout=1000 install scipy==1.7.0
pip install plinkio
pip --default-timeout=1000 install ldpred==1.0.10
pip install hickle
```

### 3.	Let us begin!

The data used in this tutorial is the test data.

#### 3.1	Step1: Coordinate data

The first step is a data synchronization step, where two or three data sets, genotypes, and summary statistics are synchronized. This generates a HDF5 file which contains the synchronized genotypes. 

**<u>!! It is very helpful to use --help for detailed options:</u>** 

```shell
ldpred coord --help
```

**<u>My code:</u>**

```shell
ldpred coord
--gf sim1_0_train ##LD Reference Genotype File.
--ssf sim1_0_ss.txt ##Summary Statistic File.
--vbim sim1_0_test.bim ##Validation SNP File. Individual Genotype data. (Optional)
--vgf sim1_0_test ##Validation genotyoe file. (Optional)
--rs SNP_ID ##SNP ID of summary data, can be rsid or other identity code.
--A1 ALT ##Effective allele of summary data
--A2 REF ##Non-Effective allele of summary data
--pos POS ##SNP base position of summary data
--chr CHR ##SNP chromosome of summary data
--pval PVAL ##Pvalue of summary data
--eff_type LOGOR ##You can choose {LINREG,LOGOR,OR}
--eff BETA ##Effect size of summary data
--se SE ##Standard error of summary data
--ncol N ##Sample size of summary data
--out test_1
```

Note: There are other codes not used here, for example, --info refers to the INFO score of summary data. Furthermore, --only-hm3 restricts the analysis to HapMap3 SNPs.  **<u>However, since the HapMap has been updated to HapMap3+, I would recommend to restrict the reference genotype to HapMap3+ variants using Plink before performing LDpred.</u>**



#### 3.2	Step2: Generate LDpred SNP weights

After generating the coordinated data file then the one can apply LDpred and run it on the synchronized dataset. This step generates two files. One is a LD file with LD information for the given LD radius, and the re-weighted effect estimates. The other file that LDpred generates contains the LDpred-adjusted effect estimates.

**<u>!! It is very helpful to use --help for detailed options:</u>**

```shell
ldpred gibbs --help
```

**<u>My code:</u>**

```shell
ldpred gibbs
--cf test_1 ##Coordinated file generated in last step
--ldr 10 ##LD radius, a value of M/3000, where M is the number of SNPs in the genome is recommended. This value will decide the running time and required memory.
--ldf LDF ##LD file
--out SNP_weight
--f 0.01 0.1 0.5 1 ##Fraction of causal variants used in Gibbs sampler
--hickle-ld ##Use hickle to save memory
--n-iter 100 ##Iteration times
```

**<u>Note: Using the test data may cause errors, which are due to the --vbim and --vgf commands at the last step.  Do not use these two commands in the last step to go through the programming based on test data.</u>**

#### 3.3	Generating Individual Risk Scores

This step calculates polygenic risk scores for the individuals in the validation data if given, otherwise it treats the LD reference genotypes as validation genotypes. 

**<u>!! It is very helpful to use --help for detailed options:</u>**

```
ldpred score --help
```

**<u>My code:</u>**

```shell
ldpred score
--gf sim1_0_train ##Validation genotype file
--rf SNP_weight ##SNP weights file
--f 0.01 0.1 0.5 1 ##causal variants fraction
--only-score ##only out put scores. I recommend construct the regression model using R instead of LDpred.
--pf-format FAM ##You can choose {LSTANDARD,STANDARD,FAM}
--pf sim1_0_train.fam ##A file with individual IDs and phenotypes
--out ind_score
```



[^1]: Vilhjálmsson BJ, Yang J, Finucane HK, et al. Modeling Linkage Disequilibrium Increases Accuracy of Polygenic Risk Scores. *Am J Hum Genet*. 2015;97(4):576-592. doi:10.1016/j.ajhg.2015.09.001
[^2]: [LDpred-2 - Basic Tutorial for Polygenic Risk Score Analyses (choishingwan.github.io)](https://choishingwan.github.io/PRS-Tutorial/ldpred/)

[^3]: https://github.com/privefl/paper-misspec/tree/main/code
[^4]: [Q and A · bvilhjal/ldpred Wiki (github.com)](https://github.com/bvilhjal/ldpred/wiki/Q-and-A)
[^5]: [2. QC of Target Data - Basic Tutorial for Polygenic Risk Score Analyses (choishingwan.github.io)](https://choishingwan.github.io/PRS-Tutorial/target/)


If you have any questions, please contact: gesangmeiduo@pku.edu.cn
