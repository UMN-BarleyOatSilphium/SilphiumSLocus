# SilphiumSLocus

This repository contains code needed to replicate the three S-locus inference and/or mapping methods described in "Development of first linkage map for Silphium integrifolium (Asteraceae) enables identification of sporophytic self-incompatibility locus", by John H. Price, Andrew R. Raduski, Yaniv Brandvain, David L. Van Tassel, and Kevin P. Smith.

Hill_climb_method is also referred to as “approach 1” for S-allele inference in the publication. This method starts with a random assignment of S-locus genotypes to each individual in the population, and then uses a hill-climbing algorithm to find the set of genotypes which best fits a user-supplied set of dominance relationships between genotypes. The results of different dominance relationships may compared to find the ones which best fit the data, and the best genotype assignments may be used as phenotypes for QTL mapping.

Single_marker_regression uses only biallelic SNP markers. This method carries out a logistic regression following the formula Si = Mij + Pij + Mij:Pij, where the success S of the ith mating is predicted by the maternal genotype M, the paternal genotype P, and the interaction of those genotypes, all for the jth marker. The calculated P-value of the maternal by paternal genotype interaction may then be used to determine association with cross incompatibility, and thus the S-locus.

The file "CrossData.csv" contains the cross data referenced in the publication and serves as a model for the data format used by the three methods.

The file “all_rules.csv” contains all of the hypothesized dominance relationships tested in the hill-climbing method and serves as a model for the data format used by the method. This format is scalable and can include as few as one dominance relationship. Pairs of genotypes with an “Odds” of 1 are predicted to be compatible, and 0 are predicted to be compatible. 

The file “three_alleles_recode.csv” contains the biallelic SNP markers used in this experiment and is used for the “single_marker_regression” method. 
