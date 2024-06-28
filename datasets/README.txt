### README ####

Here, all datasets have been collected, being:

size_data.rds - data resulting from simulations varying tree size
trees.rds - a subset of all data used for figure 5 - here, the results of 10.000 trees per model are included. The full dataset (300k trees per model) is available upon request - unfortunately it is too large to be shared this way.
balance_data.rds - data resulting from simulations varying balance (Figures 7, 8 and S5)
families.txt - data table with summary statistic information, used for plotting
normalize.txt - data table with summary statistic information, used for plotting

Phylogenies folder - within this folder are all files used to create the family-level phylogenies:
source_trees - folder that contains the original supertrees
family_lists - folder that contains per tree, a list of species with the respective family
fractions    - folder that contains per tree, a list of the fraction of species per family included in the supertree
fracced      - supertrees broken down, selecting only those subtrees with >80% species included (see fractions folder). These trees are the ones used for Figures 2 and 3.
family_level_trees - folder that contains trees with families at the tips
backbone -  the backbone tree containing all families included in the analysis

lookup.R - script to look up family associated to species
select.R - script to break down supertrees based on the fraction of species included
family_level_trees.R - script to create trees with families at the tips.
backbone/create_taxon_backbone.R script to turn distance matrix into backbone tree. Then, the family-level_trees were pasted into the resulting backbone by hand