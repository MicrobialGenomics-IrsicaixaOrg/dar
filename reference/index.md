# Package index

## Basic Functions

- [`dar`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-package.md)
  [`dar-package`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/dar-package.md)
  : dar: Differential Abundance Analysis by Consensus
- [`recipe()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/recipe.md)
  : Create a Recipe for preprocessing data
- [`prep()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/prep.md)
  : Performs all the steps defined in a Recipe
- [`bake()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/bake.md)
  : Define consensus strategies from a Recipe
- [`cool()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/cool.md)
  : Extract results from defined bake

## Step Functions - Data Preprocessing

- [`step_filter_taxa()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_taxa.md)
  : Filter taxa based on across-sample OTU abundance criteria
- [`step_filter_by_prevalence()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_prevalence.md)
  : Filter taxa by prevalence
- [`step_filter_by_abundance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_abundance.md)
  : Filter taxa by abundance
- [`step_filter_by_rarity()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_rarity.md)
  : Filter taxa by rarity
- [`step_filter_by_variance()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_filter_by_variance.md)
  : Filter taxa by variance
- [`step_subset_taxa()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_subset_taxa.md)
  : Subset taxa by taxonomic level
- [`step_rarefaction()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_rarefaction.md)
  : Resample an OTU table such that all samples have the same library
  size.

## Step Functions - Differential Abundance Methods

- [`step_aldex()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_aldex.md)
  : ALDEx2 analysis
- [`step_corncob()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_corncob.md)
  : corncob analysis
- [`step_deseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_deseq.md)
  : DESeq2 analysis
- [`step_lefse()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_lefse.md)
  : lefse analysis
- [`step_maaslin()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_maaslin.md)
  : MaAsLin2 analysis
- [`step_metagenomeseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_metagenomeseq.md)
  : MetagenomeSeq analysis
- [`step_wilcox()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/step_wilcox.md)
  : Wilcox analysis

## Exploration

- [`abundance_plt()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/abundance_plt.md)
  : Abundance boxplot

- [`corr_heatmap()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/corr_heatmap.md)
  :

  Plot otuput of the `overlap_df` function as a heatmap.

- [`exclusion_plt()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/exclusion_plt.md)
  : Plot the number of shared DA OTUs between methods.

- [`find_intersections()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/find_intersections.md)
  : Finds common OTU between method results

- [`intersection_df()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/intersection_df.md)
  : Returns data.frame with OTU intersection between methods

- [`intersection_plt()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/intersection_plt.md)
  : Plot results using UpSet plot

- [`mutual_plt()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/mutual_plt.md)
  : Mutual finding plot

- [`overlap_df()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/overlap_df.md)
  : Overlap of significant OTUs between tested methods.

## Recipe Getters and Adders

- [`add_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_tax.md)
  : Adds taxonomic level of interest in the Recipe.
- [`add_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/add_var.md)
  : Adds variable of interest to the Recipe
- [`get_phy()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_phy.md)
  : Returns phyloseq from Recipe-class object
- [`get_tax()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_tax.md)
  : Returns tax_info from Recipe-class object
- [`get_var()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/get_var.md)
  : Returns var_info from Recipe-class object
- [`otu_table()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/otu_table.md)
  : Extracts otu_table from phyloseq inside a Recipe
- [`phy_qc()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/phy_qc.md)
  : Phyloseq Quality Control Metrics
- [`sample_data()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/sample_data.md)
  : Extracts sample_data from phyloseq inside a Recipe
- [`steps_ids()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/steps_ids.md)
  : Get step_ids from recipe
- [`tax_table()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/tax_table.md)
  : Extracts tax_table from phyloseq inside a Recipe
- [`required_deps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/required_deps.md)
  : Returns required pakcages for Recipe object
- [`zero_otu()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/zero_otu.md)
  : Extract outs with all 0 values in at least on level of the variable

## Import and Export

- [`export_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/export_steps.md)
  : Export step parameters as json.
- [`import_steps()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/import_steps.md)
  : Import steps from json file
- [`read_data()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/read_data.md)
  [`validate_otu()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/read_data.md)
  [`validate_sample_data()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/read_data.md)
  [`validate_tax_table()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/read_data.md)
  [`validate_phyloseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/read_data.md)
  [`read_phyloseq()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/read_data.md)
  [`read_file()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/read_data.md)
  : Loads Phyloseq data

## Other functions

- [`contains_rarefaction()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/contains_rarefaction.md)
  : Checks if Recipe contains a rarefaction step
- [`rand_id()`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/rand_id.md)
  : Make a random identification field for steps
- [`PrepRecipe-class`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/PrepRecipe-class.md)
  : PrepRecipe-class object
- [`show(`*`<PrepRecipe>`*`)`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/Recipe-class.md)
  : Recipe-class object

## Data

- [`metaHIV_phy`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/metaHIV_phy.md)
  : Phyloseq object from metaHIV project
- [`pastry_df`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/pastry_df.md)
  : Pastery data for step id generation
- [`test_rec`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/test_rec.md)
  : Recipe for metaHIV_phy data
- [`test_prep_rec`](https://microbialgenomics-irsicaixaorg.github.io/dar/reference/test_prep_rec.md)
  : PrepRecipe for metaHIV_phy data
