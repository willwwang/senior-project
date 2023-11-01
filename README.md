This repo contains files relating to air samples taken in Upper Manhattan in early 2020.

# Sample run descriptions

Sample runs fell into four categories: 

- Acetate: Runs using an acetate + oxygenate standard. See [link 1](https://www.restek.com/en/products/reference-standards/reference-standards-by-sector/environmental--industrial-exposure/5723/) and [link 2](https://www.accustandard.com/prod0054799.html).

- CUNY: Runs using samples collected from the CUNY site in January-March 2020. Not including thin-phase EI or APCI runs.

- GS: Runs using the gas standard from Cylinder CC508270 (from 12/15). Also contains some blanks. All standard runs from days that CUNY or shutdown samples were analyzed in January-February 2021 are included.

- Shutdown: Runs using samples collected post-shutdown in East Harlem in May 2020.

Sample runs were originally outputted in a .D format. They were converted to .cdf using OpenChrom, then to .hdf using TERN.

The IGOR templates used for the analysis were based on 5-hour samples in the CUNY and shutdown periods and on an acetate/oxygenate standard fun.

# Repo files

This repo contains peak areas by compound, sample log info, R code to clean and plot the data, and resulting plots.

- benzene_toluene_ratio_analysis.R: Creates benzene_toluene_ratios.png and benzene_isooctane_ratios.png.

- compound_concentration_analysis.R: Creates concentrations.csv, sample_ptr_concs.png, and all _concs.png files.

- full_compounds.pxp: The Igor experiment containing the TERN analysis.

- GasStd_analysis.R: Creates GasStd_peakareainfo.csv, GasStd_summary.csv, and gas_standard_areas.png.

- ratio_analysis_by_compound.R: Creates butyl_acetate_ratio.png, decamethylcyclopentasiloxane_ratio.png, and selected_compound_ratios.png.

- sample_peak_area_cleaning.R: Converts the raw TERN output to the clean sample_ratios.csv file.

- sample_ratio_analysis.R: Creates compound_avg_ratios.csv, pre_post_ratios.png, and percent_change_histogram.png.

## data
Data outputs from TERN and from .R files. Also contains external PTR concentration data.

- Acetate_peakareas.csv, CUNY_peakareas.csv, GasStd_peakareas.csv, Shutdown_peakareas.csv: Raw outputs of compound peak areas from TERN.

- compound_avg_ratios.csv: The average ratios of each compound's peak area to benzene's peak area, pre- and post-shutdown.

- concentrations.csv: The calculated concentrations of selected compounds in each sample.

- GasStd_peakareainfo.csv: Compound peak areas in each gas standard run.

- GasStd_summary.csv: Average response factors for selected compounds in the gas standard runs.

- sample_ratios.csv: The peak areas, normalized peak areas, and peak area ratios for each compound in each sample. Also has associated info for each sample.

- VOC PTR data Oct2021 version_update.csv: Timestamped PTR concentration data in NYC in early 2020.

## info
Sample logs and other miscellaneous info needed for the analysis. Doesn't contain collected data.

- compound_classes.csv: Classifies compounds by "class" (more specific), and "type" (less specific).

- CUNY Sample Log.csv: The edited log for the CUNY samples. Edits to tube names to align with the sample file names, to flow rates that include ranges instead of one number, and to end times to correct errors and remove extraneous text. The first B93 sample was removed, since it was likely run as a preliminary sample and the tube was used twice.

- GasStd_filenames.csv: A list of gas standard filenames in the order they were indexed in TERN. 

- Shutdown Sample Log.csv: The edited log for the post-shutdown samples. Edits to tube names to align with the sample file names and to start times to remove extraneous text. 

- standard cylinder components 2020 2021 2022.xlsx: Concentration of compounds in the gas standard.

## plots
Plots outputted from R.

- benzene_concs.png, furfural_concs.png, toluene_concs.png: Sample concentrations of compounds over time.

- benzene_isooctane_ratios.png: A scatterplot of benzene peak areas against isooctane peak areas in the samples.

- benzene_toluene_ratios.png: A scatterplot of benzene concentrations against toluene concentrations in the samples.

- butyl_acetate_ratio.png, decamethylcyclopentasiloxane_ratio.png, selected_compound_ratios.png: Boxplots of compound:benzene peak area ratios before and after the shutdown.

- gas_standard_areas.png: Compound peak areas in the gas standards over the time that samples were run in early 2021.

- percent_change_histogram.png: The distribution of changes in peak area ratios from before to after the shutdown.

- pre_post_ratios.png: A scatterplot of pre-shutdown compound:benzene peak area ratios against post-shutdown ratios.

- sample_ptr_concs.png: Scatterplots of compound concentrations calculated from the samples against concentrations measured by PTR-MS.
