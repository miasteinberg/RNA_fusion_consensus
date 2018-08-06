# RNA_fusion_consensus

Written for use on Biowulf. Includes workflows for running Chimerascan, Mapsplice2, Ericscript, and STAR-Fusion. 

Fusions are ranked for each tool based on supporting read counts, and consensus is determined by mean rank. The rule that combines
the scores is currently set up for a specific metadata table, but this will eventually be changed to be more robust and generic.
