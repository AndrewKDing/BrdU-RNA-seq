library(readxl)
library(tidyverse)
library(ggplot)

# CHECKPOINTS
# Calculate pca and plot in ggplot
      # CHECKPOINT: calculate normalized reads // DONE
# PRIORITY PLOT: colored by date of library preparation
# PRIORITY PLOT: colored by date of RNA extraction
# PRIORITY PLOT: colored by date of sorting
# PLOT: BrdU+, colored by non-infected, 21dpi, and necropsy
# PLOT: BrdU-, colored by non-infected, 21dpi, and necropsy
# PLOT: BrdU+, colored by %BrdU
# PLOT: 21dpi only, non-infected, SIVE, SIVnoE
# plot: necropsy only, non-infected, SIVE, SIVnoE


# Principal components analysis -------------------------------------------

#going fro 500 to 20,000 genes matters a LOT. it increases principle component 1 variance from like 10% to 42%. I wonder if this means something about 
# gene signatures instead of specific differentially expressed genes?
# Using more genes looks correct, non-infected samples appear to be clustering with one another


#Note to self, data should be subsetted PRIOR to multi-dimensional scaling
