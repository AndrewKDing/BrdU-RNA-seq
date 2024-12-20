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

```{r modeling}


# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)

#Specify model 
max_form <- ~ infected + cd8_depleted + aids + sive + plus_minus + (1|animal)
interact_form <- ~ cd8_depleted + infected + timepoint*plus_minus + (1|animal)

#estimate weights using linear mixed model
v_obj_dream <- voomWithDreamWeights(dge, form, animals, BPPARAM = param)
#fitmm <- dream(v_obj_dream, form, animals)
#fitmm <- eBayes(fitmm)

topTable(fitmm, coef = "plus_minusplus", number=Inf)

v_obj_interact <- voomWithDreamWeights(dge, interact_form, animals, BPPARAM = param)
fit_interact <- dream(v_obj_interact, interact_form, animals)
fit_interact <- eBayes(fit_interact)

topTable(fit_interact, coef = "timepointnec:plus_minusplus", number=Inf) %>%
  filter(P.Value <= 0.05) %>%
  arrange(desc(logFC))
```