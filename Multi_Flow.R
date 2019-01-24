library(ezknitr)

list.folders <- c("LTEE_mCitrine_GAP1_Variants_T00", 
                  "LTEE_mCitrine_GAP1_Variants_T06", 
                  "LTEE_mCitrine_GAP1_Variants_T07", 
                  #"LTEE_mCitrine_GAP1_Variants_T08.1", 
                  "LTEE_mCitrine_GAP1_Variants_T08.3", 
                  "LTEE_mCitrine_GAP1_Variants_T11.1", 
                  "LTEE_mCitrine_GAP1_Variants_T11.2", 
                  "LTEE_mCitrine_GAP1_Variants_T13.1", 
                  "LTEE_mCitrine_GAP1_Variants_T13.2", 
                  "LTEE_mCitrine_GAP1_Variants_T14", 
                  "LTEE_mCitrine_GAP1_Variants_T15", 
                  "LTEE_mCitrine_GAP1_Variants_T18", 
                  "LTEE_mCitrine_GAP1_Variants_T22", 
                  "LTEE_mCitrine_GAP1_Variants_T25", 
                  "LTEE_mCitrine_GAP1_Variants_T27", 
                  "LTEE_mCitrine_GAP1_Variants_T29", 
                  "LTEE_mCitrine_GAP1_Variants_T34")

for(i in 1:length(list.folders)){
#fcs run sample name
name <- list.folders[i]

ezknit(file="~/Git_Synced_Projects/Flow_Cyto_Analysis/STA_LTEE_GAP1_Variants.Rmd", 
       out_dir ="~/Git_Synced_Projects/Flow_Cyto_Analysis/",  
       out_suffix=paste("_",name, sep=""), 
       params=list('name'=name), keep_html=TRUE, verbose=TRUE)
}