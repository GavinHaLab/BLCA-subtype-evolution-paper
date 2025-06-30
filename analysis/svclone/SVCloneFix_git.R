# Script by Alan Min to extract the CCF for each sample based on the SVClone output and is in the pdf ending with _ccube_sv_results.pdf
library(ggplot2)

base_dir = "/path/toSVCloneRuns/without_SVspostassign/"

samples = list.files(base_dir)
# Remove the samples that we don't want to analyze
samples = samples[samples != 'Sample1' &
                      samples != 'Sample2' &
                      samples != 'Sample3']

clonal_props = c()
for (sample in samples) {
    print(paste0(base_dir, sample, "/ccube_out/", sample,"_ccube_sv_results.RData"))
    
    # This gives you a file called doubleBreakPtsRes
    load(paste0(base_dir, sample, "/ccube_out/", sample,"_ccube_sv_results.RData"))
    
    # This is the file that SVClone gave us that should have the cluster, the number of SVs, and the proportion
    sv_csv = read.table(paste0(base_dir, sample, "/ccube_out/", sample,"_subclonal_structure.txt"), header = T)
    sv_csv$CCF = doubleBreakPtsRes$res$full.model$ccfMean
    
    # Quality control check that proportion/CCF is always the same
    tum_frac = sv_csv$proportion / sv_csv$CCF
    if (length(tum_frac) > 1) {
        first_ele_only = rep(tum_frac[1], length(tum_frac))
    } else {
        first_ele_only = tum_frac
    }
    if (!all.equal(tum_frac, first_ele_only)) {
        stop("Proportion should be the CCF times tumor fraction")
    }
    
    # Add the number of SVs that are clonal. We're defining this as CCF >= 0.9
    clonal_prop = sum(sv_csv$n_ssms[sv_csv$CCF >= 0.9])/sum(sv_csv$n_ssms)
    clonal_props = c(clonal_props, clonal_prop)
    write.table(sv_csv, file=paste0(base_dir, sample, "/ccube_out/", sample,"_subclonal_structure_with_CCF_summary.txt"), row.names = F)
}

subclonal_props = 1 - clonal_props
df1 = data.frame(samples, prop=clonal_props, classification="Clonal")
df2 = data.frame(samples, prop=subclonal_props, classification="Subclonal")
df = rbind(df1, df2)

ggplot(df) + geom_bar(aes(x=samples, y=prop, fill=classification), stat='identity') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + ylab("Proportion of clonal SVs") + xlab("Sample Name")


















