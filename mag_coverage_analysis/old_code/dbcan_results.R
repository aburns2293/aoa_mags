
# cgc output

cgc_new <- read_delim("~/aoa_mags/mag_mt_mapping/output/133_cazyme_CGC/cgc_standard.out", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

cgc_new$Culture <- 'Cul133'

cgc_full <- rbind(cgc_full, cgc_new)

write.csv(cgc_full, '~/aoa_mags/mag_coverage_analysis/cgc_output.csv', row.names = FALSE)

# hmmer output

caz_new <- read_delim("~/aoa_mags/mag_mt_mapping/output/5_cazyme_CGC/hmmer.out", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
caz_new$Culture <- 'Cul5'

caz_full <- rbind(caz_full, caz_new)

write.csv(caz_full, '~/aoa_mags/mag_coverage_analysis/caz_output.csv', row.names = FALSE)

# caz filtering based on Bei et al. parameters (eval < 1e-15, coverage > 80%)

caz_full <- caz_full %>% filter(`E Value` < 1e-15,
                                Coverage > 0.8)
