library("maftools")
library("TCGAbiolinks")
library("curl")

maf<- "/home1/vineet/pupilBio/PA220KH_somatic_variants_functotated.maf"

maffile<- read.maf(maf = maf, verbose = F)    #my original data

maf_data<- maffile@data   
maf_gene_summary<- maffile@gene.summary


#if we dont even do this i think than also the plot can ve made
GCDprojects<- getGDCprojects()
summary<- TCGAbiolinks::getProjectSummary("TCGA-CHOL")


query.maf <- GDCquery(
  project = "TCGA-CHOL",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query.maf)
maf_TCGA_data<-GDCprepare(query.maf)
maf_TCGA<- read.maf(maf_TCGA_data, isTCGA = T)

#merge Maf
merge_MAF<- merge_mafs(mafs=c(maffile,maf_TCGA),
                       useAll=T, verbose=T)


#plot maf summary

plotmafSummary(maf=maffile, rmOutlier = T, addStat = 'median',dashboard = T, titvRaw = F)
maffile@summary


#draw 
maf.titv= titv (maf=maffile, plot= F, useSyn=T)
plotTiTv(res=maf.titv)



#lollipop plot : add the name the gene  

lollipopPlot(
  maf = maffile,
  gene = 'TP53',
  AACol = 'Protein_Change',
  showMutationRate = TRUE
)

#    protein lollipop plot

plotProtein(gene="TP53", refSeqID ="NM_000546" )


#Rainfall Plots

rainfallPlot(maf=maffile, detectChangePoints = T, pointSize = 1.5)


#compare mutation load against TCGA cohort

tum.mutload=tcgaCompare(maf=maffile, cohortName="TUMOR-sample", logscale=T,capture_size=50  )

#VAF plotting: variont allele frequency
plotVaf(maf= maffile )


#:::::::::::::::::::: calculate VAF::::::::::::::::::::::::::::

maffile@data$t_alt_count
maffile@data$t_ref_count

#Calculate t_depth:
maffile@data$t_depth <- maffile@data$t_alt_count + maffile@data$t_ref_count

#calculate VAF
maffile@data$VAF <- (maffile@data$t_alt_count / maffile@data$t_depth) * 100

#Summarize VAF for Selected Genes
selected_genes <- c("TP53", "SMAD4", "APC", "ARID1A", "NUDT16L1","LRRC73",  "MED12"  )  # Replace with your genes
vaf_data <- subset(maffile@data, Hugo_Symbol %in% selected_genes)
vaf_summary <- aggregate(VAF ~ Hugo_Symbol, data = vaf_data, FUN = mean, na.rm = TRUE)
colnames(vaf_summary) <- c("Hugo_Symbol", "VAF")


#plot
ggplot(maffile@data, aes(x = Hugo_Symbol, y = VAF, size = t_depth)) +
  geom_point(aes(color = factor(Chromosome)), alpha = 0.7) +
  scale_color_brewer(palette = "Set1") +  # Contrasting color palette
  theme_minimal() +
  labs(title = "VAF Plot", x = "Gene", y = "VAF (%)") +
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),  # Font size and color
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),  # Angle x-axis labels
    plot.title = element_text(hjust = 0.5),  # Center the title
    panel.grid.major = element_line(color = "lightgrey", size = 0.2),  # Thin and light grey major grid lines
    panel.grid.minor = element_line(color = "lightgrey", size = 0.2)   # Thin and light grey minor grid lines
  )


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#..........Mutational Significance (MutSig) for your genes, you would typically need:

gene_mut_counts <- table(maffile@data$Hugo_Symbol)
gene_mut_counts <- as.data.frame(gene_mut_counts)
colnames(gene_mut_counts) <- c("Hugo_Symbol", "Mutation_Count")


background_rate <- 1e-6  # Example: 1 mutation per million bases

gene_mut_counts$Gene_Size <- 1500  # Assume 1500 bp per gene as an example
gene_mut_counts$Expected_Mutations <- gene_mut_counts$Gene_Size * background_rate
gene_mut_counts$Mutation_Count <- as.numeric(gene_mut_counts$Mutation_Count)
gene_mut_counts$Expected_Mutations <- as.numeric(gene_mut_counts$Expected_Mutations)

gene_mut_counts$p_value <- apply(
  gene_mut_counts, 1, 
  function(row) poisson.test(as.numeric(row["Mutation_Count"]), as.numeric(row["Expected_Mutations"]))$p.value
)

poisson.test(as.numeric(gene_mut_counts$Mutation_Count[1]), as.numeric(gene_mut_counts$Expected_Mutations[1]))


gene_mut_counts$p_value <- apply(
  gene_mut_counts, 1, 
  function(row) {
    mutation_count <- as.numeric(row["Mutation_Count"])
    expected_mutations <- as.numeric(row["Expected_Mutations"])
    poisson.test(mutation_count, expected_mutations)$p.value
  }
)

#>             Background Error Rate Calculation:
#......................................................................................................
library(data.table)
getSampleSummary(maffile)
getGeneSummary(maffile)

maf_data <- as.data.table(maffile@data)
# Extract relevant fields
required_columns <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "t_alt_count", "t_ref_count")
maf_subset <- maf_data[, ..required_columns]

# Calculate mutation frequency per gene: VAF variant allele frequency (Mutation freq= t_alt_count/(t_ref_count+t_alt_count))
maf_subset$mutation_frequency <- maf_subset$t_alt_count / (maf_subset$t_ref_count + maf_subset$t_alt_count)

#Estimate Background Error Rate
maf_subset$Background_Error_Rate <- maf_subset$t_ref_count / (maf_subset$t_alt_count + maf_subset$t_ref_count)

#Normalize Background Error Rate to Reads Per Million (RPM)
maf_subset$Total_Reads <- maf_subset$t_alt_count + maf_subset$t_ref_count
maf_subset$Background_Error_Rate_RPM <- (maf_subset$Background_Error_Rate / maf_subset$Total_Reads) * 1e6


# Summarize mutation frequencies by gene

gene_mutation_summary <- maf_subset[, .(Mutation_Frequency = mean(mutation_frequency, na.rm = TRUE)),
                                    by = Hugo_Symbol]
gene_mutation_summary <- gene_mutation_summary[order(-Mutation_Frequency)]  # Sort by mutation frequency

# View the top genes with highest mutation frequency
gene_mutation_summary


# Plot the top 10 most mutated genes by mutation frequency
ggplot(gene_mutation_summary[1:7], aes(x = reorder(Hugo_Symbol, -Mutation_Frequency), y = Mutation_Frequency, fill = Mutation_Frequency)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_viridis(option = "D", direction = 1) +  # Apply the viridis color scale
  labs(title = "Mutated Genes by Mutation Frequency", x = "Gene", y = "Mutation Frequency") +
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),  # Remove the major grid lines
    panel.grid.minor = element_blank(),  # Remove the minor grid lines
    axis.title.x = element_text(face = "bold", color = "black"),  # Bold and black font for x-axis title
    axis.title.y = element_text(face = "bold", color = "black", size= " 13"),  # Bold and black font for y-axis title
    axis.text.x = element_text(face = "bold", color = "black", size="13"),   # Bold and black font for x-axis labels
    axis.text.y = element_text(face = "bold", color = "black", size = "14")    # Bold and black font for y-axis labels
  )


#################################################################################################################

# Load necessary library
library(data.table)

# Step 1: Load Data
maf_data <- fread("PA221MH-N_getpileupsummaries.table")

# Step 2: Extract Relevant Fields
maf_data <- maf_data[!(contig == "contig"), .(contig, position, ref_count, alt_count, allele_frequency)]
maf_data <- maf_data[, .(AF = allele_frequency, ref_count, alt_count)]

# Step 3: Remove Zeros (Optional)
maf_data <- maf_data[AF > 0]

# Step 4: Calculate Median Background Mutation Level
median_background <- median(maf_data$AF, na.rm = TRUE)
print(paste("Median Background Mutation Level:", median_background))

# Step 5: Histogram of Allele Frequencies
hist(maf_data$AF, breaks = 50, main = "Distribution of Allele Frequencies", 
     xlab = "Allele Frequency", col = "skyblue")



max_af <- max(maf_data$AF)
print(paste("Max Allele Frequency:", max_af))

# Step 6: Check for Perfect Allele Frequencies (AF = 1)
perfect_af_count <- sum(maf_data$AF == 1)
print(paste("Number of Variants with AF = 1:", perfect_af_count))

# Step 7: Add Calculated Fields for Total Reads and Background Error Rate
maf_data[, Total_Reads := ref_count + alt_count]
maf_data[, Background_Error_Rate := ref_count / Total_Reads]

# Step 8: Calculate RPM (Reads Per Million)
maf_data[, RPM := (Background_Error_Rate / Total_Reads) * 1e6]

# Step 9: Summarize RPM
rpm_summary <- summary(maf_data$RPM)
print("RPM Summary:")
print(rpm_summary)

# Step 10: Confidence Threshold and RPM Required for Mutation Detection
confidence_threshold <- 2 * median_background
rpm_required <- (confidence_threshold / median_background) * 1e6
print(paste("Confidence Threshold:", confidence_threshold))
print(paste("Reads Per Million Required:", rpm_required))

# Output the Results as a List
results <- list(
  Median_Background = median_background,
  Max_AF = max_af,
  Perfect_AF_Count = perfect_af_count,
  RPM_Summary = rpm_summary,
  Confidence_Threshold = confidence_threshold,
  RPM_Required = rpm_required
)
print(results)

