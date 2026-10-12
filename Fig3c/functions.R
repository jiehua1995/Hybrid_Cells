# Function: Read the genomic coverage matrix
# Read the genomic coverage matrix from .bed files in a specified folder, merge them into a single table
read_genomic_coverage <- function(folder, pattern="\\.regions\\.bed$", cutoff=1){
  # Load the required packages
  require(pbapply)
  require(data.table)
  
  # List the found files
  message("The selected folder is: ", folder,"\n")
  files <- list.files(path = folder, pattern = pattern, full.names = TRUE)
  message("Number of files found: ", length(files),"\n")
  if (length(files) == 0) {
    message("No files found matching pattern: ", pattern,"\n")
  } else {
    message("Found files:\n", paste(files, collapse = "\n"),"\n")
  }
  
  # Define a function to read and preprocess each .bed file
  read_and_rename <- function(f) {
    # Read the file into a data.table with no header
    dt <- fread(f, header = FALSE)
    # Assign meaningful column names
    # Expected columns: Chromosome, Start position, End position, Gene name, Coverage value
    setnames(dt, c("Chrom", "Start", "End", "Gene", "Coverage"))
    # Extract sample name from file name by removing the suffix ".regions.bed"
    sample_name <- sub(pattern, "", basename(f))
    # Rename the 'Coverage' column to the sample name
    # This allows us to merge multiple files with different sample-specific coverage columns
    setnames(dt, "Coverage", sample_name)
    # Return the processed data.table
    return(dt)
  }
  
  # Read and process all files in parallel with progress bar
  # Each file is read, column-renamed, and returned as a data.table
  list_dt <- pblapply(files, read_and_rename)
  
  # Merge all data.tables into one, using the first four columns as keys
  # This creates a wide-format table with one column per sample
  # `all = TRUE` ensures that all unique regions across all files are included (full outer join)
  merged_dt <- Reduce(function(x, y) merge(x, y, by = c("Chrom", "Start", "End", "Gene"), all = TRUE), list_dt)
  # Replace missing values (NA) in coverage columns with 0
  # This is necessary because not every region is present in every file
  # The first 4 columns are metadata, so we loop over the remaining columns
  for (col in names(merged_dt)[-(1:4)]) {
    # Replace NA with 0
    set(merged_dt, which(is.na(merged_dt[[col]])), col, 0)
    # Replace coverage ≤ cutoff with 0
    set(merged_dt, which(merged_dt[[col]] <= cutoff), col, 0)
  }
  # Return the merged data table
  return(merged_dt)
}

# Function: Visualize the genomic coverage data
visualize_genome <- function(data_genome, selected_chromosomes, facet_by_chrom = TRUE){
  require(ggplot2)
  require(viridis)
  require(data.table)
  
  # Filter the raw data
  data_genome <- subset(data_genome, Chrom %in% selected_chromosomes)
  
  # Melt the data
  samples_dt <- melt(
    data_genome,
    id.vars = c("Chrom", "Start", "End", "Gene"),
    measure.vars = setdiff(names(data_genome), c("Chrom", "Start", "End", "Gene")),
    variable.name = "Sample",
    value.name = "Coverage"
  )
  
  # Set the Chrom factor in the specified order
  samples_dt[, Chrom := factor(Chrom, levels = selected_chromosomes)]
  
  samples_dt[, Coverage_log := log2(Coverage + 1)]
  
  # Range
  x_range <- range(samples_dt$Coverage, na.rm = TRUE)
  x_log_range <- range(samples_dt$Coverage_log, na.rm = TRUE)
  y_max <- max(
    density(samples_dt$Coverage, na.rm = TRUE)$y,
    density(samples_dt$Coverage_log, na.rm = TRUE)$y
  )
  
  # Set colors based on facet_by_chrom
  if(facet_by_chrom){
    color_aes <- aes(color = Chrom, fill = Chrom)
    scale_fill_layer <- scale_fill_viridis_d(option = "turbo", drop = FALSE)
    scale_color_layer <- scale_color_viridis_d(option = "turbo", drop = FALSE)
  } else {
    color_aes <- aes(color = "grey50", fill = "grey50")
    scale_fill_layer <- scale_fill_manual(values = "grey50", guide = "none")
    scale_color_layer <- scale_color_manual(values = "grey50", guide = "none")
  }
  
  # p1: Original distribution
  p1 <- ggplot(samples_dt) +
    geom_histogram(aes(x = Coverage, y = after_stat(density)), bins = 50, alpha = 0.4, position = "identity") +
    geom_density(aes(x = Coverage, y = after_stat(density)), linewidth = 0.8, alpha = 0.7) +
    facet_wrap(~Sample, scales = "fixed") +
    color_aes +
    scale_fill_layer +
    scale_color_layer +
    coord_cartesian(xlim = x_range, ylim = c(0, y_max)) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightblue", color = "black"),
      strip.text = element_text(color = "black")
    ) +
    labs(
      title = "Genomic Coverage Distribution per Sample (Original Scale)",
      x = "Genomic Coverage (average sequencing depth per gene)",
      y = "Density (proportion)"
    )
  print(p1)
  
  # p2: log2 distribution
  p2 <- ggplot(samples_dt) +
    geom_histogram(aes(x = Coverage_log, y = after_stat(density)), bins = 50, alpha = 0.4, position = "identity") +
    geom_density(aes(x = Coverage_log, y = after_stat(density)), linewidth = 0.8, alpha = 0.7) +
    facet_wrap(~Sample, scales = "fixed") +
    color_aes +
    scale_fill_layer +
    scale_color_layer +
    coord_cartesian(xlim = x_log_range, ylim = c(0, y_max + 1)) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightblue", color = "black"),
      strip.text = element_text(color = "black")
    ) +
    labs(
      title = "Genomic Coverage Distribution per Sample (log2 scale + 1)",
      x = "Log2(Genomic Coverage + 1)",
      y = "Density (proportion)"
    )
  print(p2)
  
  # Summary table: statistics by Sample + Chrom
  summary_dt <- samples_dt[, {
    total_rows = .N
    nonzero_rows = sum(Coverage > 0)
    nonzero_vals = Coverage[Coverage > 0]
    avg_nonzero = ifelse(nonzero_rows > 0, mean(nonzero_vals), NA_real_)
    median_nonzero = ifelse(nonzero_rows > 0, median(nonzero_vals), NA_real_)
    q1_nonzero = ifelse(nonzero_rows > 0, quantile(nonzero_vals, 0.25), NA_real_)
    q3_nonzero = ifelse(nonzero_rows > 0, quantile(nonzero_vals, 0.75), NA_real_)
    coverage = sum(Coverage)
    list(total_rows = total_rows,
         nonzero_rows = nonzero_rows,
         avg_nonzero = avg_nonzero,
         median_nonzero = median_nonzero,
         q1_nonzero = q1_nonzero,
         q3_nonzero = q3_nonzero,
         total_coverage = coverage)
  }, by = .(Sample, Chrom)]
  
  # Sort by chromosome order
  summary_dt[, Chrom := factor(Chrom, levels = selected_chromosomes)]
  summary_dt[, total_rows := as.numeric(total_rows)]
  summary_dt[, nonzero_rows := as.numeric(nonzero_rows)]
  # Replace all NA in summary_dt except Sample and Chrom with 0
  cols_to_fix <- setdiff(names(summary_dt), c("Sample", "Chrom"))
  summary_dt[, (cols_to_fix) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = cols_to_fix]
  
  
  # Plot based on facet_by_chrom
  if(facet_by_chrom){
    # Plot by chromosome
    summary_long <- melt(summary_dt, id.vars = c("Sample", "Chrom"),
                         variable.name = "Metric", value.name = "Value")
    p3 <- ggplot(summary_long, aes(x = Chrom, y = Value, fill = Chrom)) +
      geom_col() +
      facet_grid(Metric ~ Sample, scales = "free_y") +
      scale_fill_viridis_d(option = "turbo", drop = FALSE) +
      theme_bw() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "lightblue", color = "black"),
        strip.text = element_text(color = "black")
      ) +
      labs(title = "Sample × Chromosome statistics comparison",
           y = "Value",
           x = NULL)
  } else {
    # Do not distinguish chromosomes, show by sample, but facet by metric
    summary_nochrom <- summary_dt[, lapply(.SD, sum), by = Sample, 
                                  .SDcols = setdiff(names(summary_dt), c("Sample","Chrom"))]
    summary_long_nc <- melt(summary_nochrom, id.vars = "Sample",
                            variable.name = "Metric", value.name = "Value")
    p3 <- ggplot(summary_long_nc, aes(x = Sample, y = Value, fill = "grey50")) +
      geom_col() +
      facet_wrap(~Metric, scales = "free_y", ncol = 3) +  # One row per metric
      scale_fill_manual(values = "grey50", guide = "none") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_line(),
        strip.background = element_rect(fill = "lightblue", color = "black"),
        strip.text = element_text(color = "black")
      ) +
      labs(title = "Sample statistics comparison (all chromosomes combined)",
           y = "Value",
           x = "Sample")
  }
  
  print(p3)
  
  
  return(summary_dt)
}

# Function: Read totalbases for genome
read_totalbases <- function(bam_summary){
  require(data.table)
  # Read the CSV file into a data.table
  dt <- data.table::fread(bam_summary,header=TRUE,sep=",")
  return(dt)
}

# Function: normalize the genomic coverage matrix (by total bases)
normalize_genome <- function(data_genome, genome_totalbases, selected_chromosomes){
  #Filter chromosomes first
  dt <- subset(data_genome, Chrom %in% selected_chromosomes)
  # Total bases
  genome_totalbases[, Total_Bases := as.numeric(Total_Bases)]
  # Create a vector to rename
  total_bases <- setNames(genome_totalbases$Total_Bases,genome_totalbases$Sample)
  # Effect factor
  mean_total_bases <- mean(total_bases)
  # Sample columns
  sample_cols <- intersect(colnames(dt), names(total_bases))
  # Normalize each sample
  for (col in sample_cols) {
    if(!is.na(total_bases[col]) && total_bases[col] > 0){
      dt[[col]] <- dt[[col]] / total_bases[col] * mean_total_bases
    } else {
      warning(paste("Sample", col, "has missing or zero Total_Bases"))
    }
  }
  
  return(dt)
}


# Function: Read the transcriptome counts data
read_transcriptome_counts <- function(transcriptome,gtf){
  require(data.table)
  data_transcriptome <- fread(transcriptome, header = TRUE, sep = "\t")
  data_gtf <- fread(gtf, header=F, sep="\t")
  setnames(data_gtf, c("chrom", "source", "feature", "start", "end","score", "strand", "frame", "attribute"))
  gtf_gene <- data_gtf[feature == "gene"]
  gtf_gene[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', attribute)]
  gtf_gene <- gtf_gene[, .(chrom,start,end,gene_id)]
  tpm_genes_merged <- merge(data_transcriptome,gtf_gene, by="gene_id",all=FALSE)
  return(tpm_genes_merged)
}

# Function: Normalize the transcriptome counts data
normalize_transcriptome <- function(data_transcriptome, selected_chromosomes){
  # Filter to selected chromosomoes
  dt <- subset(data_transcriptome, chrom %in% selected_chromosomes)
  # 1. Find the expression columns (assuming all other columns are expression values)
  expr_cols <- setdiff(names(dt), c("gene_id", "chrom", "start", "end"))
  # 2. Compute the total expression per sample
  totals <- dt[, lapply(.SD, sum, na.rm = TRUE), .SDcols = expr_cols]
  # 3. Compute the mean of total expression across all samples
  mean_total <- mean(as.numeric(totals))
  
  # 4. Normalize each column
  for (col in expr_cols) {
    total_col <- totals[[col]]
    scale_factor <- mean_total / total_col
    dt[[col]] <- dt[[col]] * scale_factor
  }
  return(dt)
}

# Function: Visualize the transcriptome data
visualize_transcriptome <- function(data_transcriptome, selected_chromosomes, facet_by_chrom = TRUE){
  require(ggplot2)
  require(viridis)
  require(data.table)
  
  # Filter the raw data
  dt <- subset(data_transcriptome, chrom %in% selected_chromosomes)
  
  # Melt the data
  samples_dt <- melt(
    dt,
    id.vars = c("chrom", "start", "end", "gene_id"),
    measure.vars = setdiff(names(dt), c("chrom", "start", "end", "gene_id")),
    variable.name = "Sample",
    value.name = "Coverage"
  )
  
  # Set the chrom factor in the specified order
  samples_dt[, chrom := factor(chrom, levels = selected_chromosomes)]
  
  samples_dt[, Coverage_log := log2(Coverage + 1)]
  
  # Range
  x_range <- range(samples_dt$Coverage, na.rm = TRUE)
  x_log_range <- range(samples_dt$Coverage_log, na.rm = TRUE)
  y_max <- max(
    density(samples_dt$Coverage, na.rm = TRUE)$y,
    density(samples_dt$Coverage_log, na.rm = TRUE)$y
  )
  
  # Set colors based on facet_by_chrom
  if(facet_by_chrom){
    color_aes <- aes(color = chrom, fill = chrom)
    scale_fill_layer <- scale_fill_viridis_d(option = "turbo", drop = FALSE)
    scale_color_layer <- scale_color_viridis_d(option = "turbo", drop = FALSE)
  } else {
    color_aes <- aes(color = "grey50", fill = "grey50")
    scale_fill_layer <- scale_fill_manual(values = "grey50", guide = "none")
    scale_color_layer <- scale_color_manual(values = "grey50", guide = "none")
  }
  
  # p1: Original distribution
  p1 <- ggplot(samples_dt) +
    geom_histogram(aes(x = Coverage, y = after_stat(density)), bins = 50, alpha = 0.4, position = "identity") +
    geom_density(aes(x = Coverage, y = after_stat(density)), linewidth = 0.8, alpha = 0.7) +
    facet_wrap(~Sample, scales = "fixed") +
    color_aes +
    scale_fill_layer +
    scale_color_layer +
    coord_cartesian(xlim = x_range, ylim = c(0, y_max)) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightblue", color = "black"),
      strip.text = element_text(color = "black")
    ) +
    labs(
      title = "Genomic Coverage Distribution per Sample (Original Scale)",
      x = "Genomic Coverage (average sequencing depth per gene)",
      y = "Density (proportion)"
    )
  print(p1)
  
  # p2: log2 distribution
  p2 <- ggplot(samples_dt) +
    geom_histogram(aes(x = Coverage_log, y = after_stat(density)), bins = 50, alpha = 0.4, position = "identity") +
    geom_density(aes(x = Coverage_log, y = after_stat(density)), linewidth = 0.8, alpha = 0.7) +
    facet_wrap(~Sample, scales = "fixed") +
    color_aes +
    scale_fill_layer +
    scale_color_layer +
    coord_cartesian(xlim = x_log_range, ylim = c(0, y_max + 1)) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightblue", color = "black"),
      strip.text = element_text(color = "black")
    ) +
    labs(
      title = "Genomic Coverage Distribution per Sample (log2 scale + 1)",
      x = "Log2(Genomic Coverage + 1)",
      y = "Density (proportion)"
    )
  print(p2)
  
  # Summary table: statistics by Sample + Chrom
  summary_dt <- samples_dt[, {
    total_rows = .N
    nonzero_rows = sum(Coverage > 0)
    nonzero_vals = Coverage[Coverage > 0]
    avg_nonzero = ifelse(nonzero_rows > 0, mean(nonzero_vals), NA_real_)
    median_nonzero = ifelse(nonzero_rows > 0, median(nonzero_vals), NA_real_)
    q1_nonzero = ifelse(nonzero_rows > 0, quantile(nonzero_vals, 0.25), NA_real_)
    q3_nonzero = ifelse(nonzero_rows > 0, quantile(nonzero_vals, 0.75), NA_real_)
    coverage = sum(Coverage)
    list(total_rows = total_rows,
         nonzero_rows = nonzero_rows,
         avg_nonzero = avg_nonzero,
         median_nonzero = median_nonzero,
         q1_nonzero = q1_nonzero,
         q3_nonzero = q3_nonzero,
         total_coverage = coverage)
  }, by = .(Sample, chrom)]
  
  # Sort by chromosome order
  summary_dt[, chrom := factor(chrom, levels = selected_chromosomes)]
  summary_dt[, total_rows := as.numeric(total_rows)]
  summary_dt[, nonzero_rows := as.numeric(nonzero_rows)]
  # Replace all NA in summary_dt except Sample and chrom with 0
  cols_to_fix <- setdiff(names(summary_dt), c("Sample", "chrom"))
  summary_dt[, (cols_to_fix) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = cols_to_fix]
  
  
  # Plot based on facet_by_chrom
  if(facet_by_chrom){
    # Plot by chromosome
    summary_long <- melt(summary_dt, id.vars = c("Sample", "chrom"),
                         variable.name = "Metric", value.name = "Value")
    p3 <- ggplot(summary_long, aes(x = chrom, y = Value, fill = chrom)) +
      geom_col() +
      facet_grid(Metric ~ Sample, scales = "free_y") +
      scale_fill_viridis_d(option = "turbo", drop = FALSE) +
      theme_bw() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "lightblue", color = "black"),
        strip.text = element_text(color = "black")
      ) +
      labs(title = "Sample × Chromosome statistics comparison",
           y = "Value",
           x = NULL)
  } else {
    # Do not distinguish chromosomes, show by sample, but facet by metric
    summary_nochrom <- summary_dt[, lapply(.SD, sum), by = Sample, 
                                  .SDcols = setdiff(names(summary_dt), c("Sample","chrom"))]
    summary_long_nc <- melt(summary_nochrom, id.vars = "Sample",
                            variable.name = "Metric", value.name = "Value")
    p3 <- ggplot(summary_long_nc, aes(x = Sample, y = Value, fill = "grey50")) +
      geom_col() +
      facet_wrap(~Metric, scales = "free_y", ncol = 3) +  # One row per metric
      scale_fill_manual(values = "grey50", guide = "none") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_line(),
        strip.background = element_rect(fill = "lightblue", color = "black"),
        strip.text = element_text(color = "black")
      ) +
      labs(title = "Sample statistics comparison (all chromosomes combined)",
           y = "Value",
           x = "Sample")
  }
  
  print(p3)
  
  
  return(summary_dt)
}

# Function: Compare the transcriptome data to genomic data to get expression level
compare_transcriptome_genome <- function(data_transcriptome, data_genome, sample_map, genome_col_map){
  # Print the dimensions of input data.
  message("Transcriptome data rows: ", nrow(data_transcriptome), ", columns: ", ncol(data_transcriptome))
  message("Genomic data rows: ", nrow(data_genome), ", columns: ", ncol(data_genome))
  
  # Check the different of the gene ids
  missing_genes <- setdiff(data_transcriptome$gene_id, data_genome$Gene)
  if (length(missing_genes) > 0) {
    message("Warning: Some genes in transcriptome not found in genome: ", paste(head(missing_genes, 10), collapse=", "))
  }
  
  # Merge the two tables.
  merged <- merge(data_transcriptome, data_genome, by.x = "gene_id", by.y = "Gene", all.x = TRUE)
  message("Merged data rows: ", nrow(merged), ", columns: ", ncol(merged),"\n")
  
  
  # Initialize the result
  pos_cols <- c("chrom","start","end")
  result <- merged[, c(pos_cols, "gene_id"), with=FALSE]
  setnames(result, "gene_id", "gene")
  
  # Loop over each group
  for(group in names(sample_map)){
    genome_col <- genome_col_map[[group]]
    for(sample in sample_map[[group]]){
      if(!(sample %in% colnames(merged))){
        warning(sprintf("Sample %s does not exist, skipping", sample))
        next
      }
      if(!(genome_col %in% colnames(merged))){
        warning(sprintf("Genome column %s does not exist, skipping group %s", genome_col, group))
        next
      }
      
      # RNA / DNA ratio
      ratio <- ifelse(merged[[genome_col]] > 0, merged[[sample]] / merged[[genome_col]], 0)
      
      result[[sample]] <- ratio
    }
  }
  message("Expression level rows: ", nrow(result), ", columns: ", ncol(result))
  return(result)
}


# Function: visualize the comparison matrix
visualize_comparison_distribution <- function(data_comparison){
  require(data.table)
  require(ggplot2)
  require(viridis)
  require(reshape2)
  
  # Select expression columns, excluding metadata
  samples_dt <- data_comparison[, !c("chrom", "start", "end", "gene"), with = FALSE]
  
  # Transform to long format for ggplot
  samples_long <- data.table::melt(samples_dt, measure.vars = names(samples_dt),
                                   variable.name = "Sample", value.name = "Expression")
  samples_long[, Expression_log := log2(Expression + 1)]
  
  # Calculate ranges for x-axis and y-axis
  x_range <- range(samples_long$Expression, na.rm = TRUE)
  x_log_range <- range(samples_long$Expression_log, na.rm = TRUE)
  y_max <- max(
    density(samples_long$Expression, na.rm = TRUE)$y,
    density(samples_long$Expression_log, na.rm = TRUE)$y
  )
  
  # Distribution plots of original and log2 expression
  p1 <- ggplot(samples_long, aes(x = Expression)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "steelblue", alpha = 0.5) +
    geom_density(color = "darkblue", linewidth = 1) +
    facet_wrap(~Sample, scales = "fixed") +
    coord_cartesian(xlim = x_range, ylim = c(0, y_max)) +
    theme_bw() +
    labs(title = "Expression distribution (original scale)",
         y = "Density")
  print(p1)
  p2 <- ggplot(samples_long, aes(x = Expression_log)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "darkgreen", alpha = 0.5) +
    geom_density(color = "forestgreen", linewidth = 1) +
    facet_wrap(~Sample, scales = "fixed") +
    coord_cartesian(xlim = x_log_range, ylim = c(0, y_max + 0.5)) +
    theme_bw() +
    labs(title = "Expression distribution (log2 scale + 1)",
         y = "Density")
  print(p2)
  
  # Summary statistics for each sample
  summary_dt <- samples_long[, {
    total_rows = .N
    vals = Expression  # Includes zero values
    avg_all = mean(vals)
    total_coverage = sum(vals)
    list(
      total_rows = total_rows,
      avg = avg_all,
      total_expression = total_coverage
    )
  }, by = Sample]
  
  # Transform summary_dt to long format
  summary_long <- melt(summary_dt, id.vars = "Sample",
                       variable.name = "Metric", value.name = "Value")
  
  # Draw summary statistics comparison plot 
  p3 <- ggplot(summary_long, aes(x = Sample, y = Value, fill = Sample)) +
    geom_col() +
    facet_wrap(~Metric, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Summary statistics comparison",
         y = "Value") +
    scale_fill_viridis_d(option = "viridis")
  print(p3)
  
  return(list(plot_original = p1,
              plot_log = p2,
              summary_table = summary_dt,
              plot_summary = p3))
}




visualize_comparison_distribution_by_chrom <- function(data_comparison, selected_chromosomes = NULL){
  require(data.table)
  require(ggplot2)
  require(viridis)
  require(reshape2)
  
  # If selected_chromosomes is not specified, use all chromosomes in the data
  if(is.null(selected_chromosomes)){
    selected_chromosomes <- unique(data_comparison$chrom)
  }
  
  # Filter chromosomes
  data_comparison <- data_comparison[chrom %in% selected_chromosomes]
  
  # Melt the data
  expr_cols <- setdiff(colnames(data_comparison), c("chrom","start","end","gene"))
  samples_dt <- data.table::melt(
    data_comparison,
    id.vars = c("chrom","start","end","gene"),
    measure.vars = expr_cols,
    variable.name = "Sample",
    value.name = "Expression"
  )
  
  # Convert to numeric
  samples_dt[, Expression := as.numeric(Expression)]
  samples_dt[, Expression_log := log2(Expression + 1)]
  
  # By chromosome order
  samples_dt[, chrom := factor(chrom, levels = selected_chromosomes)]
  
  
  # Original expression distribution plot
  p1 <- ggplot(samples_dt) +
    geom_histogram(aes(x = Expression, fill = chrom, color = chrom),
                   bins = 50, alpha = 0.4, position = "identity") +
    geom_density(aes(x = Expression, fill = chrom), alpha = 0.3, color = NA) +
    facet_wrap(~Sample, scales = "fixed") +
    scale_fill_viridis_d(option = "turbo", drop = FALSE) +
    scale_color_viridis_d(option = "turbo", drop = FALSE) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightblue", color = "black"),
      strip.text = element_text(color = "black")
    ) +
    labs(
      title = "Expression Distribution per Sample (Original Scale, by Chromosome)",
      x = "Expression",
      y = "Density",
      fill = "Chromosome",
      color = "Chromosome"
    )
  print(p1)
  
  # log2 expression distribution plot
  p2 <- ggplot(samples_dt) +
    geom_histogram(aes(x = Expression_log, fill = chrom, color = chrom),
                   bins = 50, alpha = 0.4, position = "identity") +
    geom_density(aes(x = Expression_log, fill = chrom), alpha = 0.3, color = NA) +
    facet_wrap(~Sample, scales = "fixed") +
    scale_fill_viridis_d(option = "turbo", drop = FALSE) +
    scale_color_viridis_d(option = "turbo", drop = FALSE) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightblue", color = "black"),
      strip.text = element_text(color = "black")
    ) +
    labs(
      title = "Expression Distribution per Sample (log2 Scale + 1, by Chromosome)",
      x = "log2(Expression + 1)",
      y = "Density",
      fill = "Chromosome",
      color = "Chromosome"
    )
  print(p2)
  
  # Summary table: statistics by Sample + Chrom
  summary_dt <- samples_dt[, {
    total_rows = .N
    nonzero_rows = sum(Expression > 0)
    nonzero_vals = Expression[Expression > 0]
    avg_nonzero = ifelse(nonzero_rows>0, mean(nonzero_vals), NA_real_)
    median_nonzero = ifelse(nonzero_rows>0, median(nonzero_vals), NA_real_)
    q1_nonzero = ifelse(nonzero_rows>0, quantile(nonzero_vals, 0.25), NA_real_)
    q3_nonzero = ifelse(nonzero_rows>0, quantile(nonzero_vals, 0.75), NA_real_)
    total_expression = sum(Expression)
    list(
      total_rows = total_rows,
      nonzero_rows = nonzero_rows,
      avg_nonzero = avg_nonzero,
      median_nonzero = median_nonzero,
      q1_nonzero = q1_nonzero,
      q3_nonzero = q3_nonzero,
      total_expression = total_expression
    )
  }, by = .(Sample, chrom)]
  
  summary_dt[, chrom := factor(chrom, levels = selected_chromosomes)]
  cols_to_fix <- setdiff(names(summary_dt), c("Sample","chrom"))
  summary_dt[, (cols_to_fix) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = cols_to_fix]
  
  # Plot by Sample × Chrom
  summary_long <- melt(summary_dt, id.vars = c("Sample","chrom"),
                       variable.name = "Metric", value.name = "Value")
  p3 <- ggplot(summary_long, aes(x = chrom, y = Value, fill = chrom)) +
    geom_col() +
    facet_grid(Metric ~ Sample, scales = "free_y") +
    scale_fill_viridis_d(option = "turbo", drop = FALSE) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(fill = "lightblue", color = "black"),
      strip.text = element_text(color = "black")
    ) +
    labs(title = "Sample × Chromosome Statistics Comparison",
         y = "Value",
         x = NULL)
  print(p3)
  
  return(summary_dt)
}


# Function: Check the distribution of the expression level matrix

check_distribution <- function(data){
  # Load required packages
  require(ggplot2)
  require(fitdistrplus)
  require(patchwork)
  # For evey column in data, analyze the distribution
}


# Function: visualize the comparison data
visualize_comparison_summary <- function(data_comparison){
  
  require(data.table)
  require(ggplot2)
  require(forcats)
  
  dt <- as.data.table(data_comparison)
  
  # Keep only rows in selected_chromosomes
  dt <- dt[chrom %in% selected_chromosomes]
  
  # Select sample columns (excluding the first 4 columns)
  sample_cols <- setdiff(colnames(dt), c("chrom","start","end","gene"))
  #sample_cols <- sample_cols[!grepl("ML82|Ras3", sample_cols)]
  
  # Convert the data to long format for easier statistics and plotting
  long_dt <- data.table::melt(dt,
                              id.vars = c("chrom","start","end","gene"),
                              measure.vars = sample_cols,
                              variable.name = "sample",
                              value.name = "expression")
  
  # Compute sum, mean, and median by Sample and Chrom
  summary_dt <- long_dt[, .(
    sum_expr = sum(expression, na.rm=TRUE),
    mean_expr = mean(expression, na.rm=TRUE),
    median_expr = median(expression, na.rm=TRUE)
  ), by = .(sample, chrom)]
  
  # Set the chromosome factor order based on selected_chromosomes to ensure plotting order
  summary_dt[, chrom := factor(chrom, levels = selected_chromosomes)]
  
  # Reusable function to draw bar plots
  plot_expr <- function(dt, y_col, y_label, title){
    ggplot(dt, aes(x = chrom, y = get(y_col), fill = sample)) +
      geom_col(position = position_dodge()) +
      theme_bw() +
      labs(title = title, y = y_label, x = "Chromosome") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # Plot the sum
  p_sum <- plot_expr(summary_dt, "sum_expr", "Sum of Expression", "Total Expression by Sample and Chromosome")
  print(p_sum)
  
  # Plot the mean
  p_mean <- plot_expr(summary_dt, "mean_expr", "Mean Expression", "Mean Expression by Sample and Chromosome")
  print(p_mean)
  
  # Plot the median
  p_median <- plot_expr(summary_dt, "median_expr", "Median Expression", "Median Expression by Sample and Chromosome")
  print(p_median)
  
  return(NULL)
}

