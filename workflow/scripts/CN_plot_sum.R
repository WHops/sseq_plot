#t Load required libraries
library(ggplot2)  # For creating plots
library(dplyr)    # For data manipulation
library(tidyr)    # For data tidying
library(readr)
# Read the input file which can be .bed or .bed.gz (auto-detected.)
#' @param file_path Path to the input file
#' @return Data frame with columns "chrom", "start", "end", "name", "score", "strand", "cell"
read_input_file <- function(file_path) {
  if (grepl("\\.bed\\.gz$", file_path)) {
    # Read compressed .bed.gz file
    return(read_delim(file_path, delim = "\t", col_names = c("chrom", "start", "end", "name", "score", "strand", "cell"), col_types = cols(.default = "c"), locale = readr::default_locale(), na = c("", "NA")))
  }else {
    # Read uncompressed .bed file
    return(read.table(file_path, header = FALSE, sep = "\t",
                      col.names = c("chrom", "start", "end", "name", "score", "strand", "cell")))
  }
}

# Modify read_mappability_file
#' @param file_path Path to the input file
#' @return Data frame with columns "chr", "start", "end", "mappability", "false_reads"
read_mappability_file <- function(file_path) {
  return(read.table(file_path, header = FALSE, sep = "",
                    col.names = c("chr", "start", "end", "mappability", "false_reads")))
}

# Extract the sample name from the input file path
#' @param input_file Path to the input file
#' @return Sample name as a character vector
extract_samplename <- function(input_file) {
  # Use regular expression to extract the sample name
  match <- regexpr('BAB\\d+', input_file, perl = TRUE)
  samplename <- regmatches(input_file, match)
  
  return(samplename)
}

# Given a concatenated bed file of read locations, count reads per bin separately for + and - reads,
# and normalize by a reference region to get estimated copy numbers.
#' @param data Input data frame
#' @param bin_size Size of the bins
#' @param window_start Start position of the window
#' @param window_stop Stop position of the window
#' @return List containing the processed output data frame
process_data <- function(data, bin_size, window_start, window_stop, ref_window_start, ref_window_stop) {
  
  if(!("-" %in% data$strand)){
    data[1,'strand'] = '-'
  }
  
  if(!("+" %in% data$strand)){
    data[1,'strand'] = '+'
  }
  # Count reads per bin separately for + and - strand.
  print(dim(data))
  output <- data %>%
    mutate(bin = (start - 1) %/% bin_size * bin_size + 1) %>%
    count(strand, bin) %>%
    spread(strand, n, fill = 0) %>%
    gather(strand, read_count, c("+", "-")) %>%
    mutate(read_count = ifelse(strand == "+", read_count, -read_count)) %>%
    complete(strand, bin = seq(min(bin), max(bin), by = bin_size), fill = list(read_count = 0))
  print("hi") 
  # Get mean 'reference' read counts in a reference region (here: chrX:150Mb-152Mb)
  sum_reads_region <- output %>%
    filter(bin >= ref_window_start & bin <= ref_window_stop) %>% 
    summarize(sum_read_count_region = mean(abs(read_count)) ) %>%
    ungroup()
  
  # Normalize our read count by expected read count
  output <- output %>%
    mutate(read_count_norm = read_count / (sum_reads_region$sum_read_count_region))
  
  # Only keep bins in our window of interest
  output <- output %>%
    filter(bin >= window_start & bin <= window_stop)
  
  return(list(output = output))
}

# Custom formatting function for genomic coordinates
#' @param x Vector of genomic coordinates
#' @return Formatted genomic coordinates
format_genomic_coordinates <- function(x) {
  print(x)
  return(sapply(x, function(coord) {
    if (coord >= 1e6) {
      return(paste0(coord / 1e6, "M"))
    } else if (coord >= 1e3) {
      return(paste0(coord / 1e3, "K"))
    } else {
      return(as.character(coord))
    }
  }))
}

# Plot read density
#' @param output Output data frame
plot_read_density <- function(output, mean_mappability, bin_size, window_start, window_stop, samplename='Input sample', chromosome) {
  x_axis_interval <- pretty((window_stop - window_start)/10)[1]#2e0 * bin_size
  
  # Based on mapping simulated reads, we know which regions of the genome typically do not align well. We highlight
  # those in the plot. 
  low_mean_regions <- output %>%
    filter(mean_mappability < 50 | mean_false_reads > 5) %>%
    mutate(xmin = bin, xmax = bin + bin_size)
  
  window_start_pretty = max(ceiling(window_start/bin_size)*bin_size,0)  #(window_start)[1]
  window_end_pretty = min(max(output$bin),floor(window_stop/bin_size)*bin_size)
   
  print(window_start)
  print(window_stop)
  print(window_start_pretty)
  print(window_end_pretty)
  print(x_axis_interval) 

  print(seq(window_start_pretty, window_end_pretty, x_axis_interval))

  print(max(output$bin))
  print(min(output$bin))

  output_plus = output
  output_plus$bin = output_plus$bin+(bin_size*0.99999)
  output = rbind(output, output_plus)
  p = ggplot(output, aes(x = bin, y = read_count_norm, fill = strand)) 



  if (nrow(low_mean_regions) > 0){
    p = p + geom_rect(data = low_mean_regions, aes(xmin = xmin, xmax = xmax, ymin = -3, ymax = 3, y = Inf), fill = "lightgrey", alpha = 0.5, color = NA)
  }
    p = p +  geom_hline(yintercept = c(-3, -2, -1, 0, 1, 2, 3), color = "grey") +
    geom_area(aes(group = strand, fill=strand)) +
    scale_x_continuous(name = chromosome, breaks = seq(window_start_pretty, window_end_pretty, x_axis_interval), labels = format_genomic_coordinates) +
    scale_y_continuous(name = "Estimated copy number", breaks = c(-3, -2, -1, 0, 1, 2, 3), labels = c(3, 2, -1, 0, 1, 2, 3)) +
    theme_minimal() +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(y = 'Estimated copy number') +
    ggtitle(paste0(samplename, "Orientation-specific CN estimates in ", as.character(bin_size / 1000), " kbp bins")) +
    scale_fill_manual(values = c("+" = rgb(243 / 255, 165 / 255, 97 / 255), "-" = rgb(103 / 255, 139 / 255, 139 / 255)))
 
  print(p)
  return(p)
}

# Calculate mean mappability and false reads per bin
#' @param mappability_data Data frame containing mappability information
#' @param bin_size Size of the bins
#' @return Data frame with mean mappability and mean false reads per bin
calculate_mean_mappability <- function(mappability_data, bin_size) {
  mappability_data <- mappability_data %>%
    mutate(bin = (start - 1) %/% bin_size * bin_size + 1) %>%
    group_by(bin) %>%
    summarize(mean_mappability = mean(mappability), mean_false_reads = mean(false_reads)) %>%
    ungroup()
  
  return(mappability_data)
}

create_bedgraph <- function(data, chromosome, strand, color, file_name, bin_size = 1, append = FALSE) {
  strand_data <- data[data$strand == strand, ]
  track_line <- sprintf('track type=bedGraph name="My track (%s)" description="My custom track (%s strand)" color=%s visibility=full', strand, strand, color)
  if (!append) {
    writeLines(track_line, file_name)
  } else {
    cat(track_line, file = file_name, append = TRUE)
    cat('\n', file = file_name, append = TRUE)
  }
  chrom_size <- max(strand_data$bin)  # Get the maximum chromEnd value
  bedgraph_data <- data.frame(chromosome, strand_data$bin, pmin(strand_data$bin + bin_size, chrom_size), strand_data$read_count_norm)
  write.table(bedgraph_data, file = file_name, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, eol = "\n", na = "NA")
}


# Main script
if (sys.nframe() == 0){
  options(scipen = 999)
  # Read command-line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 8) {
    cat("Usage: Rscript script.R <input_file> <output_plotfile> <output_bedgraph> <chromosome> <window_start> <window_stop> <binsize>\n")
    quit(save="no", status=1)
  }
  
  input_file <- args[1]
  output_plotfile <- args[2]
  output_bedgraph <- args[3]
  chromosome <- args[4]
  window_start <- as.numeric(args[5])
  window_stop <- as.numeric(args[6])
  bin_size <- as.numeric(args[7])
  mappability_file <- args[8]

  ref_window_start = window_start - 2000000
  ref_window_stop = window_start

  # If one of ref window is negative, instead have the window after the region. 
  if (ref_window_start < 0){
    ref_window_start = window_stop
    ref_window_stop = window_stop + 2000000
  }

  samplename = extract_samplename(input_file)

  data <- read_input_file(input_file)
  #window_stop <- max(data$start)
  
  # Count normalized reads per bin and directionality.
  processed_data <- process_data(data, bin_size, window_start, window_stop, ref_window_start, ref_window_stop)
  output <- processed_data$output
  # Get also mappability data to highlight in the outcoming plot which regions could be un-trustworthy.
  mappability_data <- read_mappability_file(mappability_file)
  mean_mappability_data <- calculate_mean_mappability(mappability_data, bin_size)
  # Merge mean_mappability_data with output
   output <- output %>%
     left_join(mean_mappability_data, by = "bin")
  #mean_mappability_data = 'dummy'
  # Make and save the plot
  plot <- plot_read_density(output, mean_mappability_data, bin_size, window_start, window_stop, samplename, chromosome)
  ggsave(plot = plot, filename = output_plotfile, device = 'pdf', height = 3, width = 30, units = 'cm')
  
  # Filter the data frame by strand
  minus_strand <- output[output$strand == "-", ]
  plus_strand <- output[output$strand == "+", ]
  
  
  # Create BedGraph files for both strands which we can upload to UCSC.
  create_bedgraph(output, chromosome, "-", "103,139,139", output_bedgraph, bin_size = bin_size)
  create_bedgraph(output, chromosome, "+", "243,165,97", output_bedgraph, bin_size = bin_size, append = TRUE)
  

}


