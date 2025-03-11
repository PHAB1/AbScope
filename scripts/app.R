# Load necessary libraries
load_libraries <- function() {
  suppressWarnings(library(shiny))
  suppressWarnings(library(cowplot))
  suppressWarnings(library(alakazam))
  suppressWarnings(library(shazam))
  suppressWarnings(library(dplyr))
  suppressWarnings(library(scales))
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(igraph))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(Rtsne))
  suppressWarnings(library(stringdist))
  suppressWarnings(library(plotly))
  suppressWarnings(library(ggfortify))
  suppressWarnings(library(RColorBrewer))
  suppressWarnings(library(rgl))
  suppressWarnings(library(DT))
  suppressWarnings(library(factoextra))
  suppressWarnings(library(tidyverse))
  suppressWarnings(library(reshape2))
  suppressWarnings(library(htmlwidgets))
  suppressWarnings(library(svglite))
  suppressWarnings(library(plotrix))
  suppressWarnings(library(plot3D))
  suppressWarnings(library(MHTdiscrete))
  suppressWarnings(library(shinycssloaders))
  suppressWarnings(library(shinyjs))
}

# Function to summarize data
data_summary <- function(data, varname, groupnames) {
  data_sum <- data %>% 
    group_by(subject_id, gene) %>%
    summarize(sd = sd(seq_freq, na.rm = TRUE), 
              se = std.error(seq_freq, na.rm = TRUE), 
              seq_freq = mean(seq_freq, na.rm = TRUE))
  return(data_sum)
}

# Function to transform data into statistics
df2stat <- function(df, colSamples, colValues, colType, controlPattern, treatmentPattern, colSubType = FALSE) {
  df <- select(df, colSamples, colValues, colType)
  df$seq_freq <- df$seq_freq * 10000
  df <- df %>% 
    spread(key = colSamples, value = colValues) %>%
    mutate_all(~ifelse(is.na(.), 1, .))

  df$log2FoldChange <- log2((meanByPattern(df, treatmentPattern) / (meanByPattern(df, controlPattern))))
  
  df <- df %>%
    rowwise() %>%
    mutate(p_value = wilcox.test(c_across(contains(treatmentPattern)), c_across(contains(controlPattern)), paired = F)$p.value)
  
  df$p_adjusted <- p.adjust(df$p_value, method = "BH")
  df <- df %>% mutate_if(is.numeric, ~ round(., 2))

  tryCatch({ 
    df <- select(df, colType, log2FoldChange, p_value, p_adjusted, everything())
  }, error = function(e) {
    df <- select(df, colType, log2FoldChange, everything())
  })
}

# Function to calculate mean by pattern
meanByPattern <- function(df, pattern) {                          
  combined_pattern <- paste(pattern, collapse = "|")
  return(rowMeans(df[, grep(combined_pattern, colnames(df))]))
}

# Function to create colors from a column
create_colors <- function(column) {
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual', ]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  name_to_color <- data.frame(subject_id = unique(column), colors = col_vector[0:length(unique(column))])
  return(merge(data.frame(subject_id = column), name_to_color))
}

# Function to load data
load_data <- function(file_paths, perMaxSample) {
  df <- NULL
  for (f in file_paths) {
    airr_file <- airr::read_rearrangement(f)
    airr_file$seqorient <- NULL
    df <- tryCatch(
      expr = {
        df <- rbind(df, airr_file)
      },
      error = function(e) {
        df <- airr_file
      }
    )
  }
  df <- df %>%
    filter(!is.na(d_call), !is.na(c_call), !is.na(v_call), !is.na(j_call), !is.na(junction_length)) %>%
    mutate(v_call = sub(",.*", "", v_call)) %>%
    mutate(d_call = sub(",.*", "", d_call)) %>%
    mutate(j_call = sub(",.*", "", j_call)) %>%
    mutate(c_call = sub(",.*", "", c_call)) %>%
    mutate(v_call = sub("\\*.*", "", v_call)) %>%   
    mutate(d_call = sub("\\*.*", "", d_call)) %>%
    mutate(j_call = sub("\\*.*", "", j_call)) %>%
    mutate(c_call = sub("\\*.*", "", c_call)) %>%
    mutate(cprimer = sub("_.*", "", cprimer)) %>%
    mutate(across(c(duplicate_count, seq_len, clone_size_count, clone_size_freq), as.numeric))

  n_samples <- df %>%
    count(sample_id) %>%
    summarise(max_n = max(n)) %>%
    pull(max_n) * (perMaxSample / 100)
  
  n_samples <- round(n_samples)
  
  df <- df %>%
    group_by(sample_id) %>%
    group_modify(~ {
      n_group <- nrow(.x)
      .x %>% slice_sample(n = min(n_group, n_samples), weight_by = duplicate_count)
    }) %>%
    ungroup()
  
  return(df)
}

# Function to generate PCA plots
generate_pca_plot <- function(df, colorBy) {
  df <- df %>%
    select(sample_id, c_call, j_call, d_call, v_call) %>% 
    mutate(cjdv_call = paste0(c_call, j_call, d_call, v_call)) %>% 
    select(sample_id, cjdv_call) %>% 
    table(.) %>% 
    as.data.frame(.) %>% 
    mutate(Freq = Freq * 1000000 / sum(.$sample_id == sample_id)) %>%
    spread(key = "sample_id", value = "Freq")

  df <- as.data.frame(df)
  rownames(df) <- df$cjdv_call
  df$cjdv_call <- NULL
  df <- t(df)

  df <- df[, which(apply(df, 2, var) != 0)]

  pca_res <- prcomp(df, scale. = FALSE)
  eig_perc <- (pca_res$sdev^2 / sum(pca_res$sdev^2) * 100) %>% round(2)
  
  open3d(useNULL = T)
  plot3d(as.data.frame(pca_res$x[, 1:3]), size = 8, radius = 100, type = "s", col = df_tags$colors, box = FALSE, xlab = paste0("PC1 ", eig_perc[1]), ylab = paste0("PC2 ", eig_perc[2]), zlab = paste0("PC3 ", eig_perc[3]))
  if (colorBy == "subject_id") {
    text3d(as.data.frame(pca_res$x[, 1:3]), texts = c(df_tags$subject_id), cex = 0.7, pos = 3)
  } else {
    text3d(as.data.frame(pca_res$x[, 1:3]), texts = c(df_tags$sample_id), cex = 0.7, pos = 3)
  }
  rglwidget()
}

# Function to generate t-SNE plots
generate_tsne_plot <- function(df, colorBy, perplexity) {
  subsample <- function(data, column, n) {
    data %>%
      group_by(!!sym(column)) %>%
      sample_n(n, replace = TRUE) %>%
      ungroup()
  }

  bcr_subset <- subsample(df, colorBy, 1000)
  a <- bcr_subset$junction_aa
  b <- bcr_subset$junction_aa

  s_dist_hamming <- stringdistmatrix(a, b, method = "h", useNames = "string")
  s_dist_leveshtein <- stringdistmatrix(a, b, method = "lv", useNames = "string")
  s_dist <- cbind(s_dist_hamming, s_dist_leveshtein)
  s_dist[is.infinite(s_dist)] <- 0

  tsne_results <- Rtsne(as.matrix(s_dist), perplexity = perplexity, check_duplicates = FALSE, dims = 2)
  if (colorBy == "subject_id") {
    tsne_plot <- plot_ly(x = tsne_results$Y[, 1], y = tsne_results$Y[, 2], type = 'scatter', color = as.factor(bcr_subset$subject_id), text = as.factor(bcr_subset$junction_aa), mode = "markers")
  } else {
    tsne_plot <- plot_ly(x = tsne_results$Y[, 1], y = tsne_results$Y[, 2], type = 'scatter', color = as.factor(bcr_subset$sample_id), text = as.factor(bcr_subset$junction_aa), mode = "markers")
  }
  tsne_plot
}

# Function to generate network plots
generate_network_plot <- function(df, nSample, subject) {
  bcr_subset <- df %>%
    group_by(subject_id) %>%
    sample_n(size = nSample, replace = FALSE)

  a <- bcr_subset$junction_aa
  b <- bcr_subset$junction_aa

  s_dist <- stringdistmatrix(a, b, method = "h", useNames = "string")
  s_dist[s_dist < 2] <- 1
  s_dist[s_dist >= 2] <- 0

  network <- graph_from_adjacency_matrix(s_dist, mode = "undirected", weighted = T, diag = F)
  plot.igraph(network, vertex.size = 1.5, vertex.color = "grey", main = subject, edge.color = "orange", vertex.label = NA, edge.arrow.size = .1, layout = layout.fruchterman.reingold(network, niter = 1000))
}

# Function to generate abundance plots
generate_abundance_plot <- function(df, plotBy) {
  curve <- estimateAbundance(df, group = plotBy, ci = 0.95, nboot = 100, clone = "clone_id")
  abundance_plot <- plot(curve, legend_title = "", main = "Abundance Plot") + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 30), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) + theme_classic()
  abundance_plot_plotly <- ggplotly(abundance_plot)
  abundance_plot_plotly
}

# Function to generate diversity plots
generate_diversity_plot <- function(df, plotBy) {
  abund <- estimateAbundance(df, group = plotBy, nboot = 100)
  div <- alphaDiversity(abund, max_q = 6)
  diversity_plot <- plotDiversityCurve(div, legend_title = "") + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 30), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) + theme_classic()
  diversity_plot
}

# Function to generate VDJ usage plots
generate_vdj_usage_plot <- function(df, plotBy, cvdj_type_usage, mode_fg_usage) {
  df <- countGenes(df, gene = cvdj_type_usage, groups = plotBy, mode = mode_fg_usage)
  df$seq_freq <- df$seq_freq * 100
  vdj_usage_plot <- ggplot(df, aes_string(x = "gene", y = "seq_freq", fill = plotBy)) +
    theme_bw() +
    ggtitle("IGH-V-D-J Usage") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7)
  plotly::ggplotly(vdj_usage_plot)
}

# Function to generate CDR3 length plots
generate_cdr3_length_plot <- function(df, plotBy_usage_ablen, min_group_percent) {
  db_props <- aminoAcidProperties(df, seq = "junction", trim = TRUE, label = "cdr3")
  props_rm <- db_props %>% 
    count(subject_id, sample_id, c_call) %>%
    group_by(sample_id) %>%
    mutate(c_freq = n / sum(n)) %>%
    ungroup() %>%
    filter(c_freq > as.numeric(min_group_percent / 100))
  
  db_props <- semi_join(db_props, props_rm, by = c("sample_id", "c_call"))
  
  tmp_theme <- theme_bw() + theme(legend.position = "bottom")
  
  g1 <- ggplot(db_props, aes(x = c_call, y = cdr3_aa_length)) + tmp_theme +
    ggtitle("CDR3 length") + 
    xlab("Isotype") + ylab("Amino acids") +
    geom_boxplot(aes_string(fill = plotBy_usage_ablen))
  
  g2 <- ggplot(db_props, aes(x = c_call, y = cdr3_aa_gravy)) + tmp_theme + 
    ggtitle("CDR3 hydrophobicity") + 
    xlab("Isotype") + ylab("GRAVY") +
    geom_boxplot(aes_string(fill = plotBy_usage_ablen))
  
  gridPlot(g1, g2, ncol = 1)
}

# Function to generate statistics tables
generate_statistics_table <- function(df, control, treatment, cvdj_type, mode_fg) {
  treatment_samples <- unique(filter(df, subject_id == treatment)$sample_id)
  control_samples <- unique(filter(df, subject_id == control)$sample_id)
  
  family_gene <- countGenes(df, gene = cvdj_type, groups = c("sample_id"), mode = mode_fg)
  final_stat_df <- df2stat(family_gene, "sample_id", "seq_freq", "gene", control_samples, treatment_samples)
  
  final_stat_df
}

# Function to generate combinatorial statistics tables
generate_combinatorial_statistics_table <- function(df, control, treatment) {
  combinatorial_list <- list(
    list("c_call", "v_call"),
    list("c_call", "d_call"),
    list("c_call", "j_call"),
    list("j_call", "d_call"),
    list("j_call", "v_call"),
    list("d_call", "v_call")  
  )
  
  final_stat_df <- NULL
  
  treatment_samples <- unique(filter(df, subject_id == treatment)$sample_id)
  control_samples <- unique(filter(df, subject_id == control)$sample_id)
  
  for (i in 1:length(combinatorial_list)) {
    col_up <- unlist(combinatorial_list[[i]][1])
    col_down <- unlist(combinatorial_list[[i]][2])
    subset_vdj <- countGenes(df, gene = col_up, groups = c("sample_id", col_down), copy = "duplicate_count", mode = "gene")
    
    test <- df2stat(subset_vdj, "sample_id", "seq_freq", "gene", control_samples, treatment_samples, col_down)
    
    colnames(test)[1] <- "upStream"
    colnames(test)[4] <- "downStream"
    test <- select(test, "upStream", "downStream", everything())
    
    final_stat_df <- tryCatch(
      expr = {
        final_stat_df <- rbind(final_stat_df, test)
      },
      error = function(e) {
        final_stat_df <- test
      }
    )
  }
  
  final_stat_df
}

# Function to generate CDR3 intersection tables
generate_intersection_table <- function(df, db, cdr3_field) {
  df <- df %>%
    group_by(sample_id, v_call, d_call, j_call, junction_aa) %>%
    summarise(dup_count = sum(duplicate_count)) %>%
    ungroup() %>%
    arrange(desc(dup_count)) %>%
    filter(dup_count > 20) %>%
    as.data.frame(.) %>%
    mutate(junction_aa = substr(junction_aa, 2, as.integer(nchar(junction_aa)) - 1))
  
  merged_by_hamming <- data.frame()
  for (i in 1:nrow(db)) {
    for (j in 1:nrow(df)) {
      if ((nchar(db[cdr3_field][i, ]) != nchar(df["junction_aa"][j, ])) | !(grepl(df["v_call"][j, ], db["Heavy.V.Gene"][i, ]))) {
        next
      }
      hamming_dist <- sum(strsplit(db[cdr3_field][i, ], NULL)[[1]] != strsplit(df["junction_aa"][j, ], NULL)[[1]])
      if (hamming_dist < 4) {
        new_row <- cbind(db[i, ], df[j, ])
        merged_by_hamming <- rbind(merged_by_hamming, new_row)
      }
    }
  }
  
  merged_by_hamming
}
