# heatmap

ggtree_heatmap <- function(mat,
                           scale = c("row", "column", "none"),
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           show_row_tree = FALSE,
                           show_col_tree = TRUE,
                           show_rownames = FALSE,
                           show_colnames = TRUE,
                           row_distance = "euclidean",
                           col_distance = "euclidean",
                           clustering_method = "complete",
                           border_color = NA,
                           colors = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
                           legend_name = "value",
                           col_angle = 45,
                           row_fontsize = 6,
                           col_fontsize = 10,
                           tree_ratio_row = 0.15,
                           tree_ratio_col = 0.18,
                           na_value = "grey90") {
  
  scale <- match.arg(scale)
  
  # -------- 1. Input check --------
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("mat must be a matrix or data.frame")
  }
  
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("Row_", seq_len(nrow(mat)))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("Col_", seq_len(ncol(mat)))
  }
  
  # -------- 2. Scaling --------
  mat_scaled <- mat
  
  if (scale == "row") {
    mat_scaled <- t(scale(t(mat_scaled)))
  } else if (scale == "column") {
    mat_scaled <- scale(mat_scaled)
  }
  
  mat_scaled[is.nan(mat_scaled)] <- 0
  
  # -------- 3. Clustering --------
  row_hclust <- NULL
  col_hclust <- NULL
  
  row_order <- rownames(mat_scaled)
  col_order <- colnames(mat_scaled)
  
  if (cluster_rows && nrow(mat_scaled) > 1) {
    row_hclust <- hclust(
      dist(mat_scaled, method = row_distance),
      method = clustering_method
    )
    row_order <- rownames(mat_scaled)[row_hclust$order]
  }
  
  if (cluster_cols && ncol(mat_scaled) > 1) {
    col_hclust <- hclust(
      dist(t(mat_scaled), method = col_distance),
      method = clustering_method
    )
    col_order <- colnames(mat_scaled)[col_hclust$order]
  }
  
  mat_scaled <- mat_scaled[row_order, col_order, drop = FALSE]
  
  # -------- 4. Convert to long table --------
  df_long <- as.data.frame(mat_scaled) %>%
    rownames_to_column("Feature") %>%
    pivot_longer(cols = -Feature, names_to = "Sample", values_to = "Value")
  
  # Reverse the y axis so the first row is at the top
  df_long$Feature <- factor(df_long$Feature, levels = rev(row_order))
  df_long$Sample  <- factor(df_long$Sample, levels = col_order)
  
  # Symmetric color scale, closer to pheatmap(scale="row")
  max_abs <- max(abs(df_long$Value), na.rm = TRUE)
  fill_limits <- c(-max_abs, max_abs)
  
  # -------- 5. Heatmap body --------
  p_heat <- ggplot(df_long, aes(x = Sample, y = Feature, fill = Value)) +
    geom_tile(color = border_color, linewidth = 0) +
    scale_fill_gradientn(
      colours = colors,
      limits = fill_limits,
      oob = squish,
      na.value = na_value,
      name = legend_name
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.text.y = if (show_rownames) {
        element_text(size = row_fontsize, colour = "black")
      } else {
        element_blank()
      },
      axis.text.x = if (show_colnames) {
        element_text(
          angle = col_angle,
          hjust = ifelse(col_angle == 0, 0.5, 1),
          vjust = ifelse(col_angle == 0, 0.5, 1),
          size = col_fontsize,
          colour = "black"
        )
      } else {
        element_blank()
      }
    )
  
  # -------- 6. Row tree --------
  p_row_tree <- NULL
  if (show_row_tree && !is.null(row_hclust)) {
    row_phylo <- as.phylo(row_hclust)
    
    p_row_tree <- ggtree(row_phylo, size = 0.3) +
      coord_cartesian(clip = "off") +
      scale_x_reverse() +
      theme_tree2() +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)
      )
  }
  
  # -------- 7. Column tree --------
  p_col_tree <- NULL
  if (show_col_tree && !is.null(col_hclust)) {
    col_phylo <- as.phylo(col_hclust)
    
    p_col_tree <- ggtree(col_phylo, size = 0.3) +
      layout_dendrogram() +
      theme_tree2() +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)
      )
  }
  
  # -------- 8. Assemble --------
  if (is.null(p_row_tree) && is.null(p_col_tree)) {
    p_final <- p_heat
    
  } else if (!is.null(p_row_tree) && is.null(p_col_tree)) {
    p_final <- p_row_tree + p_heat +
      plot_layout(widths = c(tree_ratio_row, 1 - tree_ratio_row))
    
  } else if (is.null(p_row_tree) && !is.null(p_col_tree)) {
    p_final <- p_col_tree / p_heat +
      plot_layout(heights = c(tree_ratio_col, 1 - tree_ratio_col))
    
  } else {
    p_blank <- patchwork::plot_spacer()
    
    p_final <- (p_blank + p_col_tree) /
      (p_row_tree + p_heat) +
      plot_layout(
        widths = c(tree_ratio_row, 1 - tree_ratio_row),
        heights = c(tree_ratio_col, 1 - tree_ratio_col)
      )
  }
  
  return(list(
    plot = p_final,
    heatmap = p_heat,
    row_tree = p_row_tree,
    col_tree = p_col_tree,
    mat_scaled = mat_scaled,
    row_order = row_order,
    col_order = col_order,
    row_hclust = row_hclust,
    col_hclust = col_hclust,
    data_long = df_long
  ))
}


# BP: dot plot

# BP: dot plot

bp_dotplot <- function(df_bp) {
  
  my_family <- "sans"
  
  df_plot <- df_bp %>%
    mutate(
      GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x))),
      Description_wrap = str_wrap(Description, width = 30)
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10)
  
  # Order the y axis by Gene Ratio
  df_plot <- df_plot %>%
    mutate(
      Description_wrap = factor(
        Description_wrap,
        levels = Description_wrap[order(GeneRatio_num)]
      )
    )
  
  p <- ggplot(
    df_plot,
    aes(
      x = GeneRatio_num,
      y = Description_wrap,
      size = Count,
      color = p.adjust
    )
  ) +
    geom_point(
      alpha = 1
    ) +
    
    scale_color_viridis_c(
      option = "D",
      direction = -1,
      name = "P.adjust",
      labels = scales::label_scientific(digits = 2)
    ) +
    
    scale_size(
      range = c(1.3, 3.8),
      name = "Count"
    ) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.1, 0.15)),
      guide = guide_axis(minor.ticks = TRUE)
    ) +
    
    labs(
      x = "Orthogroup Ratio",
      y = NULL
    ) +
    
    theme_bw(base_family = my_family) +
    theme(
      # Transparent backgrounds
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      strip.background = element_rect(fill = "transparent", color = NA),
      
      # Grid
      panel.grid = element_blank(),
      
      # Border
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.1
      ),
      
      # Avoid double-thick borders
      axis.line = element_blank(),
      
      # Axis text
      axis.text.x = element_text(
        size = 6,
        color = "black",
        family = my_family
      ),
      axis.text.y = element_text(
        size = 6,
        color = "black",
        family = my_family,
        lineheight = 0.9
      ),
      
      # Axis titles
      axis.title.x = element_text(
        size = 7,
        color = "black",
        family = my_family
      ),
      axis.title.y = element_text(
        size = 7,
        color = "black",
        family = my_family
      ),
      
      # Major ticks
      axis.ticks.x = element_line(
        linewidth = 0.2,
        color = "black"
      ),
      axis.ticks.y = element_line(
        linewidth = 0.2,
        color = "black"
      ),
      axis.ticks.length = unit(1.5, "pt"),
      
      # Minor ticks
      axis.minor.ticks.x = element_line(
        linewidth = 0.1,
        color = "black"
      ),
      axis.minor.ticks.y = element_line(
        linewidth = 0.1,
        color = "black"
      ),
      axis.minor.ticks.length.x = unit(0.8, "pt"),
      axis.minor.ticks.length.y = unit(0.8, "pt"),
      
      # Compact legends
      legend.position = "right",
      legend.title = element_text(
        size = 6,
        color = "black",
        family = my_family
      ),
      legend.text = element_text(
        size = 6,
        color = "black",
        family = my_family,
        margin = margin(t = 0, r = 0, b = 0, l = 2)
      ),
      legend.key.size = unit(6, "pt"),
      legend.key.width = unit(6, "pt"),
      legend.key.height = unit(6, "pt"),
      legend.key.spacing.x = unit(1, "pt"),
      legend.key.spacing.y = unit(1, "pt"),
      
      legend.spacing.x = unit(1, "pt"),
      legend.spacing.y = unit(1, "pt"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.spacing = unit(2, "pt"),
      
      # Plot margin
      plot.margin = margin(t = 4, r = 2, b = 0, l = 2)
    ) +
    
    guides(
      # Key: P.adjust on top
      color = guide_colorbar(
        order = 1,
        barheight = unit(18, "mm"),
        barwidth  = unit(3, "mm"),
        ticks = FALSE,
        frame.colour = NA,
        title.position = "top",
        title.hjust = 0.5
      ),
      
      # Key: Count below
      size = guide_legend(
        order = 2,
        keywidth = unit(6, "pt"),
        keyheight = unit(6, "pt"),
        ncol = 1,
        bycol = TRUE,
        label.position = "right",
        label.hjust = 0,
        override.aes = list(
          alpha = 0.95
        )
      )
    )
  
  return(p)
}

# Dot plot for all ONTs
all_dotplot <- function(df_all) {
  
  my_family="sans"
  
  df_plot <- df_all %>%
    mutate(
      GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x))),
      Description_wrap = str_wrap(Description, width = 30),
      ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "CC", "MF"))
    ) %>%
    group_by(ONTOLOGY) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  p <- ggplot(
    df_plot,
    aes(
      x = GeneRatio_num,
      y = reorder_within(Description_wrap, GeneRatio_num, ONTOLOGY),
      size = Count,
      color = p.adjust
    )
  ) +
    geom_point(
      alpha = 0.95
    ) +
    
    facet_grid(
      ONTOLOGY ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    
    scale_y_reordered() +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.04, 0.15)),
      guide = guide_axis(minor.ticks = TRUE)
    ) +
    
    scale_color_viridis_c(
      option = "D",
      direction = -1,
      name = "adjusted P",
      guide = guide_colorbar(
        order = 1,
        barheight = unit(18, "mm"),
        barwidth = unit(3, "mm"),
        ticks = FALSE,
        frame.colour = NA,
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    
    scale_size(
      range = c(1.3, 3.8),
      name = "Count"
    ) +
    
    labs(
      x = "Gene Ratio",
      y = NULL
    ) +
    
    theme_bw(base_family = my_family) +
    theme(
      # Transparent backgrounds
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      strip.background = element_rect(fill = "transparent", color = NA),
      
      # Grid
      panel.grid = element_blank(),
      
      # Border
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.1
      ),
      
      # Avoid double-thick borders
      axis.line = element_blank(),
      
      # Facet labels
      strip.placement = "outside",
      strip.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        size = 6,
        family = my_family,
        color = "black",
        face = "plain",
        margin = margin(t = 1, r = 0, b = 3, l = 0)
      ),
      strip.text.y.left = element_text(
        angle = 0,
        hjust = 1,
        size = 6,
        family = my_family,
        color = "black",
        face = "plain",
        margin = margin(t = 0, r = 3, b = 0, l = 1)
      ),
      strip.clip = "off",
      strip.switch.pad.grid = unit(0, "pt"),
      
      # Axis text
      axis.text.x = element_text(
        size = 6,
        color = "black",
        family = my_family
      ),
      axis.text.y = element_text(
        size = 6,
        color = "black",
        family = my_family,
        lineheight = 0.9
      ),
      
      # Axis titles
      axis.title.x = element_text(
        size = 7,
        color = "black",
        family = my_family
      ),
      axis.title.y = element_text(
        size = 7,
        color = "black",
        family = my_family
      ),
      
      # Major ticks
      axis.ticks.x = element_line(
        linewidth = 0.2,
        color = "black"
      ),
      axis.ticks.y = element_line(
        linewidth = 0.2,
        color = "black"
      ),
      axis.ticks.length = unit(1.5, "pt"),
      
      # Minor ticks
      axis.minor.ticks.x = element_line(
        linewidth = 0.1,
        color = "black"
      ),
      axis.minor.ticks.y = element_line(
        linewidth = 0.1,
        color = "black"
      ),
      axis.minor.ticks.length.x = unit(0.8, "pt"),
      axis.minor.ticks.length.y = unit(0.8, "pt"),
      
      # Compact legends
      legend.position = "right",
      legend.title = element_text(
        size = 6,
        color = "black",
        family = my_family
      ),
      legend.text = element_text(
        size = 6,
        color = "black",
        family = my_family,
        margin = margin(t = 0, r = 0, b = 0, l = 2)
      ),
      legend.key.size = unit(6, "pt"),
      legend.key.width = unit(6, "pt"),
      legend.key.height = unit(6, "pt"),
      legend.key.spacing.x = unit(1, "pt"),
      legend.key.spacing.y = unit(1, "pt"),
      
      legend.spacing.x = unit(1, "pt"),
      legend.spacing.y = unit(1, "pt"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.spacing = unit(2, "pt"),
      
      # Facet spacing
      panel.spacing.x = unit(0.2, "lines"),
      panel.spacing.y = unit(0.4, "lines"),
      
      # Plot margin
      plot.margin = margin(t = 4, r = 2, b = 0, l = 2)
    ) +
    
    guides(
      color = guide_colorbar(
        order = 1,
        barheight = unit(18, "mm"),
        barwidth = unit(3, "mm"),
        ticks = FALSE,
        frame.colour = NA,
        title.position = "top",
        title.hjust = 0.5
      ),
      size = guide_legend(
        order = 2,
        keywidth = unit(6, "pt"),
        keyheight = unit(6, "pt"),
        ncol = 1,
        bycol = TRUE,
        label.position = "right",
        label.hjust = 0,
        override.aes = list(
          stroke = 0,
          alpha = 0.95
        )
      )
    )
  
  return(p)
}

# treeplot
modify_treeplot_pub <- function(
    p,
    # ========= Font sizes (most commonly changed) =========
    tip_label_size = 1.9,        # Tip label size; ~1.7-2.3 corresponds to 5-7 pt
    clade_label_size = 1.9,      # Group/clade label size
    legend_title_size = 6,       # Legend title size (pt)
    legend_text_size = 5.5,      # Legend text size (pt)
    
    # ========= Points and colors =========
    point_size_range = c(1.2, 3.2),   # Point size mapping range
    padj_palette = "C",               # viridis continuous palette; "C" is common (blue-green-yellow)
    padj_direction = -1,              # -1 means smaller p.adjust is more prominent
    
    # ========= Highlight background =========
    highlight_alpha = 0.35,           # Group background transparency
    
    # ========= Legend colorbar =========
    colorbar_height_mm = 18,
    colorbar_width_mm = 3.5,
    
    # ========= Legend overall =========
    legend_key_width_pt = 10,
    legend_key_height_pt = 10,
    legend_spacing_y_pt = 6,
    legend_spacing_x_pt = 2,
    
    # ========= Background / axes =========
    transparent_bg = TRUE,
    remove_panel_border = TRUE,
    remove_axes = TRUE,
    
    # ========= Enable to compress tree width =========
    compress_tree = FALSE,
    tree_width_scale = 5,
    
    # ========= Enable to display immediately with ggview =========
    show_canvas = TRUE,
    canvas_width_mm = 180,
    canvas_height_mm = 120,
    canvas_dpi = 600,
    canvas_bg = "transparent",
    
    # ========= Convert to a static object when saving (key) =========
    make_static = TRUE,
    
    # ========= Layer indices (based on your current object structure) =========
    # Your current str() looks like:
    # 3 = tip labels
    # 5 = highlight rect
    # 6 = points
    # 7 = clade/group labels
    tip_label_layer = 3,
    highlight_layer = 5,
    point_layer = 6,
    clade_label_layer = 7,
    
    verbose = TRUE
) {
  
  # -----------------------------
  # Dependencies
  # -----------------------------
  if (!inherits(p, "ggplot")) {
    stop("Input object p must be a ggplot / ggtree / treeplot object.")
  }
  
  library(ggplot2)
  library(grid)
  
  if (show_canvas && !requireNamespace("ggview", quietly = TRUE)) {
    stop("The ggview package is required: install.packages('ggview')")
  }
  
  if (make_static && !requireNamespace("ggplotify", quietly = TRUE)) {
    stop("The ggplotify package is required: install.packages('ggplotify')")
  }
  
  # -----------------------------
  # Copy the object to avoid modifying the original directly
  # -----------------------------
  p2 <- p
  
  # -----------------------------
  # Check whether there are enough layers
  # -----------------------------
  n_layers <- length(p2$layers)
  
  if (verbose) {
    message("The current object has ", n_layers, " layers.")
  }
  
  idx_needed <- c(tip_label_layer, highlight_layer, point_layer, clade_label_layer)
  if (any(idx_needed > n_layers)) {
    stop(
      "The specified layer index is out of range. Current layer count = ", n_layers,
      "; please run first:\n",
      "sapply(p$layers, function(x) class(x$geom)[1])\n",
      "to confirm the position of each layer."
    )
  }
  
  # -----------------------------
  # 1) Change the tip label font size
  # -----------------------------
  p2$layers[[tip_label_layer]]$aes_params$size <- tip_label_size
  
  # -----------------------------
  # 2) Change the clade/group label font size
  # -----------------------------
  p2$layers[[clade_label_layer]]$aes_params$size <- clade_label_size
  
  # -----------------------------
  # 3) Change the highlight background transparency
  # -----------------------------
  p2$layers[[highlight_layer]]$aes_params$alpha <- highlight_alpha
  
  # -----------------------------
  # 4) Scale the point size
  # -----------------------------
  p2 <- p2 +
    scale_size(
      range = point_size_range,
      name = "Count",
      guide = guide_legend(
        order = 2,
        override.aes = list(color = "black"),
        keyheight = unit(legend_key_height_pt, "pt"),
        keywidth  = unit(legend_key_width_pt, "pt")
      )
    )
  
  # -----------------------------
  # 5) Override p.adjust colors
  # -----------------------------
  p2 <- p2 +
    scale_color_viridis_c(
      option = padj_palette,
      direction = padj_direction,
      name = "p.adjust",
      guide = guide_colorbar(
        order = 1,
        ticks = FALSE,
        frame.colour = NA,
        barheight = unit(colorbar_height_mm, "mm"),
        barwidth  = unit(colorbar_width_mm, "mm")
      )
    )
  
  # -----------------------------
  # 6) Optional: compress tree width
  # -----------------------------
  if (compress_tree) {
    if (!is.null(p2$data$x)) {
      xmax <- max(p2$data$x, na.rm = TRUE)
      if (is.finite(xmax) && xmax > 0) {
        p2$data$x <- p2$data$x / xmax * tree_width_scale
      }
    }
  }
  
  # -----------------------------
  # 7) Unified theme
  # -----------------------------
  bg_fill <- if (transparent_bg) "transparent" else "white"
  
  p2 <- p2 +
    theme_bw(base_size = 6) +
    theme(
      panel.background = element_rect(fill = bg_fill, color = NA),
      plot.background  = element_rect(fill = bg_fill, color = NA),
      legend.background = element_rect(fill = bg_fill, color = NA),
      legend.box.background = element_rect(fill = bg_fill, color = NA),
      
      panel.grid = element_blank(),
      panel.border = if (remove_panel_border) element_blank() else element_rect(fill = NA),
      
      axis.text = if (remove_axes) element_blank() else element_text(size = 5, color = "black"),
      axis.title = if (remove_axes) element_blank() else element_text(size = 7, color = "black"),
      axis.ticks = if (remove_axes) element_blank() else element_line(linewidth = 0.2, color = "black"),
      axis.line = if (remove_axes) element_blank() else element_line(linewidth = 0.2, color = "black"),
      
      legend.position = "right",
      legend.justification = "center",
      legend.title = element_text(size = legend_title_size, color = "black"),
      legend.text  = element_text(size = legend_text_size, color = "black"),
      legend.key.width  = unit(legend_key_width_pt, "pt"),
      legend.key.height = unit(legend_key_height_pt, "pt"),
      legend.spacing.y = unit(legend_spacing_y_pt, "pt"),
      legend.spacing.x = unit(legend_spacing_x_pt, "pt"),
      
      plot.margin = margin(t = 3, r = 3, b = 2, l = 2)
    )
  
  # -----------------------------
  # 8) Display / return
  # -----------------------------
  if (show_canvas) {
    return(
      p2 + ggview::canvas(
        canvas_width_mm,
        canvas_height_mm,
        units = "mm",
        dpi = canvas_dpi,
        bg = canvas_bg
      )
    )
  }
  
  # Key: if saving, freeze into a static grob first, then convert back to a plain ggplot
  # This way ggsave won't export blank output due to interactive geoms
  if (make_static) {
    g <- grid::grid.grabExpr(print(p2))
    p2 <- ggplotify::as.ggplot(g)
  }
  
  return(p2)
}




#emap plot

modify_emapplot_pub <- function(
    p,
    
    # ========= Text =========
    label_size = 6,                 # Node label text size (pt)
    label_force = NULL,             # ggrepel repulsion force; NULL = no change
    label_lineheight = 0.9,         # Label line spacing
    font_family = "sans",           # Font for all text
    
    # ========= Points =========
    point_size_range = c(2, 7),     # Point size mapping range
    point_alpha = 1,                # Point transparency
    
    # ========= Colors =========
    color_palette = "C",            # Continuous palette (usually corresponding to p.adjust)
    color_direction = -1,           # -1 = smaller p.adjust is more prominent
    
    # ========= Edges =========
    edge_linewidth = 0.2,
    edge_alpha = 0.25,
    edge_color = "grey30",          # Edge color; NULL = no change
    
    # ========= Ellipse/group background =========
    ellipse_alpha = 0.18,
    ellipse_linewidth = 0.25,
    ellipse_show_border = TRUE,
    
    # ========= Legend =========
    legend_title_size = 6,
    legend_text_size = 5.5,
    legend_key_width_pt = 10,
    legend_key_height_pt = 10,
    legend_spacing_y_pt = 6,
    legend_spacing_x_pt = 2,
    colorbar_height_mm = 18,
    colorbar_width_mm = 3.5,
    
    # ========= Theme =========
    transparent_bg = TRUE,
    remove_axes = TRUE,
    remove_panel_border = TRUE,
    
    # ========= Layer indices (based on your current object) =========
    edge_layer = 1,
    point_layer = 2,
    ellipse_layer = 3,
    label_layer = 4,
    
    verbose = TRUE
) {
  
  if (!inherits(p, "ggplot")) {
    stop("Input object p must be a ggplot / ggtangle / emapplot object.")
  }
  
  library(ggplot2)
  library(grid)
  
  # ggplot/ggrepel size is not in pt, so convert here
  pt_to_geom_size <- function(pt) {
    (pt * 0.3527778) / .pt
  }
  
  p2 <- p
  
  # -----------------------------
  # Check the number of layers
  # -----------------------------
  n_layers <- length(p2$layers)
  if (verbose) {
    message("The current object has ", n_layers, " layers.")
  }
  
  idx_needed <- c(edge_layer, point_layer, ellipse_layer, label_layer)
  if (any(idx_needed > n_layers)) {
    stop(
      "The specified layer index is out of range. Current layer count = ", n_layers,
      "; please run first:\n",
      "sapply(p$layers, function(x) class(x$geom)[1])\n",
      "to confirm the position of each layer."
    )
  }
  
  # -----------------------------
  # 1) Modify the label layer (geom_text_repel)
  # -----------------------------
  p2$layers[[label_layer]]$aes_params$size <- pt_to_geom_size(label_size)
  p2$layers[[label_layer]]$aes_params$family <- font_family
  p2$layers[[label_layer]]$aes_params$lineheight <- label_lineheight
  
  if (!is.null(label_force)) {
    p2$layers[[label_layer]]$geom_params$force <- label_force
  }
  
  # -----------------------------
  # 2) Modify the point layer
  # Remove the edge color: stroke = 0
  # -----------------------------
  p2$layers[[point_layer]]$aes_params$alpha <- point_alpha
  p2$layers[[point_layer]]$aes_params$stroke <- 0
  
  p2 <- p2 +
    scale_size(
      range = point_size_range,
      name = "Count",
      guide = guide_legend(
        order = 2,
        override.aes = list(
          color = NA,
          fill = "black",
          stroke = 0
        ),
        keyheight = unit(legend_key_height_pt, "pt"),
        keywidth  = unit(legend_key_width_pt, "pt")
      )
    )
  
  # -----------------------------
  # 3) Modify the continuous color (usually p.adjust)
  # -----------------------------
  p2 <- p2 +
    scale_color_viridis_c(
      option = color_palette,
      direction = color_direction,
      name = "p.adjust",
      guide = guide_colorbar(
        order = 1,
        ticks = FALSE,
        frame.colour = NA,
        barheight = unit(colorbar_height_mm, "mm"),
        barwidth  = unit(colorbar_width_mm, "mm")
      )
    )
  
  # -----------------------------
  # 4) Modify the edge layer
  # -----------------------------
  p2$layers[[edge_layer]]$aes_params$linewidth <- edge_linewidth
  p2$layers[[edge_layer]]$aes_params$alpha <- edge_alpha
  
  if (!is.null(edge_color)) {
    p2$layers[[edge_layer]]$aes_params$colour <- edge_color
  }
  
  # -----------------------------
  # 5) Modify the ellipse layer
  # -----------------------------
  p2$layers[[ellipse_layer]]$aes_params$alpha <- ellipse_alpha
  p2$layers[[ellipse_layer]]$aes_params$linewidth <- ellipse_linewidth
  
  if (!ellipse_show_border) {
    p2$layers[[ellipse_layer]]$aes_params$colour <- NA
  }
  
  # -----------------------------
  # 6) Unified theme
  # -----------------------------
  bg_fill <- if (transparent_bg) "transparent" else "white"
  
  p2 <- p2 +
    theme_bw(base_size = 6) +
    theme(
      text = element_text(family = font_family, color = "black"),
      
      panel.background = element_rect(fill = bg_fill, color = NA),
      plot.background  = element_rect(fill = bg_fill, color = NA),
      legend.background = element_rect(fill = bg_fill, color = NA),
      legend.box.background = element_rect(fill = bg_fill, color = NA),
      
      panel.grid = element_blank(),
      panel.border = if (remove_panel_border) element_blank() else element_rect(fill = NA),
      
      axis.text = if (remove_axes) element_blank() else element_text(size = 5, color = "black", family = font_family),
      axis.title = if (remove_axes) element_blank() else element_text(size = 7, color = "black", family = font_family),
      axis.ticks = if (remove_axes) element_blank() else element_line(linewidth = 0.2, color = "black"),
      axis.line = if (remove_axes) element_blank() else element_line(linewidth = 0.2, color = "black"),
      
      legend.position = "right",
      legend.justification = "center",
      legend.title = element_text(size = legend_title_size, color = "black", family = font_family),
      legend.text  = element_text(size = legend_text_size, color = "black", family = font_family),
      legend.key.width  = unit(legend_key_width_pt, "pt"),
      legend.key.height = unit(legend_key_height_pt, "pt"),
      legend.spacing.y = unit(legend_spacing_y_pt, "pt"),
      legend.spacing.x = unit(legend_spacing_x_pt, "pt"),
      
      plot.margin = margin(t = 3, r = 3, b = 2, l = 2)
    )
  
  return(p2)
}
