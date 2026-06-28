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
  
  # -------- 1. 输入检查 --------
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
  
  # -------- 2. 标准化 --------
  mat_scaled <- mat
  
  if (scale == "row") {
    mat_scaled <- t(scale(t(mat_scaled)))
  } else if (scale == "column") {
    mat_scaled <- scale(mat_scaled)
  }
  
  mat_scaled[is.nan(mat_scaled)] <- 0
  
  # -------- 3. 聚类 --------
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
  
  # -------- 4. 转成长表 --------
  df_long <- as.data.frame(mat_scaled) %>%
    rownames_to_column("Feature") %>%
    pivot_longer(cols = -Feature, names_to = "Sample", values_to = "Value")
  
  # y 轴要反过来，保证第一行在最上面
  df_long$Feature <- factor(df_long$Feature, levels = rev(row_order))
  df_long$Sample  <- factor(df_long$Sample, levels = col_order)
  
  # 对称色标，更像 pheatmap(scale="row")
  max_abs <- max(abs(df_long$Value), na.rm = TRUE)
  fill_limits <- c(-max_abs, max_abs)
  
  # -------- 5. heatmap主体 --------
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
  
  # -------- 6. 行树 --------
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
  
  # -------- 7. 列树 --------
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
  
  # -------- 8. 拼图 --------
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
  
  # 按 Gene Ratio 排 y 轴顺序
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
      # 关键：P.adjust 在上面
      color = guide_colorbar(
        order = 1,
        barheight = unit(18, "mm"),
        barwidth  = unit(3, "mm"),
        ticks = FALSE,
        frame.colour = NA,
        title.position = "top",
        title.hjust = 0.5
      ),
      
      # 关键：Count 在下面
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
    # ========= 字体大小（最常改）=========
    tip_label_size = 1.9,        # tip文字大小；约 1.7~2.3 对应 5~7 pt
    clade_label_size = 1.9,      # group/clade文字大小
    legend_title_size = 6,       # 图例标题大小（pt）
    legend_text_size = 5.5,      # 图例文字大小（pt）
    
    # ========= 点和颜色 =========
    point_size_range = c(1.2, 3.2),   # 圆点映射范围
    padj_palette = "C",               # viridis连续色板，C常见为蓝绿黄
    padj_direction = -1,              # -1表示p.adjust越小颜色越醒目
    
    # ========= 高亮背景 =========
    highlight_alpha = 0.35,           # group背景透明度
    
    # ========= 图例色条 =========
    colorbar_height_mm = 18,
    colorbar_width_mm = 3.5,
    
    # ========= 图例整体 =========
    legend_key_width_pt = 10,
    legend_key_height_pt = 10,
    legend_spacing_y_pt = 6,
    legend_spacing_x_pt = 2,
    
    # ========= 背景 / 坐标轴 =========
    transparent_bg = TRUE,
    remove_panel_border = TRUE,
    remove_axes = TRUE,
    
    # ========= 如果想压缩树宽，可打开 =========
    compress_tree = FALSE,
    tree_width_scale = 5,
    
    # ========= 如果想立刻用ggview显示 =========
    show_canvas = TRUE,
    canvas_width_mm = 180,
    canvas_height_mm = 120,
    canvas_dpi = 600,
    canvas_bg = "transparent",
    
    # ========= 保存时转成静态对象（关键）=========
    make_static = TRUE,
    
    # ========= layer索引（基于你当前对象结构）=========
    # 你当前 str() 看起来是：
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
  # 依赖包
  # -----------------------------
  if (!inherits(p, "ggplot")) {
    stop("输入对象 p 必须是 ggplot / ggtree / treeplot 对象。")
  }
  
  library(ggplot2)
  library(grid)
  
  if (show_canvas && !requireNamespace("ggview", quietly = TRUE)) {
    stop("需要安装 ggview 包：install.packages('ggview')")
  }
  
  if (make_static && !requireNamespace("ggplotify", quietly = TRUE)) {
    stop("需要安装 ggplotify 包：install.packages('ggplotify')")
  }
  
  # -----------------------------
  # 复制对象，避免直接改原对象
  # -----------------------------
  p2 <- p
  
  # -----------------------------
  # 检查 layers 是否足够
  # -----------------------------
  n_layers <- length(p2$layers)
  
  if (verbose) {
    message("当前对象共有 ", n_layers, " 个 layers。")
  }
  
  idx_needed <- c(tip_label_layer, highlight_layer, point_layer, clade_label_layer)
  if (any(idx_needed > n_layers)) {
    stop(
      "指定的 layer 索引超出范围。当前 layer 数量 = ", n_layers,
      "；请先运行：\n",
      "sapply(p$layers, function(x) class(x$geom)[1])\n",
      "确认各层位置。"
    )
  }
  
  # -----------------------------
  # 1) 修改 tip label 字体大小
  # -----------------------------
  p2$layers[[tip_label_layer]]$aes_params$size <- tip_label_size
  
  # -----------------------------
  # 2) 修改 clade/group label 字体大小
  # -----------------------------
  p2$layers[[clade_label_layer]]$aes_params$size <- clade_label_size
  
  # -----------------------------
  # 3) 修改高亮背景透明度
  # -----------------------------
  p2$layers[[highlight_layer]]$aes_params$alpha <- highlight_alpha
  
  # -----------------------------
  # 4) 缩放点大小
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
  # 5) 覆盖 p.adjust 颜色
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
  # 6) 可选：压缩树宽
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
  # 7) 主题统一
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
  # 8) 显示 / 返回
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
  
  # 关键：如果要保存，先冻结成静态 grob，再转回普通 ggplot
  # 这样 ggsave 不会因为 interactive geoms 而导出空白
  if (make_static) {
    g <- grid::grid.grabExpr(print(p2))
    p2 <- ggplotify::as.ggplot(g)
  }
  
  return(p2)
}




#emap plot

modify_emapplot_pub <- function(
    p,
    
    # ========= 文字 =========
    label_size = 6,                 # 圆上标签文字大小（pt）
    label_force = NULL,             # ggrepel 排斥强度；NULL = 不改
    label_lineheight = 0.9,         # 标签行距
    font_family = "sans",           # 全部文字字体
    
    # ========= 点 =========
    point_size_range = c(2, 7),     # 点大小映射范围
    point_alpha = 1,                # 点透明度
    
    # ========= 颜色 =========
    color_palette = "C",            # 连续色板（通常对应 p.adjust）
    color_direction = -1,           # -1 = p.adjust 越小越醒目
    
    # ========= 连线 =========
    edge_linewidth = 0.2,
    edge_alpha = 0.25,
    edge_color = "grey30",          # 连线颜色；NULL = 不改
    
    # ========= ellipse/group背景 =========
    ellipse_alpha = 0.18,
    ellipse_linewidth = 0.25,
    ellipse_show_border = TRUE,
    
    # ========= 图例 =========
    legend_title_size = 6,
    legend_text_size = 5.5,
    legend_key_width_pt = 10,
    legend_key_height_pt = 10,
    legend_spacing_y_pt = 6,
    legend_spacing_x_pt = 2,
    colorbar_height_mm = 18,
    colorbar_width_mm = 3.5,
    
    # ========= 主题 =========
    transparent_bg = TRUE,
    remove_axes = TRUE,
    remove_panel_border = TRUE,
    
    # ========= layer索引（基于你当前对象）=========
    edge_layer = 1,
    point_layer = 2,
    ellipse_layer = 3,
    label_layer = 4,
    
    verbose = TRUE
) {
  
  if (!inherits(p, "ggplot")) {
    stop("输入对象 p 必须是 ggplot / ggtangle / emapplot 对象。")
  }
  
  library(ggplot2)
  library(grid)
  
  # ggplot/ggrepel 的 size 不是 pt，这里做换算
  pt_to_geom_size <- function(pt) {
    (pt * 0.3527778) / .pt
  }
  
  p2 <- p
  
  # -----------------------------
  # 检查 layer 数量
  # -----------------------------
  n_layers <- length(p2$layers)
  if (verbose) {
    message("当前对象共有 ", n_layers, " 个 layers。")
  }
  
  idx_needed <- c(edge_layer, point_layer, ellipse_layer, label_layer)
  if (any(idx_needed > n_layers)) {
    stop(
      "指定的 layer 索引超出范围。当前 layer 数量 = ", n_layers,
      "；请先运行：\n",
      "sapply(p$layers, function(x) class(x$geom)[1])\n",
      "确认各层位置。"
    )
  }
  
  # -----------------------------
  # 1) 改标签层（geom_text_repel）
  # -----------------------------
  p2$layers[[label_layer]]$aes_params$size <- pt_to_geom_size(label_size)
  p2$layers[[label_layer]]$aes_params$family <- font_family
  p2$layers[[label_layer]]$aes_params$lineheight <- label_lineheight
  
  if (!is.null(label_force)) {
    p2$layers[[label_layer]]$geom_params$force <- label_force
  }
  
  # -----------------------------
  # 2) 改点层
  # 去掉边缘颜色：stroke = 0
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
  # 3) 改连续颜色（通常是 p.adjust）
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
  # 4) 改连线层
  # -----------------------------
  p2$layers[[edge_layer]]$aes_params$linewidth <- edge_linewidth
  p2$layers[[edge_layer]]$aes_params$alpha <- edge_alpha
  
  if (!is.null(edge_color)) {
    p2$layers[[edge_layer]]$aes_params$colour <- edge_color
  }
  
  # -----------------------------
  # 5) 改 ellipse 层
  # -----------------------------
  p2$layers[[ellipse_layer]]$aes_params$alpha <- ellipse_alpha
  p2$layers[[ellipse_layer]]$aes_params$linewidth <- ellipse_linewidth
  
  if (!ellipse_show_border) {
    p2$layers[[ellipse_layer]]$aes_params$colour <- NA
  }
  
  # -----------------------------
  # 6) 统一主题
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