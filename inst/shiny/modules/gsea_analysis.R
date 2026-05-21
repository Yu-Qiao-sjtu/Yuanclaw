# =====================================================
# GSEA分析模块
# =====================================================

gsea_analysis_server <- function(input, output, session, deg_results) {
  check_runtime_deps <- function(pkgs, feature_name) {
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      showNotification(
        paste0(feature_name, " 依赖缺失: ", paste(missing_pkgs, collapse = ", "), "。请先安装后重试。"),
        type = "error",
        duration = 8
      )
      return(FALSE)
    }
    TRUE
  }

  # =====================================================
  # GSEA 模块
  # =====================================================

  # 🔧 辅助函数：转换core_enrichment为SYMBOL格式
  convert_core_enrichment_to_symbol <- function(df, deg_results) {
    if (!"core_enrichment" %in% colnames(df)) {
      return(df)
    }

    # 获取差异分析数据用于ID转换
    deg_data <- deg_results()
    res <- deg_data$deg_df
    res_clean <- res[!is.na(res$SYMBOL) & !is.na(res$ENTREZID), ]
    res_clean <- res_clean %>%
      group_by(SYMBOL) %>%
      slice(1) %>%
      ungroup()

    # 创建ENTREZID到SYMBOL的映射
    entrez_to_symbol <- setNames(res_clean$SYMBOL, res_clean$ENTREZID)

    # 转换core_enrichment列：始终转换为SYMBOL格式
    df$core_enrichment <- sapply(df$core_enrichment, function(x) {
      if (is.na(x) || !nzchar(x)) return("")

      # 检查是ENTREZID（纯数字）还是SYMBOL
      genes <- unlist(strsplit(x, "/"))

      # 如果是ENTREZID（数字），转换为SYMBOL
      if (all(grepl("^[0-9]+$", genes))) {
        symbols <- entrez_to_symbol[genes]
        symbols <- symbols[!is.na(symbols)]
        paste(symbols, collapse = "/")
      } else {
        # 已经是SYMBOL，直接返回
        x
      }
    }, USE.NAMES = FALSE)

    return(df)
  }

  gsea_results <- eventReactive(input$run_gsea, {
    req(deg_results(), input$gmt_file)
    if (!check_runtime_deps(c("clusterProfiler"), "GSEA分析")) {
      return(NULL)
    }

    showNotification("正在运行 GSEA...", type = "message")

    # 从deg_results中提取差异分析结果
    deg_data <- deg_results()
    res <- deg_data$deg_df

    id_col <- if(input$gsea_id_type == "SYMBOL") "SYMBOL" else "ENTREZID"
    res_clean <- res[!is.na(res[[id_col]]) & !is.na(res$log2FoldChange), ]

    res_clean <- res_clean %>%
      group_by(!!sym(id_col)) %>%
      filter(abs(log2FoldChange) == max(abs(log2FoldChange))) %>%
      ungroup()

    gene_list <- res_clean$log2FoldChange
    names(gene_list) <- res_clean[[id_col]]
    gene_list <- sort(gene_list, decreasing = TRUE)

    # 读取GMT文件
    gmt <- clusterProfiler::read.gmt(input$gmt_file$datapath)

    # 🔧 关键修复：如果使用SYMBOL运行但GMT是ENTREZID，转换GMT为SYMBOL
    if (input$gsea_id_type == "SYMBOL") {
      # 检查GMT中的基因是否是数字（ENTREZID）
      sample_genes <- head(gmt$gene, 100)
      if (all(grepl("^[0-9]+$", sample_genes))) {
        cat("🔄 检测到GMT使用ENTREZID，正在转换为SYMBOL...\n")
        cat(sprintf("📊 GMT文件: %d 行\n", nrow(gmt)))

        tryCatch({
          # 创建ENTREZID到SYMBOL的映射
          entrez_to_symbol <- setNames(res_clean$SYMBOL, res_clean$ENTREZID)
          cat(sprintf("📊 映射关系: %d 个ENTREZID -> %d 个SYMBOL\n",
                     length(entrez_to_symbol), sum(!is.na(entrez_to_symbol))))

          # 转换整个GMT文件
          gmt$gene_symbol <- entrez_to_symbol[as.character(gmt$gene)]

          # 统计转换结果
          n_total <- nrow(gmt)
          n_mapped <- sum(!is.na(gmt$gene_symbol))
          n_unmapped <- sum(is.na(gmt$gene_symbol))

          cat(sprintf("📊 转换结果: %d/%d 成功映射 (%.1f%%), %d 无法映射\n",
                     n_mapped, n_total, n_mapped/n_total*100, n_unmapped))

          if (n_mapped < n_total * 0.5) {
            # 如果超过50%无法映射，可能ID类型选择错误
            msg <- sprintf(
              "⚠️ GMT文件中超过50%%的基因无法映射！\n\n您的GMT文件使用ENTREZID格式，但您选择了SYMBOL。\n\n建议：\n1. 在'GMT中的ID类型'中选择'Entrez ID'\n2. 或者使用SYMBOL格式的GMT文件\n\n当前映射：%.1f%% (%d/%d)",
              n_mapped/n_total*100, n_mapped, n_total
            )
            showNotification(msg, type = "warning", duration = 10)
            cat("⚠️ 映射率过低，建议用户调整ID类型选择\n")
            # 不中断，继续使用部分映射的数据
          }

          # 过滤掉无法映射的基因
          gmt_filtered <- gmt[!is.na(gmt$gene_symbol), ]
          gmt_filtered <- gmt_filtered[, c("term", "gene_symbol")]
          colnames(gmt_filtered) <- c("term", "gene")

          if (nrow(gmt_filtered) > 0) {
            gmt <- gmt_filtered
            cat(sprintf("✅ GMT转换完成: %d 个基因集, %d 个基因\n",
                       length(unique(gmt$term)), nrow(gmt)))
          } else {
            showNotification("❌ GMT转换失败：无法将ENTREZID映射到SYMBOL\n\n请选择正确的ID类型", type = "error")
            return(NULL)
          }

        }, error = function(e) {
          error_msg <- conditionMessage(e)
          cat(sprintf("❌ GMT转换错误: %s\n", error_msg))

          # 提供用户友好的错误信息
          msg <- sprintf(
            "GMT ID类型不匹配！\n\n错误：%s\n\n您的GMT文件使用ENTREZID格式，但您选择了SYMBOL。\n\n请选择'Entrez ID'作为ID类型。",
            substr(error_msg, 1, 100)
          )
          showNotification(msg, type = "error", duration = 15)
          return(NULL)
        })
      }
    }

    tryCatch({
      gsea_res <- clusterProfiler::GSEA(gene_list,
                                        TERM2GENE = gmt,
                                        pvalueCutoff = input$gsea_pvalue,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        verbose = FALSE)

      if (nrow(gsea_res@result) == 0) {
        showNotification("GSEA 未发现显著富集通路", type = "warning")
        return(NULL)
      }
      return(gsea_res)
    }, error = function(e) {
      showNotification(paste("GSEA 运行失败:", e$message), type = "error")
      return(NULL)
    })
  })

  # === 新增：GSEA结果数据处理函数 ===
  gsea_processed_data <- reactive({
    req(gsea_results())

    gsea_obj <- gsea_results()
    df <- gsea_obj@result

    # 根据选择的导出类型处理数据
    if (input$gsea_export_type == "full") {
      # 完整结果
      return(df)
    } else if (input$gsea_export_type == "significant") {
      # 显著结果 (p < 0.05)
      sig_df <- df %>% filter(pvalue < 0.05)
      return(sig_df)
    } else if (input$gsea_export_type == "top50") {
      # Top N 结果
      top_n <- input$gsea_export_topn
      top_df <- df %>% arrange(pvalue) %>% head(top_n)
      return(top_df)
    }
  })

  # === 新增：GSEA结果下载处理器 ===

  # 下载完整结果
  output$download_gsea_full <- downloadHandler(
    filename = function() {
      paste0("GSEA_Full_Results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(gsea_results())
      df <- gsea_results()@result
      # 🔧 转换core_enrichment为SYMBOL格式
      df <- convert_core_enrichment_to_symbol(df, deg_results)
      write.csv(df, file, row.names = FALSE)
    }
  )

  # 下载显著结果
  output$download_gsea_sig <- downloadHandler(
    filename = function() {
      paste0("GSEA_Significant_Results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(gsea_results())
      df <- gsea_results()@result
      sig_df <- df %>% filter(pvalue < 0.05)
      # 🔧 转换core_enrichment为SYMBOL格式
      sig_df <- convert_core_enrichment_to_symbol(sig_df, deg_results)
      write.csv(sig_df, file, row.names = FALSE)
    }
  )

  # 下载Top N结果
  output$download_gsea_top <- downloadHandler(
    filename = function() {
      paste0("GSEA_Top", input$gsea_export_topn, "_Results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(gsea_results())
      df <- gsea_results()@result
      top_n <- input$gsea_export_topn
      top_df <- df %>% arrange(pvalue) %>% head(top_n)
      # 🔧 转换core_enrichment为SYMBOL格式
      top_df <- convert_core_enrichment_to_symbol(top_df, deg_results)
      write.csv(top_df, file, row.names = FALSE)
    }
  )

  # 下载GSEA图为SVG
  output$download_gsea_plot_svg <- downloadHandler(
    filename = function() {
      req(gsea_results())
      selected <- input$gsea_table_rows_selected
      if (length(selected)) {
        pathway_id <- gsea_results()@result$ID[selected]
      } else {
        pathway_id <- gsea_results()@result$ID[1]
      }
      paste0("GSEA_Plot_", pathway_id, "_", Sys.Date(), ".svg")
    },
    content = function(file) {
      # 重新生成plot（避免访问reactive value）
      req(gsea_results())
      gsea_obj <- gsea_results()
      selected <- input$gsea_table_rows_selected

      if (length(selected)) {
        pathway_id <- gsea_obj@result$ID[selected]
      } else {
        pathway_id <- gsea_obj@result$ID[1]
      }

      txt_col <- if(input$theme_toggle) "white" else "black"

      # 生成plot（简化版本，不带addGene以避免复杂性）
      if ("GseaVis" %in% loadedNamespaces()) {
        p <- GseaVis::gseaNb(
          object = gsea_obj,
          geneSetID = pathway_id,
          subPlot = 2,
          termWidth = 35,
          addPval = TRUE
        ) + theme(
          plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5),
          axis.title = element_text(color = txt_col, face = "bold"),
          axis.text = element_text(color = txt_col)
        )
      } else if ("enrichplot" %in% loadedNamespaces()) {
        p <- enrichplot::gseaplot2(gsea_obj, geneSetID = pathway_id) + theme(
          plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5),
          axis.title = element_text(color = txt_col, face = "bold"),
          axis.text = element_text(color = txt_col)
        )
      } else {
        p <- ggplot() + labs(title = "No plotting package available")
      }

      # 保存为SVG
      svg(file, width = 10, height = 6)
      print(p)
      dev.off()
    },
    contentType = "image/svg+xml"
  )

  # 下载山脊图为SVG
  output$download_gsea_ridge_svg <- downloadHandler(
    filename = function() {
      paste0("GSEA_Ridge_Plot_", Sys.Date(), ".svg")
    },
    content = function(file) {
      req(gsea_results())
      gsea_obj <- gsea_results()
      txt_col <- if(input$theme_toggle) "white" else "black"

      # 获取用户设置的通路数
      top_n <- suppressWarnings(as.integer(input$gsea_ridge_pathways))
      if (is.na(top_n) || top_n < 1) {
        top_n <- 10L
      }
      total_pathways <- nrow(gsea_obj@result)
      top_n <- max(1, min(top_n, total_pathways))

      # 生成ridge plot
      if ("enrichplot" %in% loadedNamespaces()) {
        p <- enrichplot::ridgeplot(gsea_obj, showCategory = top_n) +
          labs(title = sprintf("Top %d GSEA Pathways", top_n)) +
          theme(
            plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5),
            axis.title = element_text(color = txt_col, face = "bold"),
            axis.text = element_text(color = txt_col)
          )
      } else {
        p <- ggplot() + labs(title = "enrichplot package required")
      }

      # 保存为SVG
      svg(file, width = 12, height = 8)
      print(p)
      dev.off()
    },
    contentType = "image/svg+xml"
  )

  output$gsea_table <- DT::renderDataTable({
    req(gsea_results())
    df <- gsea_results()@result

    # 调试：检查数据
    cat(sprintf("📊 GSEA结果: %d 行, %d 列\n", nrow(df), ncol(df)))
    cat(sprintf("📊 列名: %s\n", paste(head(colnames(df), 10), collapse=", ")))

    # 检查是否有core_enrichment列
    has_core <- "core_enrichment" %in% colnames(df)
    cat(sprintf("📊 有core_enrichment列: %s\n", has_core))

    # 创建显示用的数据框副本
    df_show <- df

    if (has_core) {
      cat("✅ 找到core_enrichment列，正在转换为SYMBOL...\n")

      # 获取差异分析数据用于ID转换
      deg_data <- deg_results()
      res <- deg_data$deg_df

      # 创建ENTREZID到SYMBOL的映射
      res_clean <- res[!is.na(res$SYMBOL) & !is.na(res$ENTREZID), ]
      entrez_to_symbol <- setNames(res_clean$SYMBOL, res_clean$ENTREZID)

      # 转换core_enrichment列
      df_show$core_enrichment_symbol <- sapply(df_show$core_enrichment, function(core_str) {
        if (is.na(core_str) || core_str == "") {
          return("")
        }

        # 分割基因ID
        gene_ids <- unlist(strsplit(core_str, "/"))

        # 检测是否为ENTREZID（纯数字）
        if (all(grepl("^[0-9]+$", gene_ids))) {
          # 转换为SYMBOL
          gene_symbols <- entrez_to_symbol[gene_ids]
          gene_symbols <- gene_symbols[!is.na(gene_symbols)]
          return(paste(gene_symbols, collapse = "/"))
        } else {
          # 已经是SYMBOL格式
          return(core_str)
        }
      }, USE.NAMES = FALSE)

      cat(sprintf("✅ core_enrichment转换完成\n"))
      cat(sprintf("📊 示例: %s\n", df_show$core_enrichment_symbol[1]))

      # 隐藏原始的core_enrichment列，显示转换后的列
      # 重命名列以保持一致性
      df_show$core_enrichment <- df_show$core_enrichment_symbol
      df_show$core_enrichment_symbol <- NULL
    }

    # 检查df_show是否为空
    if (nrow(df_show) == 0) {
      cat("❌ 错误：df_show为空！\n")
      return(DT::datatable(data.frame(Error="No data")))
    }

    cat(sprintf("📊 准备显示: %d 行, %d 列\n", nrow(df_show), ncol(df_show)))

    # DT配置
    # 简化配置，避免DT错误
    DT::datatable(df_show,
                  selection = 'single',
                  options = list(
                    pageLength = 10,
                    scrollX = TRUE
                  ),
                  rownames = FALSE) %>%
      DT::formatRound(c("enrichmentScore", "NES", "pvalue", "p.adjust"), 4)
  })

  output$gsea_plot <- renderPlot({
    req(gsea_results())
    gsea_obj <- gsea_results()
    selected <- input$gsea_table_rows_selected

    txt_col <- if(input$theme_toggle) "white" else "black"

    if (length(selected)) {
      pathway_id <- gsea_obj@result$ID[selected]
      title_text <- pathway_id
    } else {
      pathway_id <- gsea_obj@result$ID[1]
      title_text <- paste(pathway_id, "(Default: Top 1)")
    }

    # 提取Leading Edge基因（Top N）用于在图上标记
    leading_genes <- NULL
    tryCatch({
      leading_genes_data <- extract_leading_edge_genes()

      if (!is.null(leading_genes_data) && nrow(leading_genes_data) > 0) {
        # 提取基因SYMBOL列表
        leading_genes <- leading_genes_data$gene

        # 根据排序方式确保是Top N
        top_n <- input$gsea_top_genes
        if (length(leading_genes) > top_n) {
          leading_genes <- leading_genes[1:top_n]
        }

        cat(sprintf("在GSEA图上标记 %d 个Leading Edge基因\n", length(leading_genes)))
        cat("基因列表:", paste(head(leading_genes, 10), collapse=", "), ifelse(length(leading_genes)>10, "...", ""), "\n")
      }
    }, error = function(e) {
      cat("提取Leading Edge基因失败:", e$message, "\n")
    })

    # 如果用户输入了自定义基因列表，使用自定义基因
    custom_gene_list <- input$custom_gene_list
    if (!is.null(custom_gene_list) && nzchar(trimws(custom_gene_list))) {
      # 解析自定义基因列表（支持逗号、分号、空格分隔）
      custom_genes <- unlist(strsplit(custom_gene_list, "[,;\\s]+"))
      custom_genes <- trimws(custom_genes)
      custom_genes <- custom_genes[nzchar(custom_genes)]

      if (length(custom_genes) > 0) {
        # 获取ranked gene list用于排序
        deg_data <- deg_results()
        res <- deg_data$deg_df

        # 准备ranked gene list（使用SYMBOL作为名称）
        res_clean <- res[!is.na(res$SYMBOL) & !is.na(res$log2FoldChange), ]
        res_clean <- res_clean %>%
          group_by(SYMBOL) %>%
          filter(abs(log2FoldChange) == max(abs(log2FoldChange))) %>%
          ungroup()

        gene_list <- res_clean$log2FoldChange
        names(gene_list) <- res_clean$SYMBOL
        gene_list <- sort(gene_list, decreasing = TRUE)

        # 获取通路的基因
        gmt <- clusterProfiler::read.gmt(input$gmt_file$datapath)
        pathway_genes_in_gmt <- gmt$gene[gmt$term == pathway_id]

        # 过滤出在通路中的自定义基因
        custom_genes_in_pathway <- custom_genes[custom_genes %in% pathway_genes_in_gmt]

        if (length(custom_genes_in_pathway) == 0) {
          cat("警告：自定义基因都不在通路中\n")
          leading_genes <- NULL
        } else {
          # 创建自定义基因的数据框
          custom_genes_data <- data.frame(
            gene = custom_genes_in_pathway,
            log2FoldChange = gene_list[custom_genes_in_pathway],
            stringsAsFactors = FALSE
          )

          # 移除没有log2FoldChange的基因
          custom_genes_data <- custom_genes_data[!is.na(custom_genes_data$log2FoldChange), ]

          if (nrow(custom_genes_data) == 0) {
            cat("警告：自定义基因都没有log2FoldChange值\n")
            leading_genes <- NULL
          } else {
            # 根据用户选择的方式排序自定义基因
            if (input$gsea_gene_order == "abs_logFC") {
              custom_genes_data <- custom_genes_data[order(abs(custom_genes_data$log2FoldChange), decreasing = TRUE), ]
            } else if (input$gsea_gene_order == "logFC") {
              custom_genes_data <- custom_genes_data[order(custom_genes_data$log2FoldChange, decreasing = TRUE), ]
            } else if (input$gsea_gene_order == "rank") {
              # 按在ranked list中的位置排序
              custom_genes_data$rank <- match(custom_genes_data$gene, names(gene_list))
              custom_genes_data <- custom_genes_data[order(custom_genes_data$rank), ]
            }

            # 提取排序后的基因列表
            leading_genes <- custom_genes_data$gene
            cat(sprintf("使用自定义基因列表（按%s排序）: %d 个基因\n",
                       input$gsea_gene_order, length(leading_genes)))
            cat("排序后的基因:", paste(head(leading_genes, 10), collapse=", "),
                ifelse(length(leading_genes)>10, "...", ""), "\n")
          }
        }
      }
    }

    if ("GseaVis" %in% loadedNamespaces()) {
      # ====================================================
      # 🔥 核心修复：使用 gseaNb() 的 addGene 参数
      # ====================================================
      # 原因：gseaNb() 返回 aplot 组合对象（多子图）
      # 在返回的对象上用 + 添加图层，会加到最后一个子图（下面板）
      # addGene 参数在内部直接操作 ES 曲线子图，才能正确显示
      # ====================================================

      # leading_genes 已在上方准备好（SYMBOL格式）
      # 包含：自定义基因列表 > extract_leading_edge_genes() 结果
      # 现在将 SYMBOL 列表转换为 addGene 所需的正确 ID 类型
      genes_to_label <- NULL

      if (!is.null(leading_genes) && length(leading_genes) > 0) {
        if (input$gsea_id_type == "ENTREZID") {
          # GSEA 用 ENTREZID 运行 → addGene 需要 ENTREZID
          deg_data <- deg_results()
          res_clean <- deg_data$deg_df
          res_clean <- res_clean[!is.na(res_clean$SYMBOL) & !is.na(res_clean$ENTREZID), ]
          symbol_to_entrez <- setNames(res_clean$ENTREZID, res_clean$SYMBOL)
          entrez_ids <- symbol_to_entrez[leading_genes]
          entrez_ids <- as.character(entrez_ids[!is.na(entrez_ids)])
          if (length(entrez_ids) > 0) {
            genes_to_label <- entrez_ids
            cat(sprintf("📝 addGene 使用 ENTREZID: %d 个基因\n", length(genes_to_label)))
          }
        } else {
          # GSEA 用 SYMBOL 运行 → addGene 直接用 SYMBOL
          genes_to_label <- leading_genes
          cat(sprintf("📝 addGene 使用 SYMBOL: %d 个基因\n", length(genes_to_label)))
        }
      }

      # 构建 gseaNb 参数，通过 addGene 在 ES 曲线上标注基因
      plot_args <- list(
        object    = gsea_obj,
        geneSetID = pathway_id,
        subPlot   = 2,
        termWidth = 35,
        addPval   = TRUE,
        pvalX     = input$gsea_stats_x,
        pvalY     = input$gsea_stats_y
      )

      if (!is.null(genes_to_label) && length(genes_to_label) > 0) {
        plot_args$addGene     <- genes_to_label
        plot_args$geneCol     <- if (input$theme_toggle) "#00FF00" else "#CC0000"
        plot_args$geneSize    <- 3.5
        plot_args$rmSegment   <- FALSE   # 显示从曲线到标签的连接线
        plot_args$segCol      <- if (input$theme_toggle) "#00FF00" else "#CC0000"
        plot_args$force       <- 20
        cat(sprintf("✅ 通过 addGene 参数在 ES 曲线上标注 %d 个基因\n", length(genes_to_label)))
      }

      p <- tryCatch({
        do.call(GseaVis::gseaNb, plot_args)
      }, error = function(e) {
        cat("❌ gseaNb 带 addGene 失败:", e$message, "\n", "尝试不带 addGene...\n")
        # 降级：不带 addGene 重试
        plot_args$addGene   <- NULL
        plot_args$geneCol   <- NULL
        plot_args$geneSize  <- NULL
        plot_args$rmSegment <- NULL
        plot_args$segCol    <- NULL
        plot_args$force     <- NULL
        do.call(GseaVis::gseaNb, plot_args)
      })

      p <- p + theme(
        plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5),
        axis.title = element_text(color = txt_col, face = "bold"),
        axis.text = element_text(color = txt_col),
        legend.text = element_text(color = txt_col),
        legend.title = element_text(color = txt_col),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)
      )
    } else {
      if ("enrichplot" %in% loadedNamespaces()) {
        p <- enrichplot::gseaplot2(gsea_obj, geneSetID = pathway_id, title = title_text) +
          theme(
            plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5),
            axis.title = element_text(color = txt_col, face = "bold"),
            axis.text = element_text(color = txt_col),
            legend.text = element_text(color = txt_col),
            legend.title = element_text(color = txt_col),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA)
          )
      } else {
        p <- ggplot() + labs(title = "缺少 GseaVis 或 enrichplot 包，无法绘图")
      }
    }
    print(p)
  })

  # === 新增：Leading Edge 基因提取函数 ===
  extract_leading_edge_genes <- reactive({
    req(gsea_results())

    gsea_obj <- gsea_results()
    selected <- input$gsea_table_rows_selected

    if (!length(selected)) {
      # 如果没有选择，使用第一个通路
      selected <- 1
    }

    pathway_id <- gsea_obj@result$ID[selected]

    # 提取基因列表（从原始差异分析数据）
    deg_data <- deg_results()
    res <- deg_data$deg_df

    # 根据GMT文件中的ID类型选择合适的列
    # GMT文件中的ID类型由用户在gsea_id_type中指定
    id_col_in_gmt <- if(input$gsea_id_type == "SYMBOL") "SYMBOL" else "ENTREZID"

    # 确保使用SYMBOL用于最终显示（GseaVis需要SYMBOL）
    res_clean <- res[!is.na(res$SYMBOL) & !is.na(res$log2FoldChange), ]

    # 去重并排序
    res_clean <- res_clean %>%
      group_by(SYMBOL) %>%
      filter(abs(log2FoldChange) == max(abs(log2FoldChange))) %>%
      ungroup()

    # 创建排序列表（使用SYMBOL作为名称）
    gene_list <- res_clean$log2FoldChange
    names(gene_list) <- res_clean$SYMBOL
    gene_list <- sort(gene_list, decreasing = TRUE)

    # 读取GMT文件以获取该通路的基因集
    gmt <- clusterProfiler::read.gmt(input$gmt_file$datapath)

    # 获取选中通路的基因
    pathway_genes <- gmt$gene[gmt$term == pathway_id]

    if (length(pathway_genes) == 0) {
      return(NULL)
    }

    # 如果GMT文件使用ENTREZID，需要转换为SYMBOL
    if (id_col_in_gmt == "ENTREZID") {
      # 创建ENTREZID到SYMBOL的映射
      entrez_to_symbol <- setNames(res_clean$SYMBOL, res_clean$ENTREZID)
      pathway_genes_symbol <- entrez_to_symbol[pathway_genes]
      # 移除NA值（没有映射到的基因）
      pathway_genes_symbol <- pathway_genes_symbol[!is.na(pathway_genes_symbol)]
      pathway_genes <- pathway_genes_symbol
    }

    # 🔥 新增：提取真正的Leading Edge基因
    # 首先尝试从GSEA结果的core_enrichment字段提取
    if (input$gsea_gene_order == "leading_edge") {
      tryCatch({
        # 从GSEA结果中提取core_enrichment基因
        core_enrichment_str <- gsea_obj@result$core_enrichment[selected]

        cat(sprintf("🔍 提取Leading Edge基因，selected=%d, pathway_id=%s\n", selected, pathway_id))
        cat(sprintf("🔍 core_enrichment内容: %s\n", substring(core_enrichment_str, 1, 200)))

        if (!is.na(core_enrichment_str) && nzchar(core_enrichment_str)) {
          # core_enrichment字段是用"/"分隔的基因列表
          # ⚠️ 注意：core_enrichment中的ID类型与gene_list的names类型相同
          le_genes_raw <- unlist(strsplit(core_enrichment_str, "/"))
          cat(sprintf("🔍 原始Leading Edge基因数量: %d (ID类型: %s)\n", length(le_genes_raw), input$gsea_id_type))

          # 🔧 关键修复：始终检测并转换为SYMBOL格式
          le_genes_symbol <- le_genes_raw  # 初始值

          # 检测是否为ENTREZID（纯数字）并转换为SYMBOL
          if (all(grepl("^[0-9]+$", le_genes_raw))) {
            cat("🔄 检测到ENTREZID格式，正在转换为SYMBOL...\n")
            # 创建ENTREZID到SYMBOL的映射
            entrez_to_symbol <- setNames(res_clean$SYMBOL, res_clean$ENTREZID)
            le_genes_symbol <- entrez_to_symbol[le_genes_raw]
            # 移除NA值（没有映射到的基因）
            le_genes_symbol <- le_genes_symbol[!is.na(le_genes_symbol)]
            cat(sprintf("✅ 转换后SYMBOL基因数量: %d\n", length(le_genes_symbol)))
          } else {
            cat("✅ 已经是SYMBOL格式\n")
          }

          if (length(le_genes_symbol) == 0) {
            cat("⚠️ Leading Edge基因ID转换失败，尝试使用常规方式\n")
          } else {
            # 获取这些基因的log2FoldChange
            # ⚠️ 注意：gene_list的names类型可能与le_genes不同
            # 需要使用SYMBOL作为key来查找log2FoldChange
            gene_list_symbol <- res_clean$log2FoldChange
            names(gene_list_symbol) <- res_clean$SYMBOL
            gene_list_symbol <- sort(gene_list_symbol, decreasing = TRUE)

            pathway_data <- data.frame(
              gene = le_genes_symbol,
              log2FoldChange = gene_list_symbol[le_genes_symbol],
              stringsAsFactors = FALSE
            )

            # 移除没有log2FoldChange的基因
            pathway_data <- pathway_data[!is.na(pathway_data$log2FoldChange), ]

            if (nrow(pathway_data) > 0) {
              # Leading Edge基因按在ranked list中的位置排序
              pathway_data$rank <- match(pathway_data$gene, names(gene_list_symbol))
              pathway_data <- pathway_data[order(pathway_data$rank), ]

              # 选择Top N基因
              top_n <- min(input$gsea_top_genes, nrow(pathway_data))
              pathway_data_top <- pathway_data[1:top_n, ]

              # 添加排名信息
              pathway_data_top$rank_label <- paste0("#", 1:top_n)

              cat(sprintf("✅ 提取了 %d 个真正的Leading Edge基因 (ID类型: SYMBOL)\n",
                         nrow(pathway_data_top)))
              cat("✅ Leading Edge基因示例:", paste(head(pathway_data_top$gene, 5), collapse=", "), "\n")
              return(pathway_data_top)
            }
          }
        }
      }, error = function(e) {
        cat("⚠️ 提取Leading Edge基因失败，使用常规方式:", e$message, "\n")
      })
    }

    # 如果不是leading_edge模式，或者提取失败，使用原有逻辑
    # 获取基因集中的基因及其log2FoldChange
    pathway_data <- data.frame(
      gene = pathway_genes,
      log2FoldChange = gene_list[pathway_genes],
      stringsAsFactors = FALSE
    )

    # 移除没有log2FoldChange的基因
    pathway_data <- pathway_data[!is.na(pathway_data$log2FoldChange), ]

    if (nrow(pathway_data) == 0) {
      return(NULL)
    }

    # 根据用户选择的方式排序
    if (input$gsea_gene_order == "abs_logFC") {
      pathway_data <- pathway_data[order(abs(pathway_data$log2FoldChange), decreasing = TRUE), ]
    } else if (input$gsea_gene_order == "logFC") {
      pathway_data <- pathway_data[order(pathway_data$log2FoldChange, decreasing = TRUE), ]
    } else if (input$gsea_gene_order == "rank") {
      # 按在ranked list中的位置排序
      pathway_data$rank <- match(pathway_data$gene, names(gene_list))
      pathway_data <- pathway_data[order(pathway_data$rank), ]
    }

    # 选择Top N基因
    top_n <- min(input$gsea_top_genes, nrow(pathway_data))
    pathway_data_top <- pathway_data[1:top_n, ]

    # 添加排名信息
    pathway_data_top$rank_label <- paste0("#", 1:top_n)

    return(pathway_data_top)
  })

  # === 新增：GSEA 山脊图可视化 ===
  output$gsea_ridge_plot <- renderPlot({
    req(gsea_results())
    req(input$show_gsea_ridge)

    gsea_obj <- gsea_results()
    txt_col <- if(input$theme_toggle) "white" else "black"

    # 🔧 安全转换整数，设置默认值
    top_n <- suppressWarnings(as.integer(input$gsea_ridge_pathways))
    if (is.na(top_n) || top_n < 1) {
      top_n <- 10L  # 默认显示10个通路
    }

    cat(sprintf("🎨 用户请求显示 %d 个通路的山脊图\n", top_n))

    # 使用enrichplot的ridgeplot
    if ("enrichplot" %in% loadedNamespaces()) {
      tryCatch({
        # 🔧 修复：确保top_n是有效整数
        total_pathways <- nrow(gsea_obj@result)
        top_n <- max(1, min(top_n, total_pathways))

        cat(sprintf("📊 总共有 %d 个通路，将显示前 %d 个\n", total_pathways, top_n))

        # 🔧 关键修复：使用showCategory参数限制显示数量
        # 根据enrichplot文档，showCategory接受数字或向量
        p <- enrichplot::ridgeplot(gsea_obj, showCategory = top_n) +
          labs(title = sprintf("Top %d GSEA Pathways (Total: %d)", top_n, total_pathways)) +
          theme(
            plot.title = element_text(color = txt_col, face = "bold", hjust = 0.5, size = 14),
            axis.title = element_text(color = txt_col, face = "bold"),
            axis.text = element_text(color = txt_col),
            legend.text = element_text(color = txt_col),
            legend.title = element_text(color = txt_col),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA)
          )

        print(p)
        cat("✅ 山脊图生成成功\n")
        return(NULL)
      }, error = function(e) {
        cat("❌ ridgeplot错误:", e$message, "\n")
        cat("错误详情:", conditionMessage(e), "\n")

        # 显示友好的错误消息
        showNotification(paste("ridgeplot绘图失败:", e$message), type = "warning")

        # 返回一个简单的错误图
        p <- ggplot() +
          labs(title = "山脊图生成失败") +
          geom_text(aes(x = 0.5, y = 0.5, label = "请检查GSEA结果\n或减少显示的通路数"), size = 5) +
          theme_void() +
          theme(
            plot.title = element_text(color = txt_col, hjust = 0.5, face = "bold"),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA)
          )
        print(p)
      })
    } else {
      # 如果enrichplot不可用，显示提示
      p <- ggplot() +
        labs(title = "无法生成山脊图 - 需要enrichplot包") +
        theme(
          plot.title = element_text(color = txt_col, hjust = 0.5, face = "bold"),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA)
        )
      print(p)
    }
  })


}
