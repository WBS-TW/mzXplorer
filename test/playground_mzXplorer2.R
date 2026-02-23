# =========================================================
# mzXplorer with HS detection + CCS-aware homologue search
# =========================================================

library(shiny)
library(shinyjqui)
library(shinythemes)
library(DT)
library(dplyr)
library(plotly)
library(crosstalk)
library(enviPat); data("isotopes")
library(vroom)
library(igraph)

source("playground_utils2.R")

###########  UI  ###########

ui <- shiny::navbarPage(
        "mzXplorer: interactive mass defect plots",
        theme = shinythemes::shinytheme('spacelab'),
        
        shiny::tabPanel(
                "DataAnalysis",
                shiny::fluidPage(
                        shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        shiny::fileInput(
                                                'file1', 'Choose CSV File',
                                                accept = c(
                                                        'text/csv',
                                                        'text/comma-separated-values,text/plain',
                                                        '.csv'
                                                )
                                        ),
                                        
                                        shiny::fluidRow(
                                                shiny::column(
                                                        12,
                                                        shiny::textInput("cus1", "MD formula 1", value = "CH2,O")
                                                )
                                        ),
                                        shiny::fluidRow(
                                                shiny::column(
                                                        12,
                                                        shiny::textInput("cus2", "MD formula 2", value = "Cl-H")
                                                )
                                        ),
                                        shiny::actionButton('go', 'Plot', width = "100%"),
                                        shiny::tags$br(),
                                        
                                        shiny::radioButtons(
                                                "rounding", "Rounding",
                                                choices  = c("round", "ceiling", "floor"),
                                                selected = "round",
                                                inline   = TRUE
                                        ),
                                        shiny::checkboxInput('ins', 'Show intensity/variable as size', FALSE),
                                        shiny::checkboxInput("show_leg", "Show plot legends", TRUE),
                                        
                                        # ----------------------------
                                        # Homologue search controls
                                        # ----------------------------
                                        shiny::tags$hr(), shiny::h4("Homologue search"),
                                        shiny::checkboxInput(
                                                "run_homol",
                                                "Detect homologues in filtered data",
                                                FALSE
                                        ),
                                        shiny::textInput("homol_unit", "Repeating unit (formula)", value = "CH2"),
                                        shiny::numericInput(
                                                "homol_ppm", "m/z tolerance (ppm)",
                                                value = 5, min = 1, step = 1
                                        ),
                                        
                                        # --- CCS controls (dynamic) ---
                                        shiny::uiOutput("ccs_block"),
                                        
                                        shiny::numericInput(
                                                "homol_minlen", "Minimum series length",
                                                value = 4, min = 3, step = 1
                                        ),
                                        shiny::numericInput(
                                                "homol_rttol",
                                                "Per-step RT tolerance (same units as rt); 0 = off",
                                                value = 0, min = 0, step = 0.1
                                        ),
                                        shiny::selectInput(
                                                "homol_rttrend", "RT / CCS trend",
                                                choices  = c("increasing", "decreasing", "any"),
                                                selected = "increasing"
                                        ),
                                        shiny::checkboxInput("homol_allow_gaps", "Allow single gaps (k = 2)", FALSE),
                                        shiny::numericInput(
                                                "homol_R2", "Minimum spline R² (0 = off)",
                                                value = 0.98, min = 0, max = 1, step = 0.01
                                        ),
                                        
                                        shiny::uiOutput("slide1"),
                                        shiny::uiOutput("slide2"),
                                        shiny::uiOutput("slide3"),
                                        width = 3
                                ),
                                
                                shiny::mainPanel(
                                        # Plot controls
                                        shiny::uiOutput("plotctr"),
                                        
                                        # Two plots side-by-side
                                        shiny::uiOutput("plot"),
                                        
                                        # Tables and barplots
                                        DT::DTOutput("table_selected"),
                                        shiny::tags$br(),
                                        shiny::uiOutput("intensity_ctr"),
                                        plotly::plotlyOutput("barplot"),
                                        
                                        shiny::fluidRow(
                                                shiny::column(3, shiny::downloadButton("x3", "Export Data")),
                                                shiny::column(3, shiny::actionButton("clear_homol_sel", "Clear selection"))
                                        ),
                                        
                                        shiny::tags$br(),
                                        DT::DTOutput("homol_table")
                                )
                        )
                )
        ),
        
        shiny::tabPanel(
                "Instructions",
                shiny::fluidPage(
                        shiny::titlePanel("How to use mzXplorer"),
                        shiny::br(),
                        shiny::includeMarkdown("./playground_instructions.md")
                )
        )
)

###########  Server  ###########

server <- function(input, output, session) {
        
        # Shared state
        vals <- reactiveValues(
                d = NULL, m = NULL, summ = NULL,
                series_keys = NULL, series_sel = NULL, series_cols = NULL
        )
        
        # Simple palette
        series_palette <- function(n) {
                base <- c(
                        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                        "#bcbd22", "#17becf", "#a6cee3", "#fb9a99",
                        "#33a02c", "#b15928", "#cab2d6", "#fdbf6f"
                )
                rep(base, length.out = n)
        }
        
        # ----------------------------------
        # Data + MD features
        # ----------------------------------
        MD_data <- reactive({
                req(input$file1)
                df <- vroom::vroom(input$file1$datapath, show_col_types = FALSE)
                
                # Normalize names: lowercase + trim whitespace
                nms <- trimws(tolower(names(df)))
                names(df) <- nms
                
                if (!"mz" %in% names(df) && "m/z" %in% names(df)) df$mz <- df[["m/z"]]
                if (!"intensity" %in% names(df) && "dummy" %in% names(df)) df$intensity <- df$dummy
                if (!"rt" %in% names(df) && "retention-time" %in% names(df)) df$rt <- df[["retention-time"]]
                
                stopifnot(all(c("mz","intensity","rt") %in% names(df)))
                
                # Stable ID per row to allow unambiguous joins (even if mz/rt duplicated)
                df$id <- seq_len(nrow(df))
                
                df$RMD <- round((round(df$mz) - df$mz) / df$mz * 1e6)
                df$OMD <- round((round(df$mz) - df$mz) * 1e3)
                
                mdh1 <- getmdh(df$mz, cus = input$cus1, method = input$rounding)
                mdh2 <- getmdh(df$mz, cus = input$cus2, method = input$rounding)
                name1 <- paste0("Formula1_", colnames(mdh1))
                name2 <- paste0("Formula2_", colnames(mdh2))
                mdh  <- cbind(mdh1[, -1], mdh2[, -1])
                colnames(mdh) <- c(name1[-1], name2[-1])
                
                cbind(df, mdh)
        })
        
        # Does current data have CCS?
        has_ccs <- reactive({
                "ccs" %in% names(MD_data())
        })
        
        # --------------------------
        # CCS controls block
        # --------------------------
        output$ccs_block <- renderUI({
                if (!has_ccs()) {
                        # No CCS: show information only
                        tagList(
                                tags$h5("CCS-based rules"),
                                tags$small(
                                        style = "color:#777;",
                                        "No 'ccs' column detected in the uploaded data; homologue detection will use RT only."
                                )
                        )
                } else {
                        # CCS available: show full CCS controls
                        tagList(
                                checkboxInput(
                                        "use_ccs_toggle",
                                        "Enable CCS-based homologue rules",
                                        value = FALSE
                                ),
                                uiOutput("ccs_mode_ui"),
                                uiOutput("ccs_tol_ui")
                        )
                }
        })
        
        # CCS mode + tolerance (only when toggle ON and CCS present)
        output$ccs_mode_ui <- renderUI({
                req(has_ccs())
                req(isTRUE(input$use_ccs_toggle))
                
                selectInput(
                        "ccs_mode",
                        "Use which dimension for homologue trend?",
                        choices = c(
                                "RT only"       = "rt",
                                "CCS only"      = "ccs",
                                "Both RT + CCS" = "both"
                        ),
                        selected = "ccs"
                )
        })
        
        output$ccs_tol_ui <- renderUI({
                req(has_ccs())
                req(isTRUE(input$use_ccs_toggle))
                
                numericInput(
                        "homol_ccstol",
                        "Per-step CCS tolerance (same units as ccs); 0 = off",
                        value = 0, min = 0, step = 0.1
                )
        })
        
        # -----------------
        # Sliders
        # -----------------
        output$slide1 <- renderUI({
                df <- MD_data()
                sliderInput(
                        "slide1", "Intensity range filter",
                        min   = min(df$intensity, na.rm = TRUE),
                        max   = max(df$intensity, na.rm = TRUE),
                        value = c(min(df$intensity, na.rm = TRUE),
                                  max(df$intensity, na.rm = TRUE))
                )
        })
        
        output$slide2 <- renderUI({
                df <- MD_data()
                sliderInput(
                        "slide2", "m/z range",
                        min   = min(df$mz, na.rm = TRUE),
                        max   = max(df$mz, na.rm = TRUE),
                        value = c(min(df$mz, na.rm = TRUE),
                                  max(df$mz, na.rm = TRUE))
                )
        })
        
        output$slide3 <- renderUI({
                df <- MD_data()
                sliderInput(
                        "slide3", "retention time range",
                        min   = min(df$rt, na.rm = TRUE),
                        max   = max(df$rt, na.rm = TRUE),
                        value = c(min(df$rt, na.rm = TRUE),
                                  max(df$rt, na.rm = TRUE))
                )
        })
        
        # -------------
        # Plot layout
        # -------------
        output$plot <- renderUI({
                fluidRow(
                        column(width = 6, plotlyOutput("DTPlot1")),
                        column(width = 6, plotlyOutput("DTPlot2"))
                )
        })
        
        output$plotctr <- renderUI({
                df <- MD_data()
                fluidRow(
                        h4("Plot controls"),
                        column(
                                6,
                                selectInput('xvar1', 'X variable for Plot 1',
                                            choices = names(df), selected = "rt")
                        ),
                        column(
                                6,
                                selectInput('xvar2', 'X variable for Plot 2',
                                            choices = names(df), selected = "RMD")
                        ),
                        column(
                                6,
                                selectInput('yvar1', 'Y variable for Plot 1',
                                            choices = names(df), selected = "mz")
                        ),
                        column(
                                6,
                                selectInput('yvar2', 'Y variable for Plot 2',
                                            choices = names(df), selected = "mz")
                        )
                )
        })
        
        output$intensity_ctr <- renderUI({
                df <- MD_data()
                numeric_cols <- names(dplyr::select(df, where(is.numeric)))
                
                selectInput(
                        'selectintensity',
                        'Variable for intensity (for point size & barplot)',
                        choices  = numeric_cols,
                        selected = if ('intensity' %in% numeric_cols)
                                'intensity' else numeric_cols[1]
                )
        })
        
        # --- Shared plotting function
        make_mdplot <- function(shared_data, xvar, yvar, xlabel, ylabel,
                                intensity_enabled, intensity_col,
                                series_sel, series_cols,
                                show_legend = TRUE,
                                source_id = "plot1") {
                
                m <- shared_data$data(withSelection = FALSE)
                
                # Intensity -> size mapping
                point_sizes <- if (isTRUE(intensity_enabled) &&
                                   !is.null(intensity_col) &&
                                   intensity_col %in% names(m)) {
                        vals <- suppressWarnings(as.numeric(m[[intensity_col]]))
                        rng  <- range(vals, na.rm = TRUE)
                        if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
                                rep(6, nrow(m))
                        } else {
                                3 + 7 * (vals - rng[1]) / (rng[2] - rng[1])  # 3..10 pt
                        }
                } else {
                        rep(6, nrow(m))
                }
                
                bg_opacity <- if (!is.null(series_sel) && length(series_sel) > 0) 0.05 else 0.8
                
                bg_marker <- list(
                        size    = point_sizes,
                        color   = "black",
                        opacity = bg_opacity,
                        line    = list(width = 0.3, color = "gray70")
                )
                
                x_formula <- as.formula(paste0("~", xvar))
                y_formula <- as.formula(paste0("~", yvar))
                
                fig <- plotly::plot_ly(
                        shared_data,
                        x      = x_formula,
                        y      = y_formula,
                        key    = ~.key,
                        source = source_id,
                        type   = "scattergl",
                        mode   = "markers",
                        marker = bg_marker,
                        name   = "All peaks",
                        showlegend = show_legend
                ) %>%
                        plotly::layout(
                                xaxis = list(title = xlabel),
                                yaxis = list(title = ylabel),
                                legend = list(
                                        orientation = "h",
                                        xanchor     = "center",
                                        x           = 0.5,
                                        yanchor     = "top",
                                        y           = -0.2
                                ),
                                margin = list(b = 90),
                                showlegend = show_legend
                        ) %>%
                        plotly::highlight(
                                on  = "plotly_selected",
                                off = "plotly_deselect",
                                color = I("red"),
                                selected = plotly::attrs_selected(name = "Selected")
                        )
                
                # overlay series
                if (!is.null(series_sel) && length(series_sel)) {
                        for (sid in series_sel) {
                                df_line <- m[!is.na(m$series_id) & m$series_id == sid, , drop = FALSE]
                                if (!nrow(df_line)) next
                                
                                df_line <- df_line[order(df_line$mz), ]
                                col <- series_cols[as.character(sid)]
                                
                                fig <- fig %>%
                                        plotly::add_trace(
                                                data  = df_line,
                                                x     = df_line[[xvar]],
                                                y     = df_line[[yvar]],
                                                type  = "scatter",
                                                mode  = "lines+markers",
                                                line  = list(color = col, width = 3),
                                                marker = list(
                                                        color   = col,
                                                        size    = 10,
                                                        opacity = 1,
                                                        line    = list(color = "white", width = 1)
                                                ),
                                                name       = paste0("Series ", sid),
                                                showlegend = show_legend
                                        )
                        }
                }
                
                fig
        }
        
        # ------------------------------------------
        # Unified selection: homol-table > plots
        # ------------------------------------------
        selected_keys <- reactive({
                
                # 1) If a homologue series is selected, it has priority
                if (!is.null(vals$series_keys) && length(vals$series_keys) > 0) {
                        return(vals$series_keys)
                }
                
                # 2) Otherwise, use plot-based selection
                ed1 <- plotly::event_data("plotly_selected", source = "plot1")
                ed2 <- plotly::event_data("plotly_selected", source = "plot2")
                
                keys <- c(
                        if (!is.null(ed1)) ed1$key else character(),
                        if (!is.null(ed2)) ed2$key else character()
                )
                
                unique(keys)
        })
        
        # -----------------
        # Main action with progress bar
        # -----------------
        observeEvent(input$go, {
                withProgress(message = "Running mzXplorer analysis...", value = 0, {
                        
                        # 1) Filter data
                        incProgress(0.1, detail = "Filtering data...")
                        m <- MD_data()
                        m <- m[
                                m$intensity >= input$slide1[1] &
                                        m$intensity <= input$slide1[2] &
                                        m$mz        >= input$slide2[1] &
                                        m$mz        <= input$slide2[2] &
                                        m$rt        >= input$slide3[1] &
                                        m$rt        <= input$slide3[2],
                        ]
                        
                        m_full <- m
                        
                        # 2) Homologue detection (if enabled)
                        if (isTRUE(input$run_homol)) {
                                incProgress(0.25, detail = "Detecting homologous series...")
                                
                                # safe CCS mode/tol
                                ccs_mode <- if (!has_ccs() || !isTRUE(input$use_ccs_toggle)) {
                                        "rt"
                                } else {
                                        if (is.null(input$ccs_mode)) "rt" else input$ccs_mode
                                }
                                
                                ccs_tol <- if (!has_ccs() || !isTRUE(input$use_ccs_toggle)) {
                                        0
                                } else {
                                        val <- input$homol_ccstol
                                        if (is.null(val) || is.na(val)) 0 else val
                                }
                                
                                m_homol <- find_homologues(
                                        df         = m,
                                        unit       = input$homol_unit,
                                        ppm        = input$homol_ppm,
                                        rttol      = input$homol_rttol,
                                        allow_gaps = input$homol_allow_gaps,
                                        min_length = input$homol_minlen,
                                        rt_trend   = input$homol_rttrend,
                                        R2_min     = input$homol_R2,
                                        verbose    = FALSE,
                                        ccs_mode   = ccs_mode,
                                        ccs_tol    = ccs_tol
                                )
                                
                                # Robust join via stable 'id' (avoids ambiguity for duplicate mz/rt)
                                if (!"id" %in% names(m_homol)) {
                                        stop("find_homologues must return an 'id' column for joining.")
                                }
                                
                                m <- dplyr::left_join(
                                        m_full,
                                        m_homol[, c("id", "series_id")],
                                        by = "id"
                                )
                        } else {
                                m <- m_full
                                m$series_id <- NA_integer_
                        }
                        
                        incProgress(0.2, detail = "Preparing shared data and plots...")
                        
                        # keys for crosstalk
                        m$.key <- sprintf("id%05d", seq_len(nrow(m)))
                        
                        # SharedData
                        d <- crosstalk::SharedData$new(m, key = ~.key, group = "md")
                        
                        vals$m <- m
                        vals$d <- d
                        vals$series_keys <- NULL
                        vals$series_sel  <- NULL
                        vals$series_cols <- NULL
                        
                        # 3) Plots
                        incProgress(0.15, detail = "Rendering plots...")
                        
                        output$DTPlot1 <- plotly::renderPlotly({
                                req(vals$d)
                                make_mdplot(
                                        shared_data       = vals$d,
                                        xvar              = input$xvar1,
                                        yvar              = input$yvar1,
                                        xlabel            = input$xvar1,
                                        ylabel            = input$yvar1,
                                        intensity_enabled = input$ins,
                                        intensity_col     = input$selectintensity,
                                        series_sel        = vals$series_sel,
                                        series_cols       = vals$series_cols,
                                        show_legend       = input$show_leg,
                                        source_id         = "plot1"
                                )
                        })
                        
                        output$DTPlot2 <- plotly::renderPlotly({
                                req(vals$d)
                                make_mdplot(
                                        shared_data       = vals$d,
                                        xvar              = input$xvar2,
                                        yvar              = input$yvar2,
                                        xlabel            = input$xvar2,
                                        ylabel            = input$yvar2,
                                        intensity_enabled = input$ins,
                                        intensity_col     = input$selectintensity,
                                        series_sel        = vals$series_sel,
                                        series_cols       = vals$series_cols,
                                        show_legend       = input$show_leg,
                                        source_id         = "plot2"
                                )
                        })
                        
                        # 4) Selected-data table
                        incProgress(0.15, detail = "Rendering tables...")
                        
                        output$table_selected <- DT::renderDT({
                                req(vals$m)
                                m <- vals$m
                                sel_keys <- selected_keys()
                                
                                sel_df <- if (length(sel_keys)) {
                                        m[m$.key %in% sel_keys, , drop = FALSE]
                                } else {
                                        m
                                }
                                
                                DT::datatable(
                                        sel_df,
                                        editable  = TRUE,
                                        rownames  = FALSE,
                                        selection = "none",
                                        filter    = "top",
                                        options   = list(scrollX = TRUE)
                                )
                        })
                        
                        # Barplot helper
                        nperc <- function(x) {
                                if (length(x) == 0 || all(is.na(x))) return(numeric())
                                round(x / max(x, na.rm = TRUE) * 100, 1)
                        }
                        
                        output$barplot <- plotly::renderPlotly({
                                req(vals$m)
                                m <- vals$m
                                sel_keys <- selected_keys()
                                
                                bar_out <- if (length(sel_keys)) {
                                        m[m$.key %in% sel_keys, , drop = FALSE]
                                } else {
                                        m[0, , drop = FALSE]
                                }
                                
                                if (!nrow(bar_out)) return(NULL)
                                
                                selvar    <- input$selectintensity
                                selectInt <- bar_out[[selvar]]
                                ytitle    <- paste("Relative", selvar, "(%)")
                                
                                plotly::plot_ly() %>%
                                        plotly::add_trace(
                                                data = bar_out,
                                                x    = ~mz,
                                                y    = ~nperc(selectInt),
                                                type = "bar"
                                        ) %>%
                                        plotly::layout(
                                                xaxis = list(title = "m/z"),
                                                yaxis = list(title = ytitle)
                                        )
                        })
                        
                        # 5) Export handler
                        incProgress(0.05, detail = "Setting up export...")
                        
                        output$x3 <- downloadHandler(
                                'MDplot_annotated_export.csv',
                                content = function(file) {
                                        req(vals$m)
                                        m   <- vals$m
                                        sel <- selected_keys()
                                        out <- if (length(sel)) {
                                                m[m$.key %in% sel, , drop = FALSE]
                                        } else m
                                        write.csv(out, file, row.names = FALSE)
                                }
                        )
                        
                        # 6) Homologue summary table with optional CCS summaries
                        incProgress(0.1, detail = "Summarising homologue series...")
                        
                        output$homol_table <- DT::renderDT({
                                req(isTRUE(input$run_homol))
                                req("series_id" %in% names(m))
                                ms <- m[!is.na(m$series_id), , drop = FALSE]
                                if (!nrow(ms)) return(NULL)
                                
                                base_summ <- ms %>%
                                        dplyr::group_by(series_id) %>%
                                        dplyr::summarise(
                                                n       = dplyr::n(),
                                                mz_min  = min(mz, na.rm = TRUE),
                                                mz_max  = max(mz, na.rm = TRUE),
                                                rt_min  = min(rt, na.rm = TRUE),
                                                rt_max  = max(rt, na.rm = TRUE),
                                                int_sum = sum(intensity, na.rm = TRUE),
                                                .groups = "drop"
                                        )
                                
                                if ("ccs" %in% names(ms)) {
                                        ccs_summ <- ms %>%
                                                dplyr::group_by(series_id) %>%
                                                dplyr::summarise(
                                                        ccs_min   = min(ccs, na.rm = TRUE),
                                                        ccs_max   = max(ccs, na.rm = TRUE),
                                                        ccs_range = round(ccs_max - ccs_min, 5),
                                                        .groups   = "drop"
                                                )
                                        
                                        summ <- dplyr::left_join(base_summ, ccs_summ, by = "series_id")
                                } else {
                                        summ <- base_summ
                                }
                                
                                summ <- summ %>%
                                        dplyr::arrange(dplyr::desc(n), mz_min)
                                
                                vals$summ <- summ
                                
                                DT::datatable(
                                        summ,
                                        rownames  = FALSE,
                                        selection = "multiple",
                                        options   = list(scrollX = TRUE, pageLength = 10)
                                )
                        })
                        
                        # Homologue table selection -> highlight & set priority keys
                        observeEvent(input$homol_table_rows_selected, {
                                req(vals$d, vals$m, vals$summ)
                                sel_rows <- input$homol_table_rows_selected
                                
                                if (!length(sel_rows)) {
                                        vals$d$selection(NULL)
                                        vals$series_keys <- NULL
                                        vals$series_sel  <- NULL
                                        vals$series_cols <- NULL
                                        return()
                                }
                                
                                sel_series <- vals$summ$series_id[sel_rows]
                                keys <- vals$m$.key[!is.na(vals$m$series_id) &
                                                            vals$m$series_id %in% sel_series]
                                cols <- series_palette(length(sel_series))
                                names(cols) <- as.character(sel_series)
                                
                                vals$series_keys <- keys
                                vals$series_sel  <- sel_series
                                vals$series_cols <- cols
                                
                                vals$d$selection(keys)   # for plot highlighting
                        }, ignoreInit = TRUE)
                        
                        # Clear selection
                        observeEvent(input$clear_homol_sel, {
                                req(vals$d)
                                vals$d$selection(NULL)
                                vals$series_keys <- NULL
                                vals$series_sel  <- NULL
                                vals$series_cols <- NULL
                        })
                        
                        incProgress(0.1, detail = "Done.")
                })
        })
        
        # Close the app when session completes (headless runs)
        if (!interactive()) {
                session$onSessionEnded(function() {
                        stopApp(); q("no")
                })
        }
}

shiny::shinyApp(ui, server)