# =========================================================
# mzXplorer with built-in HS detection + crosstalk-linked plots
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

source("playground_utils.R")


###########  UI  ###########

ui <- shiny::navbarPage(
  "mzXplorer: interactive mass defect plots",
  theme = shinythemes::shinytheme('spacelab'),
  shiny::tabPanel("DataAnalysis",
                  shiny::fluidPage(
                    shiny::sidebarLayout(
                      shiny::sidebarPanel(
                        shiny::fileInput('file1','Choose CSV File',
                                         accept = c('text/csv','text/comma-separated-values,text/plain','.csv')),
                        shiny::fluidRow(shiny::column(12, shiny::textInput("cus1", "MD formula 1", value = "CH2,O"))),
                        shiny::fluidRow(shiny::column(12, shiny::textInput("cus2", "MD formula 2", value = "Cl-H"))),
                        shiny::actionButton('go', 'Plot', width = "100%"),
                        shiny::tags$br(),
                        shiny::radioButtons("rounding","Rounding",
                                            choices = c("round","ceiling","floor"), selected = "round", inline = TRUE),
                        shiny::checkboxInput('ins', 'Show intensity/variable as size', FALSE),
                        shiny::checkboxInput("show_leg", "Show plot legends", TRUE),

                        # Homologue search controls
                        shiny::tags$hr(), shiny::h4("Homologue search"),
                        shiny::checkboxInput("run_homol", "Detect homologues in filtered data", FALSE),
                        shiny::textInput("homol_unit", "Repeating unit (formula)", value = "CH2"),
                        shiny::numericInput("homol_ppm", "m/z tolerance (ppm)", value = 5, min = 1, step = 1),
                        shiny::numericInput("homol_minlen", "Minimum series length", value = 4, min = 3, step = 1),
                        shiny::numericInput("homol_rttol", "Per-step RT tolerance (same units as rt); 0 = off",
                                            value = 0, min = 0, step = 0.1),
                        shiny::selectInput("homol_rttrend", "RT trend",
                                           choices = c("increasing","decreasing","any"), selected="increasing"),
                        shiny::checkboxInput("homol_allow_gaps", "Allow single gaps (k = 2)", FALSE),
                        shiny::numericInput("homol_R2", "Minimum spline R² (0 = off)", value = 0.98, min = 0, max = 1, step = 0.01),

                        shiny::uiOutput("slide1"),
                        shiny::uiOutput("slide2"),
                        shiny::uiOutput("slide3"),
                        width = 3
                      ),
                      shiny::mainPanel(
                          # --- Plot controls now in main panel ---
                          shiny::uiOutput("plotctr"),
                          #shiny::tags$br(),

                          # --- The two plots ---
                          shiny::uiOutput("plot"),
                          #shiny::tags$br(),

                          # --- Tables and barplots ---
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

  # Simple palette (no extra packages)
  series_palette <- function(n) {
    base <- c(
      "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
      "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
      "#bcbd22", "#17becf", "#a6cee3", "#fb9a99",
      "#33a02c", "#b15928", "#cab2d6", "#fdbf6f"
    )
    rep(base, length.out = n)
  }

  MD_data <- reactive({
          req(input$file1)
          df <- vroom::vroom(input$file1$datapath, show_col_types = FALSE)
          
          # Normalize names: lowercase + trim whitespace
          nms <- trimws(tolower(names(df)))
          names(df) <- nms
          
          if (!"mz" %in% names(df) && "m/z" %in% names(df)) df$mz <- df[["m/z"]]
          if (!"intensity" %in% names(df) && "dummy" %in% names(df)) df$intensity <- df$dummy
          if (!"rt" %in% names(df) && "retention-time" %in% names(df)) df$rt <- df[["retention-time"]]
          

    # Ensure required columns exist, but KEEP ALL others
    stopifnot(all(c("mz","intensity","rt") %in% names(df)))

    # Derived MD features
    df$RMD <- round((round(df$mz) - df$mz) / df$mz * 1e6)
    df$OMD <- round((round(df$mz) - df$mz) * 1e3)

    # Custom MD projections
    mdh1 <- getmdh(df$mz, cus = input$cus1, method = input$rounding)
    mdh2 <- getmdh(df$mz, cus = input$cus2, method = input$rounding)
    name1 <- paste0("Formula1_", colnames(mdh1))
    name2 <- paste0("Formula2_", colnames(mdh2))
    mdh <- cbind(mdh1[,-1], mdh2[,-1])
    colnames(mdh) <- c(name1[-1], name2[-1])

    cbind(df, mdh)
  })

  # sliders
  output$slide1 <- renderUI({
    df <- MD_data()
    sliderInput("slide1","Intensity range filter",
                min = min(df$intensity, na.rm=TRUE), max = max(df$intensity, na.rm=TRUE),
                value = c(min(df$intensity, na.rm=TRUE), max(df$intensity, na.rm=TRUE)))
  })
  output$slide2 <- renderUI({
    df <- MD_data()
    sliderInput("slide2","m/z range",
                min = min(df$mz, na.rm=TRUE), max = max(df$mz, na.rm=TRUE),
                value = c(min(df$mz, na.rm=TRUE), max(df$mz, na.rm=TRUE)))
  })
  output$slide3 <- renderUI({
    df <- MD_data()
    sliderInput("slide3","retention time range",
                min = min(df$rt, na.rm=TRUE), max = max(df$rt, na.rm=TRUE),
                value = c(min(df$rt, na.rm=TRUE), max(df$rt, na.rm=TRUE)))
  })

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
          column(6, selectInput('xvar1','X variable for Plot 1', choices = names(df),
                                selected = "rt")),
          column(6, selectInput('xvar2','X variable for Plot 2', choices = names(df),
                                selected = "RMD")),
          column(6, selectInput('yvar1','Y variable for Plot 1', choices = names(df),
                                selected = "mz")),
          column(6, selectInput('yvar2','Y variable for Plot 2', choices = names(df),
                                selected = "mz"))
      )
  })

  output$intensity_ctr <- renderUI({
      df <- MD_data()
      # Only numeric columns as choices 
      numeric_cols <- names(dplyr::select(df, where(is.numeric)))

      selectInput(
          'selectintensity',
          'Variable for intensity (for point size & barplot)',
          choices  = numeric_cols,
          selected = if ('intensity' %in% numeric_cols) 'intensity' else numeric_cols[1]
      )
  })

  # --- Shared plotting function (used by both plots)
  make_mdplot <- function(shared_data, m, xvar, yvar, xlabel, ylabel,
                          intensity_enabled, intensity_col,
                          series_sel, series_cols,
                          show_legend = TRUE) {

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

      # --- DYNAMIC BACKGROUND OPACITY ---
      # If any series are selected, dim background even more
      bg_opacity <- if (!is.null(series_sel) && length(series_sel) > 0) 0.05 else 0.8

      # Background cloud (always black)
      bg_marker <- list(
          size    = point_sizes,
          color   = "black",
          opacity = bg_opacity,                 # <--- dynamically set
          line    = list(width = 0.3, color = "gray70")
      )

      fig <- plotly::plot_ly(
          shared_data,
          x = m[[xvar]],
          y = m[[yvar]],
          type = "scattergl",
          mode = "markers",
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

      # --- OVERLAY HOMOLOGOUS SERIES (still on top, fully opaque) ---
      if (!is.null(series_sel) && length(series_sel)) {
          for (sid in series_sel) {
              # Vectorised & (not &&), keep data.frame
              df_line <- m[!is.na(m$series_id) & m$series_id == sid, , drop = FALSE]
              if (!nrow(df_line)) next

              df_line <- df_line[order(df_line$mz), ]
              col <- series_cols[as.character(sid)]

              fig <- fig %>%
                  plotly::add_trace(
                      data = df_line,
                      x    = df_line[[xvar]],
                      y    = df_line[[yvar]],
                      type = "scatter",               # SVG, draws on top
                      mode = "lines+markers",
                      line = list(color = col, width = 3),
                      marker = list(
                          color   = col,
                          size    = 10,
                          opacity = 1,                  # fully opaque
                          line    = list(color = "white", width = 1)
                      ),
                      name       = paste0("Series ", sid),
                      showlegend = show_legend
                  )
          }
      }

      fig
  }


  # Main action
  observeEvent(input$go, {
    m <- MD_data()
    m <- m[m$intensity >= input$slide1[1] &
             m$intensity <= input$slide1[2] &
             m$mz >= input$slide2[1] &
             m$mz <= input$slide2[2] &
             m$rt >= input$slide3[1] &
             m$rt <= input$slide3[2], ]

    # Keep a copy with ALL columns prior to homologue filtering
    m_full <- m

    if (isTRUE(input$run_homol)) {
      m_homol <- find_homologues(
        df = m,
        unit = input$homol_unit,
        ppm  = input$homol_ppm,
        min_length = input$homol_minlen,
        rttol = input$homol_rttol,
        rt_trend = input$homol_rttrend,
        allow_gaps = input$homol_allow_gaps,
        R2_min = input$homol_R2
      )
      # Merge series IDs back into the full dataset (key: mz + rt)
      m <- dplyr::left_join(
        m_full,
        m_homol[, c("mz","rt","series_id")],
        by = c("mz","rt")
      )
    } else {
      m <- m_full
      m$series_id <- NA_integer_
    }

    # keys for crosstalk
    m$.key <- sprintf("id%05d", seq_len(nrow(m)))

    # SharedData
    d <- crosstalk::SharedData$new(m, key = ~.key, group = "md")

    vals$m <- m
    vals$d <- d
    vals$series_keys <- NULL
    vals$series_sel  <- NULL
    vals$series_cols <- NULL

    # ---- Plot 1 ----
    output$DTPlot1 <- plotly::renderPlotly({
        req(vals$d, vals$m)
        make_mdplot(
            shared_data       = vals$d,
            m                 = vals$m,
            xvar              = input$xvar1,
            yvar              = input$yvar1,
            xlabel            = input$xvar1,
            ylabel            = input$yvar1,
            intensity_enabled = input$ins,
            intensity_col     = input$selectintensity,
            series_sel        = vals$series_sel,
            series_cols       = vals$series_cols,
            show_legend       = input$show_leg
        )
    })



    # ---- Plot 2 ----
    output$DTPlot2 <- plotly::renderPlotly({
        req(vals$d, vals$m)
        make_mdplot(
            shared_data       = vals$d,
            m                 = vals$m,
            xvar              = input$xvar2,
            yvar              = input$yvar2,
            xlabel            = input$xvar2,
            ylabel            = input$yvar2,
            intensity_enabled = input$ins,
            intensity_col     = input$selectintensity,
            series_sel        = vals$series_sel,
            series_cols       = vals$series_cols,
            show_legend       = input$show_leg
        )
    })

    # Selected rows table (shows ALL columns)
    output$table_selected <- DT::renderDT({
      T_out1 <- DT::datatable(
        m[vals$d$selection(), ],
        editable = TRUE, rownames = FALSE, selection = "none",
        filter = "top", options = list(scrollX = TRUE)
      )
      dt <- DT::datatable(
        m, editable = TRUE, rownames = FALSE, selection = "none",
        filter = "top", options = list(scrollX = TRUE)
      )
      if (NROW(T_out1) == 0) dt else T_out1
    })

    # Bar plot (on selection)
    nperc <- function(x) {
      if (length(x) == 0 || all(is.na(x))) return(numeric())
      round(x / max(x, na.rm = TRUE) * 100, 1)
    }
    output$barplot <- plotly::renderPlotly({
      bar_out  <- m[vals$d$selection(), ]
      selvar <- input$selectintensity
      selectInt <- bar_out[[selvar]]
      ytitle <- paste("Relative", selvar, "(%)")
      plotly::plot_ly() %>%
        plotly::add_trace(
          data = bar_out,
          x = bar_out$mz,
          y = nperc(selectInt),
          type = "bar"
        ) %>%
        plotly::layout(
          xaxis = list(title = "m/z"),
          yaxis = list(title = ytitle)
        )
    })

    # Export (ALL columns included)
    output$x3 <- downloadHandler(
      'MDplot_annotated_export.csv',
      content = function(file) {
        sel <- vals$d$selection()
        out <- if (length(sel)) m[m$.key %in% sel, , drop = FALSE] else m
        write.csv(out, file, row.names = FALSE)
      }
    )

    # Series summary table
    output$homol_table <- DT::renderDT({
      req(isTRUE(input$run_homol))
      req("series_id" %in% names(m))
      ms <- m[!is.na(m$series_id), , drop = FALSE]
      if (!nrow(ms)) return(NULL)
      summ <- ms %>%
        dplyr::group_by(series_id) %>%
        dplyr::summarise(
          n = dplyr::n(),
          mz_min = min(mz), mz_max = max(mz),
          rt_min = min(rt), rt_max = max(rt),
          int_sum = sum(intensity),
          .groups = "drop"
        ) %>%
        dplyr::arrange(desc(n), mz_min)

      vals$summ <- summ

      DT::datatable(
        summ,
        rownames = FALSE,
        selection = "multiple",
        options = list(scrollX = TRUE, pageLength = 10)
      )
    })

    # Table selection -> highlight & color series
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
      keys <- vals$m$.key[!is.na(vals$m$series_id) & vals$m$series_id %in% sel_series]
      cols <- series_palette(length(sel_series))
      names(cols) <- as.character(sel_series)

      vals$series_keys <- keys
      vals$series_sel  <- sel_series
      vals$series_cols <- cols

      vals$d$selection(keys)
    }, ignoreInit = TRUE)

    # Clear selection
    observeEvent(input$clear_homol_sel, {
      req(vals$d)
      vals$d$selection(NULL)
      vals$series_keys <- NULL
      vals$series_sel  <- NULL
      vals$series_cols <- NULL
    })
  })

  # Close the app when session completes (headless runs)
  if(!interactive()) {
    session$onSessionEnded(function() {
      stopApp(); q("no")
    })
  }
}

shiny::shinyApp(ui, server)
