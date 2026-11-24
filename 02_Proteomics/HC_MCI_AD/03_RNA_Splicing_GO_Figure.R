# ============================================================================

# AD11: GO BP DOWN RNA splicing figure (match supplied reference layout)

# Visualization only: NO changes to GSEA, limma, ranking or global BH FDR.

# Input: three original AD11 GSEA significant GO BP CSVs.

# ============================================================================

options(stringsAsFactors = FALSE)

packages <- c('readr','dplyr','tidyr','ggplot2','stringr','scales')

for (p in packages) if (!requireNamespace(p, quietly=TRUE))

  install.packages(p, repos='https://cloud.r-project.org')

# First look for the exact files in the AD11 output directory. If the user has

# copied the numbered CSVs to Desktop, set input_files below to their full paths.

folder <- file.path(path.expand(Sys.getenv("PROTEOME_DATA_DIR", unset = "~/Desktop")), "UDCiN_GO_BP_GSEA_MeanT_WithBatch_GlobalBH3Contrasts_nPerm10000_AD11")

input_files <- c(

  MCI_vs_HC = file.path(folder, '07_GSEA_GO_BP_MCI_vs_HC_FDR005.csv'),

  AD_vs_HC  = file.path(folder, '07_GSEA_GO_BP_AD_vs_HC_FDR005.csv'),

  AD_vs_MCI = file.path(folder, '07_GSEA_GO_BP_AD_vs_MCI_FDR005.csv')

)

if (!all(file.exists(input_files))) {

  stop('Missing AD11 input files:\n', paste(input_files[!file.exists(input_files)], collapse='\n'),

       '\nSet input_files to the exact AD11 CSV paths. Do not use older results.')

}

output_dir <- file.path(folder,'Splicing_7GO_MediumDots_TighterLayout')

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# Top-to-bottom order: identical to reference figure.

terms <- data.frame(

  ID=c('GO:0000377','GO:0000398','GO:0000375','GO:0006397',

       'GO:0008380','GO:0048024','GO:0050684'),

  Description=c(

    'RNA splicing, via transesterification reactions with bulged adenosine as nucleophile',

    'mRNA splicing, via spliceosome',

    'RNA splicing, via transesterification reactions',

    'mRNA processing',

    'RNA splicing',

    'regulation of mRNA splicing, via spliceosome',

    'regulation of mRNA processing'

  )

)

comparison_order <- c('MCI_vs_HC','AD_vs_HC','AD_vs_MCI')

comparison_labels <- c('MCI vs HC','AD vs HC','AD vs MCI')

read_comparison <- function(k) {

  z <- readr::read_csv(input_files[[k]], show_col_types=FALSE)

  required <- c('ID','NES','p.adjust','Direction')

  if (!all(required %in% names(z)))

    stop('Missing required columns in ', input_files[[k]])

  if (anyDuplicated(z$ID)) stop('Duplicated GO IDs in ', input_files[[k]])

  if ('p.adjust_global_3contrasts' %in% names(z) &&

      any(abs(z$p.adjust-z$p.adjust_global_3contrasts)>1e-10,na.rm=TRUE))

    stop('p.adjust is not the global BH value in ', input_files[[k]])

  z <- z[z$ID %in% terms$ID & z$Direction=='DOWN' &

           is.finite(z$NES) & z$NES<0 &

           is.finite(z$p.adjust) & z$p.adjust<0.05, ,drop=FALSE]

  z$Contrast <- k

  z

}

combined <- dplyr::bind_rows(lapply(comparison_order, read_comparison))

plot_data <- tidyr::expand_grid(ID=terms$ID,Contrast=comparison_order) |>

  dplyr::left_join(combined[,c('ID','Contrast','NES','p.adjust'),drop=FALSE],

                   by=c('ID','Contrast')) |>

  dplyr::left_join(terms,by='ID') |>

  dplyr::mutate(

    Comparison=factor(Contrast,levels=comparison_order,labels=comparison_labels),

    GO=factor(Description,levels=rev(terms$Description)),

    X=as.numeric(Comparison),

    AbsNES=abs(NES),

    Log10FDR=-log10(pmax(p.adjust,.Machine$double.xmin)),

    NESLabel=ifelse(is.finite(NES),sprintf('%.2f',NES),NA_character_)

  )

if (any(!is.finite(plot_data$NES)))

  stop('Not all seven selected GO terms are significant DOWN in all three comparisons.')

readr::write_csv(plot_data,file.path(output_dir,'Splicing_7GO_Plot_Data.csv'))

# Reference layout: 3 vertical columns; very light vertical gridlines;

# compact text and points; blue color gradient; right-side stacked legends.

wrapped <- stats::setNames(

  stringr::str_wrap(levels(plot_data$GO),width=37),

  levels(plot_data$GO)

)

p <- ggplot2::ggplot(plot_data,ggplot2::aes(x=X,y=GO)) +

  ggplot2::geom_point(ggplot2::aes(size=AbsNES,color=Log10FDR),alpha=0.98) +

  ggplot2::geom_text(ggplot2::aes(x=X+0.13,label=NESLabel),

                     hjust=0,size=2.8,color='#555555',na.rm=TRUE) +

  ggplot2::scale_x_continuous(

    breaks=1:3,labels=comparison_labels,limits=c(0.75,3.38),

    expand=ggplot2::expansion(mult=0)) +

  ggplot2::scale_y_discrete(labels=wrapped,drop=FALSE,

                             expand=ggplot2::expansion(add=c(0.45,0.45))) +

  ggplot2::scale_size_continuous(name='|NES|',range=c(4.6,9.5),

    breaks=c(1.5,2.0,2.5),limits=c(1.4,3.0),

    guide=ggplot2::guide_legend(order=1,override.aes=list(color='#80A7D8'))) +

  ggplot2::scale_color_gradient(name=expression(-log[10]('adjusted P')),

    low='#E2EAF8',high='#16458F',

    guide=ggplot2::guide_colorbar(order=2,barheight=grid::unit(18,'mm'),

                                  barwidth=grid::unit(5,'mm'))) +

  ggplot2::labs(x=NULL,y=NULL) +

  ggplot2::theme_minimal(base_size=10) +

  ggplot2::theme(

    panel.grid.major.y=ggplot2::element_blank(),

    panel.grid.minor=ggplot2::element_blank(),

    panel.grid.major.x=ggplot2::element_line(color='#E9EBEF',linewidth=0.38),

    axis.text.x=ggplot2::element_text(color='#202020',face='bold',size=9.5,

                                      margin=ggplot2::margin(t=8)),

    axis.text.y=ggplot2::element_text(color='#202020',size=8.0,

                                      lineheight=0.93,margin=ggplot2::margin(r=9)),

    axis.ticks=ggplot2::element_blank(),

    legend.position='right',

    legend.box='vertical',

    legend.title=ggplot2::element_text(size=8.8,color='#333333'),

    legend.text=ggplot2::element_text(size=8,color='#555555'),

    legend.key.height=grid::unit(5.2,'mm'),

    legend.spacing.y=grid::unit(3,'mm'),

    plot.margin=ggplot2::margin(12,12,9,10),

    panel.background=ggplot2::element_rect(fill='white',color=NA),

    plot.background=ggplot2::element_rect(fill='white',color=NA)

  )

print(p)

# More compact panel in both directions, with medium-sized data points.

# Smaller physical panel: 7.7 x 4.15 inches, retaining the compact spacing.

ggplot2::ggsave(file.path(output_dir,'Splicing_7GO_MediumDots_TighterLayout.png'),p,

                width=7.7,height=4.15,dpi=600,bg='white')

ggplot2::ggsave(file.path(output_dir,'Splicing_7GO_MediumDots_TighterLayout.pdf'),p,

                width=7.7,height=4.15,device='pdf',bg='white')

message('Saved to: ',output_dir)
